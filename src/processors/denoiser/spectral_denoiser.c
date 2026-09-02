/*
libspecbleach - A spectral processing library

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/denoiser_post_process.h"
#include "shared/denoiser_logic/core/denoiser_profile_core.h"
#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/denoiser_logic/processing/gain_calculator.h"
#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/denoiser_logic/processing/nlm_filter.h"
#include "shared/denoiser_logic/processing/suppression_engine.h"
#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/stft/stft_processor.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/spectral_circular_buffer.h"
#include "shared/utils/spectral_features.h"
#include "shared/utils/spectral_smoother.h"
#include "shared/utils/spectral_utils.h"
#include "shared/utils/transient_detector.h"
#include <math.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

/**
 * Unified spectral denoiser.
 *
 * Both smoothing strategies (1D temporal/spatial gain smoothing and 2D
 * Non-Local Means patch smoothing) share the same chassis: STFT analysis,
 * noise estimation/profile, transient detection, suppression, masking veto,
 * gain calculation and post-processing. Both smoothing sub-modules are fully
 * allocated at initialization so switching modes at runtime is allocation-
 * free. Both modes share a common algorithmic delay (the NLM look-ahead) so
 * reported latency never changes on a mode switch.
 */
typedef struct SbSpectralDenoiser {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t hop;

  DenoiserParameters parameters;

  float* snr_frame;                 // Current SNR frame for NLM input
  float* smoothed_snr;              // Smoothed SNR output from NLM
  float* gain_spectrum;             // Gain spectrum of the active chain
  float* gain_spectrum_b;           // Gain spectrum of the incoming chain
                                    // (transition crossfade only)
  float* noise_spectrum;            // Copy of noise profile for processing
  float* noise_spectrum_buffers[2]; // Double-buffered noise spectrum for
                                    // lock-free SPSC publication
  atomic_int
      active_noise_idx; // Index of current published noise spectrum (0 or 1)
  float* alpha;         // Oversubtraction factors (active chain)
  float* beta;          // Undersubtraction factors (active chain)
  float* alpha_b;       // Oversubtraction factors (transition chain)
  float* beta_b;        // Undersubtraction factors (transition chain)
  float* manual_noise_floor; // Manual profile floor
  TonalReducer* tonal_reducer;

  // Reusable circular buffer for aligned temporal analysis (common delay)
  SbSpectralCircularBuffer* circular_buffer;
  uint32_t layer_fft;
  uint32_t layer_noise;
  uint32_t layer_smoothed;

  SpectrumType spectrum_type;
  GainCalculationType gain_calculation_type;

  NoiseProfile* noise_profile;
  NoiseEstimator* noise_estimator;
  AdaptiveNoiseEstimator* adaptive_estimator;
  NlmFilter* nlm_filter;
  SpectralSmoother* spectrum_smoothing;
  SpectralFeatures* spectral_features;
  MaskingVeto* masking_veto;
  SuppressionEngine* suppression_engine;
  NoiseFloorManager* noise_floor_manager;
  CriticalBands* critical_bands;
  TransientDetector* transient_detector;

  float* band_energies;
  float* onset_weights;
  float* transient_mask;
  float* clean_magnitude;
  float* smoothed_magnitude; // Temporal pre-subtraction smoothed magnitude
  bool is_transient_detected;
  float transient_intensity;

  // Smoothing mode state (written by load_parameters, read by process; the
  // load/process concurrency contract forbids concurrent calls on the same
  // instance, so plain fields suffice)
  int active_mode;   // Currently rendering smoothing mode
  int pending_mode;  // Target mode during a transition
  int previous_mode; // Mode fading out during a transition
  bool in_transition;
  uint32_t transition_frames;
  uint32_t transition_pos;

  int last_adaptive_state;
  int last_noise_estimation_method;
  float aggressiveness;
  bool was_learning;
} SbSpectralDenoiser;

static int normalize_smoothing_mode(const int mode) {
  return (mode == SB_SMOOTHING_NLM_2D) ? SB_SMOOTHING_NLM_2D
                                       : SB_SMOOTHING_TEMPORAL;
}

/**
 * Bypass alignment: keeps the circular buffer, NLM history and output
 * rolling during idle/silence stretches so idle→active transitions stay
 * aligned and the reported look-ahead latency is preserved. Emits the frame
 * delayed by the configured NLM latency.
 */
static void align_bypass_frame(SbSpectralDenoiser* self, float* fft_spectrum,
                               const float* reference_spectrum) {
  spectral_circular_buffer_push(self->circular_buffer, self->layer_fft,
                                fft_spectrum);
  spectral_circular_buffer_push(self->circular_buffer, self->layer_noise,
                                self->noise_spectrum);
  spectral_circular_buffer_push(self->circular_buffer, self->layer_smoothed,
                                reference_spectrum);
  nlm_filter_calculate_snr(self->nlm_filter, reference_spectrum,
                           self->noise_spectrum, self->snr_frame);
  nlm_filter_push_frame(self->nlm_filter, self->snr_frame);
  const float* delayed_spectrum = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_fft,
      nlm_filter_get_latency_frames(self->nlm_filter));
  if (delayed_spectrum) {
    memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));
  }
  spectral_circular_buffer_advance(self->circular_buffer);
}

static bool run_nlm_chain(SbSpectralDenoiser* self, float* fft_spectrum,
                          const float* smoothed_magnitude,
                          const float* delayed_noise, float* gain_out,
                          float* alpha, float* beta);

static void run_temporal_chain(SbSpectralDenoiser* self,
                               const float* delayed_fft,
                               const float* delayed_noise, float* gain_out,
                               float* alpha, float* beta);

SpectralProcessorHandle spectral_denoiser_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, NoiseProfile* noise_profile) {

  if (!noise_profile || sample_rate == 0 || fft_size == 0 ||
      overlap_factor == 0) {
    return NULL;
  }

  SbSpectralDenoiser* self =
      (SbSpectralDenoiser*)calloc(1U, sizeof(SbSpectralDenoiser));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = (self->fft_size / 2U) + 1U;
  self->hop = self->fft_size / overlap_factor;
  self->sample_rate = sample_rate;
  self->spectrum_type = SPECTRAL_TYPE;
  self->gain_calculation_type = GAIN_ESTIMATION_TYPE;
  self->noise_profile = noise_profile;
  self->active_mode = SB_SMOOTHING_TEMPORAL;
  self->pending_mode = SB_SMOOTHING_TEMPORAL;
  self->previous_mode = SB_SMOOTHING_TEMPORAL;

  self->snr_frame = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->smoothed_snr = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->gain_spectrum = (float*)calloc(self->fft_size, sizeof(float));
  self->gain_spectrum_b = (float*)calloc(self->fft_size, sizeof(float));
  self->noise_spectrum =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noise_spectrum_buffers[0] =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noise_spectrum_buffers[1] =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  atomic_init(&self->active_noise_idx, 0);
  self->alpha = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->beta = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->alpha_b = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->beta_b = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->manual_noise_floor =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->smoothed_magnitude =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->clean_magnitude =
      (float*)calloc(self->real_spectrum_size, sizeof(float));

  if (!self->snr_frame || !self->smoothed_snr || !self->gain_spectrum ||
      !self->gain_spectrum_b || !self->noise_spectrum ||
      !self->noise_spectrum_buffers[0] || !self->noise_spectrum_buffers[1] ||
      !self->alpha || !self->beta || !self->alpha_b || !self->beta_b ||
      !self->manual_noise_floor || !self->smoothed_magnitude ||
      !self->clean_magnitude) {
    spectral_denoiser_free(self);
    return NULL;
  }

  (void)initialize_spectrum_with_value(self->gain_spectrum, self->fft_size,
                                       1.0F);
  (void)initialize_spectrum_with_value(self->gain_spectrum_b, self->fft_size,
                                       1.0F);
  (void)initialize_spectrum_with_value(self->alpha, self->real_spectrum_size,
                                       1.F);
  (void)initialize_spectrum_with_value(self->alpha_b, self->real_spectrum_size,
                                       1.F);

  // Initialize tonal reducer
  self->tonal_reducer = tonal_reducer_initialize(
      self->real_spectrum_size, self->sample_rate, self->fft_size);
  if (!self->tonal_reducer) {
    spectral_denoiser_free(self);
    return NULL;
  }

  // Circular buffer for temporal alignment (provides the common delay)
  self->circular_buffer = spectral_circular_buffer_create(DELAY_BUFFER_FRAMES);
  if (!self->circular_buffer) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->layer_fft =
      spectral_circular_buffer_add_layer(self->circular_buffer, self->fft_size);
  self->layer_noise = spectral_circular_buffer_add_layer(
      self->circular_buffer, self->real_spectrum_size);
  self->layer_smoothed = spectral_circular_buffer_add_layer(
      self->circular_buffer, self->real_spectrum_size);

  if (self->layer_fft == 0xFFFFFFFFU || self->layer_noise == 0xFFFFFFFFU ||
      self->layer_smoothed == 0xFFFFFFFFU) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->was_learning = false;
  self->aggressiveness = 0.0f;
  self->parameters.tonal_reduction = 0.0f;

  // Initialize noise estimator for learning mode
  self->noise_estimator = noise_estimation_initialize(fft_size, noise_profile);
  if (!self->noise_estimator) {
    spectral_denoiser_free(self);
    return NULL;
  }

  // NLM filter (2D smoothing strategy)
  NlmFilterConfig nlm_config = {
      .spectrum_size = self->real_spectrum_size,
      .time_buffer_size =
          NLM_SEARCH_RANGE_TIME_PAST + NLM_SEARCH_RANGE_TIME_FUTURE + 1,
      .patch_size = NLM_PATCH_SIZE,
      .paste_block_size = NLM_PASTE_BLOCK_SIZE,
      .search_range_freq = NLM_SEARCH_RANGE_FREQ,
      .search_range_time_past = NLM_SEARCH_RANGE_TIME_PAST,
      .search_range_time_future = NLM_SEARCH_RANGE_TIME_FUTURE,
      .h_parameter = NLM_DEFAULT_H_PARAMETER,
      .distance_threshold = 0.0F, // Use default (4 * h²)
  };
  self->nlm_filter = nlm_filter_initialize(nlm_config);
  if (!self->nlm_filter) {
    spectral_denoiser_free(self);
    return NULL;
  }

  // Temporal smoother (1D smoothing strategy)
  self->spectrum_smoothing = spectral_smoothing_initialize(
      self->fft_size, self->sample_rate, overlap_factor, FIXED);
  if (!self->spectrum_smoothing) {
    spectral_denoiser_free(self);
    return NULL;
  }

  // Initialize spectral features
  self->spectral_features =
      spectral_features_initialize(self->real_spectrum_size);
  if (!self->spectral_features) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->masking_veto = masking_veto_initialize(
      self->fft_size, self->sample_rate, CRITICAL_BANDS_TYPE,
      self->spectrum_type, false, USE_TEMPORAL_MASKING_DEFAULT);
  self->suppression_engine = suppression_engine_initialize(
      self->real_spectrum_size, self->sample_rate, CRITICAL_BANDS_TYPE,
      self->spectrum_type, true, USE_TEMPORAL_MASKING_DEFAULT);

  if (!self->masking_veto || !self->suppression_engine) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->noise_floor_manager = noise_floor_manager_initialize(fft_size);

  self->critical_bands = critical_bands_initialize(
      self->sample_rate, self->fft_size, CRITICAL_BANDS_TYPE);
  uint32_t num_bands = self->critical_bands
                           ? get_number_of_critical_bands(self->critical_bands)
                           : 0U;
  self->transient_detector = transient_detector_initialize(num_bands);

  self->band_energies =
      (float*)calloc(num_bands > 0U ? num_bands : 1U, sizeof(float));
  self->onset_weights =
      (float*)calloc(num_bands > 0U ? num_bands : 1U, sizeof(float));
  self->transient_mask =
      (float*)calloc(self->real_spectrum_size, sizeof(float));

  if (!self->noise_floor_manager || !self->critical_bands ||
      !self->transient_detector || !self->band_energies ||
      !self->onset_weights || !self->transient_mask) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->transition_frames = (uint32_t)fmaxf(
      1.0F, (SMOOTHING_TRANSITION_SECONDS * (float)self->sample_rate) /
                (float)self->hop);

  return self;
}

void spectral_denoiser_free(SpectralProcessorHandle instance) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  if (!self) {
    return;
  }

  if (self->noise_estimator) {
    noise_estimation_free(self->noise_estimator);
  }
  if (self->adaptive_estimator) {
    adaptive_estimator_free(self->adaptive_estimator);
  }
  if (self->nlm_filter) {
    nlm_filter_free(self->nlm_filter);
  }
  if (self->spectrum_smoothing) {
    spectral_smoothing_free(self->spectrum_smoothing);
  }
  if (self->spectral_features) {
    spectral_features_free(self->spectral_features);
  }
  if (self->masking_veto) {
    masking_veto_free(self->masking_veto);
  }
  if (self->suppression_engine) {
    suppression_engine_free(self->suppression_engine);
  }
  if (self->noise_floor_manager) {
    noise_floor_manager_free(self->noise_floor_manager);
  }

  if (self->critical_bands) {
    critical_bands_free(self->critical_bands);
  }
  if (self->transient_detector) {
    transient_detector_free(self->transient_detector);
  }
  if (self->band_energies) {
    free(self->band_energies);
  }
  if (self->onset_weights) {
    free(self->onset_weights);
  }
  if (self->transient_mask) {
    free(self->transient_mask);
  }
  if (self->clean_magnitude) {
    free(self->clean_magnitude);
  }
  if (self->smoothed_magnitude) {
    free(self->smoothed_magnitude);
  }

  free(self->snr_frame);
  free(self->smoothed_snr);
  free(self->gain_spectrum);
  free(self->gain_spectrum_b);
  free(self->noise_spectrum);
  if (self->noise_spectrum_buffers[0]) {
    free(self->noise_spectrum_buffers[0]);
  }
  if (self->noise_spectrum_buffers[1]) {
    free(self->noise_spectrum_buffers[1]);
  }
  free(self->alpha);
  free(self->beta);
  free(self->alpha_b);
  free(self->beta_b);
  if (self->manual_noise_floor) {
    free(self->manual_noise_floor);
  }

  if (self->circular_buffer) {
    spectral_circular_buffer_free(self->circular_buffer);
  }

  if (self->tonal_reducer) {
    tonal_reducer_free(self->tonal_reducer);
  }

  free(self);
}

bool load_reduction_parameters(SpectralProcessorHandle instance,
                               DenoiserParameters parameters) {
  if (!instance) {
    return false;
  }

  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  // Check if we need to initialize or re-initialize the adaptive estimator
  if (parameters.adaptive_noise) {
    AdaptiveNoiseEstimationMethod requested_method =
        (AdaptiveNoiseEstimationMethod)parameters.noise_estimation_method;

    bool needs_init = !self->adaptive_estimator ||
                      adaptive_estimator_get_method(self->adaptive_estimator) !=
                          requested_method;

    if (needs_init) {
      adaptive_estimator_free(self->adaptive_estimator);
      self->adaptive_estimator = adaptive_estimator_initialize(
          self->real_spectrum_size, self->sample_rate, self->fft_size,
          requested_method);
      self->last_adaptive_state = 0;
    }
  }

  self->parameters = parameters;

  // Runtime smoothing mode switching (allocation-free): the outgoing mode is
  // crossfaded against the incoming one over SMOOTHING_TRANSITION_SECONDS
  if (!self->in_transition) {
    const int requested = normalize_smoothing_mode(parameters.smoothing_mode);
    if (requested != self->active_mode) {
      self->previous_mode = self->active_mode;
      self->pending_mode = requested;
      self->transition_pos = 0U;
      self->in_transition = true;
    }
  } else {
    const int requested = normalize_smoothing_mode(parameters.smoothing_mode);
    if (requested != self->pending_mode && requested == self->previous_mode) {
      // Reverting to the mode that is fading out: mirror the in-progress
      // crossfade around its current blend point so the transition reverses
      // smoothly toward the original chain without a gain discontinuity
      const int outgoing = self->pending_mode;
      self->pending_mode = self->previous_mode;
      self->previous_mode = outgoing;
      self->transition_pos = self->transition_frames - self->transition_pos;
    } else {
      self->pending_mode = requested;
    }
  }

  // Update NLM h parameter based on smoothing factor
  if (self->nlm_filter) {
    nlm_filter_set_h_parameter(
        self->nlm_filter, (parameters.smoothing_factor > 0.0F)
                              ? (0.5F + (parameters.smoothing_factor * 4.5F))
                              : 0.0F);
  }

  return true;
}

bool spectral_denoiser_run(SpectralProcessorHandle instance,
                           float* fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  // 1. Preparation: Get reference spectrum and handle learning mode
  float* reference_spectrum =
      get_spectral_feature(self->spectral_features, fft_spectrum,
                           self->fft_size, self->spectrum_type);

  if (denoiser_profile_core_handle_learning_mode(
          self->noise_estimator, reference_spectrum,
          self->parameters.learn_noise, &self->was_learning)) {
    return true;
  }

  // 2. Noise Estimation: Update noise profile (Adaptive or Manual)
  DenoiserProfileCoreParams profile_params = {
      .adaptive_enabled = self->parameters.adaptive_noise,
      .spectrum_size = self->real_spectrum_size,
      .aggressiveness = &self->aggressiveness,
      .param_aggressiveness = self->parameters.aggressiveness,
      .last_adaptive_state = &self->last_adaptive_state,
      .adaptive_estimator = self->adaptive_estimator,
      .noise_profile = self->noise_profile,
      .manual_noise_floor = self->manual_noise_floor,
      .noise_spectrum = self->noise_spectrum,
      .noise_estimator = self->noise_estimator,
      .noise_profile_offset_linear =
          self->parameters.noise_profile_offset_linear,
      .tonal_noise_profile_offset_linear =
          self->parameters.tonal_noise_profile_offset_linear,
      // Previous frame's tonal mask (one-frame latency) is used to apply the
      // tonal threshold offset at detected tonal bins
      .tonal_mask = tonal_reducer_get_mask(self->tonal_reducer),
  };
  denoiser_profile_core_update(profile_params, reference_spectrum);

  // Idle bypass: no manual profile and not adaptive → skip the heavy chain.
  // Preserve latency and buffer state: push current frame, output the frame
  // delayed by the common lookahead, and advance the circular buffer so
  // idle→active transitions stay aligned.
  if (!self->parameters.adaptive_noise &&
      !is_noise_estimation_available(self->noise_profile, ROLLING_MEAN) &&
      !is_noise_estimation_available(self->noise_profile, MEDIAN) &&
      !is_noise_estimation_available(self->noise_profile, STD_DEV) &&
      !is_noise_estimation_available(self->noise_profile, CV_MASK)) {
    align_bypass_frame(self, fft_spectrum, reference_spectrum);
    return true;
  }

  // Silence bypass: input essentially silent — profile update already done
  // so aggressiveness/threshold stay responsive, but skip the heavy chain.
  {
    float max_val = 0.0f;
    for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
      if (reference_spectrum[k] > max_val) {
        max_val = reference_spectrum[k];
      }
      if (max_val > 1e-12F) {
        break;
      }
    }
    if (max_val <= 1e-12F) {
      // Keep circular buffer aligned as in idle path
      align_bypass_frame(self, fft_spectrum, reference_spectrum);

      // Safely publish noise spectrum to inactive double buffer via SPSC
      // atomic release
      int published_idx =
          atomic_load_explicit(&self->active_noise_idx, memory_order_relaxed);
      int write_noise_idx = 1 - published_idx;
      memcpy(self->noise_spectrum_buffers[write_noise_idx],
             self->noise_spectrum, self->real_spectrum_size * sizeof(float));
      atomic_store_explicit(&self->active_noise_idx, write_noise_idx,
                            memory_order_release);

      return true;
    }
  }

  // 2.1 Transient Detection via Transient Detector across Critical Bands
  // Transient detection runs on a clean signal estimate with scaled-up noise
  // subtraction to avoid false triggering from residual musical noise.
  bool transient_enabled = (self->parameters.hpss_enable != 0);
  if (transient_enabled && self->critical_bands && self->transient_detector) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      // Scale noise up using TRANSIENT_CLEAN_NOISE_SCALE to eliminate spurious
      // noise peaks
      float clean = fmaxf(reference_spectrum[k] - (TRANSIENT_CLEAN_NOISE_SCALE *
                                                   self->noise_spectrum[k]),
                          0.0F);
      self->clean_magnitude[k] = clean;
    }

    compute_critical_bands_spectrum(self->critical_bands, self->clean_magnitude,
                                    self->band_energies);
    self->is_transient_detected = transient_detector_process(
        self->transient_detector, self->band_energies, self->onset_weights,
        &self->transient_intensity);

    // Expand critical band onset weights to per-bin transient mask
    uint32_t num_bands = get_number_of_critical_bands(self->critical_bands);
    memset(self->transient_mask, 0, self->real_spectrum_size * sizeof(float));
    for (uint32_t b = 0; b < num_bands; ++b) {
      float bw = self->onset_weights[b];
      if (bw > 0.0F) {
        CriticalBandIndexes idx = get_band_indexes(self->critical_bands, b);
        uint32_t end = (idx.end_position < self->real_spectrum_size)
                           ? idx.end_position
                           : self->real_spectrum_size;
        for (uint32_t k = idx.start_position; k < end; ++k) {
          self->transient_mask[k] = fmaxf(self->transient_mask[k], bw);
        }
      }
    }
  } else {
    self->is_transient_detected = false;
    self->transient_intensity = 0.0f;
    memset(self->transient_mask, 0, self->real_spectrum_size * sizeof(float));
  }

  // 2.2 Align internal state and output to the common delayed frame
  spectral_circular_buffer_push(self->circular_buffer, self->layer_fft,
                                fft_spectrum);
  spectral_circular_buffer_push(self->circular_buffer, self->layer_noise,
                                self->noise_spectrum);

  // Compute SNR for NLM using CURRENT noise and push frame. This keeps the
  // NLM history rolling even in temporal mode so a runtime mode switch is
  // seamless and allocation-free.
  nlm_filter_calculate_snr(self->nlm_filter, reference_spectrum,
                           self->noise_spectrum, self->snr_frame);
  nlm_filter_push_frame(self->nlm_filter, self->snr_frame);

  const uint32_t nlm_delay = nlm_filter_get_latency_frames(self->nlm_filter);

  const float* nlm_intermediate_noise = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_noise, nlm_delay);
  if (!nlm_intermediate_noise) {
    nlm_intermediate_noise = self->noise_spectrum;
  }

  // NLM smoothing (runs when NLM is the active or the incoming mode). The
  // smoothed magnitude is captured explicitly so the temporal chain cannot
  // overwrite the shared alignment layer before the NLM chain consumes it.
  const bool nlm_needed =
      (self->active_mode == SB_SMOOTHING_NLM_2D) ||
      (self->in_transition && self->pending_mode == SB_SMOOTHING_NLM_2D);
  const float* nlm_smoothed = NULL;
  if (nlm_needed && nlm_filter_process(self->nlm_filter, self->smoothed_snr)) {
    nlm_filter_reconstruct_magnitude(self->nlm_filter, self->smoothed_snr,
                                     nlm_intermediate_noise, self->snr_frame);
    spectral_circular_buffer_push(self->circular_buffer, self->layer_smoothed,
                                  self->snr_frame);
    nlm_smoothed = self->snr_frame;
  }

  // Retrieve unified aligned frames at the common delay
  const float* delayed_spectrum = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_fft, nlm_delay);
  const float* delayed_noise = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_noise, nlm_delay);

  if (!delayed_spectrum) {
    delayed_spectrum = fft_spectrum;
  }
  if (!delayed_noise) {
    delayed_noise = self->noise_spectrum;
  }

  // Align output to delayed frame for post-processing
  if (delayed_spectrum != fft_spectrum) {
    memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));
  }

  // 3. Denoising Stage: dispatch the active smoothing strategy (or crossfade
  // both during a runtime mode transition)
  float* gain_a = self->gain_spectrum;
  float* gain_b = self->gain_spectrum_b;

  if (self->in_transition) {
    const float total = (float)self->transition_frames;
    const float w = (float)self->transition_pos / total; // 0 → 1

    if (self->previous_mode == SB_SMOOTHING_NLM_2D) {
      (void)run_nlm_chain(self, fft_spectrum, nlm_smoothed, delayed_noise,
                          gain_a, self->alpha, self->beta);
      run_temporal_chain(self, fft_spectrum, delayed_noise, gain_b,
                         self->alpha_b, self->beta_b);
    } else {
      run_temporal_chain(self, fft_spectrum, delayed_noise, gain_a, self->alpha,
                         self->beta);
      (void)run_nlm_chain(self, fft_spectrum, nlm_smoothed, delayed_noise,
                          gain_b, self->alpha_b, self->beta_b);
    }

    for (uint32_t k = 0U; k < self->fft_size; ++k) {
      gain_a[k] = (gain_a[k] * (1.0F - w)) + (gain_b[k] * w);
    }

    self->transition_pos++;
    if (self->transition_pos >= self->transition_frames) {
      self->active_mode = self->pending_mode;
      self->in_transition = false;
    }
  } else if (self->active_mode == SB_SMOOTHING_NLM_2D) {
    (void)run_nlm_chain(self, fft_spectrum, nlm_smoothed, delayed_noise, gain_a,
                        self->alpha, self->beta);
  } else {
    run_temporal_chain(self, fft_spectrum, delayed_noise, gain_a, self->alpha,
                       self->beta);
  }

  // 4. Post-Processing: Final gain management and mixing
  DenoiserPostProcessParams post_params = {
      .fft_size = self->fft_size,
      .real_spectrum_size = self->real_spectrum_size,
      .reduction_amount = self->parameters.reduction_amount,
      .tonal_reduction = self->parameters.tonal_reduction,
      .whitening_factor = self->parameters.whitening_factor,
      .residual_listen = self->parameters.residual_listen,
      .noise_floor_manager = self->noise_floor_manager,
      .tonal_reducer = self->tonal_reducer,
      .gain_spectrum = self->gain_spectrum,
      .noise_spectrum = delayed_noise,
      .fft_spectrum = fft_spectrum,
      .reduction_curve_bias = self->parameters.reduction_curve_bias,
  };

  denoiser_post_process_apply(post_params);

  // Finalize: Advance circular buffer write index
  spectral_circular_buffer_advance(self->circular_buffer);

  // Safely publish noise spectrum to inactive double buffer via SPSC atomic
  // release
  int published_idx =
      atomic_load_explicit(&self->active_noise_idx, memory_order_relaxed);
  int write_noise_idx = 1 - published_idx;
  memcpy(self->noise_spectrum_buffers[write_noise_idx], self->noise_spectrum,
         self->real_spectrum_size * sizeof(float));
  atomic_store_explicit(&self->active_noise_idx, write_noise_idx,
                        memory_order_release);

  return true;
}

/**
 * 2D Non-Local Means chain: NLM-smoothed magnitude feeds suppression, masking
 * veto, gain calculation. The smoothed magnitude produced by this frame's NLM
 * pass is passed in explicitly; when the NLM pass did not run or failed, the
 * last smoothed magnitude in the alignment buffer (or the delayed frame) is
 * used as an explicit fallback.
 */
static bool run_nlm_chain(SbSpectralDenoiser* self, float* fft_spectrum,
                          const float* smoothed_magnitude,
                          const float* delayed_noise, float* gain_out,
                          float* alpha, float* beta) {
  if (!smoothed_magnitude) {
    smoothed_magnitude = spectral_circular_buffer_retrieve(
        self->circular_buffer, self->layer_smoothed, 0U);
  }
  if (!smoothed_magnitude) {
    smoothed_magnitude = fft_spectrum;
  }

  // 3.1 Calculate SNR-dependent oversubtraction factors (Alpha/Beta)
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->suppression_engine, smoothed_magnitude,
                               delayed_noise, suppression_params, alpha, beta);

  // 3.2 Detect tonal components and boost alpha at tonal bins
  tonal_reducer_run(self->tonal_reducer, delayed_noise,
                    get_noise_profile(self->noise_profile, CV_MASK),
                    is_noise_estimation_available(self->noise_profile, CV_MASK),
                    alpha, self->parameters.tonal_reduction);

  // 3.3 Apply psychoacoustic veto to preserve psychoacoustically masked signal
  masking_veto_apply(self->masking_veto, smoothed_magnitude, delayed_noise,
                     fft_spectrum, alpha, self->parameters.masking_depth);

  // 3.4. Transient Protection:
  // Strictly on frequencies where transient was detected, drop alpha to
  // ALPHA_MIN (1.0)
  if (self->is_transient_detected) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_mask[k];
      if (t_weight > 0.0F) {
        float prot_factor = sqrtf(t_weight);
        alpha[k] =
            (alpha[k] * (1.0F - prot_factor)) + (ALPHA_MIN * prot_factor);
        alpha[k] = fmaxf(ALPHA_MIN, alpha[k]);
      }
    }
  }

  // 3.5. Gain Calculation
  calculate_gains(self->real_spectrum_size, self->fft_size, smoothed_magnitude,
                  delayed_noise, gain_out, alpha, beta,
                  self->gain_calculation_type);

  if (self->is_transient_detected) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_mask[k];
      if (t_weight > 0.0F) {
        gain_out[k] = fmaxf(gain_out[k], t_weight);
      }
    }
  }

  return true;
}

/**
 * 1D temporal/spatial chain: operates on the common delayed frame (uniform
 * shift of the legacy 1D output). Pre-subtraction magnitude smoothing, then
 * suppression/masking/gain, then temporal + spatial gain smoothing.
 */
static void run_temporal_chain(SbSpectralDenoiser* self,
                               const float* delayed_fft,
                               const float* delayed_noise, float* gain_out,
                               float* alpha, float* beta) {
  // Extract magnitude of the delayed frame (reuses the spectral features
  // buffer; the current-frame reference spectrum is no longer needed here)
  float* delayed_magnitude =
      get_spectral_feature(self->spectral_features, delayed_fft, self->fft_size,
                           self->spectrum_type);

  // Pre-Subtraction temporal and spatial smoothing on input magnitude to
  // suppress musical noise. Transient bins are spared bin-by-bin during gain
  // calculation and time smoothing.
  const float* effective_magnitude;

  if (self->parameters.smoothing_factor > 0.0F) {
    float smooth_alpha =
        fminf(0.95F, self->parameters.smoothing_factor * 0.90F);

    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float raw = delayed_magnitude[k];
      float prev = self->smoothed_magnitude[k];
      // Bin-by-bin adaptive smoothing: open immediately on transient bins
      // while keeping full smoothing elsewhere
      float t_w = self->transient_mask[k];
      float adapt_alpha = (1.0F - t_w) * smooth_alpha;
      self->smoothed_magnitude[k] =
          (adapt_alpha * prev) + ((1.0F - adapt_alpha) * raw);
    }

    int spatial_passes = (int)(self->parameters.smoothing_factor * 2.0F);
    for (int p = 0; p < spatial_passes; ++p) {
      spectral_smoothing_apply_spatial(self->smoothed_magnitude,
                                       self->real_spectrum_size);
    }

    effective_magnitude = self->smoothed_magnitude;
  } else {
    memcpy(self->smoothed_magnitude, delayed_magnitude,
           self->real_spectrum_size * sizeof(float));
    effective_magnitude = delayed_magnitude;
  }

  // Publish the smoothed magnitude so a pending switch to NLM starts from an
  // aligned frame
  spectral_circular_buffer_push(self->circular_buffer, self->layer_smoothed,
                                self->smoothed_magnitude);

  // Calculate SNR-dependent oversubtraction factors (Alpha/Beta)
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->suppression_engine, effective_magnitude,
                               delayed_noise, suppression_params, alpha, beta);

  // Detect tonal components and boost alpha at tonal bins
  tonal_reducer_run(self->tonal_reducer, delayed_noise,
                    get_noise_profile(self->noise_profile, CV_MASK),
                    is_noise_estimation_available(self->noise_profile, CV_MASK),
                    alpha, self->parameters.tonal_reduction);

  // Apply Structural Veto to rescue psychoacoustically masked signal
  masking_veto_apply(self->masking_veto, effective_magnitude, delayed_noise,
                     NULL, alpha, self->parameters.masking_depth);

  // When transients are detected and enabled, drop alphas firmly to ALPHA_MIN
  // (1.0) strictly on the specific frequencies where transient energy was
  // detected.
  if (self->is_transient_detected) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_mask[k];
      if (t_weight > 0.0F) {
        float prot_factor = sqrtf(t_weight);
        alpha[k] =
            (alpha[k] * (1.0F - prot_factor)) + (ALPHA_MIN * prot_factor);
        alpha[k] = fmaxf(ALPHA_MIN, alpha[k]);
      }
    }
  }

  // Gain Calculation
  calculate_gains(self->real_spectrum_size, self->fft_size, effective_magnitude,
                  delayed_noise, gain_out, alpha, beta,
                  self->gain_calculation_type);

  // Transient Protection: ensure transient bins have gain near 1.0
  if (self->is_transient_detected) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_mask[k];
      if (t_weight > 0.0F) {
        gain_out[k] = fmaxf(gain_out[k], t_weight);
      }
    }
  }

  // Temporal gain smoothing: on detected transients, the time smoother
  // attacks instantly (0 attack time) strictly on transient bins, while
  // non-transient bins receive 100% of full time smoothing.
  TimeSmoothingParameters spectral_smoothing_parameters =
      (TimeSmoothingParameters){
          .smoothing = self->parameters.smoothing_factor,
          .transient_mask = self->transient_mask,
      };
  spectral_smoothing_run(self->spectrum_smoothing,
                         spectral_smoothing_parameters, gain_out);
  if (self->parameters.smoothing_factor > 0.0f) {
    int passes = 1 + (int)(self->parameters.smoothing_factor * 2.0f);
    for (int p = 0; p < passes; ++p) {
      spectral_smoothing_apply_spatial(gain_out, self->real_spectrum_size);
    }
  }
}

const float* spectral_denoiser_get_tonal_mask(
    SpectralProcessorHandle instance) {
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  return (self && self->tonal_reducer)
             ? tonal_reducer_get_mask(self->tonal_reducer)
             : NULL;
}

uint32_t spectral_denoiser_get_peaks(SpectralProcessorHandle instance,
                                     float* peak_freqs_hz, uint32_t max_peaks) {
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  return (self && self->tonal_reducer)
             ? tonal_reducer_get_peaks(self->tonal_reducer, peak_freqs_hz,
                                       max_peaks)
             : 0;
}

const float* spectral_denoiser_get_active_noise_profile(
    SpectralProcessorHandle instance) {
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  if (!self) {
    return NULL;
  }
  int idx = atomic_load_explicit(&self->active_noise_idx, memory_order_acquire);
  return self->noise_spectrum_buffers[idx];
}

void spectral_denoiser_reset_noise_profile(SpectralProcessorHandle instance) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;
  if (!self) {
    return;
  }
  if (self->noise_estimator) {
    noise_estimation_reset(self->noise_estimator);
  }
  if (self->tonal_reducer) {
    tonal_reducer_reset(self->tonal_reducer);
  }
  if (self->manual_noise_floor) {
    memset(self->manual_noise_floor, 0,
           self->real_spectrum_size * sizeof(float));
  }
  self->was_learning = false;
  self->last_adaptive_state = 0;
}

uint32_t spectral_denoiser_get_latency_frames(
    SpectralProcessorHandle instance) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  if (!self) {
    return 0;
  }

  // Common delay: NLM look-ahead applies to both smoothing modes so the
  // reported latency never changes on a runtime mode switch
  return nlm_filter_get_latency_frames(self->nlm_filter);
}

bool spectral_denoiser_is_transient_detected(SpectralProcessorHandle instance) {
  if (!instance) {
    return false;
  }
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  return self->is_transient_detected;
}

float spectral_denoiser_get_transient_intensity(
    SpectralProcessorHandle instance) {
  if (!instance) {
    return 0.0f;
  }
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  return self->transient_intensity;
}
