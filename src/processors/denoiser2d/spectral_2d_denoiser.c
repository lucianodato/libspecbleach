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

#include "spectral_2d_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/denoiser_post_process.h"
#include "shared/denoiser_logic/core/denoiser_profile_core.h"
#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/denoiser_logic/processing/gain_calculator.h"
#include "shared/denoiser_logic/processing/hpss_filter.h"
#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/denoiser_logic/processing/nlm_filter.h"
#include "shared/denoiser_logic/processing/suppression_engine.h"
#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/utils/spectral_circular_buffer.h"
#include "shared/utils/spectral_utils.h"
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

typedef struct Spectral2DDenoiser {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t hop;

  Denoiser2DParameters parameters;

  float* snr_frame;                 // Current SNR frame for NLM input
  float* smoothed_snr;              // Smoothed SNR output from NLM
  float* gain_spectrum;             // Computed gains
  float* noise_spectrum;            // Copy of noise profile for processing
  float* noise_spectrum_buffers[2]; // Double-buffered noise spectrum for
                                    // lock-free SPSC publication
  atomic_int
      active_noise_idx; // Index of current published noise spectrum (0 or 1)
  float* alpha;         // Oversubtraction factors
  float* beta;          // Undersubtraction factors
  float* manual_noise_floor; // Manual profile floor
  TonalReducer* tonal_reducer;

  // Reusable circular buffer for aligned temporal analysis
  SbSpectralCircularBuffer* circular_buffer;
  uint32_t layer_fft;
  uint32_t layer_noise;
  uint32_t layer_nlm_smoothed;

  SpectrumType spectrum_type;
  GainCalculationType gain_calculation_type;

  NoiseProfile* noise_profile;
  NoiseEstimator* noise_estimator;
  AdaptiveNoiseEstimator* adaptive_estimator;
  NlmFilter* nlm_filter;
  SpectralFeatures* spectral_features;
  MaskingVeto* masking_veto;
  SuppressionEngine* suppression_engine;
  NoiseFloorManager* noise_floor_manager;
  HpssFilter* hpss_filter;

  float* delayed_magnitude;
  float* mask_harmonic;
  float* mask_percussive;
  float* g_h;
  float* g_p;

  int last_adaptive_state;
  int last_noise_estimation_method;
  float aggressiveness;
  bool was_learning;
} Spectral2DDenoiser;

// Header-only helpers or declarations would go here

SpectralProcessorHandle spectral_2d_denoiser_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, NoiseProfile* noise_profile) {

  if (!noise_profile || sample_rate == 0 || fft_size == 0 ||
      overlap_factor == 0) {
    return NULL;
  }

  Spectral2DDenoiser* self =
      (Spectral2DDenoiser*)calloc(1U, sizeof(Spectral2DDenoiser));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = (fft_size / 2U) + 1U;
  self->sample_rate = sample_rate;
  self->hop = fft_size / overlap_factor;
  self->spectrum_type = SPECTRAL_TYPE_2D;
  self->gain_calculation_type = GAIN_ESTIMATION_TYPE_2D;
  self->noise_profile = noise_profile;

  // Allocate buffers
  self->snr_frame = (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->snr_frame) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->smoothed_snr = (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->smoothed_snr) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  if (!self->gain_spectrum) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }
  (void)initialize_spectrum_with_value(self->gain_spectrum, fft_size, 1.0F);

  self->noise_spectrum =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->noise_spectrum) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->noise_spectrum_buffers[0] =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noise_spectrum_buffers[1] =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->noise_spectrum_buffers[0] || !self->noise_spectrum_buffers[1]) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }
  atomic_init(&self->active_noise_idx, 0);

  self->alpha = (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->alpha) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }
  (void)initialize_spectrum_with_value(self->alpha, self->real_spectrum_size,
                                       1.0F);

  self->beta = (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->beta) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->manual_noise_floor =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->manual_noise_floor) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->tonal_reducer = tonal_reducer_initialize(self->real_spectrum_size,
                                                 self->sample_rate, fft_size);
  if (!self->tonal_reducer) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  // Initialize circular buffer
  self->circular_buffer = spectral_circular_buffer_create(DELAY_BUFFER_FRAMES);
  if (!self->circular_buffer) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->layer_fft =
      spectral_circular_buffer_add_layer(self->circular_buffer, self->fft_size);
  self->layer_noise = spectral_circular_buffer_add_layer(
      self->circular_buffer, self->real_spectrum_size);
  self->layer_nlm_smoothed = spectral_circular_buffer_add_layer(
      self->circular_buffer, self->real_spectrum_size);

  if (self->layer_fft == 0xFFFFFFFFU || self->layer_noise == 0xFFFFFFFFU ||
      self->layer_nlm_smoothed == 0xFFFFFFFFU) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->was_learning = false;
  self->aggressiveness = 0.0f;
  self->parameters.tonal_reduction = 0.0f;

  // Initialize noise estimator for learning mode
  self->noise_estimator = noise_estimation_initialize(fft_size, noise_profile);
  if (!self->noise_estimator) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  // Initialize NLM filter
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
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  // Initialize spectral features
  self->spectral_features =
      spectral_features_initialize(self->real_spectrum_size);
  if (!self->spectral_features) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->masking_veto = masking_veto_initialize(
      self->fft_size, self->sample_rate, CRITICAL_BANDS_TYPE_2D,
      self->spectrum_type, false, USE_TEMPORAL_MASKING_2D_DEFAULT);
  self->suppression_engine = suppression_engine_initialize(
      self->real_spectrum_size, self->sample_rate, CRITICAL_BANDS_TYPE_2D,
      self->spectrum_type, true, USE_TEMPORAL_MASKING_2D_DEFAULT);

  if (!self->masking_veto || !self->suppression_engine) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->noise_floor_manager = noise_floor_manager_initialize(fft_size);

  HpssConfig hpss_cfg = {
      .real_spectrum_size = self->real_spectrum_size,
      .time_window_size = HPSS_TIME_WINDOW_2D_DEFAULT,
      .freq_window_size = HPSS_FREQ_WINDOW_2D_DEFAULT,
  };
  self->hpss_filter = hpss_filter_initialize(hpss_cfg);
  self->delayed_magnitude =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->mask_harmonic = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->mask_percussive =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->g_h = (float*)calloc(self->fft_size, sizeof(float));
  self->g_p = (float*)calloc(self->fft_size, sizeof(float));

  if (!self->noise_floor_manager || !self->hpss_filter ||
      !self->delayed_magnitude || !self->mask_harmonic ||
      !self->mask_percussive || !self->g_h || !self->g_p) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  return self;
}

void spectral_2d_denoiser_free(SpectralProcessorHandle instance) {
  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;

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

  free(self->snr_frame);
  free(self->smoothed_snr);
  free(self->gain_spectrum);
  free(self->noise_spectrum);
  if (self->noise_spectrum_buffers[0]) {
    free(self->noise_spectrum_buffers[0]);
  }
  if (self->noise_spectrum_buffers[1]) {
    free(self->noise_spectrum_buffers[1]);
  }
  free(self->alpha);
  free(self->beta);
  if (self->manual_noise_floor) {
    free(self->manual_noise_floor);
  }

  if (self->circular_buffer) {
    spectral_circular_buffer_free(self->circular_buffer);
  }

  if (self->tonal_reducer) {
    tonal_reducer_free(self->tonal_reducer);
  }

  if (self->hpss_filter) {
    hpss_filter_free(self->hpss_filter);
  }
  if (self->delayed_magnitude) {
    free(self->delayed_magnitude);
  }
  if (self->mask_harmonic) {
    free(self->mask_harmonic);
  }
  if (self->mask_percussive) {
    free(self->mask_percussive);
  }
  if (self->g_h) {
    free(self->g_h);
  }
  if (self->g_p) {
    free(self->g_p);
  }

  free(self);
}

bool load_2d_reduction_parameters(SpectralProcessorHandle instance,
                                  Denoiser2DParameters parameters) {
  if (!instance) {
    return false;
  }

  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;

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
  self->aggressiveness = parameters.aggressiveness;

  if (self->hpss_filter) {
    hpss_filter_set_quality_mode(self->hpss_filter,
                                 (HpssQualityMode)parameters.hpss_quality_mode);
  }

  // Update NLM h parameter based on smoothing factor (scales h up to 5.0F for
  // strong NLM patch smoothing)
  if (self->nlm_filter) {
    nlm_filter_set_h_parameter(self->nlm_filter,
                               (parameters.smoothing_factor > 0.0F)
                                   ? (0.5F + parameters.smoothing_factor * 4.5F)
                                   : 0.0F);
  }

  return true;
}

bool spectral_2d_denoiser_run(SpectralProcessorHandle instance,
                              float* fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;

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
  };
  denoiser_profile_core_update(profile_params, reference_spectrum);

  // 2.1 Align internal state and output to the delayed frame (temporal
  // plumbing)
  // 2.1.1 Push current spectra to circular buffer
  spectral_circular_buffer_push(self->circular_buffer, self->layer_fft,
                                fft_spectrum);
  spectral_circular_buffer_push(self->circular_buffer, self->layer_noise,
                                self->noise_spectrum);

  // 2.1.2 HPSS Process (operates on current reference_spectrum, outputs
  // delayed_magnitude and masks at delay L_hpss)
  if (!hpss_filter_process(self->hpss_filter, reference_spectrum,
                           self->delayed_magnitude, self->mask_harmonic,
                           self->mask_percussive)) {
    memcpy(self->delayed_magnitude, reference_spectrum,
           self->real_spectrum_size * sizeof(float));
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      self->mask_harmonic[k] = 1.0f;
      self->mask_percussive[k] = 0.0f;
    }
  }

  // 2.1.3 Compute SNR for NLM using CURRENT noise and push frame
  nlm_filter_calculate_snr(self->nlm_filter, reference_spectrum,
                           self->noise_spectrum, self->snr_frame);
  nlm_filter_push_frame(self->nlm_filter, self->snr_frame);

  const uint32_t nlm_delay = nlm_filter_get_latency_frames(self->nlm_filter);
  const uint32_t hpss_delay = hpss_filter_get_latency_frames(self->hpss_filter);
  const uint32_t total_delay =
      (hpss_delay > nlm_delay) ? hpss_delay : nlm_delay;

  float* nlm_intermediate_noise = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_noise, nlm_delay);
  if (!nlm_intermediate_noise) {
    nlm_intermediate_noise = self->noise_spectrum;
  }

  // Process NLM filter (outputs smoothed SNR at nlm_delay)
  if (nlm_filter_process(self->nlm_filter, self->smoothed_snr)) {
    nlm_filter_reconstruct_magnitude(self->nlm_filter, self->smoothed_snr,
                                     nlm_intermediate_noise, self->snr_frame);
    spectral_circular_buffer_push(self->circular_buffer,
                                  self->layer_nlm_smoothed, self->snr_frame);
  } else {
    // If NLM not ready yet, push current reference_spectrum
    spectral_circular_buffer_push(self->circular_buffer,
                                  self->layer_nlm_smoothed, reference_spectrum);
  }

  // 2.1.4 Retrieve unified aligned frames at total_delay
  float* delayed_spectrum = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_fft, total_delay);
  float* delayed_noise = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_noise, total_delay);
  float* smoothed_magnitude = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_nlm_smoothed, total_delay - nlm_delay);

  if (!delayed_spectrum) {
    delayed_spectrum = fft_spectrum;
  }
  if (!delayed_noise) {
    delayed_noise = self->noise_spectrum;
  }
  if (!smoothed_magnitude) {
    smoothed_magnitude = self->delayed_magnitude;
  }

  // Align output to delayed frame for post-processing
  if (delayed_spectrum != fft_spectrum) {
    memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));
  }

  // 3. Denoising Stage: Calculate gains and apply psychoacoustic constraints
  // 3.1 Calculate SNR-dependent oversubtraction factors (Alpha/Beta)
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->suppression_engine, smoothed_magnitude,
                               delayed_noise, suppression_params, self->alpha,
                               self->beta);

  float onset_ratio = hpss_filter_get_onset_ratio(self->hpss_filter);
  apply_onset_alpha_ducking(self->alpha, self->real_spectrum_size, onset_ratio);
  ;

  // 3.2 Detect tonal components and boost alpha at tonal bins
  tonal_reducer_run(self->tonal_reducer, delayed_noise,
                    get_noise_profile(self->noise_profile, CV_MASK),
                    is_noise_estimation_available(self->noise_profile, CV_MASK),
                    self->alpha, self->parameters.tonal_reduction);

  // 3.3 Apply psychoacoustic veto to preserve transients and moderate artifacts
  masking_veto_apply(self->masking_veto, smoothed_magnitude, delayed_noise,
                     fft_spectrum, self->alpha,
                     self->parameters.nlm_masking_protection);

  // 3.4 Dual-Path Gain Engine
  // Calculate Harmonic Gain G_H from smoothed harmonic path (NLM smoothed
  // magnitude)
  calculate_gains(self->real_spectrum_size, self->fft_size, smoothed_magnitude,
                  delayed_noise, self->g_h, self->alpha, self->beta,
                  self->gain_calculation_type);

  // Calculate Percussive Gain G_P from raw delayed magnitude (bypassing NLM
  // smoothing)
  calculate_gains(self->real_spectrum_size, self->fft_size,
                  self->delayed_magnitude, delayed_noise, self->g_p,
                  self->alpha, self->beta, self->gain_calculation_type);

  // Recombine gains: G_final = W_H * G_H + W_P * G_P
  for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
    self->gain_spectrum[k] = self->mask_harmonic[k] * self->g_h[k] +
                             self->mask_percussive[k] * self->g_p[k];
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

uint32_t spectral_2d_denoiser_get_latency_frames(
    SpectralProcessorHandle instance) {
  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;

  if (!self) {
    return 0;
  }

  uint32_t nlm_latency = nlm_filter_get_latency_frames(self->nlm_filter);
  uint32_t hpss_latency = hpss_filter_get_latency_frames(self->hpss_filter);
  return (hpss_latency > nlm_latency) ? hpss_latency : nlm_latency;
}

const float* spectral_2d_denoiser_get_tonal_mask(
    SpectralProcessorHandle instance) {
  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;
  return (self && self->tonal_reducer)
             ? tonal_reducer_get_mask(self->tonal_reducer)
             : NULL;
}

uint32_t spectral_2d_denoiser_get_peaks(SpectralProcessorHandle instance,
                                        float* peak_freqs_hz,
                                        uint32_t max_peaks) {
  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;
  return (self && self->tonal_reducer)
             ? tonal_reducer_get_peaks(self->tonal_reducer, peak_freqs_hz,
                                       max_peaks)
             : 0;
}

const float* spectral_2d_denoiser_get_active_noise_profile(
    SpectralProcessorHandle instance) {
  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;
  if (!self) {
    return NULL;
  }
  int idx = atomic_load_explicit(&self->active_noise_idx, memory_order_acquire);
  return self->noise_spectrum_buffers[idx];
}

void spectral_2d_denoiser_reset_noise_profile(
    SpectralProcessorHandle instance) {
  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;
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
