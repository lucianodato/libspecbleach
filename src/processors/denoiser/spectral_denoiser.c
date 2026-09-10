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
#include "shared/denoiser_logic/processing/bm3d_filter.h"
#include "shared/denoiser_logic/processing/dftt_filter.h"
#include "shared/denoiser_logic/processing/gain_calculator.h"
#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/denoiser_logic/processing/nlm_filter.h"
#include "shared/denoiser_logic/processing/release_shaper.h"
#include "shared/denoiser_logic/processing/suppression_engine.h"
#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/frame_rate_norm.h"
#include "shared/stft/stft_processor.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/spectral_circular_buffer.h"
#include "shared/utils/spectral_features.h"
#include "shared/utils/spectral_smoother.h"
#include "shared/utils/spectral_utils.h"
#include "shared/utils/transient_detector.h"
#include <float.h>
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
 * gain calculation and post-processing. All smoothing sub-modules are fully
 * allocated at initialization so switching modes at runtime is allocation-
 * free. Every mode shares a common algorithmic delay (the 2D look-ahead) so
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
  float* dftt_snr;                  // Refined SNR output from DFTT post-filter
  float* snr_delayed;               // Noisy SNR row aligned with NLM output
                                    // (DFTT ring input)
  float* gain_spectrum;             // Gain spectrum of the active chain
  float* gain_spectrum_b;           // Gain spectrum of the incoming chain
                                    // (transition crossfade only)
  float* noise_spectrum;            // Copy of noise profile for processing
  float* noise_spectrum_buffers[2]; // Double-buffered noise spectrum for
                                    // lock-free SPSC publication
  float* noise_bb;                  // Broadband profile (tonal peaks lifted)
  float* noise_tonal;               // Tonal residual (parallel gain path)
  float* gain_tonal;                // Tonal gain path output (shared scratch)
  atomic_int
      active_noise_idx; // Index of current published noise spectrum (0 or 1)
  float* alpha;         // Oversubtraction factors (active chain)
  float* beta;          // Undersubtraction factors (active chain)
  float* alpha_b;       // Oversubtraction factors (transition chain)
#if !TONAL_DUAL_PATH
  float* alpha_base;  // Legacy combined-path scratch (Berouti base)
  float* alpha_tonal; // Legacy combined-path scratch (tonal branch)
#endif
  float* beta_b; // Undersubtraction factors (transition chain)
#if TRANSIENT_RELIEF_PARALLEL
  float* alpha_relief; // Relief-branch alpha (all ALPHA_MIN)
  float* gain_relief;  // Relief-branch preservation gain (parallel output)
#endif
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
  Bm3dFilter* bm3d_filter;
  DfttFilter* dftt_filter;
  SpectralSmoother* spectrum_smoothing;
  ReleaseShaper* release_shaper;
  float* release_scale;
  SpectralFeatures* spectral_features;
  MaskingVeto* masking_veto;
  SuppressionEngine* suppression_engine;
  NoiseFloorManager* noise_floor_manager;
  CriticalBands* critical_bands;
  TransientDetector* transient_detector;

  // Circular buffer layers aligned at the common delay (broadband + tonal
  // profiles ride along so the parallel tonal path describes the same tile)
  uint32_t layer_noise_bb;
  uint32_t layer_noise_tonal;
  uint32_t layer_tonal_mask; // Mask aligned to the delayed noise_tonal layer
  float* band_energies;
  float* onset_weights;
  float* held_weights;        // Band weights decayed after detection (hold)
  float* transient_mask;      // Per-bin, clean-evidence gated (gain floor)
  float* transient_band_mask; // Band-level, ungated (alpha drop / smoothing)
  float* clean_magnitude;
  float* smoothed_magnitude;  // Temporal pre-subtraction smoothed magnitude
  float* knee_spectrum;       // Per-bin soft knee width (signal-dependent)
  float transient_hold_decay; // Per-hop hold decay factor
  uint32_t transient_hold_remaining; // Frames of mask hold left
  float transient_intensity;
  float hop_sec; // True hop in seconds (frame/overlap/sr); 0 = legacy derive
  bool smoothed_magnitude_seeded; // First frame seeds raw (no ramp-in)
  bool is_transient_detected;
  bool transient_protection_active; // Detected now or within hold window
  bool low_latency; // Causal 1D-only: zero look-ahead, no NLM delay

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

static bool is_nlm_family(const int mode) {
  return mode == SPECBLEACH_SMOOTHING_NLM_2D ||
         mode == SPECBLEACH_SMOOTHING_NLM_2D_DFTT;
}

static bool is_2d_family(const int mode) {
  return is_nlm_family(mode) || mode == SPECBLEACH_SMOOTHING_BM3D;
}

static int normalize_smoothing_mode(const int mode) {
  if (mode == SPECBLEACH_SMOOTHING_NLM_2D) {
    return SPECBLEACH_SMOOTHING_NLM_2D;
  }
  if (mode == SPECBLEACH_SMOOTHING_NLM_2D_DFTT) {
    return SPECBLEACH_SMOOTHING_NLM_2D_DFTT;
  }
  if (mode == SPECBLEACH_SMOOTHING_BM3D) {
    return SPECBLEACH_SMOOTHING_BM3D;
  }
  return SPECBLEACH_SMOOTHING_TEMPORAL;
}

/**
 * Bypass alignment: keeps the circular buffer, NLM history and output
 * rolling during idle/silence stretches so idle→active transitions stay
 * aligned and the reported look-ahead latency is preserved. Emits the frame
 * delayed by the configured NLM latency.
 */
static void align_bypass_frame(SbSpectralDenoiser* self, float* fft_spectrum,
                               const float* reference_spectrum) {
  if (self->low_latency) {
    return; // causal: emit current frame, no delay
  }
  // The bypass happens after the dual-path split, so the aligned layers stay
  // in the same domain as the active chain (SNR fields must not mix domains
  // across the idle→active boundary).
  spectral_circular_buffer_push(self->circular_buffer, self->layer_fft,
                                fft_spectrum);
  spectral_circular_buffer_push(self->circular_buffer, self->layer_noise,
                                self->noise_spectrum);
  spectral_circular_buffer_push(self->circular_buffer, self->layer_noise_bb,
                                self->noise_bb);
  spectral_circular_buffer_push(self->circular_buffer, self->layer_noise_tonal,
                                self->noise_tonal);
  const float* bypass_tonal_mask = tonal_reducer_get_mask(self->tonal_reducer);
  if (bypass_tonal_mask) {
    spectral_circular_buffer_push(self->circular_buffer, self->layer_tonal_mask,
                                  bypass_tonal_mask);
  }
  spectral_circular_buffer_push(self->circular_buffer, self->layer_smoothed,
                                reference_spectrum);
  nlm_filter_calculate_snr(self->nlm_filter, reference_spectrum, self->noise_bb,
                           self->snr_frame);
  nlm_filter_push_frame(self->nlm_filter, self->snr_frame);
  bm3d_filter_push_frame(self->bm3d_filter, self->snr_frame);
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
                          const float* delayed_noise_bb,
                          const float* delayed_noise_tonal,
                          const float* delayed_tonal_mask, uint32_t slot,
                          float* gain_out, float* alpha, float* beta);

static void run_temporal_chain(SbSpectralDenoiser* self,
                               const float* delayed_fft,
                               const float* delayed_noise_bb,
                               const float* delayed_noise_tonal,
                               const float* delayed_tonal_mask, uint32_t slot,
                               float* gain_out, float* alpha, float* beta);

static SpectralProcessorHandle spectral_denoiser_initialize_inner(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, const uint32_t hop_override,
    NoiseProfile* noise_profile, const bool low_latency) {

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
  self->hop =
      (hop_override > 0U) ? hop_override : (self->fft_size / overlap_factor);
  if (self->hop == 0U) {
    spectral_denoiser_free(self);
    return NULL;
  }
  self->hop_sec = sb_hop_sec(self->hop, sample_rate);
  self->sample_rate = sample_rate;
  self->spectrum_type = SPECTRAL_TYPE;
  self->gain_calculation_type = GAIN_ESTIMATION_TYPE;
  self->noise_profile = noise_profile;
  self->active_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  self->pending_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  self->previous_mode = SPECBLEACH_SMOOTHING_TEMPORAL;

  self->snr_frame = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->smoothed_snr = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->dftt_snr = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->snr_delayed = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->gain_spectrum = (float*)calloc(self->fft_size, sizeof(float));
  self->gain_spectrum_b = (float*)calloc(self->fft_size, sizeof(float));
  self->noise_spectrum =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noise_spectrum_buffers[0] =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noise_spectrum_buffers[1] =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noise_bb = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noise_tonal = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->gain_tonal = (float*)calloc(self->real_spectrum_size, sizeof(float));
  atomic_init(&self->active_noise_idx, 0);
  self->alpha = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->beta = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->alpha_b = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->beta_b = (float*)calloc(self->real_spectrum_size, sizeof(float));
#if TRANSIENT_RELIEF_PARALLEL
  self->alpha_relief = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->gain_relief = (float*)calloc(self->fft_size, sizeof(float));
#endif
#if !TONAL_DUAL_PATH
  self->alpha_base = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->alpha_tonal = (float*)calloc(self->real_spectrum_size, sizeof(float));
#endif
  self->manual_noise_floor =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->smoothed_magnitude =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->clean_magnitude =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->knee_spectrum = (float*)calloc(self->real_spectrum_size, sizeof(float));

  if (!self->snr_frame || !self->smoothed_snr || !self->dftt_snr ||
      !self->snr_delayed || !self->gain_spectrum || !self->gain_spectrum_b ||
      !self->noise_spectrum || !self->noise_spectrum_buffers[0] ||
      !self->noise_spectrum_buffers[1] || !self->noise_bb ||
      !self->noise_tonal || !self->gain_tonal || !self->alpha || !self->beta ||
      !self->alpha_b || !self->beta_b ||
#if TRANSIENT_RELIEF_PARALLEL
      !self->alpha_relief || !self->gain_relief ||
#endif
#if !TONAL_DUAL_PATH
      !self->alpha_base || !self->alpha_tonal ||
#endif
      !self->manual_noise_floor || !self->smoothed_magnitude ||
      !self->clean_magnitude || !self->knee_spectrum) {
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
#if TRANSIENT_RELIEF_PARALLEL
  // Relief branch evaluates the plain Wiener curve (no oversubtraction)
  (void)initialize_spectrum_with_value(self->alpha_relief,
                                       self->real_spectrum_size, ALPHA_MIN);
#endif

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
  self->layer_noise_bb = spectral_circular_buffer_add_layer(
      self->circular_buffer, self->real_spectrum_size);
  self->layer_noise_tonal = spectral_circular_buffer_add_layer(
      self->circular_buffer, self->real_spectrum_size);
  self->layer_tonal_mask = spectral_circular_buffer_add_layer(
      self->circular_buffer, self->real_spectrum_size);

  if (self->layer_fft == 0xFFFFFFFFU || self->layer_noise == 0xFFFFFFFFU ||
      self->layer_smoothed == 0xFFFFFFFFU ||
      self->layer_noise_bb == 0xFFFFFFFFU ||
      self->layer_noise_tonal == 0xFFFFFFFFU ||
      self->layer_tonal_mask == 0xFFFFFFFFU) {
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

  // Frame-rate normalization: fixed-ms / fixed-Hz geometry + per-hop alphas.
  const float hop_sec = self->hop_sec;
  // hop = frame/overlap_factor (exact with explicit hop, fft-derived legacy)
  const float frame_ms = hop_sec * (float)overlap_factor * 1000.0F;
  const float bin_hz = sb_bin_hz(self->sample_rate, self->fft_size);
  const SbNlmGeometry nlm_geo =
      sb_nlm_geometry_for_frame_ms(frame_ms, hop_sec, bin_hz);

  // NLM filter (2D smoothing strategy)
  NlmFilterConfig nlm_config = {
      .spectrum_size = self->real_spectrum_size,
      .time_buffer_size = nlm_geo.past + nlm_geo.future + 1,
      .patch_size = nlm_geo.patch,
      .paste_block_size = nlm_geo.paste,
      .search_range_freq = nlm_geo.search_freq,
      .search_range_time_past = nlm_geo.past,
      .search_range_time_future = nlm_geo.future,
      .h_parameter = NLM_DEFAULT_H_PARAMETER,
      .distance_threshold = 0.0F, // Use default (4 * h²)
  };
  self->nlm_filter = nlm_filter_initialize(nlm_config);
  if (!self->nlm_filter) {
    spectral_denoiser_free(self);
    return NULL;
  }

  // BM3D-lite shares the NLM geometry/latency so switches stay instant.
  Bm3dFilterConfig bm3d_config = {
      .spectrum_size = self->real_spectrum_size,
      .time_buffer_size = nlm_geo.past + nlm_geo.future + 1,
      .patch_size = nlm_geo.patch,
      .paste_block_size = nlm_geo.paste,
      .search_range_freq = nlm_geo.search_freq,
      .search_range_time_past = nlm_geo.past,
      .search_range_time_future = nlm_geo.future,
      .h_parameter = NLM_DEFAULT_H_PARAMETER,
  };
  self->bm3d_filter = bm3d_filter_initialize(bm3d_config);
  if (!self->bm3d_filter) {
    spectral_denoiser_free(self);
    return NULL;
  }

  // DFTT post-filter (paper S4.2 lite): past-only time span in ms so it adds
  // no latency on top of the NLM look-ahead; freq block in Hz so the
  // analysis stays invariant across frame sizes.
  const uint32_t dftt_span = sb_frames_for_ms(
      DFTT_TIME_MS, hop_sec, DFTT_MIN_TIME_FRAMES, DFTT_MAX_TIME_FRAMES);
  const uint32_t dftt_block = sb_bins_for_hz(
      DFTT_BLOCK_FREQ_HZ, bin_hz, DFTT_MIN_BLOCK_FREQ, DFTT_MAX_BLOCK_FREQ);
  self->dftt_filter =
      dftt_filter_initialize(self->real_spectrum_size, dftt_span, dftt_block);
  if (!self->dftt_filter) {
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
  spectral_smoothing_set_hop_samples(self->spectrum_smoothing, self->hop);

  // Adaptive release shaping (per-band closing-edge evidence)
  self->release_shaper =
      release_shaper_initialize(self->sample_rate, self->fft_size);
  if (!self->release_shaper) {
    spectral_denoiser_free(self);
    return NULL;
  }
  release_shaper_set_hop_sec(self->release_shaper, self->hop_sec);
  self->release_scale = (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->release_scale) {
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
  self->held_weights =
      (float*)calloc(num_bands > 0U ? num_bands : 1U, sizeof(float));
  self->transient_mask =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->transient_band_mask =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  // Hold decay: per-hop factor so the mask hold spans TRANSIENT_HOLD_SEC
  // regardless of frame size.
  self->transient_hold_decay =
      (self->hop_sec > 0.0F) ? expf(-self->hop_sec / TRANSIENT_HOLD_SEC) : 0.0F;
  self->transient_hold_remaining = 0U;
  self->transient_protection_active = false;

  if (!self->noise_floor_manager || !self->critical_bands ||
      !self->transient_detector || !self->band_energies ||
      !self->onset_weights || !self->transient_mask ||
      !self->transient_band_mask || !self->held_weights) {
    spectral_denoiser_free(self);
    return NULL;
  }

  if (hop_sec > 0.0F) {
    noise_estimation_set_hop_sec(self->noise_estimator, hop_sec);
    transient_detector_set_hop_sec(self->transient_detector, hop_sec);
    masking_veto_set_hop_sec(self->masking_veto, hop_sec);
    suppression_engine_set_hop_sec(self->suppression_engine, hop_sec);
  }

  const uint32_t transition_from_sec = (uint32_t)fmaxf(
      1.0F, (SMOOTHING_TRANSITION_SECONDS * (float)self->sample_rate) /
                (float)self->hop);
  self->transition_frames =
      (transition_from_sec < SMOOTHING_TRANSITION_MIN_FRAMES)
          ? SMOOTHING_TRANSITION_MIN_FRAMES
          : transition_from_sec;
  self->low_latency = low_latency;
  self->active_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  self->pending_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  self->previous_mode = SPECBLEACH_SMOOTHING_TEMPORAL;

  return self;
}

SpectralProcessorHandle spectral_denoiser_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, NoiseProfile* noise_profile) {
  return spectral_denoiser_initialize_inner(
      sample_rate, fft_size, overlap_factor, 0U, noise_profile, false);
}

SpectralProcessorHandle spectral_denoiser_initialize_with_hop(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, const uint32_t hop_samples,
    NoiseProfile* noise_profile, const bool low_latency) {
  return spectral_denoiser_initialize_inner(sample_rate, fft_size,
                                            overlap_factor, hop_samples,
                                            noise_profile, low_latency);
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
  if (self->bm3d_filter) {
    bm3d_filter_free(self->bm3d_filter);
  }
  if (self->dftt_filter) {
    dftt_filter_free(self->dftt_filter);
  }
  if (self->spectrum_smoothing) {
    spectral_smoothing_free(self->spectrum_smoothing);
  }
  if (self->release_shaper) {
    release_shaper_free(self->release_shaper);
  }
  free(self->release_scale);
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
  if (self->held_weights) {
    free(self->held_weights);
  }
  if (self->transient_mask) {
    free(self->transient_mask);
  }
  if (self->transient_band_mask) {
    free(self->transient_band_mask);
  }
  if (self->clean_magnitude) {
    free(self->clean_magnitude);
  }
  if (self->smoothed_magnitude) {
    free(self->smoothed_magnitude);
  }
  if (self->knee_spectrum) {
    free(self->knee_spectrum);
  }

  free(self->snr_frame);
  free(self->smoothed_snr);
  free(self->dftt_snr);
  free(self->snr_delayed);
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
#if TRANSIENT_RELIEF_PARALLEL
  free(self->alpha_relief);
  free(self->gain_relief);
#endif
#if !TONAL_DUAL_PATH
  free(self->alpha_base);
  free(self->alpha_tonal);
#endif
  free(self->noise_bb);
  free(self->noise_tonal);
  free(self->gain_tonal);
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
      if (self->adaptive_estimator && self->hop_sec > 0.0F) {
        adaptive_estimator_set_hop_sec(self->adaptive_estimator, self->hop_sec);
      }
      self->last_adaptive_state = 0;
    }
  }

  self->parameters = parameters;
  if (self->low_latency) {
    self->parameters.smoothing_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
    self->active_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
    self->pending_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
    self->previous_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
    self->in_transition = false;
  }

  // Runtime smoothing mode switching (allocation-free): the outgoing mode is
  // crossfaded against the incoming one over SMOOTHING_TRANSITION_SECONDS.
  // Within the 2D family (NLM <-> NLM+DFTT <-> BM3D) the switch is instant:
  // all sides share NLM history, latency and DFTT rings (pushed on every 2D
  // pass), so only the map source flips — no crossfade needed.
  if (!self->in_transition) {
    const int requested = normalize_smoothing_mode(parameters.smoothing_mode);
    if (requested != self->active_mode) {
      if (is_2d_family(requested) && is_2d_family(self->active_mode)) {
        // Crossing BM3D changes the map producer feeding the DFTT rings;
        // reset so the refinement only ever sees NLM priors (it falls back
        // to the raw NLM output until the history refills).
        if ((requested == SPECBLEACH_SMOOTHING_BM3D) !=
            (self->active_mode == SPECBLEACH_SMOOTHING_BM3D)) {
          dftt_filter_reset(self->dftt_filter);
        }
        self->active_mode = requested;
        self->pending_mode = requested;
      } else {
        self->previous_mode = self->active_mode;
        self->pending_mode = requested;
        self->transition_pos = 0U;
        self->in_transition = true;
      }
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
      // The chain-slot mapping is keyed on previous_mode, so the two tonal
      // gain histories must swap along with the mode metadata.
      tonal_reducer_swap_gain_slots(self->tonal_reducer);
    } else {
      self->pending_mode = requested;
    }
  }

  // Update NLM/BM3D h parameter based on smoothing factor
  const float h_value = (parameters.smoothing_factor > 0.0F)
                            ? (0.5F + (parameters.smoothing_factor * 4.5F))
                            : 0.0F;
  if (self->nlm_filter) {
    nlm_filter_set_h_parameter(self->nlm_filter, h_value);
  }
  if (self->bm3d_filter) {
    bm3d_filter_set_h_parameter(self->bm3d_filter, h_value);
  }

  // Update DFTT refinement strength (reduction-depth coupling, live)
  if (self->dftt_filter) {
    dftt_filter_set_strength(self->dftt_filter, parameters.dftt_strength);
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

  // 2.1 Dual-path tonal split: detect/publish the tonal mask and split the
  // noise profile into a broadband floor (bb) and a tonal residual. The
  // bb/tonal pair rides the circular buffer so the chains at the delayed
  // frame consume a split aligned to their own tile; the raw profile (this
  // frame's noise_spectrum) keeps feeding transient detection, the whitening
  // floor and the published public profile unchanged.
  // Legacy A/B build (TONAL_DUAL_PATH=0): identity split — the chains run on
  // the raw profile and the tonal gain path stays a no-op; the whole coupled
  // behavior lives in tonal_reducer_apply_alpha_boost below.
#if TONAL_DUAL_PATH
  const float split_reduction_gain = self->parameters.tonal_reduction;
#else
  const float split_reduction_gain = 1.0f;
#endif
  tonal_reducer_compute_split(
      self->tonal_reducer, self->noise_spectrum,
      get_noise_profile(self->noise_profile, CV_MASK),
      is_noise_estimation_available(self->noise_profile, CV_MASK),
      split_reduction_gain, self->noise_bb, self->noise_tonal);

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

  // 2.2 Align internal state and output to the common delayed frame
  // (skipped in low-latency mode: causal, zero look-ahead)
  const float* delayed_spectrum = fft_spectrum;
  const float* delayed_noise = self->noise_spectrum;
  const float* delayed_noise_bb = self->noise_bb;
  const float* delayed_noise_tonal = self->noise_tonal;
  const float* delayed_tonal_mask = tonal_reducer_get_mask(self->tonal_reducer);
  const float* nlm_smoothed = NULL;
  if (!self->low_latency) {
    spectral_circular_buffer_push(self->circular_buffer, self->layer_fft,
                                  fft_spectrum);
    spectral_circular_buffer_push(self->circular_buffer, self->layer_noise,
                                  self->noise_spectrum);
    spectral_circular_buffer_push(self->circular_buffer, self->layer_noise_bb,
                                  self->noise_bb);
    spectral_circular_buffer_push(self->circular_buffer,
                                  self->layer_noise_tonal, self->noise_tonal);

    // The tonal mask must ride with the residual it belongs to: the chains
    // run on the delayed noise_tonal, so the mask they evaluate must be the
    // one from that same frame, not the current one.
    const float* current_tonal_mask =
        tonal_reducer_get_mask(self->tonal_reducer);
    if (current_tonal_mask) {
      spectral_circular_buffer_push(self->circular_buffer,
                                    self->layer_tonal_mask, current_tonal_mask);
    }

    // Compute SNR for 2D filters using the broadband split (current frame)
    // and push frame. This keeps both NLM and BM3D histories rolling even in
    // temporal mode so a runtime mode switch is seamless and allocation-free.
    nlm_filter_calculate_snr(self->nlm_filter, reference_spectrum,
                             self->noise_bb, self->snr_frame);
    nlm_filter_push_frame(self->nlm_filter, self->snr_frame);
    bm3d_filter_push_frame(self->bm3d_filter, self->snr_frame);

    const uint32_t nlm_delay = nlm_filter_get_latency_frames(self->nlm_filter);

    // 2D smoothing (runs when a 2D mode is the active or the incoming mode).
    // The smoothed magnitude is captured explicitly so the temporal chain
    // cannot overwrite the shared alignment layer before the 2D chain
    // consumes it.
    const bool mode_2d_needed =
        is_2d_family(self->active_mode) ||
        (self->in_transition && is_2d_family(self->pending_mode));
    // DFTT refinement follows the DFTT mode: the active chain, or the incoming
    // side of a temporal crossfade. Rings are pushed on every 2D pass so they
    // stay warm for instant intra-family flips.
    const bool use_dftt =
        (self->active_mode == SPECBLEACH_SMOOTHING_NLM_2D_DFTT) ||
        (self->in_transition &&
         self->pending_mode == SPECBLEACH_SMOOTHING_NLM_2D_DFTT);
    // BM3D source follows the same active/incoming rule; DFTT modes always
    // read NLM so the refinement rings stay valid.
    const bool use_bm3d = (self->active_mode == SPECBLEACH_SMOOTHING_BM3D) ||
                          (self->in_transition &&
                           self->pending_mode == SPECBLEACH_SMOOTHING_BM3D);
    const bool filter_ran =
        mode_2d_needed &&
        (use_bm3d ? bm3d_filter_process(self->bm3d_filter, self->smoothed_snr)
                  : nlm_filter_process(self->nlm_filter, self->smoothed_snr));

    // Retrieve unified aligned frames at the common delay
    delayed_spectrum = spectral_circular_buffer_retrieve(
        self->circular_buffer, self->layer_fft, nlm_delay);
    delayed_noise = spectral_circular_buffer_retrieve(
        self->circular_buffer, self->layer_noise, nlm_delay);
    delayed_noise_bb = spectral_circular_buffer_retrieve(
        self->circular_buffer, self->layer_noise_bb, nlm_delay);
    delayed_noise_tonal = spectral_circular_buffer_retrieve(
        self->circular_buffer, self->layer_noise_tonal, nlm_delay);
    delayed_tonal_mask = spectral_circular_buffer_retrieve(
        self->circular_buffer, self->layer_tonal_mask, nlm_delay);

    if (!delayed_spectrum) {
      delayed_spectrum = fft_spectrum;
    }
    if (!delayed_noise) {
      delayed_noise = self->noise_spectrum;
    }
    if (!delayed_noise_bb) {
      delayed_noise_bb = self->noise_bb;
    }
    if (!delayed_noise_tonal) {
      delayed_noise_tonal = self->noise_tonal;
    }
    if (!delayed_tonal_mask) {
      delayed_tonal_mask = tonal_reducer_get_mask(self->tonal_reducer);
    }

    if (filter_ran) {
      // Noisy SNR row aligned with the 2D-emitted frame, recomputed from the
      // delayed frames so it describes the same tile the filters just
      // emitted. Feeds the DFTT rings and the confidence blend below.
      // Reuses the shared spectral_features scratch (the temporal chain
      // recomputes it on the same delayed frame anyway).
      float* delayed_reference =
          get_spectral_feature(self->spectral_features, delayed_spectrum,
                               self->fft_size, self->spectrum_type);
      nlm_filter_calculate_snr(self->nlm_filter, delayed_reference,
                               delayed_noise_bb, self->snr_delayed);
      // DFTT post-filter (paper S4.2): the noisy SNR row aligned with the
      // 2D-emitted frame — recomputed from the delayed frames so both ring
      // inputs describe the same tile — is refined while the NLM output sets
      // the suppression threshold. Falls back to the raw NLM output until the
      // DFTT history is full, or when the active mode is NLM-only.
      float* post_nlm = self->smoothed_snr;
      if (self->dftt_filter) {
        dftt_filter_push(self->dftt_filter, self->snr_delayed,
                         self->smoothed_snr);
        if (use_dftt &&
            dftt_filter_process(self->dftt_filter, self->dftt_snr)) {
          post_nlm = self->dftt_snr;
        }
      }
      // SNR-confidence blend-back (NLM/BM3D maps only): bins the smoother
      // itself rates far above the noise floor keep their raw value — map
      // smoothing drags strong harmonics — while the rest take the smoothed
      // estimate. Confidence is driven by the smoothed estimate (not the
      // instantaneous bin: noise spikes must not buy raw passthrough, or
      // gap suppression regresses). Steep 4th-power weighting
      // raw_w = e^4/(e^4+crossover^4) with crossover
      // SMOOTHING_CONFIDENCE_SNR. DFTT-refined output is exempt: its
      // flicker-suppression contract lives in the bins this blend would
      // restore to raw. Both post_nlm candidates are self-owned scratch
      // (the DFTT rings were already pushed), so blend in place.
      if (post_nlm == self->smoothed_snr) {
        const float conf_sq =
            SMOOTHING_CONFIDENCE_SNR * SMOOTHING_CONFIDENCE_SNR;
        const float conf_4th = conf_sq * conf_sq;
        const uint32_t blend_n = self->real_spectrum_size;
        for (uint32_t k = 0U; k < blend_n; k++) {
          const float e = post_nlm[k];
          const float e_sq = e * e;
          const float e_4th = e_sq * e_sq;
          const float raw_w = e_4th / (e_4th + conf_4th);
          post_nlm[k] = (raw_w * self->snr_delayed[k]) + ((1.0F - raw_w) * e);
        }
      }
      nlm_filter_reconstruct_magnitude(self->nlm_filter, post_nlm,
                                       delayed_noise_bb, self->snr_frame);
      spectral_circular_buffer_push(self->circular_buffer, self->layer_smoothed,
                                    self->snr_frame);
      nlm_smoothed = self->snr_frame;
    }
  }

  // Align output to delayed frame for post-processing
  if (!self->low_latency && delayed_spectrum != fft_spectrum) {
    memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));
  }

  // 2.3 Transient Detection via Transient Detector across Critical Bands.
  // Runs on the ALIGNED (delayed) frame so the transient mask describes the
  // same frame the gains below modify. Detecting on the current input frame
  // would fire ~46 ms before the transient is emitted and open gains on the
  // noise-only frames around it (audible noise pumping at every transient).
  // A clean signal estimate with scaled-up noise subtraction avoids false
  // triggering from residual musical noise.
  bool transient_enabled = (self->parameters.transient_protection_enable != 0);
  if (transient_enabled && self->critical_bands && self->transient_detector) {
    const float* delayed_magnitude =
        get_spectral_feature(self->spectral_features, delayed_spectrum,
                             self->fft_size, self->spectrum_type);
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      // Scale noise up using TRANSIENT_CLEAN_NOISE_SCALE to eliminate spurious
      // noise peaks
      float clean = fmaxf(delayed_magnitude[k] -
                              (TRANSIENT_CLEAN_NOISE_SCALE * delayed_noise[k]),
                          0.0F);
      self->clean_magnitude[k] = clean;
    }

    compute_critical_bands_spectrum(self->critical_bands, self->clean_magnitude,
                                    self->band_energies);
    self->is_transient_detected = transient_detector_process(
        self->transient_detector, self->band_energies, self->onset_weights,
        &self->transient_intensity);

    // Hold: keep the fired band weights alive (decayed) across the transient
    // decay tail, so oversubtraction relief outlives the detector's own
    // trigger window instead of cutting the tail off as soon as the SNR
    // drops below the trigger.
    if (self->is_transient_detected) {
      self->transient_hold_remaining =
          (self->transient_hold_decay > 0.0F)
              ? (uint32_t)((TRANSIENT_HOLD_SEC / self->hop_sec) + 0.5F)
              : 0U;
    } else if (self->transient_hold_remaining > 0U) {
      self->transient_hold_remaining--;
    }
    self->transient_protection_active =
        self->is_transient_detected || self->transient_hold_remaining > 0U;
    uint32_t num_bands = get_number_of_critical_bands(self->critical_bands);
    for (uint32_t b = 0; b < num_bands; ++b) {
      self->held_weights[b] =
          fmaxf(self->onset_weights[b],
                self->held_weights[b] * self->transient_hold_decay);
    }

    // Expand held band weights to two per-bin masks:
    // - transient_band_mask keeps the raw band weight: it only removes
    //   oversubtraction (alpha -> alpha_min) and shapes smoothing, where the
    //   Wiener curve itself keeps noise-dominant bins closed.
    // - transient_mask is additionally scaled by per-bin clean evidence and
    //   drives the hard gain floor: a band onset must not floor bins whose
    //   energy is mostly noise (clean estimate ~0 under the scaled
    //   subtraction), or every transient leaks the noise floor across the
    //   whole critical band.
    memset(self->transient_mask, 0, self->real_spectrum_size * sizeof(float));
    memset(self->transient_band_mask, 0,
           self->real_spectrum_size * sizeof(float));
    for (uint32_t b = 0; b < num_bands; ++b) {
      float bw = self->held_weights[b];
      if (bw > 0.0F) {
        CriticalBandIndexes idx = get_band_indexes(self->critical_bands, b);
        uint32_t end = (idx.end_position < self->real_spectrum_size)
                           ? idx.end_position
                           : self->real_spectrum_size;
        for (uint32_t k = idx.start_position; k < end; ++k) {
          self->transient_band_mask[k] =
              fmaxf(self->transient_band_mask[k], bw);
          float evidence = self->clean_magnitude[k] /
                           (delayed_magnitude[k] + TRANSIENT_BIN_EVIDENCE_EPS);
          self->transient_mask[k] =
              fmaxf(self->transient_mask[k], bw * evidence);
        }
      }
    }
  } else {
    self->is_transient_detected = false;
    self->transient_intensity = 0.0f;
    self->transient_protection_active = false;
    self->transient_hold_remaining = 0U;
    if (self->critical_bands) {
      uint32_t bands = get_number_of_critical_bands(self->critical_bands);
      for (uint32_t b = 0; b < bands; ++b) {
        self->held_weights[b] = 0.0F;
      }
    }
    memset(self->transient_mask, 0, self->real_spectrum_size * sizeof(float));
    memset(self->transient_band_mask, 0,
           self->real_spectrum_size * sizeof(float));
  }

  // 3. Denoising Stage: dispatch the active smoothing strategy
  // both during a runtime mode transition; low-latency is always temporal)
  float* gain_a = self->gain_spectrum;
  float* gain_b = self->gain_spectrum_b;

  if (!self->low_latency && self->in_transition) {
    const float total = (float)self->transition_frames;
    const float w = (float)self->transition_pos / total; // 0 → 1

    if (is_2d_family(self->previous_mode)) {
      (void)run_nlm_chain(self, fft_spectrum, nlm_smoothed, delayed_noise_bb,
                          delayed_noise_tonal, delayed_tonal_mask, 0U, gain_a,
                          self->alpha, self->beta);
      run_temporal_chain(self, fft_spectrum, delayed_noise_bb,
                         delayed_noise_tonal, delayed_tonal_mask, 1U, gain_b,
                         self->alpha_b, self->beta_b);
    } else {
      run_temporal_chain(self, fft_spectrum, delayed_noise_bb,
                         delayed_noise_tonal, delayed_tonal_mask, 0U, gain_a,
                         self->alpha, self->beta);
      (void)run_nlm_chain(self, fft_spectrum, nlm_smoothed, delayed_noise_bb,
                          delayed_noise_tonal, delayed_tonal_mask, 1U, gain_b,
                          self->alpha_b, self->beta_b);
    }

    for (uint32_t k = 0U; k < self->fft_size; ++k) {
      gain_a[k] = (gain_a[k] * (1.0F - w)) + (gain_b[k] * w);
    }

    self->transition_pos++;
    if (self->transition_pos >= self->transition_frames) {
      self->active_mode = self->pending_mode;
      self->in_transition = false;
      // The incoming chain built its one-pole tonal gain state in slot 1 while
      // fading in; it now becomes the active chain and reads slot 0.
      tonal_reducer_promote_gain_slot(self->tonal_reducer);
    }
  } else if (!self->low_latency && is_2d_family(self->active_mode)) {
    (void)run_nlm_chain(self, fft_spectrum, nlm_smoothed, delayed_noise_bb,
                        delayed_noise_tonal, delayed_tonal_mask, 0U, gain_a,
                        self->alpha, self->beta);
  } else {
    run_temporal_chain(self, fft_spectrum, delayed_noise_bb,
                       delayed_noise_tonal, delayed_tonal_mask, 0U, gain_a,
                       self->alpha, self->beta);
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
  if (!self->low_latency) {
    spectral_circular_buffer_advance(self->circular_buffer);
  }

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
 * 2D Non-Local Means chain: NLM-smoothed magnitude feeds broadband
 * suppression (Berouti alpha + masking veto) and a Wiener evaluation against
 * the BROADBAND noise profile only. The tonal residual is handled by the
 * parallel tonal gain path; the final bin gain is
 * min(broadband, tonal) with the transient floor re-asserted on top.
 * When the tonal path is inactive (reduction >= ~1.0 or zero mask) the tonal
 * gain is unity and the min() is a no-op.
 */
static bool run_nlm_chain(SbSpectralDenoiser* self, float* fft_spectrum,
                          const float* smoothed_magnitude,
                          const float* delayed_noise_bb,
                          const float* delayed_noise_tonal,
                          const float* delayed_tonal_mask, uint32_t slot,
                          float* gain_out, float* alpha, float* beta) {
  if (!smoothed_magnitude) {
    smoothed_magnitude = spectral_circular_buffer_retrieve(
        self->circular_buffer, self->layer_smoothed, 0U);
  }
  if (!smoothed_magnitude) {
    smoothed_magnitude = fft_spectrum;
  }

  // Calculate SNR-dependent oversubtraction factors (Alpha/Beta) on the
  // broadband profile: tonal noise no longer depresses the per-bin SNR here.
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->suppression_engine, smoothed_magnitude,
                               delayed_noise_bb, suppression_params, alpha,
                               beta);

#if TONAL_DUAL_PATH
  // Structural Veto on the broadband profile: alpha lifts where noise is
  // psychoacoustically masked. With the tonal notch decoupled into its own
  // gain path, the veto can no longer partially undo a tonal boost
  // (order-dependence removed).
  masking_veto_apply(self->masking_veto, smoothed_magnitude, delayed_noise_bb,
                     fft_spectrum, alpha, self->parameters.masking_depth);
#else
  // Legacy coupled path: parallel branches from the same Berouti base (the
  // tonal branch boosts, the veto branch preserves), combined once.
  memcpy(self->alpha_base, alpha, self->real_spectrum_size * sizeof(float));
  memcpy(self->alpha_tonal, alpha, self->real_spectrum_size * sizeof(float));
  tonal_reducer_apply_alpha_boost(self->tonal_reducer, self->alpha_tonal,
                                  self->parameters.tonal_reduction);
  masking_veto_apply(self->masking_veto, smoothed_magnitude, delayed_noise_bb,
                     fft_spectrum, alpha, self->parameters.masking_depth);
  for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
    const float boost = self->alpha_tonal[k] - self->alpha_base[k];
    const float preservation = self->alpha_base[k] - alpha[k];
    float combined = self->alpha_base[k] + boost - preservation;
    combined = fminf(combined, ALPHA_MAX_TONAL);
    alpha[k] = fmaxf(ALPHA_MIN, combined);
  }
#endif

  // Broadband gain calculation: stationary tonal peaks are gone from the
  // profile, so the 2D-smoothed Wiener gain is a smooth field.
  //
  // The NLM/DFTT chain deliberately applies NO transient relief. The 2D filter
  // already preserves onsets, and this chain has no gain smoother after
  // calculate_gains, so any per-frame alpha relief / gain blend / hard floor
  // lands in the final gain unsmoothed and injects musical noise across a
  // note's sustain (measured ~1.5x sustain MNI on pluck+noise). Transient
  // detection still runs globally for the UI; relief is applied only by the 1D
  // temporal chain, where the gain smoother owns the result.
  calculate_gains(self->real_spectrum_size, self->fft_size, smoothed_magnitude,
                  delayed_noise_bb, gain_out, alpha, beta,
                  self->gain_calculation_type, NULL);

#if TONAL_DUAL_PATH
  // Parallel tonal gain path + decision criterion. Ran on the same
  // smoothed magnitude so the min() compares like with like.
  tonal_reducer_compute_tonal_gains(
      self->tonal_reducer, slot, smoothed_magnitude, delayed_noise_tonal,
      delayed_tonal_mask, self->parameters.tonal_reduction, self->gain_tonal);
  for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
    gain_out[k] = fminf(gain_out[k], self->gain_tonal[k]);
  }
#endif

  return true;
}

/**
 * 1D temporal/spatial chain: operates on the common delayed frame (uniform
 * shift of the legacy 1D output). Pre-subtraction magnitude smoothing, then
 * suppression/masking/gain, then temporal + spatial gain smoothing.
 */
static void run_temporal_chain(SbSpectralDenoiser* self,
                               const float* delayed_fft,
                               const float* delayed_noise_bb,
                               const float* delayed_noise_tonal,
                               const float* delayed_tonal_mask, uint32_t slot,
                               float* gain_out, float* alpha, float* beta) {
  // Extract magnitude of the delayed frame (reuses the spectral features
  // buffer; the current-frame reference spectrum is no longer needed here)
  float* delayed_magnitude =
      get_spectral_feature(self->spectral_features, delayed_fft, self->fft_size,
                           self->spectrum_type);

  // Pre-Subtraction temporal stabilization on input magnitude ("time
  // smoothing of the signal spectrum"): fixed light one-pole with tau =
  // SPECTRAL_STABILIZATION_HOPS hops — frame-rate invariant by construction.
  // Gated on smoothing: slider 0 must stay a pristine raw-magnitude bypass.
  // Transient bins are spared bin-by-bin during gain calculation and time
  // smoothing.
  const float* effective_magnitude;

  if (self->parameters.smoothing_factor > 0.0F) {
    if (!self->smoothed_magnitude_seeded) {
      memcpy(self->smoothed_magnitude, delayed_magnitude,
             self->real_spectrum_size * sizeof(float));
      self->smoothed_magnitude_seeded = true;
    } else {
      const float stabilization_alpha =
          expf(-1.0F / (float)SPECTRAL_STABILIZATION_HOPS);

      for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
        float raw = delayed_magnitude[k];
        float prev = self->smoothed_magnitude[k];
        // Bin-by-bin adaptive smoothing: open immediately on transient bins
        // while keeping full smoothing elsewhere
        float t_w = self->transient_band_mask[k];
        float adapt_alpha = (1.0F - t_w) * stabilization_alpha;
        self->smoothed_magnitude[k] =
            (adapt_alpha * prev) + ((1.0F - adapt_alpha) * raw);
      }
    }

    int spatial_passes = (int)(self->parameters.smoothing_factor * 2.0F);
    for (int p = 0; p < spatial_passes; ++p) {
      spectral_smoothing_apply_spatial(self->smoothed_magnitude,
                                       self->real_spectrum_size);
    }

    effective_magnitude = self->smoothed_magnitude;

    // Adaptive release shaping: bands whose recent-energy envelope collapses
    // close fast (no spectral ghosts on the residual); the rest keep the full
    // slider release (anti-chirp protection in noise-only stretches)
    release_shaper_compute(self->release_shaper, effective_magnitude,
                           self->release_scale);
  } else {
    // Pristine bypass: raw magnitude everywhere, no smoothing state
    memcpy(self->smoothed_magnitude, delayed_magnitude,
           self->real_spectrum_size * sizeof(float));
    self->smoothed_magnitude_seeded = false;
    effective_magnitude = self->smoothed_magnitude;
  }

  // Publish the smoothed magnitude so a pending switch to NLM starts from an
  // aligned frame
  spectral_circular_buffer_push(self->circular_buffer, self->layer_smoothed,
                                self->smoothed_magnitude);

  // Calculate SNR-dependent oversubtraction factors (Alpha/Beta) on the
  // broadband profile: tonal noise no longer depresses the per-bin SNR here.
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->suppression_engine, effective_magnitude,
                               delayed_noise_bb, suppression_params, alpha,
                               beta);
#if TONAL_DUAL_PATH
  // Apply Structural Veto on the broadband profile. The tonal notch lives in
  // its own parallel gain path, so the veto can no longer partially undo a
  // tonal alpha boost (sequential order-dependence removed).
  masking_veto_apply(self->masking_veto, effective_magnitude, delayed_noise_bb,
                     NULL, alpha, self->parameters.masking_depth);
#else
  // Legacy coupled path: tonal alpha boost applied inline BEFORE the veto,
  // so the veto can partially undo the boost at masked bins (the coupled
  // behavior under measurement).
  tonal_reducer_apply_alpha_boost(self->tonal_reducer, alpha,
                                  self->parameters.tonal_reduction);
  masking_veto_apply(self->masking_veto, effective_magnitude, delayed_noise_bb,
                     NULL, alpha, self->parameters.masking_depth);
#endif
#if !TRANSIENT_RELIEF_PARALLEL
  // When transients are detected and enabled, drop alphas firmly to ALPHA_MIN
  // (1.0) strictly on the specific frequencies where transient energy was
  // detected. Band-level weight (see run_nlm_chain 3.4).
  if (self->transient_protection_active) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_band_mask[k];
      if (t_weight > 0.0F) {
        float prot_factor = sqrtf(t_weight);
        alpha[k] =
            (alpha[k] * (1.0F - prot_factor)) + (ALPHA_MIN * prot_factor);
        alpha[k] = fmaxf(ALPHA_MIN, alpha[k]);
      }
    }
  }
#endif

  // Signal-dependent knee width: bins decaying from recent signal presence
  // (stabilized energy above the current raw hop) get a wider knee so weak
  // component tails are forgiven instead of cut; steady or rising bins keep
  // the base knee. Pristine bypass keeps a zero knee everywhere.
  if (self->parameters.smoothing_factor > 0.0F) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float decay_evidence = 0.0F;
      if (delayed_noise_bb[k] > FLT_MIN) {
        decay_evidence = (effective_magnitude[k] - delayed_magnitude[k]) /
                         delayed_noise_bb[k];
      }
      self->knee_spectrum[k] =
          GAIN_WIENER_KNEE +
          (GAIN_KNEE_DECAY_BOOST *
           fminf(1.0F, fmaxf(0.0F, decay_evidence) / GAIN_KNEE_DECAY_RANGE));
    }
  } else {
    memset(self->knee_spectrum, 0, self->real_spectrum_size * sizeof(float));
  }

  // Gain Calculation on the broadband profile only: the temporal/spatial
  // smoothers below therefore never see a stationary tonal notch carved into
  // their gain field.
  calculate_gains(self->real_spectrum_size, self->fft_size, effective_magnitude,
                  delayed_noise_bb, gain_out, alpha, beta,
                  self->gain_calculation_type, self->knee_spectrum);

#if TRANSIENT_RELIEF_PARALLEL
  // Transient relief as a parallel GAIN-domain branch: blend the base gain
  // toward the plain Wiener curve (alpha = ALPHA_MIN, same knee) by the band
  // protection weight. Matches legacy at full/no protection; at partial
  // band weights it is softer than the legacy alpha lerp because the Wiener
  // curve is nonlinear in alpha, and the shared alpha is no longer mutated.
  if (self->transient_protection_active) {
    calculate_gains(self->real_spectrum_size, self->fft_size,
                    effective_magnitude, delayed_noise_bb, self->gain_relief,
                    self->alpha_relief, beta, self->gain_calculation_type,
                    self->knee_spectrum);
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      const float pf = sqrtf(self->transient_band_mask[k]);
      if (pf > 0.0F) {
        gain_out[k] = ((1.0F - pf) * gain_out[k]) + (pf * self->gain_relief[k]);
      }
    }
  }
#endif

  // Transient Protection: ensure transient bins have gain near 1.0
  if (self->transient_protection_active) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_mask[k];
      if (t_weight > 0.0F) {
        gain_out[k] = fmaxf(gain_out[k], t_weight);
      }
    }
  }

  // Temporal gain smoothing: on detected transients, the time smoother
  // attacks instantly (0 attack time) strictly on transient bins, while
  // non-transient bins receive 100% of full time smoothing. The release is
  // capped (GAIN_SMOOTHING_RELEASE_P_CAP): uncapped long releases freeze
  // gains high after offsets on gappy material, collapsing suppression
  // (measured: Temporal att 11.7 -> 4.0 dB on speech+city across the top
  // half of the slider), while the cap bounds mistracking without any fast
  // path that could chop onsets (cf. low-latency 130 ms release cap).
  TimeSmoothingParameters spectral_smoothing_parameters =
      (TimeSmoothingParameters){
          .smoothing = fminf(self->parameters.smoothing_factor,
                             GAIN_SMOOTHING_RELEASE_P_CAP),
          // Instant attack only on the detected frame; during the relief hold
          // the smoother must keep owning the gain (no held-mask override).
          .transient_mask =
              self->is_transient_detected ? self->transient_band_mask : NULL,
          .release_scale = (self->parameters.smoothing_factor > 0.0F)
                               ? self->release_scale
                               : NULL, // Unused while bypassed
      };
  spectral_smoothing_run(self->spectrum_smoothing,
                         spectral_smoothing_parameters, gain_out);
  if (self->parameters.smoothing_factor > 0.0f) {
    int passes = 1 + (int)(self->parameters.smoothing_factor * 2.0f);
    for (int p = 0; p < passes; ++p) {
      spectral_smoothing_apply_spatial(gain_out, self->real_spectrum_size);
    }
  }

  // Parallel tonal gain path + decision criterion, applied AFTER the
  // time/spatial smoothers: the notch enters the final gain at full depth
  // without being smeared into neighbors by the spatial FIR, and the
  // smoother state is never dragged around by mask flicker.
#if TONAL_DUAL_PATH
  tonal_reducer_compute_tonal_gains(
      self->tonal_reducer, slot, effective_magnitude, delayed_noise_tonal,
      delayed_tonal_mask, self->parameters.tonal_reduction, self->gain_tonal);
  for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
    gain_out[k] = fminf(gain_out[k], self->gain_tonal[k]);
  }
#endif

  // Transient Protection re-asserted after the combine: a band onset must
  // not be notched by a tonal bin underneath it (floor is now strictly >=
  // the legacy pre-smoothing behavior on transient x tonal overlap bins).
  // Keyed to the detector's current frame (the UI LED state), not the 200 ms
  // relief hold, so the smoothed gain wins again through the sustain.
  if (self->is_transient_detected) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_mask[k];
      if (t_weight > 0.0F) {
        gain_out[k] = fmaxf(gain_out[k], t_weight);
      }
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
  self->smoothed_magnitude_seeded = false;
  release_shaper_reset(self->release_shaper);
}

uint32_t spectral_denoiser_get_latency_frames(
    SpectralProcessorHandle instance) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  if (!self) {
    return 0;
  }

  // Common delay: 2D look-ahead applies to every smoothing mode so the
  // reported latency never changes on a runtime mode switch.
  // Low-latency mode is causal: zero look-ahead by construction.
  if (self->low_latency) {
    return 0;
  }
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
