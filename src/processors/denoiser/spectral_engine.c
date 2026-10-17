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

#include "spectral_engine.h"
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
#include "shared/utils/spectral_neighbors.h"
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
typedef struct DenoiserConfig {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t hop;

  SpectralEngineParameters parameters;

  SpectrumType spectrum_type;
  GainCalculationType gain_calculation_type;
  float hop_sec;    // True hop in seconds (frame/overlap/sr); 0 = legacy derive
  bool low_latency; // Causal 1D-only: zero look-ahead, no NLM delay
} DenoiserConfig;

typedef struct DenoiserSpectra {
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
  float* neighbor_power_row;        // Power row of the current input frame
                                    // (pushed into the neighbor layer, and
                                    // the only available row in low latency)
  float* neighbor_acc;              // Time-averaged neighborhood scratch
  float* blend_power;               // Neighbor-blended statistic feeding the
                                    // chains (shared scratch, consumed
                                    // synchronously inside one chain call)
  atomic_int
      active_noise_idx; // Index of current published noise spectrum (0 or 1)
  float* manual_noise_floor; // Manual profile floor
} DenoiserSpectra;

typedef struct DenoiserGains {
  float* alpha;   // Oversubtraction factors (active chain)
  float* beta;    // Undersubtraction factors (active chain)
  float* alpha_b; // Oversubtraction factors (transition chain)
  float* beta_b;  // Undersubtraction factors (transition chain)
#if TRANSIENT_RELIEF_PARALLEL
  float* alpha_relief; // Relief-branch alpha (all ALPHA_MIN)
  float* gain_relief;  // Relief-branch preservation gain (parallel output)
#endif
} DenoiserGains;

typedef struct DenoiserEngines {
  TonalReducer* tonal_reducer;

  // Reusable circular buffer for aligned temporal analysis (common delay)
  SbSpectralCircularBuffer* circular_buffer;

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
  SpectralNeighborFreqKernel* neighbor_freq_kernel; // Erb-scaled frequency
                                                    // halo (init-time build)
  MaskingVeto* masking_veto;
  SuppressionEngine* suppression_engine;
  NoiseFloorManager* noise_floor_manager;
  CriticalBands* critical_bands;
  TransientDetector* transient_detector;
} DenoiserEngines;

typedef struct DenoiserLayers {
  uint32_t layer_fft;
  uint32_t layer_noise;
  uint32_t layer_smoothed;
  uint32_t layer_power; // Power row of each frame (neighborhood statistic)
  // Circular buffer layers aligned at the common delay (broadband + tonal
  // profiles ride along so the parallel tonal path describes the same tile)
  uint32_t layer_noise_bb;
  uint32_t layer_noise_tonal;
  uint32_t layer_tonal_mask; // Mask aligned to the delayed noise_tonal layer
} DenoiserLayers;

typedef struct DenoiserTransient {
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
  bool smoothed_magnitude_seeded; // First frame seeds raw (no ramp-in)
  bool is_transient_detected;
  bool transient_protection_active; // Detected now or within hold window
} DenoiserTransient;

typedef struct DenoiserMode {
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
  float aggressiveness;
  bool was_learning;
} DenoiserMode;

typedef struct SbSpectralEngine {
  DenoiserConfig config;
  DenoiserSpectra spectra;
  DenoiserGains gains;
  DenoiserEngines engines;
  DenoiserLayers layers;
  DenoiserTransient transient;
  DenoiserMode mode;
} SbSpectralEngine;

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

static void push_aligned_frame_layers(SbSpectralEngine* self,
                                      const float* fft_spectrum) {
  // Power row of this frame rides along with the aligned layers so the
  // gain statistic can gather past/future neighborhood rows at the common
  // delay no matter which chain (or no chain) renders the frame.
  spectral_neighbors_compute_power_row(self->spectra.neighbor_power_row,
                                       fft_spectrum, self->config.fft_size);
  spectral_circular_buffer_push(self->engines.circular_buffer,
                                self->layers.layer_power,
                                self->spectra.neighbor_power_row);
  spectral_circular_buffer_push(self->engines.circular_buffer,
                                self->layers.layer_fft, fft_spectrum);
  spectral_circular_buffer_push(self->engines.circular_buffer,
                                self->layers.layer_noise,
                                self->spectra.noise_spectrum);
  spectral_circular_buffer_push(self->engines.circular_buffer,
                                self->layers.layer_noise_bb,
                                self->spectra.noise_bb);
  spectral_circular_buffer_push(self->engines.circular_buffer,
                                self->layers.layer_noise_tonal,
                                self->spectra.noise_tonal);
  // The mask must ride with this frame's tonal residual.
  const float* tonal_mask = tonal_reducer_get_mask(self->engines.tonal_reducer);
  if (tonal_mask) {
    spectral_circular_buffer_push(self->engines.circular_buffer,
                                  self->layers.layer_tonal_mask, tonal_mask);
  }
}

static void push_filter_histories(SbSpectralEngine* self,
                                  const float* reference_spectrum) {
  nlm_filter_calculate_snr(self->engines.nlm_filter, reference_spectrum,
                           self->spectra.noise_bb, self->spectra.snr_frame);
  nlm_filter_push_frame(self->engines.nlm_filter, self->spectra.snr_frame);
  bm3d_filter_push_frame(self->engines.bm3d_filter, self->spectra.snr_frame);
}

/**
 * Neighborhood statistic (reduction blend): part of the per-bin gain
 * statistic is replaced by its time-frequency neighborhood power average so
 * the suppression rule reacts to the energy of surrounding bins and frames
 * instead of single bins. Musical noise comes from independent per-bin
 * statistical fluctuations triggering sub-band gates at random; borrowing
 * neighborhood energy de-correlates those flips. High values behave like a
 * band-wide gate and keep short residual bursts after loud sections.
 *
 * At 0.0 the enforced statistic passes through untouched. In the causal
 * low-latency mode no ring history exists, so the neighborhood degrades to
 * the frequency axis of the current frame (no extra latency). Returns the
 * statistic pointer the calling chain must consume; the shared blended
 * scratch is consumed synchronously inside the chain call, so both chains
 * may run back to back during a mode transition.
 *
 * When the blend is active, the statistic feeds every downstream consumer of
 * the chain input (alpha, knee, gains, tonal floor) so the blend describes
 * one coherent reduction decision.
 */
static const float* denoiser_neighbor_statistic(SbSpectralEngine* self,
                                                const float* statistic,
                                                const float* fft_spectrum) {
  const float blend = self->config.parameters.reduction_blend;
  if (blend <= 0.0F) {
    return statistic;
  }

  const float* rows[NEIGHBOR_KERNEL_TAPS] = {NULL};
  if (self->config.low_latency) {
    // Causal mode: fill the center row from the frame being processed
    spectral_neighbors_compute_power_row(self->spectra.neighbor_power_row,
                                         fft_spectrum, self->config.fft_size);
    rows[NEIGHBOR_KERNEL_HALO] = self->spectra.neighbor_power_row;
  } else {
    // Rows ordered oldest first; future rows sit below the common delay
    // (delays stay >= 0 because the common delay is at least the kernel
    // halo). Unwritten ring slots are zero-backed, so delays beyond the
    // filled history are skipped and weights renormalize over real rows.
    const uint32_t filled = spectral_circular_buffer_get_frames_filled(
        self->engines.circular_buffer);
    for (uint32_t i = 0U; i < NEIGHBOR_KERNEL_TAPS; i++) {
      const uint32_t delay =
          nlm_filter_get_latency_frames(self->engines.nlm_filter) +
          NEIGHBOR_KERNEL_HALO - i;
      rows[i] = (delay < filled) ? spectral_circular_buffer_retrieve(
                                       self->engines.circular_buffer,
                                       self->layers.layer_power, delay)
                                 : NULL;
    }
  }

  if (!rows[NEIGHBOR_KERNEL_HALO]) {
    // Ring not filled yet: keep the enforced statistic (house fallback for
    // the first frames of the stream)
    return statistic;
  }

  const SpectralNeighborParameters neighbor_params = {
      .rows = rows,
      .rows_count = NEIGHBOR_KERNEL_TAPS,
      .spectrum_size = self->config.real_spectrum_size,
      .freq_kernel = self->engines.neighbor_freq_kernel,
  };
  if (!spectral_neighbor_average(&neighbor_params, self->spectra.neighbor_acc,
                                 self->spectra.blend_power)) {
    return statistic;
  }

  // Per-bin attack gate: bins with transient evidence keep the enforced
  // statistic (no neighborhood borrowing) so onsets do not smear from
  // neighbors; stationary bins take the full slider blend. The mask is
  // zeroed when transient protection is off (gate stays neutral).
  const float* transient_mask = self->transient.transient_band_mask;
  for (uint32_t k = 0U; k < self->config.real_spectrum_size; k++) {
    const float blend_bin = blend * (1.0F - transient_mask[k]);
    self->spectra.blend_power[k] =
        statistic[k] +
        (blend_bin * (self->spectra.blend_power[k] - statistic[k]));
  }
  return self->spectra.blend_power;
}

/**
 * Bypass alignment: keeps the circular buffer, NLM history and output
 * rolling during idle/silence stretches so idle→active transitions stay
 * aligned and the reported look-ahead latency is preserved. Emits the frame
 * delayed by the configured NLM latency.
 */
static void align_bypass_frame(SbSpectralEngine* self, float* fft_spectrum,
                               const float* reference_spectrum) {
  if (self->config.low_latency) {
    return; // causal: emit current frame, no delay
  }
  // The bypass happens after the dual-path split, so the aligned layers stay
  // in the same domain as the active chain (SNR fields must not mix domains
  // across the idle→active boundary).
  push_aligned_frame_layers(self, fft_spectrum);
  spectral_circular_buffer_push(self->engines.circular_buffer,
                                self->layers.layer_smoothed,
                                reference_spectrum);
  push_filter_histories(self, reference_spectrum);
  const float* delayed_spectrum = spectral_circular_buffer_retrieve(
      self->engines.circular_buffer, self->layers.layer_fft,
      nlm_filter_get_latency_frames(self->engines.nlm_filter));
  if (delayed_spectrum) {
    memcpy(fft_spectrum, delayed_spectrum,
           self->config.fft_size * sizeof(float));
  }
  spectral_circular_buffer_advance(self->engines.circular_buffer);
}

// Chain I/O bundles: group the per-frame spectra and outputs so chain
// signatures stay under the parameter-count limit. Passed by value (a few
// pointers); no allocation, no lifetime concerns.
typedef struct DenoiserChainFrames {
  const float* smoothed_magnitude; // NULL lets the NLM chain fall back
  const float* noise_bb;
  // cppcheck-suppress unusedStructMember ; only read when TONAL_DUAL_PATH=1
  const float* noise_tonal;
  // cppcheck-suppress unusedStructMember ; only read when TONAL_DUAL_PATH=1
  const float* tonal_mask;
} DenoiserChainFrames;

typedef struct DenoiserChainOut {
  float* gain_out;
  float* alpha;
  float* beta;
} DenoiserChainOut;

// Aligned delayed frames consumed by the chains and post stage. Filled by
// the 2D pass (or defaulted to the current frame in low-latency mode).
typedef struct DenoiserAlignedFrames {
  const float* spectrum;
  const float* noise;
  const float* noise_bb;
  const float* noise_tonal;
  const float* tonal_mask;
  const float* nlm_smoothed;
} DenoiserAlignedFrames;

static bool run_nlm_chain(SbSpectralEngine* self, const float* fft_spectrum,
                          DenoiserChainFrames frames, uint32_t slot,
                          DenoiserChainOut out);

static void run_temporal_chain(SbSpectralEngine* self, const float* delayed_fft,
                               DenoiserChainFrames frames, uint32_t slot,
                               DenoiserChainOut out);

static bool denoiser_alloc_spectra(SbSpectralEngine* self) {
  self->spectra.snr_frame =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.smoothed_snr =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.dftt_snr =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.snr_delayed =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.gain_spectrum =
      (float*)calloc(self->config.fft_size, sizeof(float));
  self->spectra.gain_spectrum_b =
      (float*)calloc(self->config.fft_size, sizeof(float));
  self->spectra.noise_spectrum =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.noise_spectrum_buffers[0] =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.noise_spectrum_buffers[1] =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.noise_bb =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.noise_tonal =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.gain_tonal =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.neighbor_power_row =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.neighbor_acc =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->spectra.blend_power =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  atomic_init(&self->spectra.active_noise_idx, 0);
  self->gains.alpha =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->gains.beta =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->gains.alpha_b =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->gains.beta_b =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
#if TRANSIENT_RELIEF_PARALLEL
  self->gains.alpha_relief =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->gains.gain_relief =
      (float*)calloc(self->config.fft_size, sizeof(float));
#endif
  self->spectra.manual_noise_floor =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->transient.smoothed_magnitude =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->transient.clean_magnitude =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->transient.knee_spectrum =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));

  if (!self->spectra.snr_frame || !self->spectra.smoothed_snr ||
      !self->spectra.dftt_snr || !self->spectra.snr_delayed ||
      !self->spectra.gain_spectrum || !self->spectra.gain_spectrum_b ||
      !self->spectra.noise_spectrum ||
      !self->spectra.noise_spectrum_buffers[0] ||
      !self->spectra.noise_spectrum_buffers[1] || !self->spectra.noise_bb ||
      !self->spectra.noise_tonal || !self->spectra.gain_tonal ||
      !self->spectra.neighbor_power_row || !self->spectra.neighbor_acc ||
      !self->spectra.blend_power || !self->gains.alpha || !self->gains.beta ||
      !self->gains.alpha_b || !self->gains.beta_b ||
#if TRANSIENT_RELIEF_PARALLEL
      !self->gains.alpha_relief || !self->gains.gain_relief ||
#endif
      !self->spectra.manual_noise_floor ||
      !self->transient.smoothed_magnitude || !self->transient.clean_magnitude ||
      !self->transient.knee_spectrum) {
    return false;
  }

  (void)initialize_spectrum_with_value(self->spectra.gain_spectrum,
                                       self->config.fft_size, 1.0F);
  (void)initialize_spectrum_with_value(self->spectra.gain_spectrum_b,
                                       self->config.fft_size, 1.0F);
  (void)initialize_spectrum_with_value(self->gains.alpha,
                                       self->config.real_spectrum_size, 1.F);
  (void)initialize_spectrum_with_value(self->gains.alpha_b,
                                       self->config.real_spectrum_size, 1.F);
#if TRANSIENT_RELIEF_PARALLEL
  // Relief branch evaluates the plain Wiener curve (no oversubtraction)
  (void)initialize_spectrum_with_value(
      self->gains.alpha_relief, self->config.real_spectrum_size, ALPHA_MIN);
#endif
  return true;
}

static bool denoiser_init_early_objects(SbSpectralEngine* self,
                                        NoiseProfile* noise_profile,
                                        uint32_t fft_size) {
  // Initialize tonal reducer
  self->engines.tonal_reducer =
      tonal_reducer_initialize(self->config.real_spectrum_size,
                               self->config.sample_rate, self->config.fft_size);
  if (!self->engines.tonal_reducer) {
    return false;
  }

  // Circular buffer for temporal alignment (provides the common delay)
  self->engines.circular_buffer =
      spectral_circular_buffer_create(DELAY_BUFFER_FRAMES);
  if (!self->engines.circular_buffer) {
    return false;
  }

  self->layers.layer_fft = spectral_circular_buffer_add_layer(
      self->engines.circular_buffer, self->config.fft_size);
  self->layers.layer_noise = spectral_circular_buffer_add_layer(
      self->engines.circular_buffer, self->config.real_spectrum_size);
  self->layers.layer_smoothed = spectral_circular_buffer_add_layer(
      self->engines.circular_buffer, self->config.real_spectrum_size);
  self->layers.layer_power = spectral_circular_buffer_add_layer(
      self->engines.circular_buffer, self->config.real_spectrum_size);
  self->layers.layer_noise_bb = spectral_circular_buffer_add_layer(
      self->engines.circular_buffer, self->config.real_spectrum_size);
  self->layers.layer_noise_tonal = spectral_circular_buffer_add_layer(
      self->engines.circular_buffer, self->config.real_spectrum_size);
  self->layers.layer_tonal_mask = spectral_circular_buffer_add_layer(
      self->engines.circular_buffer, self->config.real_spectrum_size);

  if (self->layers.layer_fft == 0xFFFFFFFFU ||
      self->layers.layer_noise == 0xFFFFFFFFU ||
      self->layers.layer_smoothed == 0xFFFFFFFFU ||
      self->layers.layer_power == 0xFFFFFFFFU ||
      self->layers.layer_noise_bb == 0xFFFFFFFFU ||
      self->layers.layer_noise_tonal == 0xFFFFFFFFU ||
      self->layers.layer_tonal_mask == 0xFFFFFFFFU) {
    return false;
  }

  self->mode.was_learning = false;
  self->mode.aggressiveness = 0.0f;
  self->config.parameters.tonal_reduction = 0.0f;

  // Initialize noise estimator for learning mode
  self->engines.noise_estimator =
      noise_estimation_initialize(fft_size, noise_profile);
  if (!self->engines.noise_estimator) {
    return false;
  }
  return true;
}

static bool denoiser_init_filters(SbSpectralEngine* self,
                                  uint32_t overlap_factor, float hop_sec,
                                  float bin_hz, SbNlmGeometry nlm_geo) {
  // NLM filter (2D smoothing strategy)
  NlmFilterConfig nlm_config = {
      .spectrum_size = self->config.real_spectrum_size,
      .time_buffer_size = nlm_geo.past + nlm_geo.future + 1,
      .patch_size = nlm_geo.patch,
      .paste_block_size = nlm_geo.paste,
      .search_range_freq = nlm_geo.search_freq,
      .search_range_time_past = nlm_geo.past,
      .search_range_time_future = nlm_geo.future,
      .h_parameter = NLM_DEFAULT_H_PARAMETER,
      .distance_threshold = 0.0F, // Use default (4 * h²)
  };
  self->engines.nlm_filter = nlm_filter_initialize(nlm_config);
  if (!self->engines.nlm_filter) {
    return false;
  }

  // BM3D-lite shares the NLM geometry/latency so switches stay instant.
  Bm3dFilterConfig bm3d_config = {
      .spectrum_size = self->config.real_spectrum_size,
      .time_buffer_size = nlm_geo.past + nlm_geo.future + 1,
      .patch_size = nlm_geo.patch,
      .paste_block_size = nlm_geo.paste,
      .search_range_freq = nlm_geo.search_freq,
      .search_range_time_past = nlm_geo.past,
      .search_range_time_future = nlm_geo.future,
      .h_parameter = NLM_DEFAULT_H_PARAMETER,
  };
  self->engines.bm3d_filter = bm3d_filter_initialize(bm3d_config);
  if (!self->engines.bm3d_filter) {
    return false;
  }

  // DFTT post-filter (paper S4.2 lite): past-only time span in ms so it adds
  // no latency on top of the NLM look-ahead; freq block in Hz so the
  // analysis stays invariant across frame sizes.
  const uint32_t dftt_span = sb_frames_for_ms(
      DFTT_TIME_MS, hop_sec, DFTT_MIN_TIME_FRAMES, DFTT_MAX_TIME_FRAMES);
  const uint32_t dftt_block = sb_bins_for_hz(
      DFTT_BLOCK_FREQ_HZ, bin_hz, DFTT_MIN_BLOCK_FREQ, DFTT_MAX_BLOCK_FREQ);
  self->engines.dftt_filter = dftt_filter_initialize(
      self->config.real_spectrum_size, dftt_span, dftt_block);
  if (!self->engines.dftt_filter) {
    return false;
  }

  // Temporal smoother (1D smoothing strategy)
  self->engines.spectrum_smoothing = spectral_smoothing_initialize(
      self->config.fft_size, self->config.sample_rate, overlap_factor, FIXED);
  if (!self->engines.spectrum_smoothing) {
    return false;
  }
  spectral_smoothing_set_hop_samples(self->engines.spectrum_smoothing,
                                     self->config.hop);

  // Adaptive release shaping (per-band closing-edge evidence)
  self->engines.release_shaper = release_shaper_initialize(
      self->config.sample_rate, self->config.fft_size);
  if (!self->engines.release_shaper) {
    return false;
  }
  release_shaper_set_hop_sec(self->engines.release_shaper,
                             self->config.hop_sec);
  self->engines.release_scale =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  if (!self->engines.release_scale) {
    return false;
  }

  // Initialize spectral features
  self->engines.spectral_features =
      spectral_features_initialize(self->config.real_spectrum_size);
  if (!self->engines.spectral_features) {
    return false;
  }

  self->engines.masking_veto = masking_veto_initialize(
      self->config.fft_size, self->config.sample_rate, CRITICAL_BANDS_TYPE,
      self->config.spectrum_type, MASKING_VETO_USE_ABSOLUTE_THRESHOLD,
      USE_TEMPORAL_MASKING_DEFAULT);
  self->engines.suppression_engine = suppression_engine_initialize(
      self->config.real_spectrum_size, self->config.sample_rate,
      CRITICAL_BANDS_TYPE);

  if (!self->engines.masking_veto || !self->engines.suppression_engine) {
    return false;
  }
  return true;
}

static bool denoiser_init_transient_tail(SbSpectralEngine* self,
                                         float hop_sec) {
  self->engines.noise_floor_manager =
      noise_floor_manager_initialize(self->config.fft_size);

  self->engines.critical_bands = critical_bands_initialize(
      self->config.sample_rate, self->config.fft_size, CRITICAL_BANDS_TYPE);
  uint32_t num_bands =
      self->engines.critical_bands
          ? get_number_of_critical_bands(self->engines.critical_bands)
          : 0U;
  self->engines.transient_detector = transient_detector_initialize(num_bands);

  self->transient.band_energies =
      (float*)calloc(num_bands > 0U ? num_bands : 1U, sizeof(float));
  self->transient.onset_weights =
      (float*)calloc(num_bands > 0U ? num_bands : 1U, sizeof(float));
  self->transient.held_weights =
      (float*)calloc(num_bands > 0U ? num_bands : 1U, sizeof(float));
  self->transient.transient_mask =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  self->transient.transient_band_mask =
      (float*)calloc(self->config.real_spectrum_size, sizeof(float));
  // Hold decay: per-hop factor so the mask hold spans TRANSIENT_HOLD_SEC
  // regardless of frame size.
  self->transient.transient_hold_decay =
      (self->config.hop_sec > 0.0F)
          ? expf(-self->config.hop_sec / TRANSIENT_HOLD_SEC)
          : 0.0F;
  self->transient.transient_hold_remaining = 0U;
  self->transient.transient_protection_active = false;

  if (!self->engines.noise_floor_manager || !self->engines.critical_bands ||
      !self->engines.transient_detector || !self->transient.band_energies ||
      !self->transient.onset_weights || !self->transient.transient_mask ||
      !self->transient.transient_band_mask || !self->transient.held_weights) {
    return false;
  }

  if (hop_sec > 0.0F) {
    noise_estimation_set_hop_sec(self->engines.noise_estimator, hop_sec);
    transient_detector_set_hop_sec(self->engines.transient_detector, hop_sec);
    masking_veto_set_hop_sec(self->engines.masking_veto, hop_sec);
    suppression_engine_set_hop_sec(self->engines.suppression_engine, hop_sec);
  }

  const uint32_t transition_from_sec = (uint32_t)fmaxf(
      1.0F, (SMOOTHING_TRANSITION_SECONDS * (float)self->config.sample_rate) /
                (float)self->config.hop);
  self->mode.transition_frames =
      (transition_from_sec < SMOOTHING_TRANSITION_MIN_FRAMES)
          ? SMOOTHING_TRANSITION_MIN_FRAMES
          : transition_from_sec;
  return true;
}

static SpectralProcessorHandle spectral_engine_initialize_inner(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, const uint32_t hop_override,
    NoiseProfile* noise_profile, const bool low_latency) {

  if (!noise_profile || sample_rate == 0 || fft_size == 0 ||
      overlap_factor == 0) {
    return NULL;
  }

  SbSpectralEngine* self =
      (SbSpectralEngine*)calloc(1U, sizeof(SbSpectralEngine));
  if (!self) {
    return NULL;
  }

  self->config.fft_size = fft_size;
  self->config.real_spectrum_size = (self->config.fft_size / 2U) + 1U;
  self->config.hop = (hop_override > 0U)
                         ? hop_override
                         : (self->config.fft_size / overlap_factor);
  if (self->config.hop == 0U) {
    goto fail;
  }
  self->config.hop_sec = sb_hop_sec(self->config.hop, sample_rate);
  self->config.sample_rate = sample_rate;
  self->config.spectrum_type = SPECTRAL_TYPE;
  self->config.gain_calculation_type = GAIN_ESTIMATION_TYPE;
  self->engines.noise_profile = noise_profile;
  self->mode.active_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  self->mode.pending_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  self->mode.previous_mode = SPECBLEACH_SMOOTHING_TEMPORAL;

  if (!denoiser_alloc_spectra(self)) {
    goto fail;
  }

  if (!denoiser_init_early_objects(self, noise_profile, fft_size)) {
    goto fail;
  }

  // Frame-rate normalization: fixed-ms / fixed-Hz geometry + per-hop alphas.
  const float hop_sec = self->config.hop_sec;
  // hop = frame/overlap_factor (exact with explicit hop, fft-derived legacy)
  const float frame_ms = hop_sec * (float)overlap_factor * 1000.0F;
  const float bin_hz =
      sb_bin_hz(self->config.sample_rate, self->config.fft_size);
  const SbNlmGeometry nlm_geo =
      sb_nlm_geometry_for_frame_ms(frame_ms, hop_sec, bin_hz);

  if (!denoiser_init_filters(self, overlap_factor, hop_sec, bin_hz, nlm_geo)) {
    goto fail;
  }

  // Erb-scaled frequency halo for the neighborhood statistic: built once
  // from the FFT bin pitch (init-time allocation, read-only afterwards)
  self->engines.neighbor_freq_kernel = spectral_neighbors_freq_kernel_create(
      bin_hz, self->config.real_spectrum_size);
  if (!self->engines.neighbor_freq_kernel) {
    goto fail;
  }

  if (!denoiser_init_transient_tail(self, hop_sec)) {
    goto fail;
  }

  self->config.low_latency = low_latency;
  self->mode.active_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  self->mode.pending_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  self->mode.previous_mode = SPECBLEACH_SMOOTHING_TEMPORAL;

  return self;

fail:
  spectral_engine_free(self);
  return NULL;
}

SpectralProcessorHandle spectral_engine_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, NoiseProfile* noise_profile) {
  return spectral_engine_initialize_inner(sample_rate, fft_size, overlap_factor,
                                          0U, noise_profile, false);
}

SpectralProcessorHandle spectral_engine_initialize_with_hop(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, const uint32_t hop_samples,
    NoiseProfile* noise_profile, const bool low_latency) {
  return spectral_engine_initialize_inner(sample_rate, fft_size, overlap_factor,
                                          hop_samples, noise_profile,
                                          low_latency);
}

static void denoiser_free_objects(SbSpectralEngine* self) {
  if (self->engines.noise_estimator) {
    noise_estimation_free(self->engines.noise_estimator);
  }
  if (self->engines.adaptive_estimator) {
    adaptive_estimator_free(self->engines.adaptive_estimator);
  }
  if (self->engines.nlm_filter) {
    nlm_filter_free(self->engines.nlm_filter);
  }
  if (self->engines.bm3d_filter) {
    bm3d_filter_free(self->engines.bm3d_filter);
  }
  if (self->engines.dftt_filter) {
    dftt_filter_free(self->engines.dftt_filter);
  }
  if (self->engines.spectrum_smoothing) {
    spectral_smoothing_free(self->engines.spectrum_smoothing);
  }
  if (self->engines.release_shaper) {
    release_shaper_free(self->engines.release_shaper);
  }
  free(self->engines.release_scale);
  if (self->engines.spectral_features) {
    spectral_features_free(self->engines.spectral_features);
  }
  if (self->engines.neighbor_freq_kernel) {
    spectral_neighbors_freq_kernel_free(self->engines.neighbor_freq_kernel);
  }
  if (self->engines.masking_veto) {
    masking_veto_free(self->engines.masking_veto);
  }
  if (self->engines.suppression_engine) {
    suppression_engine_free(self->engines.suppression_engine);
  }
  if (self->engines.noise_floor_manager) {
    noise_floor_manager_free(self->engines.noise_floor_manager);
  }

  if (self->engines.critical_bands) {
    critical_bands_free(self->engines.critical_bands);
  }
  if (self->engines.transient_detector) {
    transient_detector_free(self->engines.transient_detector);
  }
}

static void denoiser_free_buffers(SbSpectralEngine* self) {
  if (self->transient.band_energies) {
    free(self->transient.band_energies);
  }
  if (self->transient.onset_weights) {
    free(self->transient.onset_weights);
  }
  if (self->transient.held_weights) {
    free(self->transient.held_weights);
  }
  if (self->transient.transient_mask) {
    free(self->transient.transient_mask);
  }
  if (self->transient.transient_band_mask) {
    free(self->transient.transient_band_mask);
  }
  if (self->transient.clean_magnitude) {
    free(self->transient.clean_magnitude);
  }
  if (self->transient.smoothed_magnitude) {
    free(self->transient.smoothed_magnitude);
  }
  if (self->transient.knee_spectrum) {
    free(self->transient.knee_spectrum);
  }

  free(self->spectra.snr_frame);
  free(self->spectra.smoothed_snr);
  free(self->spectra.dftt_snr);
  free(self->spectra.snr_delayed);
  free(self->spectra.gain_spectrum);
  free(self->spectra.gain_spectrum_b);
  free(self->spectra.noise_spectrum);
  if (self->spectra.noise_spectrum_buffers[0]) {
    free(self->spectra.noise_spectrum_buffers[0]);
  }
  if (self->spectra.noise_spectrum_buffers[1]) {
    free(self->spectra.noise_spectrum_buffers[1]);
  }
  free(self->gains.alpha);
  free(self->gains.beta);
  free(self->gains.alpha_b);
  free(self->gains.beta_b);
#if TRANSIENT_RELIEF_PARALLEL
  free(self->gains.alpha_relief);
  free(self->gains.gain_relief);
#endif
  free(self->spectra.noise_bb);
  free(self->spectra.noise_tonal);
  free(self->spectra.gain_tonal);
  free(self->spectra.neighbor_power_row);
  free(self->spectra.neighbor_acc);
  free(self->spectra.blend_power);
  if (self->spectra.manual_noise_floor) {
    free(self->spectra.manual_noise_floor);
  }
}

void spectral_engine_free(SpectralProcessorHandle instance) {
  SbSpectralEngine* self = (SbSpectralEngine*)instance;

  if (!self) {
    return;
  }

  denoiser_free_objects(self);
  denoiser_free_buffers(self);

  if (self->engines.circular_buffer) {
    spectral_circular_buffer_free(self->engines.circular_buffer);
  }

  if (self->engines.tonal_reducer) {
    tonal_reducer_free(self->engines.tonal_reducer);
  }

  free(self);
}

static void denoiser_ensure_adaptive_estimator(
    SbSpectralEngine* self, SpectralEngineParameters parameters) {
  // Check if we need to initialize or re-initialize the adaptive estimator
  if (!parameters.adaptive_noise) {
    return;
  }
  AdaptiveNoiseEstimationMethod requested_method =
      (AdaptiveNoiseEstimationMethod)parameters.noise_estimation_method;

  bool needs_init = !self->engines.adaptive_estimator ||
                    adaptive_estimator_get_method(
                        self->engines.adaptive_estimator) != requested_method;

  if (!needs_init) {
    return;
  }
  adaptive_estimator_free(self->engines.adaptive_estimator);
  self->engines.adaptive_estimator = adaptive_estimator_initialize(
      self->config.real_spectrum_size, self->config.sample_rate,
      self->config.fft_size, requested_method);
  if (self->engines.adaptive_estimator && self->config.hop_sec > 0.0F) {
    adaptive_estimator_set_hop_sec(self->engines.adaptive_estimator,
                                   self->config.hop_sec);
  }
  self->mode.last_adaptive_state = 0;
}

static void denoiser_update_transition(SbSpectralEngine* self,
                                       SpectralEngineParameters parameters) {
  // Runtime smoothing mode switching (allocation-free): the outgoing mode is
  // crossfaded against the incoming one over SMOOTHING_TRANSITION_SECONDS.
  // Within the 2D family (NLM <-> NLM+DFTT <-> BM3D) the switch is instant:
  // all sides share NLM history, latency and DFTT rings (pushed on every 2D
  // pass), so only the map source flips — no crossfade needed.
  const int requested = normalize_smoothing_mode(parameters.smoothing_mode);
  if (!self->mode.in_transition) {
    if (requested == self->mode.active_mode) {
      return;
    }
    if (is_2d_family(requested) && is_2d_family(self->mode.active_mode)) {
      // Crossing BM3D changes the map producer feeding the DFTT rings;
      // reset so the refinement only ever sees NLM priors (it falls back
      // to the raw NLM output until the history refills).
      if ((requested == SPECBLEACH_SMOOTHING_BM3D) !=
          (self->mode.active_mode == SPECBLEACH_SMOOTHING_BM3D)) {
        dftt_filter_reset(self->engines.dftt_filter);
      }
      self->mode.active_mode = requested;
      self->mode.pending_mode = requested;
      return;
    }
    self->mode.previous_mode = self->mode.active_mode;
    self->mode.pending_mode = requested;
    self->mode.transition_pos = 0U;
    self->mode.in_transition = true;
    return;
  }
  if (requested != self->mode.pending_mode &&
      requested == self->mode.previous_mode) {
    // Reverting to the mode that is fading out: mirror the in-progress
    // crossfade around its current blend point so the transition reverses
    // smoothly toward the original chain without a gain discontinuity
    const int outgoing = self->mode.pending_mode;
    self->mode.pending_mode = self->mode.previous_mode;
    self->mode.previous_mode = outgoing;
    self->mode.transition_pos =
        self->mode.transition_frames - self->mode.transition_pos;
    // The chain-slot mapping is keyed on previous_mode, so the two tonal
    // gain histories must swap along with the mode metadata.
    tonal_reducer_swap_gain_slots(self->engines.tonal_reducer);
    return;
  }
  self->mode.pending_mode = requested;
}

static void denoiser_push_filter_params(SbSpectralEngine* self,
                                        SpectralEngineParameters parameters) {
  // Update NLM/BM3D h parameter based on smoothing factor
  const float h_value = (parameters.smoothing_factor > 0.0F)
                            ? (0.5F + (parameters.smoothing_factor * 4.5F))
                            : 0.0F;
  if (self->engines.nlm_filter) {
    nlm_filter_set_h_parameter(self->engines.nlm_filter, h_value);
  }
  if (self->engines.bm3d_filter) {
    bm3d_filter_set_h_parameter(self->engines.bm3d_filter, h_value);
  }

  // Update DFTT refinement strength (reduction-depth coupling, live)
  if (self->engines.dftt_filter) {
    dftt_filter_set_strength(self->engines.dftt_filter,
                             parameters.dftt_strength);
  }
}

bool spectral_engine_load_parameters(SpectralProcessorHandle instance,
                                     SpectralEngineParameters parameters) {
  if (!instance) {
    return false;
  }

  SbSpectralEngine* self = (SbSpectralEngine*)instance;

  denoiser_ensure_adaptive_estimator(self, parameters);

  self->config.parameters = parameters;
  if (self->config.low_latency) {
    self->config.parameters.smoothing_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
    self->mode.active_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
    self->mode.pending_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
    self->mode.previous_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
    self->mode.in_transition = false;
  }

  denoiser_update_transition(self, parameters);
  denoiser_push_filter_params(self, parameters);

  return true;
}

static bool denoiser_run_silence_bypass(SbSpectralEngine* self,
                                        float* fft_spectrum,
                                        const float* reference_spectrum) {
  // Silence bypass: input essentially silent — profile update already done
  // so aggressiveness/threshold stay responsive, but skip the heavy chain.
  float max_val = 0.0f;
  for (uint32_t k = 0U; k < self->config.real_spectrum_size; k++) {
    if (reference_spectrum[k] > max_val) {
      max_val = reference_spectrum[k];
    }
    if (max_val > 1e-12F) {
      break;
    }
  }
  if (max_val > 1e-12F) {
    return false;
  }
  // Keep circular buffer aligned as in idle path
  align_bypass_frame(self, fft_spectrum, reference_spectrum);

  // Safely publish noise spectrum to inactive double buffer via SPSC
  // atomic release
  int published_idx = atomic_load_explicit(&self->spectra.active_noise_idx,
                                           memory_order_relaxed);
  int write_noise_idx = 1 - published_idx;
  memcpy(self->spectra.noise_spectrum_buffers[write_noise_idx],
         self->spectra.noise_spectrum,
         self->config.real_spectrum_size * sizeof(float));
  atomic_store_explicit(&self->spectra.active_noise_idx, write_noise_idx,
                        memory_order_release);

  return true;
}

static void denoiser_run_2d_pass(SbSpectralEngine* self,
                                 const float* fft_spectrum,
                                 const float* reference_spectrum,
                                 DenoiserAlignedFrames* frames) {
  // 2.2 Align internal state and output to the common delayed frame
  // (skipped in low-latency mode: causal, zero look-ahead)
  frames->spectrum = fft_spectrum;
  frames->noise = self->spectra.noise_spectrum;
  frames->noise_bb = self->spectra.noise_bb;
  frames->noise_tonal = self->spectra.noise_tonal;
  frames->tonal_mask = tonal_reducer_get_mask(self->engines.tonal_reducer);
  frames->nlm_smoothed = NULL;
  if (self->config.low_latency) {
    return;
  }
  push_aligned_frame_layers(self, fft_spectrum);

  // Keep both 2D histories warm, even when temporal mode is active.
  push_filter_histories(self, reference_spectrum);

  const uint32_t nlm_delay =
      nlm_filter_get_latency_frames(self->engines.nlm_filter);

  // 2D smoothing (runs when a 2D mode is the active or the incoming mode).
  // The smoothed magnitude is captured explicitly so the temporal chain
  // cannot overwrite the shared alignment layer before the 2D chain
  // consumes it.
  const bool mode_2d_needed =
      is_2d_family(self->mode.active_mode) ||
      (self->mode.in_transition && is_2d_family(self->mode.pending_mode));
  // DFTT refinement follows the DFTT mode: the active chain, or the incoming
  // side of a temporal crossfade. Rings are pushed on every 2D pass so they
  // stay warm for instant intra-family flips.
  const bool use_dftt =
      (self->mode.active_mode == SPECBLEACH_SMOOTHING_NLM_2D_DFTT) ||
      (self->mode.in_transition &&
       self->mode.pending_mode == SPECBLEACH_SMOOTHING_NLM_2D_DFTT);
  // BM3D source follows the same active/incoming rule; DFTT modes always
  // read NLM so the refinement rings stay valid.
  const bool use_bm3d = (self->mode.active_mode == SPECBLEACH_SMOOTHING_BM3D) ||
                        (self->mode.in_transition &&
                         self->mode.pending_mode == SPECBLEACH_SMOOTHING_BM3D);
  const bool filter_ran =
      mode_2d_needed &&
      (use_bm3d ? bm3d_filter_process(self->engines.bm3d_filter,
                                      self->spectra.smoothed_snr)
                : nlm_filter_process(self->engines.nlm_filter,
                                     self->spectra.smoothed_snr));

  // Retrieve unified aligned frames at the common delay
  frames->spectrum = spectral_circular_buffer_retrieve(
      self->engines.circular_buffer, self->layers.layer_fft, nlm_delay);
  frames->noise = spectral_circular_buffer_retrieve(
      self->engines.circular_buffer, self->layers.layer_noise, nlm_delay);
  frames->noise_bb = spectral_circular_buffer_retrieve(
      self->engines.circular_buffer, self->layers.layer_noise_bb, nlm_delay);
  frames->noise_tonal = spectral_circular_buffer_retrieve(
      self->engines.circular_buffer, self->layers.layer_noise_tonal, nlm_delay);
  frames->tonal_mask = spectral_circular_buffer_retrieve(
      self->engines.circular_buffer, self->layers.layer_tonal_mask, nlm_delay);

  if (!frames->spectrum) {
    frames->spectrum = fft_spectrum;
  }
  if (!frames->noise) {
    frames->noise = self->spectra.noise_spectrum;
  }
  if (!frames->noise_bb) {
    frames->noise_bb = self->spectra.noise_bb;
  }
  if (!frames->noise_tonal) {
    frames->noise_tonal = self->spectra.noise_tonal;
  }
  if (!frames->tonal_mask) {
    frames->tonal_mask = tonal_reducer_get_mask(self->engines.tonal_reducer);
  }

  if (!filter_ran) {
    return;
  }
  // Noisy SNR row aligned with the 2D-emitted frame, recomputed from the
  // delayed frames so it describes the same tile the filters just
  // emitted. Feeds the DFTT rings and the confidence blend below.
  // Reuses the shared spectral_features scratch (the temporal chain
  // recomputes it on the same delayed frame anyway).
  const float* delayed_reference =
      get_spectral_feature(self->engines.spectral_features, frames->spectrum,
                           self->config.fft_size, self->config.spectrum_type);
  nlm_filter_calculate_snr(self->engines.nlm_filter, delayed_reference,
                           frames->noise_bb, self->spectra.snr_delayed);
  // DFTT post-filter (paper S4.2): the noisy SNR row aligned with the
  // 2D-emitted frame — recomputed from the delayed frames so both ring
  // inputs describe the same tile — is refined while the NLM output sets
  // the suppression threshold. Falls back to the raw NLM output until the
  // DFTT history is full, or when the active mode is NLM-only.
  float* post_nlm = self->spectra.smoothed_snr;
  if (self->engines.dftt_filter) {
    dftt_filter_push(self->engines.dftt_filter, self->spectra.snr_delayed,
                     self->spectra.smoothed_snr);
    if (use_dftt && dftt_filter_process(self->engines.dftt_filter,
                                        self->spectra.dftt_snr)) {
      post_nlm = self->spectra.dftt_snr;
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
  if (post_nlm == self->spectra.smoothed_snr) {
    const float conf_sq = SMOOTHING_CONFIDENCE_SNR * SMOOTHING_CONFIDENCE_SNR;
    const float conf_4th = conf_sq * conf_sq;
    const uint32_t blend_n = self->config.real_spectrum_size;
    for (uint32_t k = 0U; k < blend_n; k++) {
      const float e = post_nlm[k];
      const float e_sq = e * e;
      const float e_4th = e_sq * e_sq;
      const float raw_w = e_4th / (e_4th + conf_4th);
      post_nlm[k] =
          (raw_w * self->spectra.snr_delayed[k]) + ((1.0F - raw_w) * e);
    }
  }
  nlm_filter_reconstruct_magnitude(self->engines.nlm_filter, post_nlm,
                                   frames->noise_bb, self->spectra.snr_frame);
  spectral_circular_buffer_push(self->engines.circular_buffer,
                                self->layers.layer_smoothed,
                                self->spectra.snr_frame);
  frames->nlm_smoothed = self->spectra.snr_frame;
}

static void denoiser_run_transient(SbSpectralEngine* self,
                                   const float* delayed_spectrum,
                                   const float* delayed_noise) {
  // 2.3 Transient Detection via Transient Detector across Critical Bands.
  // Runs on the ALIGNED (delayed) frame so the transient mask describes the
  // same frame the gains below modify. Detecting on the current input frame
  // would fire ~46 ms before the transient is emitted and open gains on the
  // noise-only frames around it (audible noise pumping at every transient).
  // A clean signal estimate with scaled-up noise subtraction avoids false
  // triggering from residual musical noise.
  bool transient_enabled =
      (self->config.parameters.transient_protection_enable != 0);
  if (!transient_enabled || !self->engines.critical_bands ||
      !self->engines.transient_detector) {
    self->transient.is_transient_detected = false;
    self->transient.transient_intensity = 0.0f;
    self->transient.transient_protection_active = false;
    self->transient.transient_hold_remaining = 0U;
    if (self->engines.critical_bands) {
      uint32_t bands =
          get_number_of_critical_bands(self->engines.critical_bands);
      for (uint32_t b = 0; b < bands; ++b) {
        self->transient.held_weights[b] = 0.0F;
      }
    }
    memset(self->transient.transient_mask, 0,
           self->config.real_spectrum_size * sizeof(float));
    memset(self->transient.transient_band_mask, 0,
           self->config.real_spectrum_size * sizeof(float));
    return;
  }
  const float* delayed_magnitude =
      get_spectral_feature(self->engines.spectral_features, delayed_spectrum,
                           self->config.fft_size, self->config.spectrum_type);
  for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
    // Scale noise up using TRANSIENT_CLEAN_NOISE_SCALE to eliminate spurious
    // noise peaks
    float clean = fmaxf(
        delayed_magnitude[k] - (TRANSIENT_CLEAN_NOISE_SCALE * delayed_noise[k]),
        0.0F);
    self->transient.clean_magnitude[k] = clean;
  }

  compute_critical_bands_spectrum(self->engines.critical_bands,
                                  self->transient.clean_magnitude,
                                  self->transient.band_energies);
  self->transient.is_transient_detected = transient_detector_process(
      self->engines.transient_detector, self->transient.band_energies,
      self->transient.onset_weights, &self->transient.transient_intensity);

  // Hold: keep the fired band weights alive (decayed) across the transient
  // decay tail, so oversubtraction relief outlives the detector's own
  // trigger window instead of cutting the tail off as soon as the SNR
  // drops below the trigger.
  if (self->transient.is_transient_detected) {
    self->transient.transient_hold_remaining =
        (self->transient.transient_hold_decay > 0.0F)
            ? (uint32_t)((TRANSIENT_HOLD_SEC / self->config.hop_sec) + 0.5F)
            : 0U;
  } else if (self->transient.transient_hold_remaining > 0U) {
    self->transient.transient_hold_remaining--;
  }
  self->transient.transient_protection_active =
      self->transient.is_transient_detected ||
      self->transient.transient_hold_remaining > 0U;
  uint32_t num_bands =
      get_number_of_critical_bands(self->engines.critical_bands);
  for (uint32_t b = 0; b < num_bands; ++b) {
    self->transient.held_weights[b] = fmaxf(
        self->transient.onset_weights[b],
        self->transient.held_weights[b] * self->transient.transient_hold_decay);
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
  memset(self->transient.transient_mask, 0,
         self->config.real_spectrum_size * sizeof(float));
  memset(self->transient.transient_band_mask, 0,
         self->config.real_spectrum_size * sizeof(float));
  for (uint32_t b = 0; b < num_bands; ++b) {
    float bw = self->transient.held_weights[b];
    if (bw <= 0.0F) {
      continue;
    }
    CriticalBandIndexes idx = get_band_indexes(self->engines.critical_bands, b);
    uint32_t end = (idx.end_position < self->config.real_spectrum_size)
                       ? idx.end_position
                       : self->config.real_spectrum_size;
    for (uint32_t k = idx.start_position; k < end; ++k) {
      self->transient.transient_band_mask[k] =
          fmaxf(self->transient.transient_band_mask[k], bw);
      float evidence = self->transient.clean_magnitude[k] /
                       (delayed_magnitude[k] + TRANSIENT_BIN_EVIDENCE_EPS);
      self->transient.transient_mask[k] =
          fmaxf(self->transient.transient_mask[k], bw * evidence);
    }
  }
}

// chains write gain_b via DenoiserChainOut.out structs; clang-tidy's
// pointer-to-const hint does not see that write path.
static void denoiser_run_dispatch(
    SbSpectralEngine* self, const float* fft_spectrum,
    const DenoiserAlignedFrames* frames, float* gain_a,
    float* gain_b) { // NOLINT(readability-non-const-parameter)
  // 3. Denoising Stage: dispatch the active smoothing strategy
  // both during a runtime mode transition; low-latency is always temporal)
  if (!self->config.low_latency && self->mode.in_transition) {
    const float total = (float)self->mode.transition_frames;
    const float w = (float)self->mode.transition_pos / total; // 0 → 1

    DenoiserChainFrames frames_in = {
        .smoothed_magnitude = frames->nlm_smoothed,
        .noise_bb = frames->noise_bb,
        .noise_tonal = frames->noise_tonal,
        .tonal_mask = frames->tonal_mask,
    };
    if (is_2d_family(self->mode.previous_mode)) {
      DenoiserChainOut out_a = {.gain_out = gain_a,
                                .alpha = self->gains.alpha,
                                .beta = self->gains.beta};
      DenoiserChainOut out_b = {.gain_out = gain_b,
                                .alpha = self->gains.alpha_b,
                                .beta = self->gains.beta_b};
      (void)run_nlm_chain(self, fft_spectrum, frames_in, 0U, out_a);
      run_temporal_chain(self, fft_spectrum, frames_in, 1U, out_b);
    } else {
      DenoiserChainOut out_a = {.gain_out = gain_a,
                                .alpha = self->gains.alpha,
                                .beta = self->gains.beta};
      DenoiserChainOut out_b = {.gain_out = gain_b,
                                .alpha = self->gains.alpha_b,
                                .beta = self->gains.beta_b};
      run_temporal_chain(self, fft_spectrum, frames_in, 0U, out_a);
      (void)run_nlm_chain(self, fft_spectrum, frames_in, 1U, out_b);
    }

    for (uint32_t k = 0U; k < self->config.fft_size; ++k) {
      gain_a[k] = (gain_a[k] * (1.0F - w)) + (gain_b[k] * w);
    }

    self->mode.transition_pos++;
    if (self->mode.transition_pos >= self->mode.transition_frames) {
      self->mode.active_mode = self->mode.pending_mode;
      self->mode.in_transition = false;
      // The incoming chain built its one-pole tonal gain state in slot 1 while
      // fading in; it now becomes the active chain and reads slot 0.
      tonal_reducer_promote_gain_slot(self->engines.tonal_reducer);
    }
    return;
  }
  if (!self->config.low_latency && is_2d_family(self->mode.active_mode)) {
    DenoiserChainFrames frames_in = {
        .smoothed_magnitude = frames->nlm_smoothed,
        .noise_bb = frames->noise_bb,
        .noise_tonal = frames->noise_tonal,
        .tonal_mask = frames->tonal_mask,
    };
    DenoiserChainOut out_a = {.gain_out = gain_a,
                              .alpha = self->gains.alpha,
                              .beta = self->gains.beta};
    (void)run_nlm_chain(self, fft_spectrum, frames_in, 0U, out_a);
    return;
  }
  DenoiserChainFrames frames_in = {
      .smoothed_magnitude = frames->nlm_smoothed,
      .noise_bb = frames->noise_bb,
      .noise_tonal = frames->noise_tonal,
      .tonal_mask = frames->tonal_mask,
  };
  DenoiserChainOut out_a = {
      .gain_out = gain_a, .alpha = self->gains.alpha, .beta = self->gains.beta};
  run_temporal_chain(self, fft_spectrum, frames_in, 0U, out_a);
}

bool spectral_engine_run(SpectralProcessorHandle instance,
                         float* fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  SbSpectralEngine* self = (SbSpectralEngine*)instance;

  // 1. Preparation: Get reference spectrum and handle learning mode
  float* reference_spectrum =
      get_spectral_feature(self->engines.spectral_features, fft_spectrum,
                           self->config.fft_size, self->config.spectrum_type);

  if (denoiser_profile_core_handle_learning_mode(
          self->engines.noise_estimator, reference_spectrum,
          self->config.parameters.learn_noise, &self->mode.was_learning)) {
    return true;
  }

  // 2. Noise Estimation: Update noise profile (Adaptive or Manual)
  DenoiserProfileCoreParams profile_params = {
      .adaptive_enabled = self->config.parameters.adaptive_noise,
      .spectrum_size = self->config.real_spectrum_size,
      .aggressiveness = &self->mode.aggressiveness,
      .param_aggressiveness = self->config.parameters.aggressiveness,
      .last_adaptive_state = &self->mode.last_adaptive_state,
      .adaptive_estimator = self->engines.adaptive_estimator,
      .noise_profile = self->engines.noise_profile,
      .manual_noise_floor = self->spectra.manual_noise_floor,
      .noise_spectrum = self->spectra.noise_spectrum,
      .noise_estimator = self->engines.noise_estimator,
      .noise_profile_offset_linear =
          self->config.parameters.noise_profile_offset_linear,
      .tonal_noise_profile_offset_linear =
          self->config.parameters.tonal_noise_profile_offset_linear,
      // Previous frame's tonal mask (one-frame latency) is used to apply the
      // tonal threshold offset at detected tonal bins
      .tonal_mask = tonal_reducer_get_mask(self->engines.tonal_reducer),
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
  const float split_reduction_gain = self->config.parameters.tonal_reduction;
#else
  const float split_reduction_gain = 1.0f;
#endif
  tonal_reducer_compute_split(
      self->engines.tonal_reducer, self->spectra.noise_spectrum,
      get_noise_profile(self->engines.noise_profile, CV_MASK),
      is_noise_estimation_available(self->engines.noise_profile, CV_MASK),
      split_reduction_gain, self->spectra.noise_bb, self->spectra.noise_tonal);

  // Idle bypass: no manual profile and not adaptive → skip the heavy chain.
  // Preserve latency and buffer state: push current frame, output the frame
  // delayed by the common lookahead, and advance the circular buffer so
  // idle→active transitions stay aligned.
  if (!self->config.parameters.adaptive_noise &&
      !is_noise_estimation_available(self->engines.noise_profile,
                                     ROLLING_MEAN) &&
      !is_noise_estimation_available(self->engines.noise_profile, MEDIAN) &&
      !is_noise_estimation_available(self->engines.noise_profile, STD_DEV) &&
      !is_noise_estimation_available(self->engines.noise_profile, CV_MASK)) {
    align_bypass_frame(self, fft_spectrum, reference_spectrum);
    return true;
  }

  // Silence bypass: input essentially silent — profile update already done
  // so aggressiveness/threshold stay responsive, but skip the heavy chain.
  if (denoiser_run_silence_bypass(self, fft_spectrum, reference_spectrum)) {
    return true;
  }

  // 2.2 Align internal state and output to the common delayed frame
  // (skipped in low-latency mode: causal, zero look-ahead)
  DenoiserAlignedFrames frames;
  denoiser_run_2d_pass(self, fft_spectrum, reference_spectrum, &frames);
  const float* delayed_spectrum = frames.spectrum;
  const float* delayed_noise = frames.noise;
  const float* nlm_smoothed = frames.nlm_smoothed;
  (void)nlm_smoothed;

  // Align output to delayed frame for post-processing
  if (!self->config.low_latency && delayed_spectrum != fft_spectrum) {
    memcpy(fft_spectrum, delayed_spectrum,
           self->config.fft_size * sizeof(float));
  }

  denoiser_run_transient(self, delayed_spectrum, delayed_noise);

  // 3. Denoising Stage: dispatch the active smoothing strategy
  // both during a runtime mode transition; low-latency is always temporal)
  float* gain_a = self->spectra.gain_spectrum;
  float* gain_b = self->spectra.gain_spectrum_b;

  denoiser_run_dispatch(self, fft_spectrum, &frames, gain_a, gain_b);

  // 4. Post-Processing: Final gain management and mixing
  DenoiserPostProcessParams post_params = {
      .fft_size = self->config.fft_size,
      .real_spectrum_size = self->config.real_spectrum_size,
      .reduction_amount = self->config.parameters.reduction_amount,
      .tonal_reduction = self->config.parameters.tonal_reduction,
      .whitening_factor = self->config.parameters.whitening_factor,
      .residual_listen = self->config.parameters.residual_listen,
      .noise_floor_manager = self->engines.noise_floor_manager,
      .tonal_reducer = self->engines.tonal_reducer,
      .gain_spectrum = self->spectra.gain_spectrum,
      .noise_spectrum = delayed_noise,
      .fft_spectrum = fft_spectrum,
      .reduction_curve_bias = self->config.parameters.reduction_curve_bias,
  };

  denoiser_post_process_apply(post_params);

  // Finalize: Advance circular buffer write index
  if (!self->config.low_latency) {
    spectral_circular_buffer_advance(self->engines.circular_buffer);
  }

  // Safely publish noise spectrum to inactive double buffer via SPSC atomic
  // release
  int published_idx = atomic_load_explicit(&self->spectra.active_noise_idx,
                                           memory_order_relaxed);
  int write_noise_idx = 1 - published_idx;
  memcpy(self->spectra.noise_spectrum_buffers[write_noise_idx],
         self->spectra.noise_spectrum,
         self->config.real_spectrum_size * sizeof(float));
  atomic_store_explicit(&self->spectra.active_noise_idx, write_noise_idx,
                        memory_order_release);

  return true;
}

/**
 * 2D Non-Local Means chain: NLM-smoothed magnitude feeds broadband
 * suppression (Berouti alpha) and a Wiener evaluation against the
 * BROADBAND noise profile only. The tonal residual is handled by the
 * parallel tonal gain path; the final bin gain is
 * min(broadband, tonal) with the transient floor re-asserted on top and
 * the masking audibility floor relieving it last.
 * When the tonal path is inactive (reduction >= ~1.0 or zero mask) the tonal
 * gain is unity and the min() is a no-op.
 */
static bool run_nlm_chain(SbSpectralEngine* self, const float* fft_spectrum,
                          DenoiserChainFrames frames, uint32_t slot,
                          DenoiserChainOut out) {
  const float* smoothed_magnitude = frames.smoothed_magnitude;
  float* gain_out = out.gain_out;
  float* alpha = out.alpha;
  float* beta = out.beta;
  if (!smoothed_magnitude) {
    smoothed_magnitude = spectral_circular_buffer_retrieve(
        self->engines.circular_buffer, self->layers.layer_smoothed, 0U);
  }
  if (!smoothed_magnitude) {
    smoothed_magnitude = fft_spectrum;
  }

  // Neighborhood blend: alpha, veto, broadband and tonal gains all consume
  // the neighbor-informed statistic so the blend stays one coherent
  // reduction decision in the 2D chains too.
  smoothed_magnitude =
      denoiser_neighbor_statistic(self, smoothed_magnitude, fft_spectrum);

  // Backward masking needs the frame one hop AFTER the one being processed:
  // the delayed frame (nlm_delay hops back) plus one hop, not the current
  // frame, which is desynchronized by the whole NLM look-ahead.
  const uint32_t nlm_delay =
      nlm_filter_get_latency_frames(self->engines.nlm_filter);
  // The veto consumes power spectra: fetch the future frame's power layer
  // (aligned with the NLM look-ahead) or convert the current input frame.
  const float* future_spectrum = NULL;
  if (nlm_delay > 0U) {
    future_spectrum = spectral_circular_buffer_retrieve(
        self->engines.circular_buffer, self->layers.layer_power,
        nlm_delay - 1U);
  }
  if (!future_spectrum) {
    spectral_neighbors_compute_power_row(self->spectra.neighbor_power_row,
                                         fft_spectrum, self->config.fft_size);
    future_spectrum = self->spectra.neighbor_power_row;
  }

  // Calculate SNR-dependent oversubtraction factors (Alpha/Beta) on the
  // broadband profile: tonal noise no longer depresses the per-bin SNR here.
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->config.parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->engines.suppression_engine,
                               smoothed_magnitude, frames.noise_bb,
                               suppression_params, alpha, beta);

#if !TONAL_DUAL_PATH
  // Legacy coupled path: tonal alpha boost applied inline (the veto moved
  // to the gain domain — audibility floor applied below, after gains).
  tonal_reducer_apply_alpha_boost(self->engines.tonal_reducer, alpha,
                                  self->config.parameters.tonal_reduction);
  for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
    alpha[k] = fminf(alpha[k], ALPHA_MAX_TONAL);
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
  calculate_gains(
      (GainCalcArgs){.real_spectrum_size = self->config.real_spectrum_size,
                     .fft_size = self->config.fft_size,
                     .spectrum = smoothed_magnitude,
                     .noise_spectrum = frames.noise_bb,
                     .gain_spectrum = gain_out,
                     .alpha = alpha,
                     .beta = beta,
                     .type = self->config.gain_calculation_type,
                     .knee = NULL});

#if TONAL_DUAL_PATH
  // Parallel tonal gain path + decision criterion. Ran on the same
  // smoothed magnitude so the min() compares like with like.
  tonal_reducer_compute_tonal_gains(
      self->engines.tonal_reducer, slot, smoothed_magnitude, frames.noise_tonal,
      frames.tonal_mask, self->config.parameters.tonal_reduction,
      self->spectra.gain_tonal);
  for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
    gain_out[k] = fminf(gain_out[k], self->spectra.gain_tonal[k]);
  }
#endif

  // Masking audibility floor applied last on the combined gain: relief is
  // depth-lerped and bounded by the masked-audibility level, so it cannot
  // excavate unsmoothed flicker the way the old alpha modulation did in
  // this chain (which has no gain smoother after calculate_gains).
  masking_veto_apply_floor(self->engines.masking_veto, smoothed_magnitude,
                           frames.noise_bb, future_spectrum, gain_out,
                           self->config.parameters.masking_depth);

  return true;
}

/**
 * 1D temporal/spatial chain: operates on the common delayed frame (uniform
 * shift of the legacy 1D output). Pre-subtraction magnitude smoothing, then
 * suppression/masking/gain, then temporal + spatial gain smoothing.
 */
static void temporal_prepare_magnitude(SbSpectralEngine* self,
                                       const float* delayed_fft,
                                       const float** effective_magnitude,
                                       const float** delayed_magnitude) {
  // Extract magnitude of the delayed frame (reuses the spectral features
  // buffer; the current-frame reference spectrum is no longer needed here)
  *delayed_magnitude =
      get_spectral_feature(self->engines.spectral_features, delayed_fft,
                           self->config.fft_size, self->config.spectrum_type);

  // Pre-Subtraction temporal stabilization on input magnitude ("time
  // smoothing of the signal spectrum"): fixed light one-pole with tau =
  // SPECTRAL_STABILIZATION_HOPS hops — frame-rate invariant by construction.
  // Gated on smoothing: slider 0 must stay a pristine raw-magnitude bypass.
  // Transient bins are spared bin-by-bin during gain calculation and time
  // smoothing.
  if (self->config.parameters.smoothing_factor <= 0.0F) {
    // Pristine bypass: raw magnitude everywhere, no smoothing state
    memcpy(self->transient.smoothed_magnitude, *delayed_magnitude,
           self->config.real_spectrum_size * sizeof(float));
    self->transient.smoothed_magnitude_seeded = false;
    *effective_magnitude = self->transient.smoothed_magnitude;
    return;
  }
  if (!self->transient.smoothed_magnitude_seeded) {
    memcpy(self->transient.smoothed_magnitude, *delayed_magnitude,
           self->config.real_spectrum_size * sizeof(float));
    self->transient.smoothed_magnitude_seeded = true;
  } else {
    const float stabilization_alpha =
        expf(-1.0F / (float)SPECTRAL_STABILIZATION_HOPS);

    for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
      float raw = (*delayed_magnitude)[k];
      float prev = self->transient.smoothed_magnitude[k];
      // Bin-by-bin adaptive smoothing: open immediately on transient bins
      // while keeping full smoothing elsewhere
      float t_w = self->transient.transient_band_mask[k];
      float adapt_alpha = (1.0F - t_w) * stabilization_alpha;
      self->transient.smoothed_magnitude[k] =
          (adapt_alpha * prev) + ((1.0F - adapt_alpha) * raw);
    }
  }

  int spatial_passes = (int)(self->config.parameters.smoothing_factor * 2.0F);
  for (int p = 0; p < spatial_passes; ++p) {
    spectral_smoothing_apply_spatial(self->transient.smoothed_magnitude,
                                     self->config.real_spectrum_size);
  }

  *effective_magnitude = self->transient.smoothed_magnitude;

  // Adaptive release shaping: bands whose recent-energy envelope collapses
  // close fast (no spectral ghosts on the residual); the rest keep the full
  // slider release (anti-chirp protection in noise-only stretches)
  release_shaper_compute(self->engines.release_shaper, *effective_magnitude,
                         self->engines.release_scale);
}

static void temporal_apply_suppression(SbSpectralEngine* self,
                                       const float* effective_magnitude,
                                       DenoiserChainFrames frames,
                                       DenoiserChainOut out) {
  // Calculate SNR-dependent oversubtraction factors (Alpha/Beta) on the
  // broadband profile: tonal noise no longer depresses the per-bin SNR here.
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->config.parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->engines.suppression_engine,
                               effective_magnitude, frames.noise_bb,
                               suppression_params, out.alpha, out.beta);
#if !TONAL_DUAL_PATH
  // Legacy coupled path: tonal alpha boost applied inline (the veto moved
  // to the gain domain — audibility floor applied in
  // temporal_compute_gains, before the gain smoothers).
  tonal_reducer_apply_alpha_boost(self->engines.tonal_reducer, out.alpha,
                                  self->config.parameters.tonal_reduction);
  for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
    out.alpha[k] = fminf(out.alpha[k], ALPHA_MAX_TONAL);
  }
#endif
#if !TRANSIENT_RELIEF_PARALLEL
  // When transients are detected and enabled, drop alphas firmly to ALPHA_MIN
  // (1.0) strictly on the specific frequencies where transient energy was
  // detected. Band-level weight (see run_nlm_chain 3.4).
  if (self->transient.transient_protection_active) {
    for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
      float t_weight = self->transient.transient_band_mask[k];
      if (t_weight > 0.0F) {
        float prot_factor = sqrtf(t_weight);
        out.alpha[k] =
            (out.alpha[k] * (1.0F - prot_factor)) + (ALPHA_MIN * prot_factor);
        out.alpha[k] = fmaxf(ALPHA_MIN, out.alpha[k]);
      }
    }
  }
#endif
}

static void temporal_compute_knee(SbSpectralEngine* self,
                                  const float* effective_magnitude,
                                  const float* delayed_magnitude,
                                  const float* delayed_noise_bb) {
  // Signal-dependent knee width: bins decaying from recent signal presence
  // (stabilized energy above the current raw hop) get a wider knee so weak
  // component tails are forgiven instead of cut; steady or rising bins keep
  // the base knee. Pristine bypass keeps a zero knee everywhere.
  if (self->config.parameters.smoothing_factor <= 0.0F) {
    memset(self->transient.knee_spectrum, 0,
           self->config.real_spectrum_size * sizeof(float));
    return;
  }
  for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
    float decay_evidence = 0.0F;
    if (delayed_noise_bb[k] > FLT_MIN) {
      decay_evidence =
          (effective_magnitude[k] - delayed_magnitude[k]) / delayed_noise_bb[k];
    }
    self->transient.knee_spectrum[k] =
        GAIN_WIENER_KNEE +
        (GAIN_KNEE_DECAY_BOOST *
         fminf(1.0F, fmaxf(0.0F, decay_evidence) / GAIN_KNEE_DECAY_RANGE));
  }
}

static void temporal_compute_gains(SbSpectralEngine* self,
                                   const float* effective_magnitude,
                                   DenoiserChainFrames frames,
                                   DenoiserChainOut out) {
  // Gain Calculation on the broadband profile only: the temporal/spatial
  // smoothers below therefore never see a stationary tonal notch carved into
  // their gain field.
  calculate_gains(
      (GainCalcArgs){.real_spectrum_size = self->config.real_spectrum_size,
                     .fft_size = self->config.fft_size,
                     .spectrum = effective_magnitude,
                     .noise_spectrum = frames.noise_bb,
                     .gain_spectrum = out.gain_out,
                     .alpha = out.alpha,
                     .beta = out.beta,
                     .type = self->config.gain_calculation_type,
                     .knee = self->transient.knee_spectrum});

#if TRANSIENT_RELIEF_PARALLEL
  // Transient relief as a parallel GAIN-domain branch: blend the base gain
  // toward the plain Wiener curve (alpha = ALPHA_MIN, same knee) by the band
  // protection weight. Matches legacy at full/no protection; at partial
  // band weights it is softer than the legacy alpha lerp because the Wiener
  // curve is nonlinear in alpha, and the shared alpha is no longer mutated.
  if (self->transient.transient_protection_active) {
    calculate_gains(
        (GainCalcArgs){.real_spectrum_size = self->config.real_spectrum_size,
                       .fft_size = self->config.fft_size,
                       .spectrum = effective_magnitude,
                       .noise_spectrum = frames.noise_bb,
                       .gain_spectrum = self->gains.gain_relief,
                       .alpha = self->gains.alpha_relief,
                       .beta = out.beta,
                       .type = self->config.gain_calculation_type,
                       .knee = self->transient.knee_spectrum});
    for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
      const float pf = sqrtf(self->transient.transient_band_mask[k]);
      if (pf > 0.0F) {
        out.gain_out[k] =
            ((1.0F - pf) * out.gain_out[k]) + (pf * self->gains.gain_relief[k]);
      }
    }
  }
#endif

  // Transient Protection: ensure transient bins have gain near 1.0
  if (self->transient.transient_protection_active) {
    for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
      float t_weight = self->transient.transient_mask[k];
      if (t_weight > 0.0F) {
        out.gain_out[k] = fmaxf(out.gain_out[k], t_weight);
      }
    }
  }

  // Masking audibility floor: relief joins the gain field BEFORE the
  // time/spatial smoothers so the smoother owns the transition. The floor
  // reports the model's masking threshold (its temporal masking layer owns
  // all time dynamics) and replaces the old alpha modulation that
  // flickered with SFM classification.
  masking_veto_apply_floor(self->engines.masking_veto, effective_magnitude,
                           frames.noise_bb, NULL, out.gain_out,
                           self->config.parameters.masking_depth);
}

static void temporal_smooth_gains(SbSpectralEngine* self,
                                  DenoiserChainOut out) {
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
          .smoothing = fminf(self->config.parameters.smoothing_factor,
                             GAIN_SMOOTHING_RELEASE_P_CAP),
          // Instant attack only on the detected frame; during the relief hold
          // the smoother must keep owning the gain (no held-mask override).
          .transient_mask = self->transient.is_transient_detected
                                ? self->transient.transient_band_mask
                                : NULL,
          .release_scale = (self->config.parameters.smoothing_factor > 0.0F)
                               ? self->engines.release_scale
                               : NULL, // Unused while bypassed
      };
  spectral_smoothing_run(self->engines.spectrum_smoothing,
                         spectral_smoothing_parameters, out.gain_out);
  if (self->config.parameters.smoothing_factor > 0.0f) {
    int passes = 1 + (int)(self->config.parameters.smoothing_factor * 2.0f);
    for (int p = 0; p < passes; ++p) {
      spectral_smoothing_apply_spatial(out.gain_out,
                                       self->config.real_spectrum_size);
    }
  }
}

static void temporal_apply_tonal_min(SbSpectralEngine* self, uint32_t slot,
                                     const float* effective_magnitude,
                                     DenoiserChainFrames frames,
                                     DenoiserChainOut out) {
  // Parallel tonal gain path + decision criterion, applied AFTER the
  // time/spatial smoothers: the notch enters the final gain at full depth
  // without being smeared into neighbors by the spatial FIR, and the
  // smoother state is never dragged around by mask flicker.
#if TONAL_DUAL_PATH
  tonal_reducer_compute_tonal_gains(
      self->engines.tonal_reducer, slot, effective_magnitude,
      frames.noise_tonal, frames.tonal_mask,
      self->config.parameters.tonal_reduction, self->spectra.gain_tonal);
  for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
    out.gain_out[k] = fminf(out.gain_out[k], self->spectra.gain_tonal[k]);
  }
#else
  (void)slot;
  (void)effective_magnitude;
  (void)frames;
  (void)out;
#endif

  // Transient Protection re-asserted after the combine: a band onset must
  // not be notched by a tonal bin underneath it (floor is now strictly >=
  // the legacy pre-smoothing behavior on transient x tonal overlap bins).
  // Keyed to the detector's current frame (the UI LED state), not the 200 ms
  // relief hold, so the smoothed gain wins again through the sustain.
  if (self->transient.is_transient_detected) {
    for (uint32_t k = 0U; k < self->config.real_spectrum_size; ++k) {
      float t_weight = self->transient.transient_mask[k];
      if (t_weight > 0.0F) {
        out.gain_out[k] = fmaxf(out.gain_out[k], t_weight);
      }
    }
  }
}
static void run_temporal_chain(SbSpectralEngine* self, const float* delayed_fft,
                               DenoiserChainFrames frames, uint32_t slot,
                               DenoiserChainOut out) {
  const float* delayed_noise_bb = frames.noise_bb;

  const float* effective_magnitude = NULL;
  const float* delayed_magnitude = NULL;
  temporal_prepare_magnitude(self, delayed_fft, &effective_magnitude,
                             &delayed_magnitude);

  // Publish the smoothed magnitude so a pending switch to NLM starts from an
  // aligned frame
  spectral_circular_buffer_push(self->engines.circular_buffer,
                                self->layers.layer_smoothed,
                                self->transient.smoothed_magnitude);

  // Neighborhood blend: the gain decision statistic (alpha, knee, gains,
  // tonal floor) consumes the neighbor-informed result below; the smoothing
  // state published on layer_smoothed above stays on the unblended statistic.
  effective_magnitude =
      denoiser_neighbor_statistic(self, effective_magnitude, delayed_fft);

  temporal_apply_suppression(self, effective_magnitude, frames, out);
  temporal_compute_knee(self, effective_magnitude, delayed_magnitude,
                        delayed_noise_bb);
  temporal_compute_gains(self, effective_magnitude, frames, out);
  temporal_smooth_gains(self, out);
  temporal_apply_tonal_min(self, slot, effective_magnitude, frames, out);
}

const float* spectral_engine_get_tonal_mask(SpectralProcessorHandle instance) {
  const SbSpectralEngine* self = (const SbSpectralEngine*)instance;
  return (self && self->engines.tonal_reducer)
             ? tonal_reducer_get_mask(self->engines.tonal_reducer)
             : NULL;
}

uint32_t spectral_engine_get_peaks(SpectralProcessorHandle instance,
                                   float* peak_freqs_hz, uint32_t max_peaks) {
  const SbSpectralEngine* self = (const SbSpectralEngine*)instance;
  return (self && self->engines.tonal_reducer)
             ? tonal_reducer_get_peaks(self->engines.tonal_reducer,
                                       peak_freqs_hz, max_peaks)
             : 0;
}

const float* spectral_engine_get_active_noise_profile(
    SpectralProcessorHandle instance) {
  const SbSpectralEngine* self = (const SbSpectralEngine*)instance;
  if (!self) {
    return NULL;
  }
  int idx = atomic_load_explicit(&self->spectra.active_noise_idx,
                                 memory_order_acquire);
  return self->spectra.noise_spectrum_buffers[idx];
}

void spectral_engine_reset_noise_profile(SpectralProcessorHandle instance) {
  SbSpectralEngine* self = (SbSpectralEngine*)instance;
  if (!self) {
    return;
  }
  if (self->engines.noise_estimator) {
    noise_estimation_reset(self->engines.noise_estimator);
  }
  if (self->engines.tonal_reducer) {
    tonal_reducer_reset(self->engines.tonal_reducer);
  }
  if (self->spectra.manual_noise_floor) {
    memset(self->spectra.manual_noise_floor, 0,
           self->config.real_spectrum_size * sizeof(float));
  }
  if (self->spectra.noise_spectrum) {
    memset(self->spectra.noise_spectrum, 0,
           self->config.real_spectrum_size * sizeof(float));
  }
  if (self->spectra.noise_spectrum_buffers[0]) {
    memset(self->spectra.noise_spectrum_buffers[0], 0,
           self->config.real_spectrum_size * sizeof(float));
  }
  if (self->spectra.noise_spectrum_buffers[1]) {
    memset(self->spectra.noise_spectrum_buffers[1], 0,
           self->config.real_spectrum_size * sizeof(float));
  }
  if (self->spectra.noise_bb) {
    memset(self->spectra.noise_bb, 0,
           self->config.real_spectrum_size * sizeof(float));
  }
  if (self->spectra.noise_tonal) {
    memset(self->spectra.noise_tonal, 0,
           self->config.real_spectrum_size * sizeof(float));
  }
  self->mode.was_learning = false;
  self->mode.last_adaptive_state = 0;
  self->transient.smoothed_magnitude_seeded = false;
  release_shaper_reset(self->engines.release_shaper);
}

uint32_t spectral_engine_get_latency_frames(SpectralProcessorHandle instance) {
  const SbSpectralEngine* self = (const SbSpectralEngine*)instance;

  if (!self) {
    return 0;
  }

  // Common delay: 2D look-ahead applies to every smoothing mode so the
  // reported latency never changes on a runtime mode switch.
  // Low-latency mode is causal: zero look-ahead by construction.
  if (self->config.low_latency) {
    return 0;
  }
  return nlm_filter_get_latency_frames(self->engines.nlm_filter);
}

bool spectral_engine_is_transient_detected(SpectralProcessorHandle instance) {
  if (!instance) {
    return false;
  }
  const SbSpectralEngine* self = (const SbSpectralEngine*)instance;
  return self->transient.is_transient_detected;
}

float spectral_engine_get_transient_intensity(
    SpectralProcessorHandle instance) {
  if (!instance) {
    return 0.0f;
  }
  const SbSpectralEngine* self = (const SbSpectralEngine*)instance;
  return self->transient.transient_intensity;
}
