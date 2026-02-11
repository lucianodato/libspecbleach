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
#include "shared/denoiser_logic/core/denoise_mixer.h"
#include "shared/denoiser_logic/core/denoiser_core.h"
#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/denoiser_logic/processing/gain_calculator.h"
#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/denoiser_logic/processing/nlm_filter.h"
#include "shared/denoiser_logic/processing/suppression_engine.h"
#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/utils/spectral_utils.h"
#include <float.h>
#include <stdlib.h>
#include <string.h>

typedef struct Spectral2DDenoiser {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t hop;

  Denoiser2DParameters parameters;

  float* snr_frame;          // Current SNR frame for NLM input
  float* smoothed_snr;       // Smoothed SNR output from NLM
  float* gain_spectrum;      // Computed gains
  float* noise_spectrum;     // Copy of noise profile for processing
  float* alpha;              // Oversubtraction factors
  float* beta;               // Undersubtraction factors
  float* manual_noise_floor; // Manual profile floor
  TonalReducer* tonal_reducer;

  // Delay buffer for audio alignment
  float* spectral_delay_buffer;
  float* noise_delay_buffer;
  float* magnitude_delay_buffer;
  uint32_t delay_buffer_write_index;

  SpectrumType spectrum_type;
  GainCalculationType gain_calculation_type;

  NoiseProfile* noise_profile;
  NoiseEstimator* noise_estimator;
  AdaptiveNoiseEstimator* adaptive_estimator;
  NlmFilter* nlm_filter;
  SpectralFeatures* spectral_features;
  MaskingVeto* masking_veto;
  SuppressionEngine* suppression_engine;
  DenoiseMixer* mixer;
  NoiseFloorManager* noise_floor_manager;

  int last_adaptive_state;
  int last_noise_estimation_method;
  float aggressiveness;
  bool was_learning;
} Spectral2DDenoiser;

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
  self->spectrum_type = SPECTRAL_TYPE_GENERAL;
  self->gain_calculation_type = GAIN_ESTIMATION_TYPE;
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

  // Allocate spectral delay buffer
  // We need to store full FFT frames (packed complex)
  // Assuming fft_size floats is sufficient for the packed format used by STFT
  // processor
  self->spectral_delay_buffer =
      (float*)calloc((size_t)DELAY_BUFFER_FRAMES * fft_size, sizeof(float));
  if (!self->spectral_delay_buffer) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->noise_delay_buffer = (float*)calloc(
      (size_t)DELAY_BUFFER_FRAMES * self->real_spectrum_size, sizeof(float));
  if (!self->noise_delay_buffer) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->magnitude_delay_buffer = (float*)calloc(
      (size_t)DELAY_BUFFER_FRAMES * self->real_spectrum_size, sizeof(float));
  if (!self->magnitude_delay_buffer) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  self->delay_buffer_write_index = 0;
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
      self->fft_size, self->sample_rate, self->spectrum_type);
  self->suppression_engine =
      suppression_engine_initialize(self->real_spectrum_size, self->sample_rate,
                                    CRITICAL_BANDS_TYPE, self->spectrum_type);

  if (!self->masking_veto || !self->suppression_engine) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  // Initialize mixer
  self->mixer = denoise_mixer_initialize(fft_size, sample_rate, self->hop);
  if (!self->mixer) {
    spectral_2d_denoiser_free(self);
    return NULL;
  }

  // Initialize noise floor manager
  self->noise_floor_manager =
      noise_floor_manager_initialize(fft_size, sample_rate, self->hop);
  if (!self->noise_floor_manager) {
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
  if (self->mixer) {
    denoise_mixer_free(self->mixer);
  }
  if (self->noise_floor_manager) {
    noise_floor_manager_free(self->noise_floor_manager);
  }

  free(self->snr_frame);
  free(self->smoothed_snr);
  free(self->gain_spectrum);
  free(self->noise_spectrum);
  free(self->alpha);
  free(self->beta);
  free(self->spectral_delay_buffer);
  free(self->noise_delay_buffer);
  free(self->magnitude_delay_buffer);
  free(self->manual_noise_floor);
  tonal_reducer_free(self->tonal_reducer);

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

  // Update NLM h parameter based on smoothing factor
  if (self->nlm_filter && parameters.smoothing_factor > 0.0F) {
    nlm_filter_set_h_parameter(self->nlm_filter, parameters.smoothing_factor);
  }

  return true;
}

static void push_to_delay_buffers(Spectral2DDenoiser* self,
                                  const float* fft_spectrum,
                                  const float* magnitude);
static void process_nlm_input(Spectral2DDenoiser* self,
                              const float* reference_spectrum);
static void calculate_reduction_gains(Spectral2DDenoiser* self,
                                      const float* smoothed_snr,
                                      const float* delayed_noise,
                                      const float* delayed_magnitude);

bool spectral_2d_denoiser_run(SpectralProcessorHandle instance,
                              float* fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;

  // Get reference spectrum (power or magnitude)
  float* reference_spectrum =
      get_spectral_feature(self->spectral_features, fft_spectrum,
                           self->fft_size, self->spectrum_type);

  if (denoiser_core_handle_learning_mode(
          self->noise_estimator, reference_spectrum,
          self->parameters.learn_noise, &self->was_learning)) {
    return true;
  }

  // --- Denoising mode: use NLM for 2D smoothing ---

  // 1. Store current spectral frame in delay buffer
  push_to_delay_buffers(self, fft_spectrum, reference_spectrum);

  // 2. Determine which noise profile to use
  DenoiserCoreProfileParams profile_params = {
      .adaptive_enabled = self->parameters.adaptive_noise,
      .spectrum_size = self->real_spectrum_size,
      .aggressiveness = self->aggressiveness,
      .param_aggressiveness = self->parameters.aggressiveness,
      .last_adaptive_state = &self->last_adaptive_state,
      .adaptive_estimator = self->adaptive_estimator,
      .noise_profile = self->noise_profile,
      .manual_noise_floor = self->manual_noise_floor,
      .noise_spectrum = self->noise_spectrum,
  };
  denoiser_core_update_noise_profile(profile_params, reference_spectrum);

  // Update internal aggressiveness state if it changed during update
  self->aggressiveness = profile_params.aggressiveness;
  if (*profile_params.last_adaptive_state) {
    self->aggressiveness = self->parameters.aggressiveness;
  }

  // 3. Store the noise spectrum in the delay buffer
  // 4. Compute SNR for NLM using the CURRENT adaptive noise (Whitening)
  // 5. Push frame to NLM filter
  process_nlm_input(self, reference_spectrum);

  // 6. Process if NLM buffer is ready
  if (nlm_filter_is_ready(self->nlm_filter)) {
    // Determine retrieve index for aligned data
    uint32_t read_index = (self->delay_buffer_write_index +
                           DELAY_BUFFER_FRAMES - NLM_SEARCH_RANGE_TIME_FUTURE) %
                          DELAY_BUFFER_FRAMES;

    float* delayed_spectrum =
        &self->spectral_delay_buffer[(size_t)read_index * self->fft_size];
    float* delayed_noise = &self->noise_delay_buffer[(size_t)read_index *
                                                     self->real_spectrum_size];
    float* delayed_magnitude_spectrum =
        &self->magnitude_delay_buffer[(size_t)read_index *
                                      self->real_spectrum_size];

    // Moderating the NLM reduction via Masking Veto
    if (nlm_filter_process(self->nlm_filter, self->smoothed_snr)) {
      calculate_reduction_gains(self, self->smoothed_snr, delayed_noise,
                                delayed_magnitude_spectrum);

      DenoiserCorePostProcessParams post_params = {
          .fft_size = self->fft_size,
          .real_spectrum_size = self->real_spectrum_size,
          .reduction_amount = self->parameters.reduction_amount,
          .tonal_reduction = self->parameters.tonal_reduction,
          .whitening_factor = self->parameters.whitening_factor,
          .mixer_whitening_factor = self->parameters.whitening_factor,
          .residual_listen = self->parameters.residual_listen,
          .noise_floor_manager = self->noise_floor_manager,
          .tonal_reducer = self->tonal_reducer,
          .mixer = self->mixer,
          .gain_spectrum = self->gain_spectrum,
          .noise_spectrum = self->noise_spectrum,
          .fft_spectrum =
              fft_spectrum, // Will be overwritten by delayed_spectrum first
      };

      // Copy delayed spectrum to output first as we process the delayed frame
      memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));

      // We need to pass the *delayed* noise spectrum to post processing,
      // but denoiser_core_apply_post_processing uses params.noise_spectrum.
      // In 2D, self->noise_spectrum is the CURRENT frame's noise, but we are
      // processing a DELAYED frame. We must use delayed_noise. However,
      // calculating_reduction_gains used delayed_noise (via pointer mod or
      // copy?) Wait, calculate_reduction_gains puts result in
      // self->gain_spectrum. But NFM needs noise spectrum to fill gaps.
      // *Critical*: We must use delayed_noise in post_params.
      post_params.noise_spectrum = delayed_noise;

      denoiser_core_apply_post_processing(post_params);
    } else {
      // NLM process failed unexpectedly
      memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));
    }
  } else {
    // If NLM not ready yet (startup latency), output the delayed spectrum
    uint32_t read_index = (self->delay_buffer_write_index +
                           DELAY_BUFFER_FRAMES - NLM_SEARCH_RANGE_TIME_FUTURE) %
                          DELAY_BUFFER_FRAMES;
    float* delayed_spectrum =
        &self->spectral_delay_buffer[(size_t)read_index * self->fft_size];
    memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));
  }

  // Advance delay buffer write index
  self->delay_buffer_write_index =
      (self->delay_buffer_write_index + 1) % DELAY_BUFFER_FRAMES;

  return true;
}

static void push_to_delay_buffers(Spectral2DDenoiser* self,
                                  const float* fft_spectrum,
                                  const float* magnitude) {
  memcpy(&self->spectral_delay_buffer[(size_t)self->delay_buffer_write_index *
                                      self->fft_size],
         fft_spectrum, self->fft_size * sizeof(float));

  memcpy(&self->magnitude_delay_buffer[(size_t)self->delay_buffer_write_index *
                                       self->real_spectrum_size],
         magnitude, self->real_spectrum_size * sizeof(float));
}

static void process_nlm_input(Spectral2DDenoiser* self,
                              const float* reference_spectrum) {
  // 3. Store the noise spectrum in the delay buffer
  memcpy(&self->noise_delay_buffer[(size_t)self->delay_buffer_write_index *
                                   self->real_spectrum_size],
         self->noise_spectrum, self->real_spectrum_size * sizeof(float));

  // 4. Compute SNR for NLM using the CURRENT adaptive noise
  for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
    float denom = self->noise_spectrum[k] > FLT_MIN ? self->noise_spectrum[k]
                                                    : SPECTRAL_EPSILON;
    self->snr_frame[k] = reference_spectrum[k] / denom;
  }

  // 5. Push frame to NLM filter
  nlm_filter_push_frame(self->nlm_filter, self->snr_frame);
}

static void calculate_reduction_gains(Spectral2DDenoiser* self,
                                      const float* smoothed_snr,
                                      const float* delayed_noise,
                                      const float* delayed_magnitude) {
  // 1. Convert smoothed SNR back to spectral domain
  // We reuse self->snr_frame as a temp buffer for smoothed_magnitude
  float* smoothed_magnitude = self->snr_frame;
  for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
    float denom =
        delayed_noise[k] > FLT_MIN ? delayed_noise[k] : SPECTRAL_EPSILON;
    smoothed_magnitude[k] = smoothed_snr[k] * denom;
  }

  // 2. Calculate SNR-dependent oversubtraction factors
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->suppression_engine, smoothed_magnitude,
                               delayed_noise, suppression_params, self->alpha,
                               self->beta);

  // 3. Detect tonal components
  tonal_reducer_run(self->tonal_reducer, delayed_noise,
                    get_noise_profile(self->noise_profile, MAX),
                    get_noise_profile(self->noise_profile, MEDIAN), self->alpha,
                    self->parameters.tonal_reduction);

  // 4. Apply psychoacoustic veto
  masking_veto_apply(self->masking_veto, smoothed_magnitude, delayed_magnitude,
                     delayed_noise, self->alpha, 1.0F,
                     self->parameters.nlm_masking_protection,
                     self->parameters.masking_elasticity);

  calculate_gains(self->real_spectrum_size, self->fft_size, smoothed_magnitude,
                  delayed_noise, self->gain_spectrum, self->alpha, self->beta,
                  self->gain_calculation_type);
}

uint32_t spectral_2d_denoiser_get_latency_frames(
    SpectralProcessorHandle instance) {
  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;

  if (!self || !self->nlm_filter) {
    return 0;
  }

  return nlm_filter_get_latency_frames(self->nlm_filter);
}
