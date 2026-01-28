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
#include "shared/gain_estimation/gain_estimators.h"
#include "shared/noise_estimation/noise_estimator.h"
#include "shared/post_estimation/nlm_filter.h"
#include "shared/post_estimation/noise_floor_manager.h"
#include "shared/utils/denoise_mixer.h"
#include "shared/utils/spectral_features.h"
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

  float* snr_frame;      // Current SNR frame for NLM input
  float* smoothed_snr;   // Smoothed SNR output from NLM
  float* gain_spectrum;  // Computed gains
  float* noise_spectrum; // Copy of noise profile for processing
  float* alpha;          // Oversubtraction factors
  float* beta;           // Undersubtraction factors

  // Delay buffer for audio alignment
  float* spectral_delay_buffer;
  uint32_t delay_buffer_write_index;

  SpectrumType spectrum_type;
  GainEstimationType gain_estimation_type;

  NoiseProfile* noise_profile;
  NoiseEstimator* noise_estimator;
  NlmFilter* nlm_filter;
  SpectralFeatures* spectral_features;
  DenoiseMixer* mixer;
  NoiseFloorManager* noise_floor_manager;
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
  self->gain_estimation_type = GAIN_ESTIMATION_TYPE;
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
  self->delay_buffer_write_index = 0;

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
  if (self->nlm_filter) {
    nlm_filter_free(self->nlm_filter);
  }
  if (self->spectral_features) {
    spectral_features_free(self->spectral_features);
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

  free(self);
}

bool load_2d_reduction_parameters(SpectralProcessorHandle instance,
                                  Denoiser2DParameters parameters) {
  if (!instance) {
    return false;
  }

  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;
  self->parameters = parameters;

  // Update NLM h parameter based on smoothing factor
  if (self->nlm_filter && parameters.smoothing_factor > 0.0F) {
    nlm_filter_set_h_parameter(self->nlm_filter, parameters.smoothing_factor);
  }

  return true;
}

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

  if (self->parameters.learn_noise > 0) {
    // Learning mode: update noise profile for all modes
    // No latency compensation needed here as we are just learning (output is
    // passthrough or silenced by host?) Usually plugins bypass processing
    // during learn. But if we passthrough, we should also delay audio to match
    // switch behavior? For simplicity, no delay in learn mode, but usually
    // learning overrides processing.

    for (int mode = ROLLING_MEAN; mode <= MAX; mode++) {
      noise_estimation_run(self->noise_estimator, (NoiseEstimatorType)mode,
                           reference_spectrum);
    }
  } else if (is_noise_estimation_available(
                 self->noise_profile, self->parameters.noise_reduction_mode)) {
    // Denoising mode: use NLM for 2D smoothing

    // 1. Store current spectral frame in delay buffer
    memcpy(&self->spectral_delay_buffer[(size_t)self->delay_buffer_write_index *
                                        self->fft_size],
           fft_spectrum, self->fft_size * sizeof(float));

    // Get noise profile
    memcpy(self->noise_spectrum,
           get_noise_profile(self->noise_profile,
                             self->parameters.noise_reduction_mode),
           self->real_spectrum_size * sizeof(float));

    // Compute SNR frame: snr[k] = spectrum[k] / noise[k]
    for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
      if (self->noise_spectrum[k] > FLT_MIN) {
        self->snr_frame[k] = reference_spectrum[k] / self->noise_spectrum[k];
      } else {
        self->snr_frame[k] = reference_spectrum[k];
      }
    }

    // Push frame to NLM filter
    nlm_filter_push_frame(self->nlm_filter, self->snr_frame);

    // Process if NLM buffer is ready
    if (nlm_filter_process(self->nlm_filter, self->smoothed_snr)) {
      // NLM output corresponds to frame at (current -
      // NLM_SEARCH_RANGE_TIME_FUTURE) because we wait for future frames to be
      // available before processing the target frame.

      // 2. Retrieve delayed spectral frame matching the NLM latency
      uint32_t read_index =
          (self->delay_buffer_write_index + DELAY_BUFFER_FRAMES -
           NLM_SEARCH_RANGE_TIME_FUTURE) %
          DELAY_BUFFER_FRAMES;
      float* delayed_spectrum =
          &self->spectral_delay_buffer[(size_t)read_index * self->fft_size];

      // Convert smoothed SNR back to spectral domain for gain calculation
      // smoothed_magnitude[k] = smoothed_snr[k] * noise[k]
      float* smoothed_magnitude = self->snr_frame; // Reuse buffer
      for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
        smoothed_magnitude[k] = self->smoothed_snr[k] * self->noise_spectrum[k];
      }

      // Initialize alpha/beta for simple gain estimation
      for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
        self->alpha[k] = DEFAULT_OVERSUBTRACTION;
        self->beta[k] = self->parameters.reduction_amount;
      }

      // Estimate gains using smoothed spectrum
      estimate_gains(self->real_spectrum_size, self->fft_size,
                     smoothed_magnitude, self->noise_spectrum,
                     self->gain_spectrum, self->alpha, self->beta,
                     self->gain_estimation_type);

      // Apply noise floor management
      noise_floor_manager_apply(
          self->noise_floor_manager, self->real_spectrum_size, self->fft_size,
          self->gain_spectrum, self->noise_spectrum,
          self->parameters.reduction_amount, self->parameters.whitening_factor);

      // Apply gains via mixer to the DELAYED spectrum
      DenoiseMixerParameters mixer_params = {
          .noise_level = self->parameters.reduction_amount,
          .residual_listen = self->parameters.residual_listen,
          .whitening_amount = 0.0F,
      };

      // Note: denoise_mixer_run modifies the input spectrum array in-place
      // We process the delayed spectrum and copy it back to fft_spectrum output
      // Or we can just process delayed_spectrum and memcpy result to
      // fft_spectrum? denoise_mixer_run(mixer, input_and_output_spectrum,
      // gains, ...)

      // We must not modify the delay buffer in place if we were to read it
      // again (we won't). But let's copy delayed spectrum to fft_spectrum
      // first, then apply mix in place on fft_spectrum.
      memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));

      denoise_mixer_run(self->mixer, fft_spectrum, self->gain_spectrum,
                        mixer_params);
    } else {
      // If NLM not ready yet (startup latency)
      // We should ideally output silence or 0-gain to avoid glitches?
      // Or just pass through current frame (which will be mixed with future
      // ones later... bad). Silence is safer for latency filling.
      memset(fft_spectrum, 0, self->fft_size * sizeof(float));
    }

    // Advance delay buffer write index
    self->delay_buffer_write_index =
        (self->delay_buffer_write_index + 1) % DELAY_BUFFER_FRAMES;
  }

  return true;
}

uint32_t spectral_2d_denoiser_get_latency_frames(
    SpectralProcessorHandle instance) {
  Spectral2DDenoiser* self = (Spectral2DDenoiser*)instance;

  if (!self || !self->nlm_filter) {
    return 0;
  }

  return nlm_filter_get_latency_frames(self->nlm_filter);
}
