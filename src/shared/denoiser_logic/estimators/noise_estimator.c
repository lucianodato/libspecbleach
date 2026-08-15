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

#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/utils/spectral_circular_buffer.h"
#include "shared/utils/spectral_smoother.h"
#include "shared/utils/spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseEstimator {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  SbSpectralCircularBuffer* median_buffer;
  uint32_t layer_median;

  float* welford_mean;
  float* welford_m2;
  uint32_t welford_count;

  NoiseProfile* noise_profile;
};

NoiseEstimator* noise_estimation_initialize(const uint32_t fft_size,
                                            NoiseProfile* noise_profile) {
  NoiseEstimator* self = (NoiseEstimator*)calloc(1U, sizeof(NoiseEstimator));

  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = (self->fft_size / 2U) + 1U;

  self->noise_profile = noise_profile;
  self->median_buffer =
      spectral_circular_buffer_create(NUMBER_OF_MEDIAN_SPECTRUM);

  if (!self->median_buffer) {
    noise_estimation_free(self);
    return NULL;
  }

  self->layer_median = spectral_circular_buffer_add_layer(
      self->median_buffer, self->real_spectrum_size);

  self->welford_mean = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->welford_m2 = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->welford_count = 0U;

  if (!self->welford_mean || !self->welford_m2) {
    noise_estimation_free(self);
    return NULL;
  }

  return self;
}

void noise_estimation_free(NoiseEstimator* self) {
  if (!self) {
    return;
  }

  // Don't free noise profile used as reference here

  if (self->median_buffer) {
    spectral_circular_buffer_free(self->median_buffer);
  }

  free(self->welford_mean);
  free(self->welford_m2);

  free(self);
}

void noise_estimation_reset(NoiseEstimator* self) {
  if (!self) {
    return;
  }
  self->welford_count = 0U;
  if (self->welford_mean) {
    memset(self->welford_mean, 0, self->real_spectrum_size * sizeof(float));
  }
  if (self->welford_m2) {
    memset(self->welford_m2, 0, self->real_spectrum_size * sizeof(float));
  }
  if (self->median_buffer) {
    spectral_circular_buffer_clear(self->median_buffer);
  }
}

static void update_rolling_mean(NoiseEstimator* self, float* noise_profile,
                                const float* signal_spectrum,
                                NoiseEstimatorType type) {
  get_rolling_mean_spectrum(
      noise_profile, signal_spectrum,
      get_noise_profile_block_count(self->noise_profile, type),
      self->real_spectrum_size);
  increment_block_count(self->noise_profile, type);
}

static void update_median(NoiseEstimator* self, float* noise_profile,
                          const float* signal_spectrum,
                          NoiseEstimatorType type) {
  spectral_circular_buffer_push(self->median_buffer, self->layer_median,
                                signal_spectrum);
  spectral_circular_buffer_advance(self->median_buffer);

  const uint32_t blocks = NUMBER_OF_MEDIAN_SPECTRUM;
  const float* history_frames[NUMBER_OF_MEDIAN_SPECTRUM];

  for (uint32_t i = 0; i < blocks; i++) {
    history_frames[i] = spectral_circular_buffer_retrieve(
        self->median_buffer, self->layer_median, i + 1);
  }

  if (get_rolling_median_spectrum(noise_profile, history_frames, blocks,
                                  self->real_spectrum_size)) {
    set_noise_profile_available(self->noise_profile, type);
  }
}

static void update_welford(NoiseEstimator* self, const float* signal_spectrum) {
  self->welford_count++;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    float delta = signal_spectrum[k] - self->welford_mean[k];
    self->welford_mean[k] += delta / (float)self->welford_count;
    float delta2 = signal_spectrum[k] - self->welford_mean[k];
    self->welford_m2[k] += delta * delta2;
  }
}

bool noise_estimation_run(NoiseEstimator* self,
                          const NoiseEstimatorType noise_estimator_type,
                          float* signal_spectrum) {
  if (!self || !signal_spectrum) {
    return false;
  }

  float* noise_profile =
      get_noise_profile(self->noise_profile, noise_estimator_type);

  switch (noise_estimator_type) {
    case ROLLING_MEAN:
      update_rolling_mean(self, noise_profile, signal_spectrum,
                          noise_estimator_type);
      break;
    case MEDIAN:
      update_median(self, noise_profile, signal_spectrum, noise_estimator_type);
      break;
    case STD_DEV:
      update_welford(self, signal_spectrum);
      break;
    case CV_MASK:
      break;
    default:
      break;
  }

  return true;
}

void noise_estimation_finalize(NoiseEstimator* self,
                               NoiseEstimatorType noise_estimator_type) {
  if (!self) {
    return;
  }

  if (noise_estimator_type == STD_DEV) {
    float* std_dev_profile = get_noise_profile(self->noise_profile, STD_DEV);
    for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
      float variance =
          self->welford_m2[k] / fmaxf(1.0f, (float)(self->welford_count - 1));
      std_dev_profile[k] = sqrtf(variance);
    }
    set_noise_profile_available(self->noise_profile, STD_DEV);
    return;
  }

  if (noise_estimator_type == CV_MASK) {
    float* cv_mask = get_noise_profile(self->noise_profile, CV_MASK);
    float* std_dev_profile = get_noise_profile(self->noise_profile, STD_DEV);
    for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
      float cv = std_dev_profile[k] / (self->welford_mean[k] + 1e-12f);
      // Map CV to a 0.0 - 1.0 mask. Rayleigh noise CV is ~0.523.
      float mask = 1.0f - (cv / 0.52f);
      cv_mask[k] = fmaxf(0.0f, fminf(1.0f, mask));
    }
    // Apply a light 3-bin spatial smooth so the mask catches the "skirts" of
    // the hum
    spectral_smoothing_apply_spatial(cv_mask, self->real_spectrum_size);
    set_noise_profile_available(self->noise_profile, CV_MASK);
    return;
  }

  float* noise_profile =
      get_noise_profile(self->noise_profile, noise_estimator_type);

  if (noise_profile && is_noise_estimation_available(self->noise_profile,
                                                     noise_estimator_type)) {
    // Basic refinement
    interpolate_spectrum_gaps(noise_profile, self->real_spectrum_size,
                              NOISE_ESTIMATION_INTERPOLATION_THRESHOLD);
    smooth_spectrum(noise_profile, self->real_spectrum_size,
                    NOISE_ESTIMATION_SMOOTHING_FACTOR);
  }
}
