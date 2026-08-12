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

#include "spectral_smoother.h"
#include "shared/configurations.h"
#include "shared/utils/spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct SpectralSmoother {
  uint32_t fft_size;
  uint32_t sample_rate;
  uint32_t overlap_factor;
  uint32_t real_spectrum_size;
  TimeSmoothingType type;

  float* smoothed_spectrum_previous;
  bool initialized;
};

SpectralSmoother* spectral_smoothing_initialize(uint32_t fft_size,
                                                uint32_t sample_rate,
                                                uint32_t overlap_factor,
                                                TimeSmoothingType type) {
  if (fft_size == 0U) {
    return NULL;
  }

  SpectralSmoother* self =
      (SpectralSmoother*)calloc(1U, sizeof(SpectralSmoother));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->sample_rate = (sample_rate > 0U) ? sample_rate : 44100U;
  self->overlap_factor =
      (overlap_factor > 0U) ? overlap_factor : OVERLAP_FACTOR_1D;
  self->real_spectrum_size = (fft_size / 2U) + 1U;
  self->type = type;

  self->smoothed_spectrum_previous =
      (float*)calloc(self->real_spectrum_size, sizeof(float));

  if (!self->smoothed_spectrum_previous) {
    spectral_smoothing_free(self);
    return NULL;
  }

  self->initialized = false;
  return self;
}

void spectral_smoothing_free(SpectralSmoother* self) {
  if (!self) {
    return;
  }
  if (self->smoothed_spectrum_previous) {
    free(self->smoothed_spectrum_previous);
  }
  free(self);
}

bool spectral_smoothing_run(SpectralSmoother* self,
                            TimeSmoothingParameters parameters, float* gains) {
  if (!self || !gains) {
    return false;
  }

  float smoothing = parameters.smoothing;
  if (smoothing <= 0.0F) {
    memcpy(self->smoothed_spectrum_previous, gains,
           sizeof(float) * self->real_spectrum_size);
    self->initialized = true;
    return true;
  }

  if (!self->initialized) {
    memcpy(self->smoothed_spectrum_previous, gains,
           sizeof(float) * self->real_spectrum_size);
    self->initialized = true;
    return true;
  }

  // Map smoothing factor (0.0 to 1.0) to release time (10ms to 150ms)
  float p = fmaxf(0.0f, fminf(1.0f, smoothing));
  float tau_sec =
      GAIN_SMOOTHING_MIN_RELEASE_SEC +
      (p * (GAIN_SMOOTHING_MAX_RELEASE_SEC - GAIN_SMOOTHING_MIN_RELEASE_SEC));
  float dt = ((float)self->fft_size / (float)self->overlap_factor) /
             (float)self->sample_rate;
  float alpha = expf(-dt / tau_sec);

  uint32_t k = 0U;
  for (k = 0U; k < self->real_spectrum_size; k++) {
    gains[k] = (alpha * self->smoothed_spectrum_previous[k]) +
               ((1.0F - alpha) * gains[k]);
    self->smoothed_spectrum_previous[k] = gains[k];
  }

  return true;
}

void spectral_smoothing_apply_spatial(float* data, uint32_t size) {
  if (!data || size < 2U) {
    return;
  }

  float prev = data[0];
  uint32_t i = 0U;
  for (i = 1U; i < size; i++) {
    float curr = data[i];
    data[i] = 0.25F * prev + 0.5F * curr +
              0.25F * (i + 1 < size ? data[i + 1] : curr);
    prev = curr;
  }
}

void spectral_smoothing_apply_simple_temporal(float* current, float* memory,
                                              uint32_t size, float smoothing) {
  if (!current || !memory || size == 0U) {
    return;
  }

  uint32_t i = 0U;
  for (i = 0U; i < size; i++) {
    current[i] = smoothing * memory[i] + (1.0F - smoothing) * current[i];
    memory[i] = current[i];
  }
}
