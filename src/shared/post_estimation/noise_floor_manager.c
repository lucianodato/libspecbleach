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

#include "noise_floor_manager.h"
#include "spectral_whitening.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseFloorManager {
  SpectralWhitening* whitening;
  float* whitening_weights;
  uint32_t real_spectrum_size;
};

NoiseFloorManager* noise_floor_manager_initialize(const uint32_t fft_size,
                                                  const uint32_t sample_rate,
                                                  const uint32_t hop) {
  NoiseFloorManager* self =
      (NoiseFloorManager*)calloc(1U, sizeof(NoiseFloorManager));
  if (!self) {
    return NULL;
  }

  self->real_spectrum_size = fft_size / 2U + 1U;

  self->whitening = spectral_whitening_initialize(fft_size);
  if (!self->whitening) {
    free(self);
    return NULL;
  }

  self->whitening_weights =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->whitening_weights) {
    spectral_whitening_free(self->whitening);
    free(self);
    return NULL;
  }

  return self;
}

void noise_floor_manager_free(NoiseFloorManager* self) {
  if (!self) {
    return;
  }
  if (self->whitening) {
    spectral_whitening_free(self->whitening);
  }
  if (self->whitening_weights) {
    free(self->whitening_weights);
  }
  free(self);
}

void noise_floor_manager_apply(NoiseFloorManager* self,
                               uint32_t real_spectrum_size, uint32_t fft_size,
                               float* gain_spectrum, const float* noise_profile,
                               float reduction_amount, float whitening_factor) {
  if (!self || !gain_spectrum || !noise_profile) {
    return;
  }

  // 1. Calculate whitening weights (including tapering)
  spectral_whitening_get_weights(self->whitening, whitening_factor,
                                 noise_profile, self->whitening_weights);

  // 2. Apply biasing + frequency-dependent floor
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float floor = reduction_amount * self->whitening_weights[k];
    if (floor > 1.0f)
      floor = 1.0f;

    float range = 1.0f - floor;
    gain_spectrum[k] = floor + (range * gain_spectrum[k]);
  }

  // 3. Symmetric copy
  for (uint32_t k = 1U; k < fft_size - k; k++) {
    gain_spectrum[fft_size - k] = gain_spectrum[k];
  }
}
