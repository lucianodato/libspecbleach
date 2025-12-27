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

#include "spectral_whitening.h"
#include "../configurations.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct SpectralWhitening {
  float* tapering_window;
  uint32_t fft_size;
  uint32_t real_spectrum_size;
};

SpectralWhitening* spectral_whitening_initialize(const uint32_t fft_size) {
  SpectralWhitening* self =
      (SpectralWhitening*)calloc(1U, sizeof(SpectralWhitening));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;

  self->tapering_window =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->tapering_window) {
    free(self);
    return NULL;
  }

  // Pre-calculate Right half of Hamming window for HF tapering
  uint32_t N = self->real_spectrum_size * 2 - 1;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    uint32_t n = k + self->real_spectrum_size - 1;
    self->tapering_window[k] =
        0.54f - 0.46f * cosf(2.0f * (float)M_PI * (float)n / (float)(N - 1));
  }

  return self;
}

void spectral_whitening_free(SpectralWhitening* self) {
  if (!self)
    return;
  if (self->tapering_window)
    free(self->tapering_window);
  free(self);
}

void spectral_whitening_get_weights(SpectralWhitening* self,
                                    float whitening_factor,
                                    const float* noise_profile,
                                    float* weights_out) {
  if (!self || !weights_out || !noise_profile)
    return;

  float noise_peak = 1e-12f;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    if (noise_profile[k] > noise_peak) {
      noise_peak = noise_profile[k];
    }
  }

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    float weight = 1.0f;
    if (whitening_factor > 0.0f && noise_profile[k] > 1e-12f) {
      // Power-law valley filling
      weight = powf(noise_peak / noise_profile[k], whitening_factor);
    }
    // Weights include tapering
    weights_out[k] = weight * self->tapering_window[k];
  }
}
