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

#include "shared/denoiser_logic/processing/spectral_whitening.h"
#include "shared/configurations.h"
#include "shared/utils/spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct SpectralWhitening {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  float* smoothed_profile;
};

SpectralWhitening* spectral_whitening_initialize(const uint32_t fft_size) {
  if (fft_size == 0U) {
    return NULL;
  }
  SpectralWhitening* self =
      (SpectralWhitening*)calloc(1U, sizeof(SpectralWhitening));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = (self->fft_size / 2U) + 1U;

  self->smoothed_profile =
      (float*)calloc(self->real_spectrum_size, sizeof(float));

  if (!self->smoothed_profile) {
    spectral_whitening_free(self);
    return NULL;
  }

  return self;
}

void spectral_whitening_free(SpectralWhitening* self) {
  if (!self) {
    return;
  }
  if (self->smoothed_profile) {
    free(self->smoothed_profile);
  }
  free(self);
}

void spectral_whitening_get_ideal_reduction_db(SpectralWhitening* self,
                                               float reduction_amount,
                                               const float* noise_profile,
                                               float* ideal_reduction_db_out) {
  if (!self || !ideal_reduction_db_out || !noise_profile) {
    return;
  }

  // Calculate broadband reduction in dB
  float r_db = -20.0f * log10f(reduction_amount + 1e-12f);
  if (r_db < 0.0f) {
    r_db = 0.0f;
  }

  // 1. Create a smoothed copy of the noise profile in pre-allocated buffer
  memcpy(self->smoothed_profile, noise_profile,
         self->real_spectrum_size * sizeof(float));
  smooth_spectrum(self->smoothed_profile, self->real_spectrum_size, 0.5f);

  // 2. Compute Peak and Reference in the Logarithmic (dB) Domain
  float pmax_db = -200.0f;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    float db = 10.0f * log10f(self->smoothed_profile[k] + SPECTRAL_EPSILON);
    if (db > pmax_db) {
      pmax_db = db;
    }
  }

  // To find a stable visual reference line, take the mean dB of the valid
  // spectrum.
  float sum_db = 0.0f;
  uint32_t count = 0U;
  float threshold_db = pmax_db - 60.0f;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    float db = 10.0f * log10f(self->smoothed_profile[k] + SPECTRAL_EPSILON);
    if (db > threshold_db) {
      sum_db += db;
      count++;
    }
  }
  float pref_db = (count > 0U) ? (sum_db / (float)count) : -120.0f;

  float delta_max = fmaxf(pmax_db - pref_db, 0.0f);

  float excess = 0.0f;
  if (r_db > delta_max && r_db > 1.0f) {
    excess = r_db - delta_max;
  }

  // 3. Compute Ideal Reduction in dB for 100% Whitening
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    float db = 10.0f * log10f(self->smoothed_profile[k] + SPECTRAL_EPSILON);
    float delta = (db > pref_db) ? (db - pref_db) : 0.0f;

    // The ideal reduction perfectly flattens the noise profile and shifts it
    // down by any excess
    ideal_reduction_db_out[k] = delta + excess;
  }
}
