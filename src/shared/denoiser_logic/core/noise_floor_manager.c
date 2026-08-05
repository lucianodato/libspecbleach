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

#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/processing/spectral_whitening.h"
#include "shared/utils/spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseFloorManager {
  SpectralWhitening* whitening;
  float* whitening_weights;
  uint32_t real_spectrum_size;
  uint32_t fft_size;
};

NoiseFloorManager* noise_floor_manager_initialize(const uint32_t fft_size) {
  if (fft_size == 0U) {
    return NULL;
  }
  NoiseFloorManager* self =
      (NoiseFloorManager*)calloc(1U, sizeof(NoiseFloorManager));
  if (!self) {
    return NULL;
  }

  self->real_spectrum_size = (fft_size / 2U) + 1U;
  self->fft_size = fft_size;

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
                               float reduction_amount,
                               float tonal_reduction_amount,
                               const float* tonal_mask,
                               float whitening_factor) {
  if (!self || !gain_spectrum || !noise_profile) {
    return;
  }

  if (real_spectrum_size != self->real_spectrum_size ||
      fft_size != self->fft_size) {
    real_spectrum_size = self->real_spectrum_size;
    fft_size = self->fft_size;
  }

  if (reduction_amount >= 0.999f && tonal_reduction_amount >= 0.999f) {
    // Transparency Guard: If both reduction paths are at 0dB (1.0f),
    // we force unity gain to ensure BIT TRANSPARENCY and skip processing.
    for (uint32_t k = 0U; k < real_spectrum_size; k++) {
      gain_spectrum[k] = 1.0f;
    }
    sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
    return;
  }

  // 1. Calculate ideal whitening reduction in dB
  spectral_whitening_get_ideal_reduction_db(self->whitening, reduction_amount,
                                            noise_profile,
                                            self->whitening_weights);

  // 2. Apply biasing + frequency-dependent floor
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float mask = (tonal_mask) ? tonal_mask[k] : 0.0f;
    if (mask > 0.0f) {
      mask = sqrtf(sqrtf(mask));
    }

    // Calculate the available reduction budget (in linear amplitude)
    float dual_path_reduction =
        (reduction_amount * (1.0f - mask)) + (tonal_reduction_amount * mask);

    // Convert budget to dB
    float r_dp_db = -20.0f * log10f(dual_path_reduction + 1e-12f);
    if (r_dp_db < 0.0f) {
      r_dp_db = 0.0f;
    }

    // Ideal reduction needed to perfectly flatten the profile
    float ideal_db = self->whitening_weights[k];

    // Cap the ideal reduction by the available budget to prevent digging holes
    float flattened_db = ideal_db;
    if (flattened_db > r_dp_db) {
      flattened_db = r_dp_db;
    }

    // Interpolate between deep notches (0% whitening) and perfectly flat (100%
    // whitening)
    float target_db = ((1.0f - whitening_factor) * r_dp_db) +
                      (whitening_factor * flattened_db);

    // Convert back to linear floor
    float whitened_floor = powf(10.0f, -target_db / 20.0f);

    // Safety clamp (Floor cannot exceed Unity Gain)
    if (whitened_floor > 1.0f) {
      whitened_floor = 1.0f;
    }

    /*
     * CLAMPING FLOOR:
     * We use fmaxf to ensure the gain never drops below the whitened floor.
     * This "replaces" the suppressed noise with the floor level,
     * effectively shaving off the spectral shape of the noise residual.
     */
    gain_spectrum[k] = fmaxf(whitened_floor, gain_spectrum[k]);
  }

  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}
