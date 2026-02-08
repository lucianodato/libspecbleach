/*
libspecbleach - A spectral processing library

Copyright 2022-2026 Luciano Dato <lucianodato@gmail.com>

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

#include "suppression_engine.h"
#include "../configurations.h"
#include <math.h>
#include <stdlib.h>

struct SuppressionEngine {
  uint32_t real_spectrum_size;
};

SuppressionEngine* suppression_engine_initialize(uint32_t real_spectrum_size) {
  SuppressionEngine* self =
      (SuppressionEngine*)calloc(1U, sizeof(SuppressionEngine));
  if (!self) {
    return NULL;
  }

  self->real_spectrum_size = real_spectrum_size;
  return self;
}

void suppression_engine_free(SuppressionEngine* self) {
  if (self) {
    free(self);
  }
}

void suppression_engine_calculate(SuppressionEngine* self,
                                  const float* reference_spectrum,
                                  const float* noise_spectrum,
                                  float suppression_strength, float* alpha,
                                  float* beta) {
  if (!self || !reference_spectrum || !noise_spectrum || !alpha || !beta) {
    return;
  }

  /*
   * Berouti-style SNR-dependent oversubtraction.
   * Alpha (Oversubtraction) decreases as SNR increases.
   * suppression_strength scales the maximum Alpha.
   */

  const float alpha_max_user =
      ALPHA_MIN + (suppression_strength * (ALPHA_MAX - ALPHA_MIN));
  const float snr_range_db =
      SUPPRESSION_HIGHER_SNR_DB - SUPPRESSION_LOWER_SNR_DB;

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float snr_db = 10.0F * log10f(reference_spectrum[k] /
                                        (noise_spectrum[k] + SPECTRAL_EPSILON));

    if (snr_db <= SUPPRESSION_LOWER_SNR_DB) {
      alpha[k] = alpha_max_user;
    } else if (snr_db >= SUPPRESSION_HIGHER_SNR_DB) {
      alpha[k] = ALPHA_MIN;
    } else {
      // Linear interpolation between alpha_max_user and ALPHA_MIN
      const float normalized_snr =
          (snr_db - SUPPRESSION_LOWER_SNR_DB) / snr_range_db;
      alpha[k] =
          alpha_max_user - (normalized_snr * (alpha_max_user - ALPHA_MIN));
    }

    /* Beta (Undersubtraction) is kept minimal as we use a dedicated noise floor
     * manager */
    beta[k] = 0.0F;
  }
}
