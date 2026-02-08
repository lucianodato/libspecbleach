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

#include "masking_veto.h"
#include "../configurations.h"
#include "../pre_estimation/masking_estimator.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

struct MaskingVeto {
  uint32_t real_spectrum_size;
  MaskingEstimator* masking_estimator;
  float* clean_signal_estimation;
  float* masking_thresholds;
};

MaskingVeto* masking_veto_initialize(uint32_t fft_size, uint32_t sample_rate,
                                     SpectrumType spectrum_type) {
  MaskingVeto* self = (MaskingVeto*)calloc(1U, sizeof(MaskingVeto));
  if (!self) {
    return NULL;
  }

  self->real_spectrum_size = (fft_size / 2U) + 1U;
  self->masking_estimator =
      masking_estimation_initialize(fft_size, sample_rate, spectrum_type);

  if (!self->masking_estimator) {
    masking_veto_free(self);
    return NULL;
  }

  self->clean_signal_estimation =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->masking_thresholds =
      (float*)calloc(self->real_spectrum_size, sizeof(float));

  if (!self->clean_signal_estimation || !self->masking_thresholds) {
    masking_veto_free(self);
    return NULL;
  }

  // Relative-only Masking: Disable the absolute threshold floor.
  // This ensures that the Veto only triggers when actual SIGNAL (music)
  // is providing masking, and doesn't shut down in silence or pure noise.
  masking_estimation_set_use_absolute_threshold(self->masking_estimator, false);

  return self;
}

void masking_veto_free(MaskingVeto* self) {
  if (!self) {
    return;
  }

  masking_estimation_free(self->masking_estimator);
  free(self->clean_signal_estimation);
  free(self->masking_thresholds);
  free(self);
}

void masking_veto_apply(MaskingVeto* self, const float* smoothed_spectrum,
                        const float* noisy_spectrum,
                        const float* noise_spectrum, float* alpha,
                        float floor_alpha, float depth, float elasticity) {
  if (!self || !smoothed_spectrum || !noisy_spectrum || !noise_spectrum ||
      !alpha || depth < 0.0F) {
    return;
  }

  // 1. Estimate clean signal magnitude from SMOOTHED signal for masking
  // thresholds. We use the smoothed signal as the "Masker" because it defines
  // the steady psychoacoustic floor.
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    self->clean_signal_estimation[k] =
        fmaxf(smoothed_spectrum[k] - noise_spectrum[k], 0.0F);
  }

  if (!compute_masking_thresholds(self->masking_estimator,
                                  self->clean_signal_estimation,
                                  self->masking_thresholds)) {
    return;
  }

  /**
   * 2. Compare NOISY signal to masking thresholds to determine audibility.
   * If noise is louder than thresholds, it's AUDIBLE (alpha stays high).
   * If noise is quieter than thresholds, it's MASKED (alpha drops to floor).
   */
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float nmr_db =
        10.0F * log10f(noise_spectrum[k] /
                       (self->masking_thresholds[k] + SPECTRAL_EPSILON));

    float audibility = 1.0F;
    if (nmr_db <= 0.0F) {
      audibility = 0.0F;
    } else if (nmr_db >= 20.0F) {
      audibility = 1.0F;
    } else {
      audibility = nmr_db / 20.0F;
    }

    // 3. Ratio-based Conjunction Logic with Elasticity:
    // Smoothing Ratio: how much energy NLM preserved relative to input.
    const float smoothing_ratio = fminf(
        smoothed_spectrum[k] / (noisy_spectrum[k] + SPECTRAL_EPSILON), 1.0F);

    const float initial_alpha = alpha[k];

    // Protection is moderated by elasticity.
    // At elasticity 0.0, we have pure psychoacoustic protection.
    // At elasticity 1.0, we have no protection (oversubtraction leaks 100%).
    const float protection_amount = (1.0F - audibility) * depth;
    const float protection_effective = protection_amount * (1.0F - elasticity);

    // Target moderated alpha based on psychoacoustic masking (Legacy floor)
    const float moderated_alpha = floor_alpha + ((initial_alpha - floor_alpha) *
                                                 (1.0F - protection_effective));

    // Final Alpha: Scale the moderated alpha by the smoothing ratio ONLY
    // in masked regions to rescue signal energy (transients and masked detail)
    // blurred by the smoothing process. Moderated by elasticity.
    alpha[k] = moderated_alpha *
               (1.0F - (protection_effective * (1.0F - smoothing_ratio)));
  }
}
