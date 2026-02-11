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

#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/configurations.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/masking_estimator.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct MaskingVeto {
  uint32_t real_spectrum_size;
  MaskingEstimator* masking_estimator;
  CriticalBands* critical_bands_helper;
  float* clean_signal_estimation;
  float* masking_thresholds;
  float* band_audibility;
  float* band_centers; // Store center bin for each band
  uint32_t sample_rate;
};

MaskingVeto* masking_veto_initialize(uint32_t fft_size, uint32_t sample_rate,
                                     CriticalBandType critical_band_type,
                                     SpectrumType spectrum_type) {
  MaskingVeto* self = (MaskingVeto*)calloc(1U, sizeof(MaskingVeto));
  if (!self) {
    return NULL;
  }

  self->real_spectrum_size = (fft_size / 2U) + 1U;
  self->masking_estimator = masking_estimation_initialize(
      fft_size, sample_rate, critical_band_type, spectrum_type);

  // Helper to access band indexes for classification
  self->critical_bands_helper =
      critical_bands_initialize(sample_rate, fft_size, critical_band_type);

  if (!self->masking_estimator || !self->critical_bands_helper) {
    masking_veto_free(self);
    return NULL;
  }

  const uint32_t num_bands =
      get_number_of_critical_bands(self->critical_bands_helper);

  self->clean_signal_estimation =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->masking_thresholds =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->band_audibility = (float*)calloc(num_bands, sizeof(float));
  self->band_centers = (float*)calloc(num_bands, sizeof(float));
  self->sample_rate = sample_rate;

  if (!self->clean_signal_estimation || !self->masking_thresholds ||
      !self->band_audibility || !self->band_centers) {
    masking_veto_free(self);
    return NULL;
  }

  // Pre-calculate band center bins for interpolation
  for (uint32_t j = 0; j < num_bands; j++) {
    const CriticalBandIndexes indexes =
        get_band_indexes(self->critical_bands_helper, j);
    self->band_centers[j] =
        (float)indexes.start_position +
        ((float)(indexes.end_position - indexes.start_position) / 2.0F);
  }

  // Relative-only Masking for the Veto: Disable the absolute threshold floor.
  // The Veto should only trigger when SIGNAL is providing masking.
  masking_estimation_set_use_absolute_threshold(self->masking_estimator, false);

  return self;
}

void masking_veto_free(MaskingVeto* self) {
  if (!self) {
    return;
  }

  masking_estimation_free(self->masking_estimator);
  critical_bands_free(self->critical_bands_helper);
  free(self->clean_signal_estimation);
  free(self->masking_thresholds);
  free(self->band_audibility);
  free(self->band_centers);
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

  // 1. Estimate clean signal magnitude from SMOOTHED signal
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    self->clean_signal_estimation[k] =
        fmaxf(smoothed_spectrum[k] - noise_spectrum[k], 0.0F);
  }

  // 2. Compute psychoacoustic masking thresholds
  if (!compute_masking_thresholds(self->masking_estimator,
                                  self->clean_signal_estimation,
                                  self->masking_thresholds)) {
    return;
  }

  /**
   * 3. Calculate Band-wise Audibility
   * We calculate the Noise-to-Mask Ratio (NMR) for each critical band.
   */
  const uint32_t num_bands =
      get_number_of_critical_bands(self->critical_bands_helper);

  for (uint32_t j = 0U; j < num_bands; j++) {
    const CriticalBandIndexes indexes =
        get_band_indexes(self->critical_bands_helper, j);

    float band_noise_energy = 0.0F;
    float band_threshold = 0.0F;

    for (uint32_t k = indexes.start_position; k < indexes.end_position; k++) {
      band_noise_energy += noise_spectrum[k];
      band_threshold += self->masking_thresholds[k];
    }

    const uint32_t bins_in_band = indexes.end_position - indexes.start_position;
    if (bins_in_band > 0) {
      const float nmr_db =
          10.0F * log10f((band_noise_energy + SPECTRAL_EPSILON) /
                         (band_threshold + SPECTRAL_EPSILON));

      // Map NMR (Noise-to-Mask Ratio) to 0.0-1.0 audibility range
      if (nmr_db <= 0.0F) {
        self->band_audibility[j] = 0.0F;
      } else if (nmr_db >= 20.0F) {
        self->band_audibility[j] = 1.0F;
      } else {
        self->band_audibility[j] = nmr_db / 20.0F;
      }
    } else {
      self->band_audibility[j] = 0.0F;
    }
  }

  /**
   * 4. Apply Interpolated Veto
   * Instead of applying in blocks, we interpolate audibility for each bin
   * between the centers of the critical bands. This structurally eliminates
   * sharp gain transitions and Bark band boundaries.
   */
  uint32_t current_band = 0;
  for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
    // Find the two bands to interpolate between for bin k
    while (current_band < num_bands - 1 &&
           (float)k > self->band_centers[current_band + 1]) {
      current_band++;
    }

    float bin_audibility;
    if (current_band < num_bands - 1) {
      // Linear interpolation between band centers
      const float x0 = self->band_centers[current_band];
      const float x1 = self->band_centers[current_band + 1];
      const float y0 = self->band_audibility[current_band];
      const float y1 = self->band_audibility[current_band + 1];

      if (x1 > x0) {
        const float t = ((float)k - x0) / (x1 - x0);
        bin_audibility = (y0 * (1.0F - fmaxf(0.0F, fminf(1.0F, t)))) +
                         (y1 * fmaxf(0.0F, fminf(1.0F, t)));
      } else {
        bin_audibility = y0;
      }
    } else {
      bin_audibility = self->band_audibility[num_bands - 1];
    }

    const float initial_alpha = alpha[k];

    // Protection is moderated by elasticity
    const float protection_amount = (1.0F - bin_audibility) * depth;
    const float protection_effective = protection_amount * (1.0F - elasticity);

    // Apply the Veto transformation to the gain factor (alpha).
    // In fully masked areas, alpha is reduced towards floor_alpha.
    alpha[k] = floor_alpha +
               ((initial_alpha - floor_alpha) * (1.0F - protection_effective));
  }
}
