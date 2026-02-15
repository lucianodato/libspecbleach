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
#include "shared/utils/transient_detector.h"
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
  float* stable_clean_signal;
  float* masking_thresholds;
  float* band_audibility;
  float* band_centers; // Store center bin for each band
  float* onset_weights;
  TransientDetector* transient_detector;
  uint32_t sample_rate;
  float* future_clean_estimation_buf;
  float* band_energies_buf;
};

MaskingVeto* masking_veto_initialize(uint32_t fft_size, uint32_t sample_rate,
                                     CriticalBandType critical_band_type,
                                     SpectrumType spectrum_type,
                                     bool use_absolute_threshold,
                                     bool use_temporal_masking) {

  MaskingVeto* self = (MaskingVeto*)calloc(1U, sizeof(MaskingVeto));
  if (!self) {
    return NULL;
  }

  self->real_spectrum_size = (fft_size / 2U) + 1U;
  self->masking_estimator = masking_estimation_initialize(
      fft_size, sample_rate, critical_band_type, spectrum_type,
      use_absolute_threshold, use_temporal_masking);

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
  self->stable_clean_signal =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->masking_thresholds =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->band_audibility = (float*)calloc(num_bands, sizeof(float));
  self->band_centers = (float*)calloc(num_bands, sizeof(float));
  self->onset_weights = (float*)calloc(num_bands, sizeof(float));
  self->transient_detector = transient_detector_initialize(num_bands);
  self->future_clean_estimation_buf =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->band_energies_buf = (float*)calloc(num_bands, sizeof(float));
  self->sample_rate = sample_rate;

  if (!self->clean_signal_estimation || !self->stable_clean_signal ||
      !self->masking_thresholds || !self->band_audibility ||
      !self->band_centers || !self->onset_weights ||
      !self->transient_detector || !self->future_clean_estimation_buf ||
      !self->band_energies_buf) {
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

  return self;
}

void masking_veto_free(MaskingVeto* self) {
  if (!self) {
    return;
  }

  masking_estimation_free(self->masking_estimator);
  critical_bands_free(self->critical_bands_helper);
  free(self->clean_signal_estimation);
  free(self->stable_clean_signal);
  free(self->masking_thresholds);
  free(self->band_audibility);
  free(self->band_centers);
  free(self->onset_weights);
  transient_detector_free(self->transient_detector);
  free(self->future_clean_estimation_buf);
  free(self->band_energies_buf);
  free(self);
}

void masking_veto_apply(MaskingVeto* self, const float* smoothed_spectrum,
                        const float* noisy_spectrum,
                        const float* noise_spectrum,
                        const float* future_spectrum, float* alpha,
                        float floor_alpha, float depth, float elasticity) {
  if (!self || !smoothed_spectrum || !noisy_spectrum || !noise_spectrum ||
      !alpha || depth < 0.0F || self->real_spectrum_size == 0U) {
    return;
  }

  // 1. Estimate clean signal magnitude from SMOOTHED signal
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float current_clean =
        fmaxf(smoothed_spectrum[k] - noise_spectrum[k], 0.0F);

    // Stabilize the clean signal estimation for the denoiser (prevents gargle
    // in the spreading skirts)
    self->stable_clean_signal[k] =
        (MASKING_VETO_SMOOTHING * current_clean) +
        ((1.0F - MASKING_VETO_SMOOTHING) * self->stable_clean_signal[k]);
    self->clean_signal_estimation[k] = self->stable_clean_signal[k];
  }

  // 1.1 Estimate future clean signal magnitude for backward masking
  float* future_clean_estimation = NULL;
  if (future_spectrum) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
      self->future_clean_estimation_buf[k] =
          fmaxf(future_spectrum[k] - noise_spectrum[k], 0.0F);
    }
    future_clean_estimation = self->future_clean_estimation_buf;
  }

  // 2. Compute psychoacoustic masking thresholds
  if (!compute_masking_thresholds(
          self->masking_estimator, self->clean_signal_estimation,
          future_clean_estimation, self->masking_thresholds)) {
    return;
  }

  /**
   * 3. Calculate Band-wise Audibility
   * We calculate the Noise-to-Mask Ratio (NMR) for each critical band.
   */
  const uint32_t num_bands =
      get_number_of_critical_bands(self->critical_bands_helper);

  // 2.1 Calculate band energies for the transient detector
  for (uint32_t j = 0U; j < num_bands; j++) {
    const CriticalBandIndexes indexes =
        get_band_indexes(self->critical_bands_helper, j);
    float energy = 0.0F;
    for (uint32_t k = indexes.start_position; k < indexes.end_position; k++) {
      energy += self->clean_signal_estimation[k];
    }
    self->band_energies_buf[j] = energy;
  }

  // 2.2 Update onset weights using the innovation ratio
  (void)transient_detector_process(
      self->transient_detector, self->band_energies_buf, self->onset_weights);

  for (uint32_t j = 0U; j < num_bands; j++) {
    const CriticalBandIndexes indexes =
        get_band_indexes(self->critical_bands_helper, j);

    const uint32_t bins_in_band = indexes.end_position - indexes.start_position;
    if (bins_in_band > 0) {
      float band_signal_energy = self->band_energies_buf[j];
      float band_noise_energy = 0.0F;
      float band_threshold = 0.0F;

      for (uint32_t k = indexes.start_position; k < indexes.end_position; k++) {
        band_noise_energy += noise_spectrum[k];
        band_threshold += self->masking_thresholds[k];
      }

      // 5. SNR-based selectivity
      const float snr_db =
          10.0F * log10f((band_signal_energy + SPECTRAL_EPSILON) /
                         (band_noise_energy + SPECTRAL_EPSILON));
      // Flat SNR threshold for simple selectivity
      const float snr_threshold = MASKING_VETO_SNR_THRESHOLD;
      const float snr_weight = fminf(
          fmaxf((snr_db - snr_threshold) / MASKING_VETO_SNR_RANGE, 0.0F), 1.0F);

      // Final selectivity is strictly SNR + Transient protection
      const float selectivity_weight =
          fmaxf(snr_weight, self->onset_weights[j]);

      const float nmr_db =
          10.0F * log10f((band_noise_energy + SPECTRAL_EPSILON) /
                         (band_threshold + SPECTRAL_EPSILON));

      // Map NMR (Noise-to-Mask Ratio) to 0.0-1.0 audibility range
      float audibility;
      if (nmr_db <= 0.0F) {
        audibility = 0.0F;
      } else if (nmr_db >= MASKING_VETO_NMR_RANGE) {
        audibility = 1.0F;
      } else {
        audibility = nmr_db / MASKING_VETO_NMR_RANGE;
      }

      // If selectivity is low (selectivity_weight -> 0), audibility -> 1.0
      // (bin is treated as noise/unprotected).
      self->band_audibility[j] = fmaxf(audibility, 1.0F - selectivity_weight);
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
