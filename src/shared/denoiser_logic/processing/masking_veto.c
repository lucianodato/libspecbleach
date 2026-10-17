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
#include "shared/frame_rate_norm.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/masking_estimator.h"
#include "shared/utils/spectral_smoother.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

struct MaskingVeto {
  uint32_t real_spectrum_size;
  SpectrumType spectrum_type; // Gain floor exponent (power vs magnitude)
  MaskingEstimator* masking_estimator;
  CriticalBands* critical_bands_helper;
  float smooth;
  float floor_smooth_memory; // Mapper regularization memory (tau at hop)
  const float* temporal_persistence; // Model's per-band temporal claim
  float* clean_signal_estimation;
  float* stable_clean_signal;
  float* masking_thresholds;
  float* band_floor_state;          // Regularized floor state per band
  float* band_floor_smooth;         // Per-band floor (spatially smoothed plate)
  float* band_masker_mask;          // Per-band validated-masker gate (0/1)
  float* band_centers;              // Store center bin for each band
  const float* absolute_thresholds; // Per-bin AT (owned by the estimator)
  uint32_t sample_rate;
  float* future_clean_estimation_buf;
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
  self->spectrum_type = spectrum_type;
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
  self->band_floor_state = (float*)calloc(num_bands, sizeof(float));
  self->band_floor_smooth = (float*)calloc(num_bands, sizeof(float));
  self->band_masker_mask = (float*)calloc(num_bands, sizeof(float));
  self->band_centers = (float*)calloc(num_bands, sizeof(float));
  self->future_clean_estimation_buf =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->sample_rate = sample_rate;
  self->absolute_thresholds =
      masking_estimation_get_absolute_thresholds(self->masking_estimator);
  self->temporal_persistence =
      masking_estimation_get_temporal_persistence(self->masking_estimator);
  self->smooth = MASKING_VETO_SMOOTHING;
  self->floor_smooth_memory =
      sb_alpha_from_tau(MASKING_VETO_FLOOR_SMOOTH_SEC, LEGACY_REF_HOP_SEC);

  if (!self->clean_signal_estimation || !self->stable_clean_signal ||
      !self->masking_thresholds || !self->band_floor_state ||
      !self->band_floor_smooth || !self->band_masker_mask ||
      !self->band_centers || !self->future_clean_estimation_buf ||
      !self->temporal_persistence) {
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
  free(self->band_floor_state);
  free(self->band_floor_smooth);
  free(self->band_masker_mask);
  free(self->band_centers);
  free(self->future_clean_estimation_buf);
  free(self);
}

static const float* mv_update_clean_est(MaskingVeto* self,
                                        const float* smoothed_spectrum,
                                        const float* noise_spectrum,
                                        const float* future_spectrum) {
  // 1. Estimate clean signal magnitude from SMOOTHED signal
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float current_clean = fmaxf(
        smoothed_spectrum[k] - (MASKING_VETO_NOISE_GATE * noise_spectrum[k]),
        0.0F);

    // Stabilize the clean signal estimation (prevents gargle in spreading
    // skirts). self->smooth is the memory (recurrent) weight, retuned per
    // hop to preserve the time constant -- same convention as the temporal
    // smoother below.
    self->stable_clean_signal[k] =
        ((1.0F - self->smooth) * current_clean) +
        (self->smooth * self->stable_clean_signal[k]);
    self->clean_signal_estimation[k] = self->stable_clean_signal[k];
  }

  // 1.1 Estimate future clean signal magnitude for backward masking
  if (!future_spectrum) {
    return NULL;
  }
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    self->future_clean_estimation_buf[k] = fmaxf(
        future_spectrum[k] - (MASKING_VETO_NOISE_GATE * noise_spectrum[k]),
        0.0F);
  }
  return self->future_clean_estimation_buf;
}

/**
 * Audibility closure per band: the gain at which the suppressed residual
 * lands exactly on the masking threshold. Power spectra carry gain
 * quadratically (residual power = gain^2 * noise power), so the floor is
 * sqrt(threshold/noise); magnitude spectra are linear. Noise already below
 * the threshold floors at 1.0 (no processing in that band, iZotope
 * "Masking" semantics); very audible noise floors near 0.0 (no relief).
 * @param relief_ceiling Maximum lift for this claim; reduction must stay
 *   the owner of the output even where the model says inaudible (see
 *   MASKING_VETO_RELIEF_CEILING vs the persistence-only ceiling).
 */
static float mv_floor_target(float band_noise_energy, float band_threshold,
                             SpectrumType spectrum_type, float relief_ceiling) {
  const float ratio = (band_threshold + SPECTRAL_EPSILON) /
                      (band_noise_energy + SPECTRAL_EPSILON);
  const float floor = (spectrum_type == POWER_SPECTRUM) ? sqrtf(ratio) : ratio;
  // Relief fades to zero as the noise becomes clearly audible above the
  // mask (legacy protection ramp over MASKING_VETO_RELIEF_RANGE_DB).
  const float nmr_db = -10.0F * log10f(ratio + SPECTRAL_EPSILON);
  const float protection = 1.0F - (nmr_db / MASKING_VETO_RELIEF_RANGE_DB);
  const float relief =
      fminf(floor, relief_ceiling) * fminf(fmaxf(protection, 0.0F), 1.0F);
  return fminf(fmaxf(relief, 0.0F), 1.0F);
}

static bool mv_compute_band_masker_validation(MaskingVeto* self,
                                              const float* noise_spectrum,
                                              const float* future_spectrum) {
  const uint32_t num_bands =
      get_number_of_critical_bands(self->critical_bands_helper);
  if (num_bands == 0U) {
    return false;
  }

  for (uint32_t j = 0U; j < num_bands; j++) {
    const CriticalBandIndexes indexes =
        get_band_indexes(self->critical_bands_helper, j);
    float band_noise_energy = 0.0F;
    float band_clean_energy = 0.0F;
    float band_future_energy = 0.0F;

    for (uint32_t k = indexes.start_position; k < indexes.end_position; k++) {
      band_noise_energy += noise_spectrum[k];
      band_clean_energy += self->clean_signal_estimation[k];
      if (future_spectrum) {
        band_future_energy += future_spectrum[k];
      }
    }

    self->band_masker_mask[j] =
        (band_clean_energy >
             MASKING_VETO_MASKER_MIN_RATIO * band_noise_energy ||
         band_future_energy > MASKING_VETO_MASKER_MIN_RATIO * band_noise_energy)
            ? 1.0F
            : 0.0F;
  }
  return true;
}

static void mv_compute_band_floor(MaskingVeto* self,
                                  const float* noise_spectrum,
                                  const float* future_clean_estimation) {
  /**
   * 3. Band-wise audibility floor: report the masking threshold as the
   * model computes it — frequency masking per frame in the estimator,
   * temporal persistence (forward masking decays, backward masking off the
   * lookahead frame) inside the model's temporal layer. No attack/release
   * envelope of its own: this layer owns audibility (what relief is valid),
   * not time. The hard instant-reset of the alpha veto that chopped
   * decaying tails stays gone because the threshold itself decays with the
   * note. Relief requires a masker: the current signal (note sustain), an
   * incoming one (the future frame pre-masks the current frame, backward
   * masking), or the model's persisted temporal claim (forward masking
   * still in force after a recent masker). In bands with none of these,
   * there is nothing to mask, so suppression runs at full depth.
   */
  const uint32_t num_bands =
      get_number_of_critical_bands(self->critical_bands_helper);

  for (uint32_t j = 0U; j < num_bands; j++) {
    const CriticalBandIndexes indexes =
        get_band_indexes(self->critical_bands_helper, j);

    float band_noise_energy = 0.0F;
    float band_threshold = 0.0F;
    float band_clean_energy = 0.0F;
    float band_future_energy = 0.0F;
    float band_at_energy = 0.0F;

    for (uint32_t k = indexes.start_position; k < indexes.end_position; k++) {
      band_noise_energy += noise_spectrum[k];
      band_threshold += self->masking_thresholds[k];
      band_clean_energy += self->clean_signal_estimation[k];
      if (future_clean_estimation) {
        band_future_energy += future_clean_estimation[k];
      }
      if (self->absolute_thresholds) {
        band_at_energy += self->absolute_thresholds[k];
      }
    }

    // The model emits band-level masking power, but it is stored per bin;
    // averaging it back keeps the threshold on the same scale as the summed
    // band noise (a per-bin comparison inflates the ratio by bins-in-band
    // and made noise-only sections read as "masked").
    band_threshold /=
        fmaxf((float)(indexes.end_position - indexes.start_position), 1.0F);

    // Relief requires either a real masker or inaudible noise that needs
    // no digging: noise below the absolute-hearing threshold is already
    // inaudible, so suppression runs at full depth instead of holding it.
    const bool audacious_noise =
        self->absolute_thresholds == NULL ||
        (band_at_energy > 0.0F && band_noise_energy > band_at_energy);

    // A real masker: current signal, or an incoming one (backward masking).
    // The model's own temporal claim counts as well: when its temporal
    // layer says a recent masker's forward masking is still in force in
    // this band, the residual there is still borderline-inaudible and
    // relief must not reset to full depth at the instant the masker
    // estimate drops (this is what preserved floors through word gaps).
    // No veto-side envelope: the claim decays exactly per the model's
    // forward/backward masking law. Reduction precedence still holds: a
    // claim with no live masker behind it only preserves low-level detail,
    // so it is capped at the stricter persistence ceiling.
    const bool live_masker =
        band_clean_energy > MASKING_VETO_MASKER_MIN_RATIO * band_noise_energy ||
        band_future_energy > MASKING_VETO_MASKER_MIN_RATIO * band_noise_energy;
    const bool masker_present =
        live_masker || (self->temporal_persistence &&
                        self->temporal_persistence[j] >
                            MASKING_VETO_PERSISTENCE_RATIO * band_noise_energy);

    const bool valid = masker_present && audacious_noise;
    // The two relief ceilings are the evidence-class hierarchy and may
    // diverge if tuning separates them again; the split stays even while
    // both map to 0.30.
    const float relief_ceiling = live_masker
                                     // NOLINTNEXTLINE(bugprone-branch-clone)
                                     ? MASKING_VETO_RELIEF_CEILING
                                     : MASKING_VETO_PERSISTENCE_RELIEF_CEILING;
    const float target =
        valid ? mv_floor_target(band_noise_energy, band_threshold,
                                self->spectrum_type, relief_ceiling)
              : 0.0F;

    // Mapper regularization: symmetric one-pole (processing-time tau, not a
    // psychoacoustic decay) keeps the model-informed floor from flickering
    // when the lookahead/clean estimates cross the noise. The model's
    // temporal layer alone owns all psychoacoustic persistence in time.
    self->band_floor_state[j] =
        ((1.0F - self->floor_smooth_memory) * target) +
        (self->floor_smooth_memory * self->band_floor_state[j]);
    self->band_floor_smooth[j] = self->band_floor_state[j];
  }

  // Spectral stabilization: smooth the floor ACROSS bands to prevent sharp
  // spectral edges at Bark boundaries.
  spectral_smoothing_apply_spatial(self->band_floor_smooth, num_bands);
}

static void mv_apply_floor(const MaskingVeto* self, float* gain, float depth,
                           uint32_t num_bands) {
  /**
   * 4. Interpolated audibility relief, applied per bin between band centers
   * (avoids sharp gain transitions at Bark band boundaries):
   *   gain = max(gain, gain + depth * (floor - gain))
   * Depth-lerped upward only: the worst-case residual per band is bounded
   * by the masked-audibility level, so relief flicker cannot excavate.
   */
  uint32_t current_band = 0;
  for (uint32_t k = 0; k < self->real_spectrum_size; k++) {
    while (current_band < num_bands - 1 &&
           (float)k > self->band_centers[current_band + 1]) {
      current_band++;
    }

    float bin_floor;
    if (current_band < num_bands - 1) {
      const float x0 = self->band_centers[current_band];
      const float x1 = self->band_centers[current_band + 1];
      const float y0 = self->band_floor_smooth[current_band];
      const float y1 = self->band_floor_smooth[current_band + 1];

      if (x1 > x0) {
        const float t = ((float)k - x0) / (x1 - x0);
        const float t_clamped = fmaxf(0.0F, fminf(1.0F, t));
        bin_floor = (y0 * (1.0F - t_clamped)) + (y1 * t_clamped);
      } else {
        bin_floor = y0;
      }
    } else {
      bin_floor = self->band_floor_smooth[num_bands - 1];
    }

    if (bin_floor > gain[k]) {
      gain[k] += depth * (bin_floor - gain[k]);
    }
  }
}

void masking_veto_apply_floor(MaskingVeto* self, const float* spectrum,
                              const float* noise_spectrum,
                              const float* future_spectrum, float* gain,
                              float depth) {
  if (!self || !spectrum || !noise_spectrum || !gain || !(depth > 0.0F) ||
      self->real_spectrum_size == 0U) {
    return;
  }

  // 1. Estimate clean signal magnitude from the smoothed spectrum
  const float* future_clean_estimation =
      mv_update_clean_est(self, spectrum, noise_spectrum, future_spectrum);

  // 2. Validate per-band masker evidence: bands whose clean estimate does
  // not exceed the noise are untrusted this frame, so the model's
  // frequency threshold and temporal memory are fed only by real maskers
  // (biased noise estimates otherwise accumulate phantom persistence in
  // noise-only sections right through any persistence-based relief).
  mv_compute_band_masker_validation(self, noise_spectrum,
                                    future_clean_estimation);

  // 3. Compute psychoacoustic masking thresholds
  if (!compute_masking_thresholds(
          self->masking_estimator, self->clean_signal_estimation,
          future_clean_estimation, self->band_masker_mask,
          self->masking_thresholds)) {
    return;
  }

  // 4 + 5. Band floors (stabilized) interpolated into a relief on the gain
  mv_compute_band_floor(self, noise_spectrum, future_clean_estimation);
  mv_apply_floor(self, gain, fminf(depth, 1.0F),
                 get_number_of_critical_bands(self->critical_bands_helper));
}

void masking_veto_set_hop_sec(MaskingVeto* self, float hop_sec) {
  if (!self || !(hop_sec > 0.0F)) {
    return;
  }
  self->smooth =
      sb_alpha_retuned(MASKING_VETO_SMOOTHING, LEGACY_REF_HOP_SEC, hop_sec);
  self->floor_smooth_memory =
      sb_alpha_from_tau(MASKING_VETO_FLOOR_SMOOTH_SEC, hop_sec);
  masking_estimation_set_hop_sec(self->masking_estimator, hop_sec);
}
