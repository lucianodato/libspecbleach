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

#include "masking_estimator.h"
#include "../configurations.h"
#include "../frame_rate_norm.h"
#include "shared/utils/absolute_hearing_thresholds.h"
#include "shared/utils/critical_bands.h"
#include <math.h>
#include <stdlib.h>

// Note: Psychoacoustic constants are now imported from configurations.h

typedef struct SpreadingParams {
  float s_up;
  float s_total;
  float s_offset;
  float norm_factor;
  float y_shift;
} SpreadingParams;

static SpreadingParams compute_spreading_params(float level_db);
static float evaluate_spreading_gain(float dz, const SpreadingParams* params);

/**
 * compute_tonality_factor: SFM-based NMT/TMN classification.
 * Uses Spectral Flatness Measure to differentiate between Tone-masking-Noise
 * and Noise-masking-Tone, adjusting the masking offset accordingly.
 */
static float compute_tonality_factor(MaskingEstimator* self,
                                     const float* spectrum, uint32_t band);

typedef struct MaskingConfig {
  uint32_t fft_size;
  uint32_t sample_rate;
  uint32_t number_critical_bands;
  uint32_t real_spectrum_size;
  float spectral_additivity_exponent;
  bool use_temporal_masking;
  float backward_decay;
  bool absolute_threshold_enabled;
} MaskingConfig;

typedef struct MaskingOwned {
  AbsoluteHearingThresholds* reference_spectrum;
  CriticalBands* critical_bands;
  CriticalBandIndexes band_indexes;
} MaskingOwned;

typedef struct MaskingState {
  float* critical_bands_spectrum;
  float* critical_bands_reference_spectrum;
  float* spreading_matrix; // Matrix for simultaneous masking
  float* masking_offset;
  float* previous_thresholds;
  float* forward_decays;
  float* future_thresholds;
  float* absolute_threshold_cb;
} MaskingState;

typedef struct MaskingScratch {
  // Temporary buffers to avoid VLAs and stack overflows
  float* future_cb_spectrum_buf;
  float* bark_levels_buf;
  float* spreaded_future_buf;
  float* spreaded_current_buf;
  SpreadingParams* spreading_params_buf;
} MaskingScratch;

struct MaskingEstimator {
  MaskingConfig config;
  MaskingOwned owned;
  MaskingState state;
  MaskingScratch scratch;
};

MaskingEstimator* masking_estimation_initialize(
    const uint32_t fft_size, const uint32_t sample_rate,
    CriticalBandType critical_band_type, SpectrumType spectrum_type,
    bool use_absolute_threshold, bool use_temporal_masking) {

  MaskingEstimator* self =
      (MaskingEstimator*)calloc(1U, sizeof(MaskingEstimator));

  if (!self) {
    return NULL;
  }

  self->config.fft_size = fft_size;
  self->config.real_spectrum_size = (self->config.fft_size / 2U) + 1U;
  self->config.sample_rate = sample_rate;

  self->owned.critical_bands = critical_bands_initialize(
      self->config.sample_rate, self->config.fft_size, critical_band_type);
  if (!self->owned.critical_bands) {
    masking_estimation_free(self);
    return NULL;
  }
  self->config.number_critical_bands =
      get_number_of_critical_bands(self->owned.critical_bands);

  self->state.critical_bands_spectrum =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->state.critical_bands_reference_spectrum =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->state.spreading_matrix =
      (float*)calloc((size_t)self->config.number_critical_bands *
                         (size_t)self->config.number_critical_bands,
                     sizeof(float));
  self->state.masking_offset =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->state.previous_thresholds =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->state.future_thresholds =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->state.forward_decays =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->state.absolute_threshold_cb =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));

  self->scratch.future_cb_spectrum_buf =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->scratch.bark_levels_buf =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->scratch.spreaded_future_buf =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->scratch.spreaded_current_buf =
      (float*)calloc(self->config.number_critical_bands, sizeof(float));
  self->scratch.spreading_params_buf = (SpreadingParams*)calloc(
      self->config.number_critical_bands, sizeof(SpreadingParams));

  self->owned.reference_spectrum = absolute_hearing_thresholds_initialize(
      self->config.sample_rate, self->config.fft_size, spectrum_type);

  self->config.spectral_additivity_exponent = SPECTRAL_ADDITIVITY_EXPONENT_PEAQ;
  self->config.use_temporal_masking = use_temporal_masking;
  self->config.absolute_threshold_enabled = use_absolute_threshold;

  if (!self->state.critical_bands_spectrum ||
      !self->state.critical_bands_reference_spectrum ||
      !self->state.spreading_matrix || !self->state.masking_offset ||
      !self->state.previous_thresholds || !self->state.future_thresholds ||
      !self->state.forward_decays || !self->state.absolute_threshold_cb ||
      !self->scratch.future_cb_spectrum_buf || !self->scratch.bark_levels_buf ||
      !self->scratch.spreaded_future_buf ||
      !self->scratch.spreaded_current_buf ||
      !self->scratch.spreading_params_buf || !self->owned.reference_spectrum) {
    masking_estimation_free(self);
    return NULL;
  }

  // Temporal masking decay constants. hop_time must be the TRUE hop
  // (frame/overlap/sr); call masking_estimation_set_hop_sec() after init when
  // the true hop is known (FFT-padded sizes overestimate it slightly).
  const float hop_time = (float)fft_size / (4.0F * (float)sample_rate);

  // Frequency-dependent forward masking (Low: 100ms, High: 25ms)
  for (uint32_t j = 0U; j < self->config.number_critical_bands; j++) {
    const float bark = fminf((float)j, 24.0F);
    const float weight = bark / 24.0F; // 0 to 1
    const float tau = ((1.0F - weight) * FORWARD_MASKING_TAU_LOW_SEC) +
                      (weight * FORWARD_MASKING_TAU_HIGH_SEC);
    self->state.forward_decays[j] = expf(-hop_time / tau);
  }

  // Backward masking (10ms) remains constant across frequency
  self->config.backward_decay = expf(-hop_time / BACKWARD_MASKING_TAU_SEC);

  return self;
}

void masking_estimation_set_hop_sec(MaskingEstimator* self, float hop_sec) {
  if (!self || !(hop_sec > 0.0F) || !self->state.forward_decays) {
    return;
  }
  for (uint32_t j = 0U; j < self->config.number_critical_bands; j++) {
    const float bark = fminf((float)j, 24.0F);
    const float weight = bark / 24.0F;
    const float tau = ((1.0F - weight) * FORWARD_MASKING_TAU_LOW_SEC) +
                      (weight * FORWARD_MASKING_TAU_HIGH_SEC);
    self->state.forward_decays[j] = expf(-hop_sec / tau);
  }
  self->config.backward_decay = expf(-hop_sec / BACKWARD_MASKING_TAU_SEC);
}

void masking_estimation_free(MaskingEstimator* self) {
  if (!self) {
    return;
  }
  absolute_hearing_thresholds_free(self->owned.reference_spectrum);
  critical_bands_free(self->owned.critical_bands);

  free(self->state.critical_bands_spectrum);
  free(self->state.critical_bands_reference_spectrum);
  free(self->state.spreading_matrix);
  free(self->state.masking_offset);
  free(self->state.previous_thresholds);
  free(self->state.future_thresholds);
  free(self->state.forward_decays);
  free(self->state.absolute_threshold_cb);
  free(self->scratch.future_cb_spectrum_buf);
  free(self->scratch.bark_levels_buf);
  free(self->scratch.spreaded_future_buf);
  free(self->scratch.spreaded_current_buf);
  free(self->scratch.spreading_params_buf);

  free(self);
}

bool compute_masking_thresholds(MaskingEstimator* self, const float* spectrum,
                                const float* future_spectrum,
                                float* masking_thresholds) {
  if (!self || !spectrum || !masking_thresholds) {
    return false;
  }

  compute_critical_bands_spectrum(self->owned.critical_bands, spectrum,
                                  self->state.critical_bands_spectrum);

  const float spectral_p = self->config.spectral_additivity_exponent;
  const float spectral_inv_p = 1.0F / spectral_p;

  // 1. Calculate spreaded future spectrum (Frequency Masking only)
  if (future_spectrum) {
    compute_critical_bands_spectrum(self->owned.critical_bands, future_spectrum,
                                    self->scratch.future_cb_spectrum_buf);

    for (uint32_t j = 0U; j < self->config.number_critical_bands; j++) {
      self->scratch.bark_levels_buf[j] =
          (10.F *
           log10f(self->scratch.future_cb_spectrum_buf[j] + SPECTRAL_EPSILON)) +
          DB_FS_TO_SPL_REF;
      self->scratch.spreading_params_buf[j] =
          compute_spreading_params(self->scratch.bark_levels_buf[j]);
    }

    for (uint32_t i = 0U; i < self->config.number_critical_bands; i++) {
      float spreaded_p = 0.F;
      for (uint32_t j = 0U; j < self->config.number_critical_bands; j++) {
        const float dz = (float)i - (float)j;
        const float gain =
            evaluate_spreading_gain(dz, &self->scratch.spreading_params_buf[j]);
        spreaded_p +=
            powf(self->scratch.future_cb_spectrum_buf[j] * gain, spectral_p);
      }
      self->scratch.spreaded_future_buf[i] = powf(spreaded_p, spectral_inv_p);
    }

    for (uint32_t j = 0U; j < self->config.number_critical_bands; j++) {
      const float tonality_factor =
          compute_tonality_factor(self, future_spectrum, j);
      const float bark_idx = fminf((float)(j + 1), 25.0F);
      // Tone-masking-noise (TMN) has higher offset (~16-41dB depending on Bark)
      // Noise-masking-tone (NMT) has lower offset (~6dB)
      const float offset = (tonality_factor * (TMN_OFFSET_BASE + bark_idx)) +
                           (NMT_OFFSET_DB * (1.F - tonality_factor));

      self->state.future_thresholds[j] = powf(
          10.F,
          (log10f(self->scratch.spreaded_future_buf[j] + SPECTRAL_EPSILON) -
           (offset / 10.F)));
    }
  }

  for (uint32_t j = 0U; j < self->config.number_critical_bands; j++) {
    self->scratch.bark_levels_buf[j] =
        (10.F *
         log10f(self->state.critical_bands_spectrum[j] + SPECTRAL_EPSILON)) +
        DB_FS_TO_SPL_REF;
    self->scratch.spreading_params_buf[j] =
        compute_spreading_params(self->scratch.bark_levels_buf[j]);
  }

  for (uint32_t i = 0U; i < self->config.number_critical_bands; i++) {
    float spreaded_p = 0.F;
    for (uint32_t j = 0U; j < self->config.number_critical_bands; j++) {
      const float dz = (float)i - (float)j;
      const float gain =
          evaluate_spreading_gain(dz, &self->scratch.spreading_params_buf[j]);
      spreaded_p +=
          powf(self->state.critical_bands_spectrum[j] * gain, spectral_p);
    }
    self->scratch.spreaded_current_buf[i] = powf(spreaded_p, spectral_inv_p);
  }

  for (uint32_t j = 0U; j < self->config.number_critical_bands; j++) {

    const float tonality_factor = compute_tonality_factor(self, spectrum, j);
    const float bark_idx = fminf((float)(j + 1), 25.0F);

    self->state.masking_offset[j] =
        (tonality_factor * (TMN_OFFSET_BASE + bark_idx)) +
        (NMT_OFFSET_DB * (1.F - tonality_factor));

    // 1. Calculate frequency masking threshold for current frame
    float threshold =
        powf(10.F,
             (log10f(self->scratch.spreaded_current_buf[j] + SPECTRAL_EPSILON) -
              (self->state.masking_offset[j] / 10.F)));

    // 2. Combine with temporal masking using Power Law (p=0.6)
    // Total_T = (T_freq^p + T_forward^p + T_backward^p)^(1/p)
    // This model (Johnston, 1988) better reflects the non-linear summation
    // of multiple maskers compared to simple linear addition.
    if (self->config.use_temporal_masking) {
      float threshold_p = powf(threshold, POWER_LAW_EXPONENT);

      // Add forward masking contribution
      float forward_threshold =
          self->state.previous_thresholds[j] * self->state.forward_decays[j];
      threshold_p += powf(forward_threshold, POWER_LAW_EXPONENT);

      // Add backward masking contribution if available
      if (future_spectrum) {
        float backward_threshold =
            self->state.future_thresholds[j] * self->config.backward_decay;
        threshold_p += powf(backward_threshold, POWER_LAW_EXPONENT);
      }

      threshold = powf(threshold_p, 1.0F / POWER_LAW_EXPONENT);
    }

    self->state.previous_thresholds[j] =
        threshold; // Update state for next frame

    self->owned.band_indexes = get_band_indexes(self->owned.critical_bands, j);

    for (uint32_t k = self->owned.band_indexes.start_position;
         k < self->owned.band_indexes.end_position; k++) {
      masking_thresholds[k] = threshold;
    }
  }

  if (self->config.absolute_threshold_enabled) {
    apply_thresholds_as_floor(self->owned.reference_spectrum,
                              masking_thresholds);
  }

  return true;
}

static SpreadingParams compute_spreading_params(float level_db) {
  const float s_up =
      fminf(fmaxf(S_MAX_UPWARD - ((level_db - S_LEVEL_REF_DB) * S_SLOPE_FACTOR),
                  S_MIN_UPWARD),
            S_MAX_UPWARD);
  const float s_total = (S_DOWNWARD + s_up) * 0.5F;
  const float s_offset = (S_DOWNWARD - s_up) * 0.5F;

  // By algebra, norm equals sqrt(S_DOWNWARD*s_up)
  const float norm_factor = sqrtf(S_DOWNWARD * s_up);
  const float inv_norm_factor = 1.0F / norm_factor;
  const float y_shift = s_offset * inv_norm_factor;

  SpreadingParams params = {
      .s_up = s_up,
      .s_total = s_total,
      .s_offset = s_offset,
      .norm_factor = norm_factor,
      .y_shift = y_shift,
  };
  return params;
}

static float evaluate_spreading_gain(float dz, const SpreadingParams* params) {
  const float y = dz + params->y_shift;
  const float sf_db = params->norm_factor + (params->s_offset * y) -
                      (params->s_total * sqrtf(1.0F + (y * y)));

  return powf(10.0F, sf_db * 0.1F);
}

static float compute_tonality_factor(MaskingEstimator* self,
                                     const float* spectrum, uint32_t band) {
  float sum_bins = 0.F;
  float sum_log_bins = 0.F;

  self->owned.band_indexes = get_band_indexes(self->owned.critical_bands, band);

  for (uint32_t k = self->owned.band_indexes.start_position;
       k < self->owned.band_indexes.end_position; k++) {
    const float val = fmaxf(spectrum[k], SPECTRAL_EPSILON);
    sum_bins += val;
    sum_log_bins += log10f(val);
  }

  float bins_in_band = (float)self->owned.band_indexes.end_position -
                       (float)self->owned.band_indexes.start_position;

  if (bins_in_band <= 1.0F) {
    return 1.0F;
  }

  const float sfm =
      10.F * ((sum_log_bins / bins_in_band) - log10f(sum_bins / bins_in_band));

  // SFM range is mapped to [0, 1] tonality factor
  const float tonality_factor =
      fminf(fmaxf((sfm - SFM_MAX_DB) / (SFM_MIN_DB - SFM_MAX_DB), 0.0F), 1.0F);

  return tonality_factor;
}
