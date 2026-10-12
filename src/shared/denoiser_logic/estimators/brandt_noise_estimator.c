/*
libspecbleach - A spectral processing library

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

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

#include "shared/denoiser_logic/estimators/brandt_noise_estimator.h"
#include "shared/configurations.h"
#include "shared/frame_rate_norm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct BrandtNoiseEstimator {
  uint32_t spectrum_size;
  uint32_t history_size;
  uint32_t history_index; // Circular buffer write head
  uint32_t stats_interval;
  float history_duration_ms; // For true-hop history rebuild (init-time only)

  float* history_buffer;      // Size: spectrum_size * history_size
  float* sorted_history;      // Sorted copy of each bin's history
  float* last_noise_spectrum; // Persisted noise to keep during rejection

  float percentile;
  uint32_t trim_count; // Number of items to average (history_size * percentile)
  float correction_factor;
  float correction_factors[5]; // Pre-calculated for search
  uint32_t active_frame_count;
  bool is_first_frame;
};

// Helper: Calculate bias correction factor for trimmed mean of exponential dist
// Formula: Factor = 1 / (1 + ( (1-P)/P * ln(1-P) ))
static float calculate_correction_factor(float p) {
  if (p <= 0.0f || p >= 1.0f) {
    return 1.0f; // Invalid P, no correction
  }
  float term = (1.0f - p) / p * logf(1.0f - p);
  float denominator = 1.0f + term;
  if (fabsf(denominator) < BRANDT_ESTIMATOR_BIAS_EPSILON) {
    return 1.0f; // Avoid division by zero
  }
  return 1.0f / denominator;
}

BrandtNoiseEstimator* brandt_noise_estimator_initialize(
    uint32_t spectrum_size, float history_duration_ms, uint32_t sample_rate,
    uint32_t fft_size) {
  float percentile = BRANDT_DEFAULT_PERCENTILE; // Advised for music restoration
  BrandtNoiseEstimator* self =
      (BrandtNoiseEstimator*)calloc(1, sizeof(BrandtNoiseEstimator));
  if (!self) {
    return NULL;
  }

  self->spectrum_size = spectrum_size;
  self->percentile = percentile;
  self->history_duration_ms = history_duration_ms;

  // Calculate history size from duration. hop == frame/OVERLAP_FACTOR;
  // the old *0.5 (50% overlap) assumption is fixed to /OVERLAP_FACTOR here.
  float ms_per_frame = (float)fft_size * 1000.0f / (float)sample_rate;
  // Overlap consideration: Usually frame step is hop_size.
  // Assuming hop = fft_size/2 or similar? The caller usually provides
  // parameters. If exact duration needed, we might need hop_size.
  // Fixed 4x-overlap conversion (was 50% approximation):
  float frame_duration = ms_per_frame / (float)OVERLAP_FACTOR;
  if (frame_duration < BRANDT_ESTIMATOR_MIN_DURATION_MS) {
    frame_duration = BRANDT_ESTIMATOR_MIN_DURATION_MS;
  }

  self->history_size = (uint32_t)(history_duration_ms / frame_duration);
  if (self->history_size < BRANDT_ESTIMATOR_MIN_HISTORY_FRAMES) {
    self->history_size = BRANDT_ESTIMATOR_MIN_HISTORY_FRAMES; // Minimum history
  }

  self->trim_count = (uint32_t)((float)self->history_size * percentile);
  if (self->trim_count < 1) {
    self->trim_count = 1; // At least min
  }
  if (self->trim_count > self->history_size) {
    self->trim_count = self->history_size;
  }

  self->correction_factor = calculate_correction_factor(percentile);

  static const float p_candidates[] = {0.1f, 0.25f, 0.5f, 0.75f, 1.0f};
  for (int i = 0; i < 5; i++) {
    self->correction_factors[i] = calculate_correction_factor(p_candidates[i]);
  }

  // Allocate buffers
  self->history_buffer = (float*)calloc(
      (size_t)self->spectrum_size * self->history_size, sizeof(float));
  self->sorted_history = (float*)calloc(
      (size_t)self->spectrum_size * self->history_size, sizeof(float));
  self->last_noise_spectrum =
      (float*)calloc(self->spectrum_size, sizeof(float));

  if (!self->history_buffer || !self->sorted_history ||
      !self->last_noise_spectrum) {
    brandt_noise_estimator_free(self);
    return NULL;
  }

  self->is_first_frame = true;
  self->active_frame_count = 0;
  self->stats_interval = BRANDT_ESTIMATOR_STATS_UPDATE_INTERVAL_FRAMES;
  return self;
}

void brandt_noise_estimator_free(BrandtNoiseEstimator* self) {
  if (self) {
    free(self->history_buffer);
    free(self->sorted_history);
    free(self->last_noise_spectrum);
    free(self);
  }
}

// Helper: Fast inline Shell sort for float arrays
static inline void sb_sort_floats_inline(float* arr, uint32_t n) {
  if (n < 2) {
    return;
  }
  uint32_t h = 1;
  while (h < n / 3) {
    h = (3 * h) + 1;
  }
  while (h >= 1) {
    for (uint32_t i = h; i < n; i++) {
      float temp = arr[i];
      uint32_t j = i;
      while (j >= h && arr[j - h] > temp) {
        arr[j] = arr[j - h];
        j -= h;
      }
      arr[j] = temp;
    }
    h /= 3;
  }
}

// Replace one value in an already sorted history without sorting the full
// window. The old value is guaranteed to be present because it is read from
// the corresponding circular-history slot immediately before that slot is
// overwritten.
static inline bool sb_update_sorted_history(float* sorted, uint32_t n,
                                            float old_value, float new_value) {
  if (n == 0U) {
    return false;
  }

  uint32_t low = 0U;
  uint32_t high = n;
  while (low < high) {
    uint32_t middle = low + ((high - low) / 2U);
    if (sorted[middle] < old_value) {
      low = middle + 1U;
    } else {
      high = middle;
    }
  }
  if (low == n || sorted[low] != old_value) {
    return false;
  }

  memmove(&sorted[low], &sorted[low + 1U],
          (size_t)(n - low - 1U) * sizeof(float));

  // Insert after equal values. This leaves the sequence sorted and avoids
  // searching through duplicate runs on subsequent replacements.
  uint32_t remaining = n - 1U;
  low = 0U;
  high = remaining;
  while (low < high) {
    uint32_t middle = low + ((high - low) / 2U);
    if (sorted[middle] <= new_value) {
      low = middle + 1U;
    } else {
      high = middle;
    }
  }
  memmove(&sorted[low + 1U], &sorted[low],
          (size_t)(remaining - low) * sizeof(float));
  sorted[low] = new_value;
  return true;
}

static void sb_rebuild_sorted_history(BrandtNoiseEstimator* self,
                                      uint32_t bin) {
  const float* history =
      &self->history_buffer[(size_t)bin * self->history_size];
  float* sorted = &self->sorted_history[(size_t)bin * self->history_size];
  memcpy(sorted, history, self->history_size * sizeof(float));
  sb_sort_floats_inline(sorted, self->history_size);
}

static void sb_fill_jittered_history(BrandtNoiseEstimator* self, uint32_t bin,
                                     float value) {
  uint32_t jitter_counts[11] = {0U};
  size_t bin_offset = (size_t)bin * self->history_size;
  float* history = &self->history_buffer[bin_offset];
  float* sorted = &self->sorted_history[bin_offset];

  for (uint32_t t = 0U; t < self->history_size; t++) {
    uint32_t jitter_index = (t + bin) % 11U;
    float jitter = 1.0f + (0.01f * ((float)jitter_index - 5.0f) / 5.0f);
    history[t] = value * jitter;
    jitter_counts[jitter_index]++;
  }

  uint32_t sorted_index = 0U;
  for (uint32_t i = 0U; i < 11U; i++) {
    uint32_t jitter_index = (value < 0.0f) ? (10U - i) : i;
    float jitter = 1.0f + (0.01f * ((float)jitter_index - 5.0f) / 5.0f);
    float jittered_value = value * jitter;
    for (uint32_t count = 0U; count < jitter_counts[jitter_index]; count++) {
      sorted[sorted_index++] = jittered_value;
    }
  }
}

static float calculate_ad_norm(const float* sorted, uint32_t q, float mu,
                               float b) {
  if (mu < 1e-15f) {
    return 1.0f;
  }
  float mu_inv = 1.0f / mu;
  float exp_b_mu = (b * mu_inv > 20.0f) ? 0.0f : expf(-b * mu_inv);
  float denom = 1.0f - exp_b_mu;
  if (fabsf(denom) < SPECTRAL_EPSILON) {
    return 1.0f;
  }

  float abs_diff_sum = 0.0f;
  float q_inv = 1.0f / (float)q;
  float denom_inv = 1.0f / denom;

  for (uint32_t i = 0; i < q; i++) {
    float val = sorted[i] * mu_inv;
    float f_te;
    if (val > 20.0f) {
      // Since sorted array is ascending, all remaining elements will also have
      // val > 20.0f
      for (uint32_t j = i; j < q; j++) {
        float f_emp = (float)(j + 1) * q_inv;
        abs_diff_sum += fabsf(f_emp - denom_inv);
      }
      break;
    }
    f_te = (1.0f - expf(-val)) * denom_inv;
    float f_empirical = (float)(i + 1) * q_inv;
    abs_diff_sum += fabsf(f_empirical - f_te);
  }
  return abs_diff_sum * (2.0f * q_inv);
}

bool brandt_noise_estimator_run(BrandtNoiseEstimator* self,
                                const float* spectrum, float* noise_spectrum) {
  if (!self || !spectrum || !noise_spectrum || self->history_size == 0) {
    return false;
  }

  float frame_energy = 0.F;
  for (uint32_t k = 0U; k < self->spectrum_size; k++) {
    frame_energy += spectrum[k];
  }
  frame_energy /= (float)self->spectrum_size;

  if (self->is_first_frame && frame_energy > ESTIMATOR_SILENCE_THRESHOLD) {
    float inv_factor = 1.0f / calculate_correction_factor(0.5f);
    for (uint32_t k = 0; k < self->spectrum_size; k++) {
      float val = spectrum[k] * inv_factor;
      self->last_noise_spectrum[k] = spectrum[k];
      sb_fill_jittered_history(self, k, val);
    }
    self->is_first_frame = false;
  }

  if (frame_energy < ESTIMATOR_SILENCE_THRESHOLD) {
    memcpy(noise_spectrum, self->last_noise_spectrum,
           self->spectrum_size * sizeof(float));
    return true;
  }

  // Record current frame into circular buffer
  uint32_t current_idx = self->history_index;
  for (uint32_t k = 0; k < self->spectrum_size; k++) {
    size_t bin_offset = (size_t)k * self->history_size;
    float old_value = self->history_buffer[bin_offset + current_idx];
    self->history_buffer[bin_offset + current_idx] = spectrum[k];
    if (!sb_update_sorted_history(&self->sorted_history[bin_offset],
                                  self->history_size, old_value, spectrum[k])) {
      sb_rebuild_sorted_history(self, k);
    }
  }
  self->history_index = (current_idx + 1) % self->history_size;

  // Subsample expensive statistical update every N frames (or on first frame)
  const uint32_t stats_interval =
      (self->stats_interval > 0U)
          ? self->stats_interval
          : BRANDT_ESTIMATOR_STATS_UPDATE_INTERVAL_FRAMES;
  bool update_stats = self->is_first_frame ||
                      ((self->active_frame_count % stats_interval) == 0);
  self->active_frame_count++;

  if (update_stats) {
    static const float p_candidates[] = {0.1f, 0.25f, 0.5f, 0.75f, 1.0f};
    uint32_t q_candidates[5];
    for (int i = 0; i < 5; i++) {
      q_candidates[i] = (uint32_t)(p_candidates[i] * (float)self->history_size);
    }

    for (uint32_t k = 0; k < self->spectrum_size; k++) {
      const float* sorted =
          &self->sorted_history[(size_t)k * self->history_size];

      float min_ad_norm = 2.0f;
      float best_mu = self->last_noise_spectrum[k];

      float prefix_sum = 0.0f;
      uint32_t candidate_index = 0U;
      for (uint32_t j = 0; j < self->history_size; j++) {
        prefix_sum += sorted[j];
        while (candidate_index < 5U &&
               q_candidates[candidate_index] <= j + 1U) {
          uint32_t q = q_candidates[candidate_index];
          if (q >= 10U) {
            float mu_trunc = prefix_sum / (float)q;
            if (mu_trunc > ESTIMATOR_SILENCE_THRESHOLD) {
              float b = sorted[q - 1U];
              float mu_full =
                  mu_trunc * self->correction_factors[candidate_index];
              float ad_norm = calculate_ad_norm(sorted, q, mu_full, b);
              if (ad_norm < min_ad_norm) {
                min_ad_norm = ad_norm;
                best_mu = mu_full;
              }
            }
          }
          candidate_index++;
        }
      }

      if (1.0f - min_ad_norm >= BRANDT_MIN_CONFIDENCE) {
        self->last_noise_spectrum[k] = best_mu;
      }
      noise_spectrum[k] = self->last_noise_spectrum[k];
    }
  } else {
    // Copy persistent noise estimate
    memcpy(noise_spectrum, self->last_noise_spectrum,
           self->spectrum_size * sizeof(float));
  }

  return true;
}

void brandt_noise_estimator_set_state(BrandtNoiseEstimator* self,
                                      const float* initial_profile) {
  if (!self || !initial_profile) {
    return;
  }

  // Start with a known state: Fill history with this profile
  float inverse_factor = 1.0f / self->correction_factor;

  for (uint32_t k = 0; k < self->spectrum_size; k++) {
    float val = initial_profile[k] * inverse_factor;
    self->last_noise_spectrum[k] = initial_profile[k];
    sb_fill_jittered_history(self, k, val);
  }
  self->is_first_frame = false;
}

void brandt_noise_estimator_update_seed(BrandtNoiseEstimator* self,
                                        const float* seed_profile) {
  if (!self || !seed_profile) {
    return;
  }
  // Similar to set_state but maybe only updates part of history?
  // For now, treat same as set_state to ensure quick convergence/reset.
  brandt_noise_estimator_set_state(self, seed_profile);
  self->is_first_frame = false;
}

void brandt_noise_estimator_apply_floor(BrandtNoiseEstimator* self,
                                        const float* floor_profile) {
  if (!self || !floor_profile) {
    return;
  }
  float inverse_factor = 1.0f / self->correction_factor;
  for (uint32_t k = 0; k < self->spectrum_size; k++) {
    float floor_val = floor_profile[k] * inverse_factor;
    size_t bin_offset = (size_t)k * self->history_size;
    for (uint32_t t = 0; t < self->history_size; t++) {
      if (self->history_buffer[bin_offset + t] < floor_val) {
        self->history_buffer[bin_offset + t] = floor_val;
      }
      if (self->sorted_history[bin_offset + t] < floor_val) {
        self->sorted_history[bin_offset + t] = floor_val;
      }
    }
  }
}

void brandt_noise_estimator_set_hop_sec(BrandtNoiseEstimator* self,
                                        float hop_sec) {
  if (!self || !(hop_sec > 0.0F)) {
    return;
  }
  uint32_t interval =
      sb_frames_for_ms(BRANDT_ESTIMATOR_STATS_UPDATE_MS, hop_sec, 1U, 16U);
  self->stats_interval = interval;

  // Rebuild history storage for the true hop when it differs from the
  // fft-derived approximation used at init. Init-time only: never called
  // from the audio thread (spectral_denoiser calls it during setup).
  float hop_ms = hop_sec * 1000.0F;
  if (hop_ms < BRANDT_ESTIMATOR_MIN_DURATION_MS) {
    hop_ms = BRANDT_ESTIMATOR_MIN_DURATION_MS;
  }
  uint32_t history_size = (uint32_t)(self->history_duration_ms / hop_ms);
  if (history_size < BRANDT_ESTIMATOR_MIN_HISTORY_FRAMES) {
    history_size = BRANDT_ESTIMATOR_MIN_HISTORY_FRAMES;
  }
  if (history_size == self->history_size) {
    return;
  }
  float* history_buffer =
      (float*)calloc((size_t)self->spectrum_size * history_size, sizeof(float));
  float* sorted_history =
      (float*)calloc((size_t)self->spectrum_size * history_size, sizeof(float));
  if (!history_buffer || !sorted_history) {
    free(history_buffer);
    free(sorted_history);
    return; // keep existing storage on allocation failure
  }
  free(self->history_buffer);
  free(self->sorted_history);
  self->history_buffer = history_buffer;
  self->sorted_history = sorted_history;
  self->history_size = history_size;
  self->history_index = 0U;
  self->trim_count = (uint32_t)((float)history_size * self->percentile);
  if (self->trim_count < 1U) {
    self->trim_count = 1U;
  }
  if (self->trim_count > history_size) {
    self->trim_count = history_size;
  }
}
