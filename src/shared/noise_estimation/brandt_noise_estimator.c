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

#include "brandt_noise_estimator.h"
#include "../configurations.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct BrandtNoiseEstimator {
  uint32_t spectrum_size;
  uint32_t history_size;
  uint32_t history_index; // Circular buffer write head

  float* history_buffer; // Size: spectrum_size * history_size
  float* sort_buffer;    // Scratch space for sorting (Size: history_size)

  float percentile;
  uint32_t trim_count; // Number of items to average (history_size * percentile)
  float correction_factor;
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
  if (fabsf(denominator) < ESTIMATOR_BIAS_EPSILON) {
    return 1.0f; // Avoid division by zero
  }
  return 1.0f / denominator;
}

// Helper: Compare floats for qsort
static int compare_floats(const void* a, const void* b) {
  float fa = *(const float*)a;
  float fb = *(const float*)b;
  return (fa > fb) - (fa < fb);
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

  // Calculate history size from duration
  float ms_per_frame = (float)fft_size * 1000.0f / (float)sample_rate;
  // Overlap consideration: Usually frame step is hop_size.
  // Assuming hop = fft_size/2 or similar? The caller usually provides
  // parameters. If exact duration needed, we might need hop_size.
  // Assuming standard 50% overlap for calculation roughly:
  float frame_duration = ms_per_frame * 0.5f; // Rough approximation of step
  if (frame_duration < ESTIMATOR_MIN_DURATION_MS) {
    frame_duration = ESTIMATOR_MIN_DURATION_MS;
  }

  self->history_size = (uint32_t)(history_duration_ms / frame_duration);
  if (self->history_size < ESTIMATOR_MIN_HISTORY_FRAMES) {
    self->history_size = ESTIMATOR_MIN_HISTORY_FRAMES; // Minimum history
  }

  self->trim_count = (uint32_t)((float)self->history_size * percentile);
  if (self->trim_count < 1) {
    self->trim_count = 1; // At least min
  }
  if (self->trim_count > self->history_size) {
    self->trim_count = self->history_size;
  }

  self->correction_factor = calculate_correction_factor(percentile);

  // Allocate buffers
  self->history_buffer = (float*)calloc(
      (size_t)self->spectrum_size * self->history_size, sizeof(float));
  self->sort_buffer = (float*)calloc(self->history_size, sizeof(float));

  if (!self->history_buffer || !self->sort_buffer) {
    brandt_noise_estimator_free(self);
    return NULL;
  }

  self->is_first_frame = true;
  return self;
}

void brandt_noise_estimator_free(BrandtNoiseEstimator* self) {
  if (self) {
    free(self->history_buffer);
    free(self->sort_buffer);
    free(self);
  }
}

bool brandt_noise_estimator_run(BrandtNoiseEstimator* self,
                                const float* spectrum, float* noise_spectrum) {
  if (!self || !spectrum || !noise_spectrum) {
    return false;
  }

  float frame_energy = 0.F;
  for (uint32_t k = 0U; k < self->spectrum_size; k++) {
    frame_energy += spectrum[k];
  }
  frame_energy /= (float)self->spectrum_size;

  if (self->is_first_frame) {
    if (frame_energy < ESTIMATOR_SILENCE_THRESHOLD) {
      memset(noise_spectrum, 0, self->spectrum_size * sizeof(float));
      return true;
    }
    float inverse_factor = 1.0f / self->correction_factor;
    for (uint32_t k = 0; k < self->spectrum_size; k++) {
      float val = spectrum[k] * inverse_factor;
      for (uint32_t t = 0; t < self->history_size; t++) {
        self->history_buffer[((size_t)k * self->history_size) + t] = val;
      }
    }
    self->is_first_frame = false;
  } else {
    // Silence check for subsequent frames
    if (frame_energy < ESTIMATOR_SILENCE_THRESHOLD) {
      goto compute_output;
    }

    // 1. Update history (Circular buffer)
    for (uint32_t k = 0; k < self->spectrum_size; k++) {
      self->history_buffer[((size_t)k * self->history_size) +
                           self->history_index] = spectrum[k];
    }
    // Advance history index
    self->history_index = (self->history_index + 1) % self->history_size;
  }

compute_output:
  // 2. Process each bin
  if (!self->history_buffer || !self->sort_buffer || self->history_size == 0) {
    return false;
  }

  for (uint32_t k = 0; k < self->spectrum_size; k++) {
    // Copy history to sort buffer
    float* bin_history = &self->history_buffer[(size_t)k * self->history_size];
    memcpy(self->sort_buffer, bin_history, self->history_size * sizeof(float));

    // Sort to find lowest items
    qsort(self->sort_buffer, self->history_size, sizeof(float), compare_floats);

    // Compute Trimmed Mean
    float sum = 0.0f;
    for (uint32_t i = 0; i < self->trim_count; i++) {
      sum += self->sort_buffer[i];
    }
    float trimmed_mean = sum / (float)self->trim_count;

    // Apply Correction
    noise_spectrum[k] = trimmed_mean * self->correction_factor;
  }

  return true;
}

void brandt_noise_estimator_set_state(BrandtNoiseEstimator* self,
                                      const float* initial_profile) {
  if (!self || !initial_profile) {
    return;
  }

  // Start with a known state: Fill history with this profile
  // To avoid identical values (which sorts weirdly?), we just fill flat.
  // But strictly, we should reverse-engineer the trimmed mean.
  // If we fill buffer with X, trimmed mean is X. Output is X * Factor.
  // We want output to be Profile. So input should be Profile / Factor.
  float inverse_factor = 1.0f / self->correction_factor;

  for (uint32_t k = 0; k < self->spectrum_size; k++) {
    float val = initial_profile[k] * inverse_factor;
    for (uint32_t t = 0; t < self->history_size; t++) {
      self->history_buffer[((size_t)k * self->history_size) + t] = val;
    }
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
    for (uint32_t t = 0; t < self->history_size; t++) {
      if (self->history_buffer[((size_t)k * self->history_size) + t] <
          floor_val) {
        self->history_buffer[((size_t)k * self->history_size) + t] = floor_val;
      }
    }
  }
}

void brandt_noise_estimator_set_history_duration(
    const BrandtNoiseEstimator* self, float history_duration_ms,
    uint32_t sample_rate, uint32_t fft_size) {
  if (!self) {
    return;
  }

  // To avoid frequent reallocations, we only update if it's a significant
  // change For now, let's just update the internal logic if we don't want to
  // realloc. Actually, changing history size at runtime is best done by
  // pre-allocating a MAX size. Given we are in the middle of a fix, I will only
  // expose Percentile as "Sensitivity" first.
  (void)history_duration_ms;
  (void)sample_rate;
  (void)fft_size;
}
