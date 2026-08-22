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

#include "hpss_filter.h"
#include "shared/configurations.h"
#include "shared/utils/simd_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct HpssFilter {
  HpssConfig config;
  bool is_enabled;
  bool is_initialized;

  float* prev_h; // Previous frame harmonic magnitude estimate H_{t-1, k}
  float* h;      // Current frame harmonic magnitude estimate H_{t, k}
  float* p;      // Current frame percussive magnitude estimate P_{t, k}

  float percussive_ratio;
  float energy_ratio_threshold;
  float margin_factor;
  bool is_transient_detected;
};

HpssFilter* hpss_filter_initialize(HpssConfig config) {
  if (config.real_spectrum_size == 0U) {
    return NULL;
  }

  HpssFilter* self = (HpssFilter*)calloc(1U, sizeof(HpssFilter));
  if (!self) {
    return NULL;
  }

  self->config = config;
  self->is_enabled = true;
  self->is_initialized = false;
  self->energy_ratio_threshold = HPSS_TRANSIENT_ENERGY_RATIO_THRESHOLD;
  self->margin_factor = HPSS_TRANSIENT_MARGIN_FACTOR;
  self->percussive_ratio = 0.0f;
  self->is_transient_detected = false;

  self->prev_h = (float*)calloc(config.real_spectrum_size, sizeof(float));
  self->h = (float*)calloc(config.real_spectrum_size, sizeof(float));
  self->p = (float*)calloc(config.real_spectrum_size, sizeof(float));

  if (!self->prev_h || !self->h || !self->p) {
    hpss_filter_free(self);
    return NULL;
  }

  return self;
}

void hpss_filter_free(HpssFilter* self) {
  if (!self) {
    return;
  }

  if (self->prev_h) {
    free(self->prev_h);
  }
  if (self->h) {
    free(self->h);
  }
  if (self->p) {
    free(self->p);
  }
  free(self);
}

void hpss_filter_set_enabled(HpssFilter* self, bool enabled) {
  if (!self) {
    return;
  }
  self->is_enabled = enabled;
}

uint32_t hpss_filter_get_latency_frames(const HpssFilter* self) {
  (void)self;
  return 0U; // Sliding HPSS is strictly causal with zero lookahead frames
}

float hpss_filter_get_onset_ratio(const HpssFilter* self) {
  return self ? self->percussive_ratio : 0.0f;
}

bool hpss_filter_is_transient_detected(const HpssFilter* self) {
  return (self && self->is_enabled) ? self->is_transient_detected : false;
}

bool hpss_filter_process(HpssFilter* self, const float* current_magnitude,
                         const float* noise_floor, float* current_magnitude_out,
                         float* mask_harmonic_out, float* mask_percussive_out) {
  if (!self || !current_magnitude) {
    return false;
  }

  const uint32_t spectrum_size = self->config.real_spectrum_size;

  if (current_magnitude_out) {
    memcpy(current_magnitude_out, current_magnitude,
           spectrum_size * sizeof(float));
  }

  if (!self->is_enabled) {
    if (mask_harmonic_out) {
      for (uint32_t k = 0U; k < spectrum_size; ++k) {
        mask_harmonic_out[k] = 1.0f;
      }
    }
    if (mask_percussive_out) {
      memset(mask_percussive_out, 0, spectrum_size * sizeof(float));
    }
    self->percussive_ratio = 0.0f;
    self->is_transient_detected = false;
    return true;
  }

  if (!self->is_initialized) {
    for (uint32_t k = 0U; k < spectrum_size; ++k) {
      self->prev_h[k] = current_magnitude[k];
      self->h[k] = current_magnitude[k];
      self->p[k] = 0.0f;
    }
    self->is_initialized = true;
  }

  sb_simd_state_t simd_state = sb_simd_enable_ftz_daz();

  // 1. Initialize H and P for current frame
  for (uint32_t k = 0U; k < spectrum_size; ++k) {
    float mag = current_magnitude[k];
    self->h[k] = HPSS_SLIDING_SMOOTH_FACTOR * mag;
    self->p[k] = HPSS_SLIDING_SMOOTH_FACTOR * mag;
  }

  // 2. Sliding HPSS Iterations (Ono / Tachibana ISMIR 2008)
  // Harmonic continuity: H_{t,k} aligns with H_{t-1,k} (temporal continuity)
  // Percussive continuity: P_{t,k} aligns with (P_{t,k-1} + P_{t,k+1})/2
  // (frequency continuity)
  for (uint32_t iter = 0U; iter < HPSS_SLIDING_ITERATIONS; ++iter) {
    for (uint32_t k = 0U; k < spectrum_size; ++k) {
      float mag = current_magnitude[k];
      if (mag <= SPECTRAL_EPSILON) {
        self->h[k] = 0.0f;
        self->p[k] = 0.0f;
        continue;
      }

      // Harmonic reference from previous frame
      float h_ref = self->prev_h[k];

      // Percussive reference from adjacent frequency bins
      float p_prev = (k > 0U) ? self->p[k - 1U] : self->p[k];
      float p_next = (k + 1U < spectrum_size) ? self->p[k + 1U] : self->p[k];
      float p_ref = HPSS_SLIDING_SMOOTH_FACTOR * (p_prev + p_next);

      // Auxiliary function update
      float h_sq = h_ref * h_ref;
      float p_sq = p_ref * p_ref;
      float denom = h_sq + p_sq;

      float w_h = (denom > SPECTRAL_EPSILON) ? (h_sq / denom)
                                             : HPSS_SLIDING_SMOOTH_FACTOR;
      float w_p = 1.0f - w_h;

      self->h[k] = w_h * mag;
      self->p[k] = w_p * mag;
    }
  }

  // 3. Compute Soft Masks with Noise Floor Adaptive Gating & Excess Margin
  float total_mag_sq = 0.0f;
  float percussive_mag_sq = 0.0f;
  const float current_margin = self->margin_factor;

  for (uint32_t k = 0U; k < spectrum_size; ++k) {
    float h_val = self->h[k];
    float p_val = self->p[k];

    // Margin test against horizontal sustain
    float p_diff = p_val - (current_margin * h_val);

    // If noise floor is provided, ensure percussive energy also emerges above
    // noise floor
    if (noise_floor != NULL && p_diff > 0.0f) {
      float noise_level = noise_floor[k];
      if (p_val < HPSS_NOISE_FLOOR_GATE_MULTIPLIER * noise_level) {
        p_diff = 0.0f; // Sub-noise fluctuation: treat as stationary noise
      }
    }

    float p_excess =
        (k >= (uint32_t)HPSS_BASS_CUTOFF_BINS && p_diff > 0.0f) ? p_diff : 0.0f;

    float h_sq = h_val * h_val;
    float p_excess_sq = p_excess * p_excess;
    float sum_sq = h_sq + p_excess_sq;

    float w_p = 0.0f;
    float w_h = 1.0f;
    if (sum_sq > SPECTRAL_EPSILON) {
      float raw_w_p = p_excess_sq / sum_sq;
      // Apply smooth power compression to prevent pre-echo ghosting
      w_p = raw_w_p * raw_w_p;
      w_h = 1.0f - w_p;
    }

    if (mask_harmonic_out) {
      mask_harmonic_out[k] = w_h;
    }
    if (mask_percussive_out) {
      mask_percussive_out[k] = w_p;
    }

    float mag_sq = current_magnitude[k] * current_magnitude[k];
    if (k >= (uint32_t)HPSS_BASS_CUTOFF_BINS) {
      total_mag_sq += mag_sq;
      percussive_mag_sq += w_p * mag_sq;
    }

    // Update state for next frame with exponential temporal tracking
    self->prev_h[k] =
        HPSS_SLIDING_SMOOTH_FACTOR * (self->prev_h[k] + current_magnitude[k]);
  }

  sb_simd_restore_state(simd_state);

  self->percussive_ratio = (total_mag_sq > HPSS_RATIO_DENOMINATOR_FLOOR)
                               ? (percussive_mag_sq / total_mag_sq)
                               : 0.0f;

  self->is_transient_detected =
      (self->percussive_ratio >= self->energy_ratio_threshold);

  return true;
}
