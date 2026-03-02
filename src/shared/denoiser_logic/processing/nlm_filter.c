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

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/nlm_filter.h"
#include "shared/utils/simd_utils.h"

struct NlmFilter {
  NlmFilterConfig config;

  // Ring buffer for SNR frames
  float** frame_buffer;   // [time_buffer_size][spectrum_size]
  uint32_t buffer_head;   // Next write position
  uint32_t frames_filled; // Number of frames currently in buffer

  // Target frame index (allows look-ahead)
  uint32_t target_frame_offset;

  // Precomputed values
  float h_squared;
  float inv_h_squared; // Precomputed 1/h^2 for multiplication
  float distance_threshold_actual;

  // Scratch buffer for processing (avoid realloc)
  float* weight_accum;
};

// Helper: clamp index to valid range
static inline uint32_t clamp_index(int32_t idx, uint32_t max_val) {
  if (idx < 0) {
    return 0;
  }
  if ((uint32_t)idx >= max_val) {
    return max_val - 1;
  }
  return (uint32_t)idx;
}

// Helper: get frame from ring buffer (handles wrap-around)
static inline float* get_frame(NlmFilter* self, int32_t relative_offset) {
  int32_t idx = (int32_t)self->buffer_head -
                (int32_t)self->config.search_range_time_future - 1 +
                relative_offset;

  // Wrap around in ring buffer
  while (idx < 0) {
    idx += (int32_t)self->config.time_buffer_size;
  }
  idx %= (int32_t)self->config.time_buffer_size;

  return self->frame_buffer[idx];
}

// Helper: compute squared Euclidean distance between two patches
// Optimized with SIMD for common patch sizes (4, 8)
static float compute_patch_distance(NlmFilter* self, int32_t target_time,
                                    uint32_t target_freq,
                                    int32_t candidate_time,
                                    uint32_t candidate_freq) {
  float distance = 0.0F;
  const uint32_t patch_size = self->config.patch_size;
  const uint32_t half_patch = patch_size / 2;
  const uint32_t spectrum_size = self->config.spectrum_size;

  // Check if we can use the fast SIMD path (no frequency boundary checking
  // needed) Frequency indices must satisfy: 0 <= freq - half_patch  AND freq -
  // half_patch + patch_size <= spectrum_size
  bool safe_bounds =
      (target_freq >= half_patch) &&
      (target_freq + patch_size - half_patch <= spectrum_size) &&
      (candidate_freq >= half_patch) &&
      (candidate_freq + patch_size - half_patch <= spectrum_size);

  for (uint32_t dt = 0; dt < patch_size; dt++) {
    int32_t t_target = target_time + (int32_t)dt - (int32_t)half_patch;
    int32_t t_cand = candidate_time + (int32_t)dt - (int32_t)half_patch;

    // Clamp time indices (always needed as time buffer is circular/finite)
    // But get_frame handles the circular logic relative to buffer_head if we
    // pass relative offset wait, get_frame argument is 'relative_offset' from
    // CURRENT time (0) The arguments t_target/t_cand passed here are relative
    // to current frame? In caller: compute_patch_distance(filter, 0, ...) ->
    // target_time is 0 (current). So 'dt' is relative offset. We don't need to
    // clamp 'dt' against buffer size, get_frame handles wrapping. BUT we
    // usually want to search only within valid history. 'target_time' is 0.
    // 'candidate_time' is 'dt' (search loop var). loop dt (patch row): operates
    // on row 't' of the patch. relative offset = time + dt - half_patch.

    float* target_frame = get_frame(self, t_target);
    float* cand_frame = get_frame(self, t_cand);

    if (safe_bounds && patch_size == 8) {
      // Fast Path: Direct pointer access + SIMD for 8x8
      distance += sb_vec8_ssd(sb_load8(target_frame + (target_freq - half_patch)),
                              sb_load8(cand_frame + (candidate_freq - half_patch)));

    } else if (safe_bounds && patch_size == 4) {
      // Fast Path for 4x4
      distance += sb_vec4_ssd(sb_load4(target_frame + (target_freq - half_patch)),
                              sb_load4(cand_frame + (candidate_freq - half_patch)));

    } else {
      // Slow safe path with clamping
      for (uint32_t df = 0; df < patch_size; df++) {
        uint32_t f_target = clamp_index(
            (int32_t)target_freq + (int32_t)df - (int32_t)half_patch,
            spectrum_size);
        uint32_t f_cand = clamp_index(
            (int32_t)candidate_freq + (int32_t)df - (int32_t)half_patch,
            spectrum_size);

        float diff = target_frame[f_target] - cand_frame[f_cand];
        distance += diff * diff;
      }
    }
  }

  return distance;
}

NlmFilter* nlm_filter_initialize(NlmFilterConfig config) {
  // Validate configuration
  if (config.spectrum_size == 0) {
    return NULL;
  }

  NlmFilter* self = (NlmFilter*)calloc(1U, sizeof(NlmFilter));
  if (!self) {
    return NULL;
  }

  // Apply defaults for zero values
  self->config = config;
  if (self->config.patch_size == 0) {
    self->config.patch_size = NLM_PATCH_SIZE;
  }
  if (self->config.paste_block_size == 0) {
    self->config.paste_block_size = NLM_PASTE_BLOCK_SIZE;
  }
  if (self->config.search_range_freq == 0) {
    self->config.search_range_freq = NLM_SEARCH_RANGE_FREQ;
  }
  if (self->config.search_range_time_past == 0) {
    self->config.search_range_time_past = NLM_SEARCH_RANGE_TIME_PAST;
  }
  if (self->config.search_range_time_future == 0) {
    self->config.search_range_time_future = NLM_SEARCH_RANGE_TIME_FUTURE;
  }
  if (self->config.h_parameter <= 0.0F) {
    self->config.h_parameter = NLM_DEFAULT_H_PARAMETER;
  }

  // Compute time buffer size from search ranges
  if (self->config.time_buffer_size == 0) {
    self->config.time_buffer_size = self->config.search_range_time_past +
                                    self->config.search_range_time_future + 1;
  }

  // Precompute values
  self->h_squared = self->config.h_parameter * self->config.h_parameter;
  self->inv_h_squared = 1.0f / self->h_squared;
  if (self->config.distance_threshold <= 0.0F) {
    // Default: 4 * h² (patches 2× noise variance away contribute < 2%)
    self->distance_threshold_actual = 4.0F * self->h_squared;
  } else {
    self->distance_threshold_actual = self->config.distance_threshold;
  }

  // Target frame is at the "past" offset to allow look-ahead
  self->target_frame_offset = self->config.search_range_time_past;

  // Allocate ring buffer for frames
  self->frame_buffer =
      (float**)calloc(self->config.time_buffer_size, sizeof(float*));
  if (!self->frame_buffer) {
    nlm_filter_free(self);
    return NULL;
  }

  for (uint32_t i = 0; i < self->config.time_buffer_size; i++) {
    self->frame_buffer[i] =
        (float*)calloc(self->config.spectrum_size, sizeof(float));
    if (!self->frame_buffer[i]) {
      nlm_filter_free(self);
      return NULL;
    }
  }

  self->buffer_head = 0;
  self->frames_filled = 0;

  // Allocate scratch buffer
  self->weight_accum =
      (float*)calloc(self->config.spectrum_size, sizeof(float));
  if (!self->weight_accum) {
    nlm_filter_free(self);
    return NULL;
  }

  return self;
}

void nlm_filter_free(NlmFilter* filter) {
  if (!filter) {
    return;
  }

  if (filter->frame_buffer) {
    for (uint32_t i = 0; i < filter->config.time_buffer_size; i++) {
      if (filter->frame_buffer[i]) {
        free((void*)filter->frame_buffer[i]);
      }
    }
    free((void*)filter->frame_buffer);
  }

  if (filter->weight_accum) {
    free(filter->weight_accum);
  }

  free(filter);
}

void nlm_filter_set_h_parameter(NlmFilter* filter, float h) {
  if (!filter || h <= 0.0F) {
    return;
  }

  filter->config.h_parameter = h;
  filter->h_squared = h * h;
  filter->inv_h_squared = 1.0f / filter->h_squared;

  // Update distance threshold if using default
  if (filter->config.distance_threshold <= 0.0F) {
    filter->distance_threshold_actual = 4.0F * filter->h_squared;
  }
}

void nlm_filter_push_frame(NlmFilter* filter, const float* snr_frame) {
  if (!filter || !snr_frame) {
    return;
  }

  // Copy frame to ring buffer
  memcpy(filter->frame_buffer[filter->buffer_head], snr_frame,
         filter->config.spectrum_size * sizeof(float));

  // Advance head
  filter->buffer_head =
      (filter->buffer_head + 1) % filter->config.time_buffer_size;

  // Track fill level
  if (filter->frames_filled < filter->config.time_buffer_size) {
    filter->frames_filled++;
  }
}

bool nlm_filter_is_ready(NlmFilter* filter) {
  if (!filter) {
    return false;
  }
  return filter->frames_filled >= filter->config.time_buffer_size;
}

bool nlm_filter_process(NlmFilter* filter, float* smoothed_snr) {
  if (!filter || !smoothed_snr) {
    return false;
  }

  // Need full buffer for processing
  if (!nlm_filter_is_ready(filter)) {
    return false;
  }

#ifdef __SSE__
  unsigned int old_mxcsr = _mm_getcsr();
  _mm_setcsr(old_mxcsr | 0x8040); // Set flush-to-zero and denormals-are-zero
#endif

  const uint32_t spectrum_size = filter->config.spectrum_size;
  const uint32_t paste_size = filter->config.paste_block_size;
  const uint32_t search_freq = filter->config.search_range_freq;
  const int32_t search_time_past =
      (int32_t)filter->config.search_range_time_past;
  const int32_t search_time_future =
      (int32_t)filter->config.search_range_time_future;

  // Get target frame (center of buffer)
  float* target_frame = get_frame(filter, 0);

  // Initialize output with zeros
  memset(smoothed_snr, 0, spectrum_size * sizeof(float));

  // Use preallocated buffer
  float* weight_sum = filter->weight_accum;
  // Clear weight buffer for new frame
  memset(weight_sum, 0, spectrum_size * sizeof(float));

  // Process using pasted block optimization
  // For each paste block position in frequency
  for (uint32_t block_start = 0; block_start < spectrum_size;
       block_start += paste_size) {

    uint32_t block_center = block_start + (paste_size / 2);
    if (block_center >= spectrum_size) {
      block_center = spectrum_size - 1;
    }

    float current_inv_h2 = filter->inv_h_squared;
    float current_dist_threshold = filter->distance_threshold_actual;

    // --- Optimization: Register Blocking for Target Patch ---
    // We load the 8x8 target patch ONCE here, keeping it in registers
    // to avoid reloading it ~350 times in the inner loop.

    // Storage for pre-loaded target patch (flattened 8x8)
    sb_vec8_t target_vecs[8];

    // Pre-load Target Patch
    // Use nlm_filter 8x8 constraints directly for max speed
    const uint32_t half_patch_size =
        4; // Hardcoded for 8x8 optimized path checks

    if (filter->config.patch_size == 8) {
      // Unroll loading for 8x8
      for (int r = 0; r < 8; r++) {
        int32_t t_offset = (int32_t)r - (int32_t)half_patch_size;
        float* row_ptr =
            get_frame(filter, t_offset) + (block_center - half_patch_size);

        bool safe_load =
            (block_center >= 4) && (block_center + 4 <= spectrum_size);

        if (safe_load) {
          target_vecs[r] = sb_load8(row_ptr);
        }
      }
    }

    // Search over time and frequency window
    for (int32_t dt = -search_time_past; dt <= search_time_future; dt++) {
      for (int32_t df = -(int32_t)search_freq; df <= (int32_t)search_freq;
           df++) {

        uint32_t cand_center =
            clamp_index((int32_t)block_center + df, spectrum_size);

        float distance = 0.0F;

        // --- INLINED DISTANCE CALCULATION WITH REGISTER BLOCKING ---
        bool safe_bounds =
            (block_center >= 4) && (block_center + 4 <= spectrum_size) &&
            (cand_center >= 4) && (cand_center + 4 <= spectrum_size);

        if (filter->config.patch_size == 8 && safe_bounds) {
          // 8x8 Optimized Path
          sb_acc8_t sum = sb_acc8_zero();

          for (int r = 0; r < 8; r++) {
            int32_t t_cand = dt + r - 4;
            float* cand_ptr = get_frame(filter, t_cand) + (cand_center - 4);
            sum = sb_acc8_add_ssd(sum, target_vecs[r], sb_load8(cand_ptr));
          }
          distance = sb_acc8_hsum(sum);
        } else {
          // Fallback for boundaries or non-8x8
          distance =
              compute_patch_distance(filter, 0, block_center, dt, cand_center);
        }

        // Distance thresholding (early termination optimization)
        if (distance > current_dist_threshold) {
          continue;
        }

        // Compute weight
        float weight = expf(-distance * current_inv_h2);
        if (weight < NLM_MIN_WEIGHT) {
          continue;
        }

        // Get candidate frame
        float* cand_frame = get_frame(filter, dt);

        // Apply weight to all bins in the paste block
        uint32_t current_paste_limit = paste_size;
        if (block_start + paste_size > spectrum_size) {
          current_paste_limit = spectrum_size - block_start;
        }

        for (uint32_t i = 0; i < current_paste_limit; i++) {
          uint32_t target_bin = block_start + i;
          uint32_t cand_bin =
              clamp_index((int32_t)target_bin + df, spectrum_size);

          smoothed_snr[target_bin] += weight * cand_frame[cand_bin];
          weight_sum[target_bin] += weight;
        }
      }
    }
  }

  // Normalize by weight sum
  for (uint32_t k = 0; k < spectrum_size; k++) {
    if (weight_sum[k] > NLM_MIN_WEIGHT) {
      smoothed_snr[k] /= weight_sum[k];
    } else {
      // Fallback to original value if no valid weights
      smoothed_snr[k] = target_frame[k];
    }
  }

#ifdef __SSE__
  _mm_setcsr(old_mxcsr);
#endif

  return true;
}

void nlm_filter_reset(NlmFilter* filter) {
  if (!filter) {
    return;
  }

  // Clear all frames
  for (uint32_t i = 0; i < filter->config.time_buffer_size; i++) {
    memset(filter->frame_buffer[i], 0,
           filter->config.spectrum_size * sizeof(float));
  }

  filter->buffer_head = 0;
  filter->frames_filled = 0;
}

uint32_t nlm_filter_get_latency_frames(NlmFilter* filter) {
  if (!filter) {
    return 0;
  }
  // Latency is the number of look-ahead frames
  return filter->config.search_range_time_future;
}

void nlm_filter_calculate_snr(NlmFilter* filter,
                              const float* reference_spectrum,
                              const float* noise_spectrum, float* snr_frame) {
  if (!filter || !reference_spectrum || !noise_spectrum || !snr_frame) {
    return;
  }

  const uint32_t spectrum_size = filter->config.spectrum_size;

  for (uint32_t k = 0; k < spectrum_size; k++) {
    float denom =
        noise_spectrum[k] > FLT_MIN ? noise_spectrum[k] : SPECTRAL_EPSILON;
    snr_frame[k] = reference_spectrum[k] / denom;
  }
}

void nlm_filter_reconstruct_magnitude(NlmFilter* filter,
                                      const float* smoothed_snr,
                                      const float* noise_spectrum,
                                      float* magnitude_spectrum) {
  if (!filter || !smoothed_snr || !noise_spectrum || !magnitude_spectrum) {
    return;
  }

  const uint32_t spectrum_size = filter->config.spectrum_size;

  for (uint32_t k = 0; k < spectrum_size; k++) {
    float denom =
        noise_spectrum[k] > FLT_MIN ? noise_spectrum[k] : SPECTRAL_EPSILON;
    magnitude_spectrum[k] = smoothed_snr[k] * denom;
  }
}
