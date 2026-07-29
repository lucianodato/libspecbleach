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

#ifndef NLM_FILTER_INTERNAL_H
#define NLM_FILTER_INTERNAL_H

#include "shared/denoiser_logic/processing/nlm_filter.h"
#include <stdint.h>

typedef bool (*nlm_process_impl_fn)(NlmFilter* filter, float* smoothed_snr);

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

  // Pre-computed frame pointer cache (populated once per process call)
  // Indexed by [time_past + 4 + dt] where dt ranges from
  // -search_range_time_past - 4 to +search_range_time_future + 4
  float** frame_ptrs;
  uint32_t total_time_span; // search_time_past + search_time_future + 1 + 8

  // Scratch buffer for processing (avoid realloc)
  float* weight_accum;

  // Function pointer for runtime architecture dispatch
  nlm_process_impl_fn process_fn;

  // Number of threads for parallel processing
  uint32_t num_threads;
};

// Helper: clamp index to valid range
static inline __attribute__((unused)) uint32_t clamp_index(int32_t idx,
                                                           uint32_t max_val) {
  if (idx < 0) {
    return 0;
  }
  if ((uint32_t)idx >= max_val) {
    return max_val - 1;
  }
  return (uint32_t)idx;
}

// Helper: get frame from ring buffer (handles wrap-around)
static inline __attribute__((unused)) float* get_frame(
    NlmFilter* self, int32_t relative_offset) {
  const int32_t size = (int32_t)self->config.time_buffer_size;
  // Center is at (head - future - 1)
  int32_t idx = (int32_t)self->buffer_head -
                (int32_t)self->config.search_range_time_future - 1 +
                relative_offset;

  // Standard mathematical modulo for negative numbers
  idx = ((idx % size) + size) % size;

  return self->frame_buffer[idx];
}

// Pre-compute all frame pointers for the current processing window.
// After this call, frame_ptrs[search_time_past + 4 + dt] gives the frame
// at relative offset dt (where dt ranges from -past-4 to +future+4).
static inline __attribute__((unused)) void populate_frame_ptrs(
    NlmFilter* self) {
  const int32_t past = (int32_t)self->config.search_range_time_past;
  const int32_t future = (int32_t)self->config.search_range_time_future;

  // We need a halo of 4 frames for the 8x8 patch comparison
  for (int32_t dt = -past - 4; dt <= future + 4; dt++) {
    self->frame_ptrs[past + 4 + dt] = get_frame(self, dt);
  }
}

// O(1) frame lookup using pre-computed pointer cache.
// dt ranges from -search_time_past - 4 to +search_time_future + 4.
static inline __attribute__((unused)) float* cached_get_frame(NlmFilter* self,
                                                              int32_t dt) {
  return self
      ->frame_ptrs[(int32_t)self->config.search_range_time_past + 4 + dt];
}

#include "shared/configurations.h"
#include "shared/utils/simd_utils.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Helper: compute squared Euclidean distance between two patches
static inline float compute_patch_distance(NlmFilter* self, int32_t target_time,
                                           uint32_t target_freq,
                                           int32_t candidate_time,
                                           uint32_t candidate_freq) {
  float distance = 0.0F;
  const uint32_t patch_size = self->config.patch_size;
  const uint32_t half_patch = patch_size / 2;
  const uint32_t spectrum_size = self->config.spectrum_size;

  bool safe_bounds =
      (target_freq >= half_patch) &&
      (target_freq + patch_size - half_patch <= spectrum_size) &&
      (candidate_freq >= half_patch) &&
      (candidate_freq + patch_size - half_patch <= spectrum_size);

  for (uint32_t dt = 0; dt < patch_size; dt++) {
    int32_t t_target = target_time + (int32_t)dt - (int32_t)half_patch;
    int32_t t_cand = candidate_time + (int32_t)dt - (int32_t)half_patch;

    const float* target_frame = get_frame(self, t_target);
    const float* cand_frame = get_frame(self, t_cand);

    if (safe_bounds && patch_size == 8) {
      distance +=
          sb_vec8_ssd(sb_load8(target_frame + (target_freq - half_patch)),
                      sb_load8(cand_frame + (candidate_freq - half_patch)));

    } else if (safe_bounds && patch_size == 4) {
      distance +=
          sb_vec4_ssd(sb_load4(target_frame + (target_freq - half_patch)),
                      sb_load4(cand_frame + (candidate_freq - half_patch)));

    } else {
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

static inline bool nlm_filter_process_core(NlmFilter* filter,
                                           float* smoothed_snr) {
  if (!filter || !smoothed_snr) {
    return false;
  }

  if (!nlm_filter_is_ready(filter)) {
    return false;
  }

  sb_simd_state_t old_simd_state = sb_simd_enable_ftz_daz();

  const uint32_t spectrum_size = filter->config.spectrum_size;
  const uint32_t paste_size = filter->config.paste_block_size;
  const uint32_t search_freq = filter->config.search_range_freq;
  const int32_t search_time_past =
      (int32_t)filter->config.search_range_time_past;
  const int32_t search_time_future =
      (int32_t)filter->config.search_range_time_future;

  populate_frame_ptrs(filter);

  float* target_frame = cached_get_frame(filter, 0);

  memset(smoothed_snr, 0, spectrum_size * sizeof(float));

  float* weight_sum = filter->weight_accum;
  memset(weight_sum, 0, spectrum_size * sizeof(float));

#if SB_HAS_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(filter->num_threads)
#endif
  for (uint32_t block_start = 0; block_start < spectrum_size;
       block_start += paste_size) {

    uint32_t block_center = block_start + (paste_size / 2);
    if (block_center >= spectrum_size) {
      block_center = spectrum_size - 1;
    }

    uint32_t current_paste_limit = paste_size;
    if (block_start + paste_size > spectrum_size) {
      current_paste_limit = spectrum_size - block_start;
    }

    float target_snr_sum = 0.0F;
    for (uint32_t i = 0; i < current_paste_limit; i++) {
      target_snr_sum += target_frame[block_start + i];
    }
    if (target_snr_sum < 1e-6F) {
      continue;
    }

    float current_inv_h2 = filter->inv_h_squared;
    float current_dist_threshold = filter->distance_threshold_actual;

    if (spectrum_size > 1) {
      float normalized_freq = (float)block_center / (float)(spectrum_size - 1);
      float freq_scale = 1.0F + NLM_FREQ_DEPENDENT_SMOOTHING_SCALE *
                                    normalized_freq * normalized_freq;
      float h_val = filter->config.h_parameter * freq_scale;
      float h_squared = h_val * h_val;
      current_inv_h2 = 1.0F / h_squared;
      if (filter->config.distance_threshold <= 0.0F) {
        current_dist_threshold = NLM_DISTANCE_THRESHOLD_MULTIPLIER * h_squared;
      } else {
        current_dist_threshold =
            filter->config.distance_threshold * freq_scale * freq_scale;
      }
    }

    sb_vec8_t target_vecs[8];
    const uint32_t half_patch_size = 4;

    bool safe_block = (block_center >= half_patch_size) &&
                      (block_center + half_patch_size <= spectrum_size);

    if (filter->config.patch_size == 8) {
      for (int r = 0; r < 8; r++) {
        if (safe_block) {
          int32_t t_offset = (int32_t)r - (int32_t)half_patch_size;
          float* row_ptr = cached_get_frame(filter, t_offset) +
                           (block_center - half_patch_size);
          target_vecs[r] = sb_load8(row_ptr);
        } else {
          target_vecs[r] = sb_set8(0.0f);
        }
      }
    }

    for (int32_t dt = -search_time_past; dt <= search_time_future; dt++) {
      float* cand_rows[8];
      if (filter->config.patch_size == 8) {
        for (int r = 0; r < 8; r++) {
          cand_rows[r] = cached_get_frame(filter, dt + r - 4);
        }
      }

      for (int32_t df = -(int32_t)search_freq; df <= (int32_t)search_freq;
           df++) {

        uint32_t cand_center =
            clamp_index((int32_t)block_center + df, spectrum_size);

        float distance = 0.0F;

        bool safe_bounds = safe_block && (cand_center >= half_patch_size) &&
                           (cand_center + half_patch_size <= spectrum_size);

        if (filter->config.patch_size == 8 && safe_bounds) {
          uint32_t cand_f_start = cand_center - 4;
          float* cand_row_ptrs[8] = {
              cand_rows[0] + cand_f_start, cand_rows[1] + cand_f_start,
              cand_rows[2] + cand_f_start, cand_rows[3] + cand_f_start,
              cand_rows[4] + cand_f_start, cand_rows[5] + cand_f_start,
              cand_rows[6] + cand_f_start, cand_rows[7] + cand_f_start,
          };
          distance = sb_vec8_patch_ssd(target_vecs, cand_row_ptrs);
        } else {
          distance =
              compute_patch_distance(filter, 0, block_center, dt, cand_center);
        }

        if (distance > current_dist_threshold) {
          continue;
        }

        float weight = sb_fast_expf(-distance * current_inv_h2);
        if (weight < NLM_MIN_WEIGHT) {
          continue;
        }

        float* cand_frame = cached_get_frame(filter, dt);

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

  for (uint32_t k = 0; k < spectrum_size; k++) {
    if (weight_sum[k] > NLM_MIN_WEIGHT) {
      smoothed_snr[k] /= weight_sum[k];
    } else {
      smoothed_snr[k] = target_frame[k];
    }
  }

  sb_simd_restore_state(old_simd_state);

  return true;
}

// Generic implementation (SSE/NEON/Scalar)
bool nlm_filter_process_generic(NlmFilter* filter, float* smoothed_snr);

// AVX implementation
bool nlm_filter_process_avx(NlmFilter* filter, float* smoothed_snr);

#endif /* NLM_FILTER_INTERNAL_H */
