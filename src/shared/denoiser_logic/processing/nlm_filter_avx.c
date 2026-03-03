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

// Compile this file only when AVX is enabled globally, or when dynamic dispatch
// compilation flags enable AVX specifically for this translation unit.

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/nlm_filter_internal.h"
#include "shared/utils/simd_utils.h"

static float compute_patch_distance_avx(NlmFilter* self, int32_t target_time,
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

    float* target_frame = get_frame(self, t_target);
    float* cand_frame = get_frame(self, t_cand);

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

bool nlm_filter_process_avx(NlmFilter* filter, float* smoothed_snr) {
  sb_simd_state_t old_simd_state = sb_simd_enable_ftz_daz();

  const uint32_t spectrum_size = filter->config.spectrum_size;
  const uint32_t paste_size = filter->config.paste_block_size;
  const uint32_t search_freq = filter->config.search_range_freq;
  const int32_t search_time_past =
      (int32_t)filter->config.search_range_time_past;
  const int32_t search_time_future =
      (int32_t)filter->config.search_range_time_future;

  float* target_frame = get_frame(filter, 0);

  memset(smoothed_snr, 0, spectrum_size * sizeof(float));

  float* weight_sum = filter->weight_accum;
  memset(weight_sum, 0, spectrum_size * sizeof(float));

  for (uint32_t block_start = 0; block_start < spectrum_size;
       block_start += paste_size) {

    uint32_t block_center = block_start + (paste_size / 2);
    if (block_center >= spectrum_size) {
      block_center = spectrum_size - 1;
    }

    float current_inv_h2 = filter->inv_h_squared;
    float current_dist_threshold = filter->distance_threshold_actual;

    sb_vec8_t target_vecs[8];

    const uint32_t half_patch_size =
        4; // Hardcoded for 8x8 optimized path checks

    bool safe_block = (block_center >= half_patch_size) &&
                      (block_center + half_patch_size <= spectrum_size);

    if (filter->config.patch_size == 8) {
      for (int r = 0; r < 8; r++) {
        if (safe_block) {
          int32_t t_offset = (int32_t)r - (int32_t)half_patch_size;
          float* row_ptr =
              get_frame(filter, t_offset) + (block_center - half_patch_size);
          target_vecs[r] = sb_load8(row_ptr);
        } else {
          target_vecs[r] = sb_set8(0.0f);
        }
      }
    }

    for (int32_t dt = -search_time_past; dt <= search_time_future; dt++) {
      for (int32_t df = -(int32_t)search_freq; df <= (int32_t)search_freq;
           df++) {

        uint32_t cand_center =
            clamp_index((int32_t)block_center + df, spectrum_size);

        float distance = 0.0F;

        bool safe_bounds = safe_block && (cand_center >= half_patch_size) &&
                           (cand_center + half_patch_size <= spectrum_size);

        if (filter->config.patch_size == 8 && safe_bounds) {
          float* cand_row_ptrs[8];
          for (int r = 0; r < 8; r++) {
            cand_row_ptrs[r] =
                get_frame(filter, dt + r - 4) + (cand_center - 4);
          }

          distance = sb_vec8_patch_ssd(target_vecs, cand_row_ptrs);
        } else {
          distance = compute_patch_distance_avx(filter, 0, block_center, dt,
                                                cand_center);
        }

        if (distance > current_dist_threshold) {
          continue;
        }

        float weight = expf(-distance * current_inv_h2);
        if (weight < NLM_MIN_WEIGHT) {
          continue;
        }

        float* cand_frame = get_frame(filter, dt);

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
