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

#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/nlm_filter.h"
#include "shared/denoiser_logic/processing/patch_filter_context.h"
#include "shared/utils/general_utils.h"
#include "shared/utils/simd_utils.h"
#include "shared/utils/thread_pool.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef bool (*nlm_process_impl_fn)(NlmFilter* filter, float* smoothed_snr);

struct NlmFilter {
  NlmFilterConfig config;
  PatchFilterContext context;

  // Target frame index (allows look-ahead)
  uint32_t target_frame_offset;

  // Precomputed values
  float h_squared;
  float inv_h_squared; // Precomputed 1/h^2 for multiplication
  float distance_threshold_actual;

  // Scratch buffer for processing (avoid realloc)
  float* weight_accum;

  // Function pointer for runtime architecture dispatch
  nlm_process_impl_fn process_fn;
};

// Chunked single-row SSD over n contiguous bins (8-wide + 4-wide + scalar
// tail). Safe (in-bounds) segments only; edges use the clamped scalar path.
static inline SB_UNUSED float sb_row_ssd_n(const float* a, const float* b,
                                           uint32_t n) {
  float ssd = 0.0F;
  uint32_t i = 0U;
  if (n >= 8U) {
    sb_acc8_t acc = sb_acc8_zero();
    for (; i + 8U <= n; i += 8U) {
      acc = sb_acc8_add_ssd(acc, sb_load8(a + i), sb_load8(b + i));
    }
    ssd += sb_acc8_hsum(acc);
  }
  if (i + 4U <= n) {
    ssd += sb_vec4_ssd(sb_load4(a + i), sb_load4(b + i));
    i += 4U;
  }
  for (; i < n; i++) {
    float diff = a[i] - b[i];
    ssd += diff * diff;
  }
  return ssd;
}

// Full-patch SSD between preloaded target rows and candidate frame rows.
// tgt_rows/cand_rows hold patch_size row pointers; each row spans
// patch_size contiguous bins. Chunked 8/4/1 per row; for patch_size == 4
// this matches the old per-row vec4 path bit-for-bit.
static inline SB_UNUSED float sb_patch_ssd_n(float* const* tgt_rows,
                                             float* const* cand_rows,
                                             uint32_t patch_size) {
  float ssd = 0.0F;
  for (uint32_t r = 0U; r < patch_size; r++) {
    ssd += sb_row_ssd_n(tgt_rows[r], cand_rows[r], patch_size);
  }
  return ssd;
}

// Helper: compute squared Euclidean distance between two patches
static inline SB_UNUSED float compute_patch_distance(NlmFilter* self,
                                                     int32_t target_time,
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

  if (safe_bounds && patch_size == 8) {
    // Legacy fast path (kept bit-identical for the default 46ms option).
    for (uint32_t dt = 0; dt < patch_size; dt++) {
      int32_t t_target = target_time + (int32_t)dt - (int32_t)half_patch;
      int32_t t_cand = candidate_time + (int32_t)dt - (int32_t)half_patch;

      const float* target_frame =
          patch_filter_context_get_frame(&self->context, t_target);
      const float* cand_frame =
          patch_filter_context_get_frame(&self->context, t_cand);
      distance +=
          sb_vec8_ssd(sb_load8(target_frame + (target_freq - half_patch)),
                      sb_load8(cand_frame + (candidate_freq - half_patch)));
    }
    return distance;
  }

  if (safe_bounds) {
    for (uint32_t dt = 0; dt < patch_size; dt++) {
      int32_t t_target = target_time + (int32_t)dt - (int32_t)half_patch;
      int32_t t_cand = candidate_time + (int32_t)dt - (int32_t)half_patch;

      const float* target_frame =
          patch_filter_context_get_frame(&self->context, t_target);
      const float* cand_frame =
          patch_filter_context_get_frame(&self->context, t_cand);
      distance +=
          sb_row_ssd_n(target_frame + (target_freq - half_patch),
                       cand_frame + (candidate_freq - half_patch), patch_size);
    }
    return distance;
  }

  for (uint32_t dt = 0; dt < patch_size; dt++) {
    int32_t t_target = target_time + (int32_t)dt - (int32_t)half_patch;
    int32_t t_cand = candidate_time + (int32_t)dt - (int32_t)half_patch;

    const float* target_frame =
        patch_filter_context_get_frame(&self->context, t_target);
    const float* cand_frame =
        patch_filter_context_get_frame(&self->context, t_cand);

    distance = patch_filter_accumulate_patch_row_ssd(
        distance, target_frame, cand_frame, target_freq, candidate_freq,
        patch_size, half_patch, spectrum_size);
  }

  return distance;
}

// Context for the block-smoothing loop. Every block accumulates into a
// disjoint range of smoothed_snr/weight_sum, so contiguous static partitioning
// is race-free and produces the same per-bin accumulation order (and therefore
// the same output) as sequential execution.
typedef struct {
  NlmFilter* filter;
  float* smoothed_snr;
  float* weight_sum;
  float* target_frame;
} NlmBlockTask;

// Per-block working state for nlm_process_block_range. Groups the scalars
// and preloaded target data shared by the stage helpers so each helper stays
// under the parameter-count limit. All arrays live on the caller's stack
// (same footprint as the previous inline version); no allocation.
typedef struct {
  NlmFilter* filter;
  float* smoothed_snr;
  float* weight_sum;
  uint32_t block_start;
  uint32_t current_paste_limit;
  uint32_t block_center;
  uint32_t patch_size;
  uint32_t half_patch_size;
  uint32_t spectrum_size;
  bool safe_block;
  sb_vec8_t target_vecs[8];
  // Flat preloaded target patch for the generalized (non-8) path:
  // patch_size rows x NLM_MAX_PATCH_FRAMES stride, 1KB max on stack.
  float target_patch[NLM_MAX_PATCH_FRAMES * NLM_MAX_PATCH_FRAMES];
  float* tgt_rows[NLM_MAX_PATCH_FRAMES];
} NlmBlockCtx;

static inline SB_UNUSED bool nlm_preload_target(NlmBlockCtx* b) {
  NlmFilter* filter = b->filter;
  const uint32_t block_center = b->block_center;
  const uint32_t patch_size = b->patch_size;
  const uint32_t half_patch_size = b->half_patch_size;
  const uint32_t spectrum_size = b->spectrum_size;

  // Patch rows span [center - half, center + (patch - half) - 1]; for odd
  // patch sizes the upper reach is half + 1 bins, so the bound must use
  // (patch_size - half_patch_size), not half_patch_size.
  b->safe_block =
      (block_center >= half_patch_size) &&
      (block_center + (patch_size - half_patch_size) <= spectrum_size);

  if (patch_size == 8) {
    for (int r = 0; r < 8; r++) {
      if (b->safe_block) {
        int32_t t_offset = (int32_t)r - (int32_t)half_patch_size;
        const float* row_ptr =
            patch_filter_context_cached_get_frame(&filter->context, t_offset) +
            (block_center - half_patch_size);
        b->target_vecs[r] = sb_load8(row_ptr);
      } else {
        b->target_vecs[r] = sb_set8(0.0f);
      }
    }
    return b->safe_block;
  }
  if (b->safe_block && patch_size <= NLM_MAX_PATCH_FRAMES) {
    for (uint32_t r = 0; r < patch_size; r++) {
      int32_t t_offset = (int32_t)r - (int32_t)half_patch_size;
      const float* row_ptr =
          patch_filter_context_cached_get_frame(&filter->context, t_offset) +
          (block_center - half_patch_size);
      float* dst = &b->target_patch[(size_t)r * NLM_MAX_PATCH_FRAMES];
      memcpy(dst, row_ptr, patch_size * sizeof(float));
      b->tgt_rows[r] = dst;
    }
  }
  return b->safe_block;
}

static inline SB_UNUSED void nlm_load_candidate_rows(NlmBlockCtx* b, int32_t dt,
                                                     float** cand_rows,
                                                     float** cand_rows_n) {
  NlmFilter* filter = b->filter;
  const uint32_t patch_size = b->patch_size;
  const uint32_t half_patch_size = b->half_patch_size;
  if (patch_size == 8) {
    for (int r = 0; r < 8; r++) {
      cand_rows[r] =
          patch_filter_context_cached_get_frame(&filter->context, dt + r - 4);
    }
    return;
  }
  if (patch_size <= NLM_MAX_PATCH_FRAMES) {
    for (uint32_t r = 0; r < patch_size; r++) {
      cand_rows_n[r] = patch_filter_context_cached_get_frame(
          &filter->context, dt + (int32_t)r - (int32_t)half_patch_size);
    }
  }
}

static inline SB_UNUSED float nlm_dispatch_distance(NlmBlockCtx* b, int32_t dt,
                                                    uint32_t cand_center,
                                                    float** cand_rows,
                                                    float** cand_rows_n) {
  NlmFilter* filter = b->filter;
  const uint32_t block_center = b->block_center;
  const uint32_t patch_size = b->patch_size;
  const uint32_t half_patch_size = b->half_patch_size;
  const uint32_t spectrum_size = b->spectrum_size;

  bool safe_bounds =
      b->safe_block && (cand_center >= half_patch_size) &&
      (cand_center + (patch_size - half_patch_size) <= spectrum_size);

  if (patch_size == 8 && safe_bounds && cand_rows[0]) {
    uint32_t cand_f_start = cand_center - 4;
    float* cand_row_ptrs[8] = {
        cand_rows[0] + cand_f_start, cand_rows[1] + cand_f_start,
        cand_rows[2] + cand_f_start, cand_rows[3] + cand_f_start,
        cand_rows[4] + cand_f_start, cand_rows[5] + cand_f_start,
        cand_rows[6] + cand_f_start, cand_rows[7] + cand_f_start,
    };
    return sb_vec8_patch_ssd(b->target_vecs, cand_row_ptrs);
  }
  if (patch_size != 8 && safe_bounds && b->safe_block &&
      patch_size <= NLM_MAX_PATCH_FRAMES && cand_rows_n[0]) {
    const uint32_t cand_f_start = cand_center - half_patch_size;
    float* cand_ptrs[NLM_MAX_PATCH_FRAMES];
    float* tgt_ptrs[NLM_MAX_PATCH_FRAMES];
    for (uint32_t r = 0; r < patch_size; r++) {
      cand_ptrs[r] = cand_rows_n[r] + cand_f_start;
      tgt_ptrs[r] = b->tgt_rows[r];
    }
    return sb_patch_ssd_n(tgt_ptrs, cand_ptrs, patch_size);
  }
  return compute_patch_distance(filter, 0, block_center, dt, cand_center);
}

static inline SB_UNUSED void nlm_accumulate_paste(NlmBlockCtx* b, int32_t df,
                                                  float weight,
                                                  const float* cand_frame) {
  for (uint32_t i = 0; i < b->current_paste_limit; i++) {
    uint32_t target_bin = b->block_start + i;
    uint32_t cand_bin =
        patch_filter_clamp_index((int32_t)target_bin + df, b->spectrum_size);

    b->smoothed_snr[target_bin] += weight * cand_frame[cand_bin];
    b->weight_sum[target_bin] += weight;
  }
}

static inline SB_UNUSED void nlm_process_block_range(void* raw_ctx,
                                                     uint32_t first_block,
                                                     uint32_t block_count) {
  const NlmBlockTask* ctx = (const NlmBlockTask*)raw_ctx;
  NlmFilter* filter = ctx->filter;
  float* smoothed_snr = ctx->smoothed_snr;
  float* weight_sum = ctx->weight_sum;
  const float* target_frame = ctx->target_frame;

  const uint32_t spectrum_size = filter->config.spectrum_size;
  const uint32_t paste_size = filter->config.paste_block_size;
  const uint32_t search_freq = filter->config.search_range_freq;
  const int32_t search_time_past =
      (int32_t)filter->config.search_range_time_past;
  const int32_t search_time_future =
      (int32_t)filter->config.search_range_time_future;
  const float current_inv_h2 = filter->inv_h_squared;
  const float current_dist_threshold = filter->distance_threshold_actual;
  const uint32_t patch_size = filter->config.patch_size;
  const uint32_t half_patch_size = patch_size / 2;

  for (uint32_t block = 0; block < block_count; block++) {
    const uint32_t block_start = (first_block + block) * paste_size;

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

    NlmBlockCtx b = {
        .filter = filter,
        .smoothed_snr = smoothed_snr,
        .weight_sum = weight_sum,
        .block_start = block_start,
        .current_paste_limit = current_paste_limit,
        .block_center = block_center,
        .patch_size = patch_size,
        .half_patch_size = half_patch_size,
        .spectrum_size = spectrum_size,
        .safe_block = false,
    };
    nlm_preload_target(&b);

    for (int32_t dt = -search_time_past; dt <= search_time_future; dt++) {
      float* cand_rows[8] = {NULL};
      // Candidate row pointers for the generalized path (any patch <= 16).
      float* cand_rows_n[NLM_MAX_PATCH_FRAMES] = {NULL};
      nlm_load_candidate_rows(&b, dt, cand_rows, cand_rows_n);

      for (int32_t df = -(int32_t)search_freq; df <= (int32_t)search_freq;
           df++) {

        uint32_t cand_center =
            patch_filter_clamp_index((int32_t)block_center + df, spectrum_size);

        float distance =
            nlm_dispatch_distance(&b, dt, cand_center, cand_rows, cand_rows_n);

        if (distance > current_dist_threshold) {
          continue;
        }

        float weight = sb_fast_expf(-distance * current_inv_h2);
        if (weight < NLM_MIN_WEIGHT) {
          continue;
        }

        const float* cand_frame =
            patch_filter_context_cached_get_frame(&filter->context, dt);

        nlm_accumulate_paste(&b, df, weight, cand_frame);
      }
    }
  }
}

static inline SB_UNUSED bool nlm_filter_process_core(NlmFilter* filter,
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

  patch_filter_context_populate_frame_ptrs(&filter->context);

  float* target_frame =
      patch_filter_context_cached_get_frame(&filter->context, 0);

  if (filter->config.h_parameter <= 0.0F || filter->h_squared <= 0.0F) {
    memcpy(smoothed_snr, target_frame, spectrum_size * sizeof(float));
    sb_simd_restore_state(old_simd_state);
    return true;
  }

  memset(smoothed_snr, 0, spectrum_size * sizeof(float));

  float* weight_sum = filter->weight_accum;
  memset(weight_sum, 0, spectrum_size * sizeof(float));

  NlmBlockTask task_ctx = {filter, smoothed_snr, weight_sum, target_frame};
  const uint32_t num_blocks = (spectrum_size + paste_size - 1) / paste_size;

  if (filter->context.pool) {
    sb_thread_pool_parallel_for(filter->context.pool, num_blocks,
                                nlm_process_block_range, &task_ctx);
  } else {
    nlm_process_block_range(&task_ctx, 0, num_blocks);
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
