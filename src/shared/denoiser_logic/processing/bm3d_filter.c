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

#include "shared/denoiser_logic/processing/bm3d_filter.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/patch_filter_context.h"
#include "shared/utils/simd_utils.h"
#include "shared/utils/thread_pool.h"

struct Bm3dFilter {
  Bm3dFilterConfig config;
  PatchFilterContext context;
  float h_squared;
  float inv_h_squared;
};

// Scalar clamped patch SSD (edges / non-8 patch sizes).
static float bm3d_patch_ssd_scalar(Bm3dFilter* self, uint32_t target_freq,
                                   int32_t cand_dt, uint32_t cand_freq) {
  const uint32_t patch = self->config.patch_size;
  const uint32_t half = patch / 2;
  const uint32_t n = self->config.spectrum_size;
  float ssd = 0.0F;
  for (uint32_t dt = 0; dt < patch; dt++) {
    const float* t_row = patch_filter_context_cached_get_frame(
        &self->context, (int32_t)dt - (int32_t)half);
    const float* c_row = patch_filter_context_cached_get_frame(
        &self->context, cand_dt + (int32_t)dt - (int32_t)half);
    ssd = patch_filter_accumulate_patch_row_ssd(
        ssd, (PatchRowArgs){.target_row = t_row,
                            .candidate_row = c_row,
                            .target_freq = target_freq,
                            .candidate_freq = cand_freq,
                            .patch_size = patch,
                            .half_patch_size = half,
                            .spectrum_size = n});
  }
  return ssd;
}

// One paste block = one shared match list amortized over paste bins, then a
// cheap per-bin collaborative hard-threshold shrink. Blocks write disjoint
// bins: race-free. (50%-overlap Hann-tapered dual grids measured identical
// SD/MNI/SISDR at 2x cost: the residual distortion is engine-bound, not
// seam-bound — single grid kept.)
typedef struct {
  Bm3dFilter* filter;
  float* out;
  const float* target_frame;
} Bm3dBlockTask;

// Per-block working state for bm3d_process_block_range. Groups the scalars,
// preloaded target vectors and match lists shared by the stage helpers so
// each helper stays under the parameter-count limit. Same stack footprint
// as the previous inline version; no allocation.
typedef struct {
  Bm3dFilter* filter;
  const float* target;
  float* out;
  uint32_t center;
  uint32_t limit;
  uint32_t block_start;
  uint32_t patch;
  uint32_t half;
  uint32_t n;
  int32_t past;
  int32_t future;
  int32_t search_f;
  float dist_thresh;
  float thr;
  bool use_vec8;
  sb_vec8_t target_vecs[8];
  int32_t m_dt[BM3D_STACK_MAX - 1];
  int32_t m_df[BM3D_STACK_MAX - 1];
  float m_w[BM3D_STACK_MAX - 1];
  float m_d[BM3D_STACK_MAX - 1];
} Bm3dBlockCtx;

static inline SB_UNUSED bool bm3d_preload_target(Bm3dBlockCtx* b) {
  Bm3dFilter* filter = b->filter;
  // Reference patch at the block center; SIMD fast path for patch == 8.
  b->use_vec8 = (b->patch == 8) && (b->center >= b->half) &&
                (b->center + (b->patch - b->half) <= b->n);
  if (!b->use_vec8) {
    return false;
  }
  for (int r = 0; r < 8; r++) {
    b->target_vecs[r] = sb_load8(patch_filter_context_cached_get_frame(
                                     &filter->context, r - (int32_t)b->half) +
                                 (b->center - b->half));
  }
  return true;
}

static inline SB_UNUSED uint32_t bm3d_worst_slot(const float* m_d) {
  uint32_t worst = 0;
  for (uint32_t i = 1; i < BM3D_STACK_MAX - 1; i++) {
    if (m_d[i] > m_d[worst]) {
      worst = i;
    }
  }
  return worst;
}

static inline SB_UNUSED void bm3d_search_block(Bm3dBlockCtx* b) {
  Bm3dFilter* filter = b->filter;
  for (uint32_t i = 0; i < BM3D_STACK_MAX - 1; i++) {
    b->m_d[i] = b->dist_thresh;
    b->m_w[i] = 0.0F;
    b->m_dt[i] = 0;
    b->m_df[i] = 0;
  }
  for (int32_t dt = -b->past; dt <= b->future; dt++) {
    float* cand_rows[8] = {NULL};
    if (b->use_vec8) {
      for (int r = 0; r < 8; r++) {
        cand_rows[r] = patch_filter_context_cached_get_frame(
            &filter->context, dt + r - (int32_t)b->half);
      }
    }
    for (int32_t df = -b->search_f; df <= b->search_f; df++) {
      if (dt == 0 && df == 0) {
        continue;
      }
      uint32_t cf = patch_filter_clamp_index((int32_t)b->center + df, b->n);
      float d;
      if (b->use_vec8 && cf >= b->half && cf + (b->patch - b->half) <= b->n) {
        const uint32_t fs = cf - b->half;
        float* ptrs[8] = {cand_rows[0] + fs, cand_rows[1] + fs,
                          cand_rows[2] + fs, cand_rows[3] + fs,
                          cand_rows[4] + fs, cand_rows[5] + fs,
                          cand_rows[6] + fs, cand_rows[7] + fs};
        d = sb_vec8_patch_ssd(b->target_vecs, ptrs);
      } else {
        d = bm3d_patch_ssd_scalar(filter, b->center, dt, cf);
      }
      if (d >= b->dist_thresh) {
        continue;
      }
      float w = sb_fast_expf(-d * filter->inv_h_squared);
      if (w < NLM_MIN_WEIGHT) {
        continue;
      }
      uint32_t worst = bm3d_worst_slot(b->m_d);
      if (d < b->m_d[worst]) {
        b->m_d[worst] = d;
        b->m_w[worst] = w;
        b->m_dt[worst] = dt;
        b->m_df[worst] = df;
      }
    }
  }
}

static inline SB_UNUSED void bm3d_collab_bin(Bm3dBlockCtx* b, uint32_t k) {
  Bm3dFilter* filter = b->filter;
  // Weighted mean of (mean + kept residuals): matches whose residual
  // survives the hard threshold keep their value, the rest collapse to
  // the stack mean. Flat regions denoise harder than NLM averaging while
  // strong structure survives via residuals. Measured better than a
  // canonical stack-DCT hard-threshold on the low-SNR white-noise case:
  // SSD-matched stacks have low variance, so the residual-keep preserves
  // detail the transform-domain cut averages away.
  float mean = b->target[k]; // anchor weight 1
  float w_sum = 1.0F;
  for (uint32_t m = 0; m < BM3D_STACK_MAX - 1; m++) {
    if (b->m_w[m] <= 0.0F) {
      continue;
    }
    mean +=
        b->m_w[m] * patch_filter_context_cached_get_frame(
                        &filter->context, b->m_dt[m])[patch_filter_clamp_index(
                        (int32_t)k + b->m_df[m], b->n)];
    w_sum += b->m_w[m];
  }
  mean /= w_sum;
  float num = mean; // anchor: mean * 1
  float den = 1.0F;
  for (uint32_t m = 0; m < BM3D_STACK_MAX - 1; m++) {
    if (b->m_w[m] <= 0.0F) {
      continue;
    }
    float v = patch_filter_context_cached_get_frame(
        &filter->context,
        b->m_dt[m])[patch_filter_clamp_index((int32_t)k + b->m_df[m], b->n)];
    float r = v - mean;
    if (fabsf(r) > b->thr) {
      num += b->m_w[m] * (mean + r);
      den += b->m_w[m];
    } else {
      num += b->m_w[m] * mean;
      den += b->m_w[m];
    }
  }
  // confidence blend-back against the raw bin (trust raw at high SNR)
  // now happens once at engine level for all 2D modes, on the aligned
  // delayed map.
  b->out[k] = num / den;
}

static void bm3d_process_block_range(void* raw_ctx, uint32_t first_block,
                                     uint32_t block_count) {
  const Bm3dBlockTask* ctx = (const Bm3dBlockTask*)raw_ctx;
  Bm3dFilter* filter = ctx->filter;
  float* out = ctx->out;
  const float* target = ctx->target_frame;

  const uint32_t n = filter->config.spectrum_size;
  const uint32_t paste = filter->config.paste_block_size;
  const uint32_t patch = filter->config.patch_size;
  const uint32_t half = patch / 2;
  const int32_t past = (int32_t)filter->config.search_range_time_past;
  const int32_t future = (int32_t)filter->config.search_range_time_future;
  const int32_t search_f = (int32_t)filter->config.search_range_freq;
  const float dist_thresh =
      BM3D_DISTANCE_THRESHOLD_MULTIPLIER * filter->h_squared;
  const float thr = BM3D_THRESHOLD_K * filter->config.h_parameter;

  for (uint32_t b = 0; b < block_count; b++) {
    const uint32_t block_start = (first_block + b) * paste;
    uint32_t center = block_start + (paste / 2U);
    if (center >= n) {
      center = n - 1;
    }
    uint32_t limit = paste;
    if (block_start + paste > n) {
      limit = n - block_start;
    }

    float block_sum = 0.0F;
    for (uint32_t i = 0; i < limit; i++) {
      block_sum += target[block_start + i];
    }
    if (block_sum < 1e-6F) {
      memcpy(out + block_start, target + block_start, limit * sizeof(float));
      continue;
    }

    Bm3dBlockCtx blk = {
        .filter = filter,
        .target = target,
        .out = out,
        .center = center,
        .limit = limit,
        .block_start = block_start,
        .patch = patch,
        .half = half,
        .n = n,
        .past = past,
        .future = future,
        .search_f = search_f,
        .dist_thresh = dist_thresh,
        .thr = thr,
        .use_vec8 = false,
    };

    // Reference patch at the block center; SIMD fast path for patch == 8.
    bm3d_preload_target(&blk);

    // Top-(STACK_MAX-1) matches; the target anchors the stack at combine.
    bm3d_search_block(&blk);

    // Per-bin collaborative shrink with the shared match list.
    for (uint32_t i = 0; i < limit; i++) {
      bm3d_collab_bin(&blk, block_start + i);
    }
  }
}

Bm3dFilter* bm3d_filter_initialize(Bm3dFilterConfig config) {
  if (config.spectrum_size == 0) {
    return NULL;
  }
  Bm3dFilter* self = (Bm3dFilter*)calloc(1U, sizeof(Bm3dFilter));
  if (!self) {
    return NULL;
  }
  self->config = config;
  sb_resolve_patch_geometry_defaults(
      &self->config.patch_size, &self->config.paste_block_size,
      &self->config.search_range_freq, &self->config.search_range_time_past,
      &self->config.search_range_time_future);
  if (self->config.time_buffer_size == 0) {
    self->config.time_buffer_size = self->config.search_range_time_past +
                                    self->config.search_range_time_future + 1;
  }
  // Geometry validation: the frame ring needs one slot per relative frame
  // (-past..+future) and the pointer cache assumes the patch half-width fits
  // the halo. Anything else reads aliased or out-of-bounds frames.
  if (self->config.time_buffer_size <
          self->config.search_range_time_past +
              self->config.search_range_time_future + 1U ||
      self->config.patch_size / 2U > NLM_HALO_FRAMES) {
    bm3d_filter_free(self);
    return NULL;
  }
  bm3d_filter_set_h_parameter(self, config.h_parameter <= 0.0F
                                        ? BM3D_DEFAULT_H_PARAMETER
                                        : config.h_parameter);

  if (!patch_filter_context_initialize(
          &self->context, self->config.spectrum_size,
          self->config.time_buffer_size, self->config.search_range_time_past,
          self->config.search_range_time_future, self->config.num_threads)) {
    bm3d_filter_free(self);
    return NULL;
  }
  return self;
}

void bm3d_filter_free(Bm3dFilter* filter) {
  if (!filter) {
    return;
  }
  patch_filter_context_free(&filter->context);
  free(filter);
}

void bm3d_filter_set_h_parameter(Bm3dFilter* filter, float h) {
  if (!filter) {
    return;
  }
  if (h <= 0.0F) {
    filter->config.h_parameter = 0.0F;
    filter->h_squared = 0.0F;
    filter->inv_h_squared = 0.0F;
    return;
  }
  float target_h = h > BM3D_MAX_H_PARAMETER ? BM3D_MAX_H_PARAMETER : h;
  filter->config.h_parameter = target_h;
  filter->h_squared = target_h * target_h;
  filter->inv_h_squared = 1.0F / filter->h_squared;
}

void bm3d_filter_push_frame(Bm3dFilter* filter, const float* snr_frame) {
  if (!filter || !snr_frame) {
    return;
  }
  patch_filter_context_push_frame(&filter->context, snr_frame);
}

bool bm3d_filter_is_ready(const Bm3dFilter* filter) {
  return filter && patch_filter_context_is_ready(&filter->context);
}

void bm3d_filter_reset(Bm3dFilter* filter) {
  if (!filter) {
    return;
  }
  patch_filter_context_reset(&filter->context);
}

uint32_t bm3d_filter_get_latency_frames(const Bm3dFilter* filter) {
  if (!filter) {
    return 0;
  }
  return patch_filter_context_get_latency_frames(&filter->context);
}

bool bm3d_filter_process(Bm3dFilter* filter, float* smoothed_snr) {
  if (!filter || !smoothed_snr || !bm3d_filter_is_ready(filter)) {
    return false;
  }
  sb_simd_state_t old_state = sb_simd_enable_ftz_daz();
  patch_filter_context_populate_frame_ptrs(&filter->context);
  const float* target =
      patch_filter_context_cached_get_frame(&filter->context, 0);

  if (filter->config.h_parameter <= 0.0F) {
    memcpy(smoothed_snr, target, filter->config.spectrum_size * sizeof(float));
    sb_simd_restore_state(old_state);
    return true;
  }

  const uint32_t n = filter->config.spectrum_size;
  const uint32_t paste = filter->config.paste_block_size;
  const uint32_t num_blocks = (n + paste - 1U) / paste;
  Bm3dBlockTask task = {filter, smoothed_snr, target};
  if (filter->context.pool) {
    sb_thread_pool_parallel_for(filter->context.pool, num_blocks,
                                bm3d_process_block_range, &task);
  } else {
    bm3d_process_block_range(&task, 0, num_blocks);
  }
  sb_simd_restore_state(old_state);
  return true;
}

void bm3d_filter_calculate_snr(const Bm3dFilter* filter,
                               const float* reference_spectrum,
                               const float* noise_spectrum, float* snr_frame) {
  if (!filter || !reference_spectrum || !noise_spectrum || !snr_frame) {
    return;
  }
  const uint32_t n = filter->config.spectrum_size;
  for (uint32_t k = 0; k < n; k++) {
    float denom = noise_spectrum[k] > NLM_SNR_NOISE_FLOOR_MIN
                      ? noise_spectrum[k]
                      : NLM_SNR_NOISE_FLOOR_MIN;
    snr_frame[k] = sqrtf(reference_spectrum[k] / denom);
  }
}

void bm3d_filter_reconstruct_magnitude(const Bm3dFilter* filter,
                                       const float* smoothed_snr,
                                       const float* noise_spectrum,
                                       float* magnitude_spectrum) {
  if (!filter || !smoothed_snr || !noise_spectrum || !magnitude_spectrum) {
    return;
  }
  const uint32_t n = filter->config.spectrum_size;
  for (uint32_t k = 0; k < n; k++) {
    float denom = noise_spectrum[k] > NLM_SNR_NOISE_FLOOR_MIN
                      ? noise_spectrum[k]
                      : NLM_SNR_NOISE_FLOOR_MIN;
    float s = smoothed_snr[k];
    magnitude_spectrum[k] = s * s * denom;
  }
}
