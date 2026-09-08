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
#include "shared/utils/simd_utils.h"
#include "shared/utils/thread_pool.h"

struct Bm3dFilter {
  Bm3dFilterConfig config;
  float** frame_buffer;
  uint32_t buffer_head;
  uint32_t frames_filled;
  float h_squared;
  float inv_h_squared;
  // Pre-computed frame pointer cache (same halo scheme as the NLM filter).
  float** frame_ptrs;
  uint32_t total_time_span;
  SbThreadPool* pool;
  uint32_t num_threads;
};

static inline uint32_t bm3d_clamp(int32_t idx, uint32_t max_val) {
  if (idx < 0) {
    return 0;
  }
  if ((uint32_t)idx >= max_val) {
    return max_val - 1;
  }
  return (uint32_t)idx;
}

static inline float* bm3d_get_frame(Bm3dFilter* self, int32_t rel) {
  const int32_t size = (int32_t)self->config.time_buffer_size;
  int32_t idx = (int32_t)self->buffer_head -
                (int32_t)self->config.search_range_time_future - 1 + rel;
  idx = ((idx % size) + size) % size;
  return self->frame_buffer[idx];
}

static inline void bm3d_populate_ptrs(Bm3dFilter* self) {
  const int32_t past = (int32_t)self->config.search_range_time_past;
  const int32_t future = (int32_t)self->config.search_range_time_future;
  for (int32_t dt = -past - (int32_t)NLM_HALO_FRAMES;
       dt <= future + (int32_t)NLM_HALO_FRAMES; dt++) {
    self->frame_ptrs[past + (int32_t)NLM_HALO_FRAMES + dt] =
        bm3d_get_frame(self, dt);
  }
}

static inline float* bm3d_cached_frame(Bm3dFilter* self, int32_t dt) {
  return self->frame_ptrs[(int32_t)self->config.search_range_time_past +
                          (int32_t)NLM_HALO_FRAMES + dt];
}

// Scalar clamped patch SSD (edges / non-8 patch sizes).
static float bm3d_patch_ssd_scalar(Bm3dFilter* self, uint32_t target_freq,
                                   int32_t cand_dt, uint32_t cand_freq) {
  const uint32_t patch = self->config.patch_size;
  const uint32_t half = patch / 2;
  const uint32_t n = self->config.spectrum_size;
  float ssd = 0.0F;
  for (uint32_t dt = 0; dt < patch; dt++) {
    const float* t_row = bm3d_cached_frame(self, (int32_t)dt - (int32_t)half);
    const float* c_row =
        bm3d_cached_frame(self, cand_dt + (int32_t)dt - (int32_t)half);
    for (uint32_t df = 0; df < patch; df++) {
      uint32_t tf =
          bm3d_clamp((int32_t)target_freq + (int32_t)df - (int32_t)half, n);
      uint32_t cf =
          bm3d_clamp((int32_t)cand_freq + (int32_t)df - (int32_t)half, n);
      float d = t_row[tf] - c_row[cf];
      ssd += d * d;
    }
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
    uint32_t center = block_start + paste / 2;
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

    // Reference patch at the block center; SIMD fast path for patch == 8.
    bool use_vec8 =
        (patch == 8) && (center >= half) && (center + (patch - half) <= n);
    sb_vec8_t target_vecs[8];
    if (use_vec8) {
      for (int r = 0; r < 8; r++) {
        target_vecs[r] = sb_load8(bm3d_cached_frame(filter, r - (int32_t)half) +
                                  (center - half));
      }
    }

    // Top-(STACK_MAX-1) matches; the target anchors the stack at combine.
    int32_t m_dt[BM3D_STACK_MAX - 1];
    int32_t m_df[BM3D_STACK_MAX - 1];
    float m_w[BM3D_STACK_MAX - 1];
    float m_d[BM3D_STACK_MAX - 1];
    for (uint32_t i = 0; i < BM3D_STACK_MAX - 1; i++) {
      m_d[i] = dist_thresh;
      m_w[i] = 0.0F;
      m_dt[i] = 0;
      m_df[i] = 0;
    }
    for (int32_t dt = -past; dt <= future; dt++) {
      float* cand_rows[8] = {NULL};
      if (use_vec8) {
        for (int r = 0; r < 8; r++) {
          cand_rows[r] = bm3d_cached_frame(filter, dt + r - (int32_t)half);
        }
      }
      for (int32_t df = -search_f; df <= search_f; df++) {
        if (dt == 0 && df == 0) {
          continue;
        }
        uint32_t cf = bm3d_clamp((int32_t)center + df, n);
        float d;
        if (use_vec8 && cf >= half && cf + (patch - half) <= n) {
          const uint32_t fs = cf - half;
          float* ptrs[8] = {cand_rows[0] + fs, cand_rows[1] + fs,
                            cand_rows[2] + fs, cand_rows[3] + fs,
                            cand_rows[4] + fs, cand_rows[5] + fs,
                            cand_rows[6] + fs, cand_rows[7] + fs};
          d = sb_vec8_patch_ssd(target_vecs, ptrs);
        } else {
          d = bm3d_patch_ssd_scalar(filter, center, dt, cf);
        }
        if (d >= dist_thresh) {
          continue;
        }
        float w = sb_fast_expf(-d * filter->inv_h_squared);
        if (w < NLM_MIN_WEIGHT) {
          continue;
        }
        uint32_t worst = 0;
        for (uint32_t i = 1; i < BM3D_STACK_MAX - 1; i++) {
          if (m_d[i] > m_d[worst]) {
            worst = i;
          }
        }
        if (d < m_d[worst]) {
          m_d[worst] = d;
          m_w[worst] = w;
          m_dt[worst] = dt;
          m_df[worst] = df;
        }
      }
    }

    // Per-bin collaborative shrink with the shared match list.
    for (uint32_t i = 0; i < limit; i++) {
      const uint32_t k = block_start + i;
      // Weighted mean of (mean + kept residuals): matches whose residual
      // survives the hard threshold keep their value, the rest collapse to
      // the stack mean. Flat regions denoise harder than NLM averaging while
      // strong structure survives via residuals. Measured better than a
      // canonical stack-DCT hard-threshold on the low-SNR white-noise case:
      // SSD-matched stacks have low variance, so the residual-keep preserves
      // detail the transform-domain cut averages away.
      {
        float mean = target[k]; // anchor weight 1
        float w_sum = 1.0F;
        for (uint32_t m = 0; m < BM3D_STACK_MAX - 1; m++) {
          if (m_w[m] <= 0.0F) {
            continue;
          }
          mean += m_w[m] *
                  bm3d_cached_frame(
                      filter, m_dt[m])[bm3d_clamp((int32_t)k + m_df[m], n)];
          w_sum += m_w[m];
        }
        mean /= w_sum;
        float num = mean; // anchor: mean * 1
        float den = 1.0F;
        for (uint32_t m = 0; m < BM3D_STACK_MAX - 1; m++) {
          if (m_w[m] <= 0.0F) {
            continue;
          }
          float v = bm3d_cached_frame(
              filter, m_dt[m])[bm3d_clamp((int32_t)k + m_df[m], n)];
          float r = v - mean;
          if (fabsf(r) > thr) {
            num += m_w[m] * (mean + r);
            den += m_w[m];
          } else {
            num += m_w[m] * mean;
            den += m_w[m];
          }
        }
        // confidence blend-back against the raw bin (trust raw at high SNR)
        // now happens once at engine level for all 2D modes, on the aligned
        // delayed map.
        out[k] = num / den;
      }
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
  if (self->config.time_buffer_size == 0) {
    self->config.time_buffer_size = self->config.search_range_time_past +
                                    self->config.search_range_time_future + 1;
  }
  bm3d_filter_set_h_parameter(self, config.h_parameter <= 0.0F
                                        ? BM3D_DEFAULT_H_PARAMETER
                                        : config.h_parameter);

  self->num_threads = self->config.num_threads > 0U ? self->config.num_threads
                                                    : NLM_NUM_THREADS_DEFAULT;
  if (self->num_threads > NLM_MAX_THREADS) {
    self->num_threads = NLM_MAX_THREADS;
  }
  if (self->num_threads > 1U) {
    self->pool = sb_thread_pool_create(self->num_threads - 1U);
    if (!self->pool) {
      self->num_threads = 1U;
    }
  } else {
    self->num_threads = 1U;
  }

  self->frame_buffer =
      (float**)calloc(self->config.time_buffer_size, sizeof(float*));
  if (!self->frame_buffer) {
    bm3d_filter_free(self);
    return NULL;
  }
  for (uint32_t i = 0; i < self->config.time_buffer_size; i++) {
    self->frame_buffer[i] =
        (float*)calloc(self->config.spectrum_size, sizeof(float));
    if (!self->frame_buffer[i]) {
      bm3d_filter_free(self);
      return NULL;
    }
  }
  self->total_time_span = self->config.search_range_time_past +
                          self->config.search_range_time_future + 1U +
                          2U * NLM_HALO_FRAMES;
  self->frame_ptrs = (float**)calloc(self->total_time_span, sizeof(float*));
  if (!self->frame_ptrs) {
    bm3d_filter_free(self);
    return NULL;
  }
  return self;
}

void bm3d_filter_free(Bm3dFilter* filter) {
  if (!filter) {
    return;
  }
  if (filter->frame_buffer) {
    for (uint32_t i = 0; i < filter->config.time_buffer_size; i++) {
      free(filter->frame_buffer[i]);
    }
    free(filter->frame_buffer);
  }
  free(filter->frame_ptrs);
  sb_thread_pool_free(filter->pool);
  filter->pool = NULL;
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
  memcpy(filter->frame_buffer[filter->buffer_head], snr_frame,
         filter->config.spectrum_size * sizeof(float));
  filter->buffer_head =
      (filter->buffer_head + 1) % filter->config.time_buffer_size;
  if (filter->frames_filled < filter->config.time_buffer_size) {
    filter->frames_filled++;
  }
}

bool bm3d_filter_is_ready(Bm3dFilter* filter) {
  return filter && filter->frames_filled >= filter->config.time_buffer_size;
}

void bm3d_filter_reset(Bm3dFilter* filter) {
  if (!filter) {
    return;
  }
  for (uint32_t i = 0; i < filter->config.time_buffer_size; i++) {
    memset(filter->frame_buffer[i], 0,
           filter->config.spectrum_size * sizeof(float));
  }
  filter->buffer_head = 0;
  filter->frames_filled = 0;
}

uint32_t bm3d_filter_get_latency_frames(Bm3dFilter* filter) {
  if (!filter) {
    return 0;
  }
  return filter->config.search_range_time_future;
}

bool bm3d_filter_process(Bm3dFilter* filter, float* smoothed_snr) {
  if (!filter || !smoothed_snr || !bm3d_filter_is_ready(filter)) {
    return false;
  }
  sb_simd_state_t old_state = sb_simd_enable_ftz_daz();
  bm3d_populate_ptrs(filter);
  float* target = bm3d_cached_frame(filter, 0);

  if (filter->config.h_parameter <= 0.0F) {
    memcpy(smoothed_snr, target, filter->config.spectrum_size * sizeof(float));
    sb_simd_restore_state(old_state);
    return true;
  }

  const uint32_t n = filter->config.spectrum_size;
  const uint32_t paste = filter->config.paste_block_size;
  const uint32_t num_blocks = (n + paste - 1U) / paste;
  Bm3dBlockTask task = {filter, smoothed_snr, target};
  if (filter->pool) {
    sb_thread_pool_parallel_for(filter->pool, num_blocks,
                                bm3d_process_block_range, &task);
  } else {
    bm3d_process_block_range(&task, 0, num_blocks);
  }
  sb_simd_restore_state(old_state);
  return true;
}

void bm3d_filter_calculate_snr(Bm3dFilter* filter,
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

void bm3d_filter_reconstruct_magnitude(Bm3dFilter* filter,
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
