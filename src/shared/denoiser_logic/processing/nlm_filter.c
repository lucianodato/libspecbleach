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
#include "shared/denoiser_logic/processing/nlm_filter_internal.h"
#include "shared/utils/simd_utils.h"

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

NlmFilter* nlm_filter_initialize(NlmFilterConfig config) {
  if (config.spectrum_size == 0) {
    return NULL;
  }

  NlmFilter* self = (NlmFilter*)calloc(1U, sizeof(NlmFilter));
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
  if (self->config.h_parameter <= 0.0F) {
    self->config.h_parameter = NLM_DEFAULT_H_PARAMETER;
  }

  if (self->config.time_buffer_size == 0) {
    self->config.time_buffer_size = self->config.search_range_time_past +
                                    self->config.search_range_time_future + 1;
  }

  self->h_squared = self->config.h_parameter * self->config.h_parameter;
  self->inv_h_squared = 1.0f / self->h_squared;
  if (self->config.distance_threshold <= 0.0F) {
    self->distance_threshold_actual = 4.0F * self->h_squared;
  } else {
    self->distance_threshold_actual = self->config.distance_threshold;
  }

  self->target_frame_offset = self->config.search_range_time_past;

  self->num_threads = 4;
#ifdef _OPENMP
  const char* omp_env = getenv("OMP_NUM_THREADS");
  if (omp_env) {
    int val = atoi(omp_env);
    if (val > 0) {
      self->num_threads = (uint32_t)val;
    }
  }
#endif

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

  self->weight_accum =
      (float*)calloc(self->config.spectrum_size, sizeof(float));
  if (!self->weight_accum) {
    nlm_filter_free(self);
    return NULL;
  }

  self->total_time_span = self->config.search_range_time_past +
                          self->config.search_range_time_future + 1 + 8;
  self->frame_ptrs = (float**)calloc(self->total_time_span, sizeof(float*));
  if (!self->frame_ptrs) {
    nlm_filter_free(self);
    return NULL;
  }

  self->process_fn = nlm_filter_process_generic;

#if defined(__x86_64__) || defined(__i386__)
  if (__builtin_cpu_supports("avx")) {
    self->process_fn = nlm_filter_process_avx;
  }
#endif

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

  if (filter->frame_ptrs) {
    free(filter->frame_ptrs);
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

  if (filter->config.distance_threshold <= 0.0F) {
    filter->distance_threshold_actual = 4.0F * filter->h_squared;
  }
}

void nlm_filter_push_frame(NlmFilter* filter, const float* snr_frame) {
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

  if (!nlm_filter_is_ready(filter)) {
    return false;
  }

  return filter->process_fn(filter, smoothed_snr);
}

#ifdef _OPENMP
#include <omp.h>
#endif

bool nlm_filter_process_generic(NlmFilter* filter, float* smoothed_snr) {
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

  const float current_inv_h2 = filter->inv_h_squared;
  const float current_dist_threshold =
      filter->distance_threshold_actual;

#pragma omp parallel for schedule(dynamic) num_threads(filter->num_threads)
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

    for (int32_t dt = -search_time_past; dt <= search_time_future;
         dt++) {
      float* cand_rows[8];
      if (filter->config.patch_size == 8) {
        for (int r = 0; r < 8; r++) {
          cand_rows[r] = cached_get_frame(filter, dt + r - 4);
        }
      }

      for (int32_t df = -(int32_t)search_freq;
           df <= (int32_t)search_freq; df++) {

        uint32_t cand_center =
            clamp_index((int32_t)block_center + df, spectrum_size);

        float distance = 0.0F;

        bool safe_bounds =
            safe_block && (cand_center >= half_patch_size) &&
            (cand_center + half_patch_size <= spectrum_size);

        if (filter->config.patch_size == 8 && safe_bounds) {
          sb_acc8_t sum = sb_acc8_zero();
          uint32_t cand_f_start = cand_center - 4;

          sum = sb_acc8_add_ssd(sum, target_vecs[0],
                                sb_load8(cand_rows[0] + cand_f_start));
          sum = sb_acc8_add_ssd(sum, target_vecs[1],
                                sb_load8(cand_rows[1] + cand_f_start));
          sum = sb_acc8_add_ssd(sum, target_vecs[2],
                                sb_load8(cand_rows[2] + cand_f_start));
          sum = sb_acc8_add_ssd(sum, target_vecs[3],
                                sb_load8(cand_rows[3] + cand_f_start));
          sum = sb_acc8_add_ssd(sum, target_vecs[4],
                                sb_load8(cand_rows[4] + cand_f_start));
          sum = sb_acc8_add_ssd(sum, target_vecs[5],
                                sb_load8(cand_rows[5] + cand_f_start));
          sum = sb_acc8_add_ssd(sum, target_vecs[6],
                                sb_load8(cand_rows[6] + cand_f_start));
          sum = sb_acc8_add_ssd(sum, target_vecs[7],
                                sb_load8(cand_rows[7] + cand_f_start));

          distance = sb_acc8_hsum(sum);
        } else {
          distance = compute_patch_distance(filter, 0, block_center,
                                            dt, cand_center);
        }

        if (distance > current_dist_threshold) {
          continue;
        }

        float weight = sb_fast_expf(-distance * current_inv_h2);
        if (weight < NLM_MIN_WEIGHT) {
          continue;
        }

        float* cand_frame = cached_get_frame(filter, dt);

        if (current_paste_limit == 8) {
          uint32_t t0 = block_start, t1 = block_start + 1,
                   t2 = block_start + 2, t3 = block_start + 3,
                   t4 = block_start + 4, t5 = block_start + 5,
                   t6 = block_start + 6, t7 = block_start + 7;

          uint32_t c0 = clamp_index((int32_t)t0 + df, spectrum_size),
                   c1 = clamp_index((int32_t)t1 + df, spectrum_size),
                   c2 = clamp_index((int32_t)t2 + df, spectrum_size),
                   c3 = clamp_index((int32_t)t3 + df, spectrum_size),
                   c4 = clamp_index((int32_t)t4 + df, spectrum_size),
                   c5 = clamp_index((int32_t)t5 + df, spectrum_size),
                   c6 = clamp_index((int32_t)t6 + df, spectrum_size),
                   c7 = clamp_index((int32_t)t7 + df, spectrum_size);

          smoothed_snr[t0] += weight * cand_frame[c0];
          weight_sum[t0] += weight;
          smoothed_snr[t1] += weight * cand_frame[c1];
          weight_sum[t1] += weight;
          smoothed_snr[t2] += weight * cand_frame[c2];
          weight_sum[t2] += weight;
          smoothed_snr[t3] += weight * cand_frame[c3];
          weight_sum[t3] += weight;
          smoothed_snr[t4] += weight * cand_frame[c4];
          weight_sum[t4] += weight;
          smoothed_snr[t5] += weight * cand_frame[c5];
          weight_sum[t5] += weight;
          smoothed_snr[t6] += weight * cand_frame[c6];
          weight_sum[t6] += weight;
          smoothed_snr[t7] += weight * cand_frame[c7];
          weight_sum[t7] += weight;
        } else if (current_paste_limit == 4) {
          uint32_t t0 = block_start, t1 = block_start + 1,
                   t2 = block_start + 2, t3 = block_start + 3;

          uint32_t c0 = clamp_index((int32_t)t0 + df, spectrum_size),
                   c1 = clamp_index((int32_t)t1 + df, spectrum_size),
                   c2 = clamp_index((int32_t)t2 + df, spectrum_size),
                   c3 = clamp_index((int32_t)t3 + df, spectrum_size);

          smoothed_snr[t0] += weight * cand_frame[c0];
          weight_sum[t0] += weight;
          smoothed_snr[t1] += weight * cand_frame[c1];
          weight_sum[t1] += weight;
          smoothed_snr[t2] += weight * cand_frame[c2];
          weight_sum[t2] += weight;
          smoothed_snr[t3] += weight * cand_frame[c3];
          weight_sum[t3] += weight;
        } else {
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

void nlm_filter_reset(NlmFilter* filter) {
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

uint32_t nlm_filter_get_latency_frames(NlmFilter* filter) {
  if (!filter) {
    return 0;
  }
  return filter->config.search_range_time_future;
}

void nlm_filter_calculate_snr(NlmFilter* filter,
                              const float* reference_spectrum,
                              const float* noise_spectrum, float* snr_frame) {
  if (!filter || !reference_spectrum || !noise_spectrum || !snr_frame) {
    return;
  }

  const uint32_t spectrum_size = filter->config.spectrum_size;
  uint32_t k = 0;

  sb_vec8_t eps = sb_set8(SPECTRAL_EPSILON);
  sb_vec8_t flt_min = sb_set8(FLT_MIN);

  for (; k + 7 < spectrum_size; k += 8) {
    sb_vec8_t noise = sb_load8(noise_spectrum + k);
    sb_vec8_t mask = sb_gt8(noise, flt_min);
    sb_vec8_t denom = sb_sel8(mask, noise, eps);
    sb_vec8_t snr = sb_div8(sb_load8(reference_spectrum + k), denom);
    sb_store8(snr_frame + k, snr);
  }

  for (; k < spectrum_size; k++) {
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
  uint32_t k = 0;

  sb_vec8_t eps = sb_set8(SPECTRAL_EPSILON);
  sb_vec8_t flt_min = sb_set8(FLT_MIN);

  for (; k + 7 < spectrum_size; k += 8) {
    sb_vec8_t noise = sb_load8(noise_spectrum + k);
    sb_vec8_t mask = sb_gt8(noise, flt_min);
    sb_vec8_t denom = sb_sel8(mask, noise, eps);
    sb_vec8_t mag = sb_mul8(sb_load8(smoothed_snr + k), denom);
    sb_store8(magnitude_spectrum + k, mag);
  }

  for (; k < spectrum_size; k++) {
    float denom =
        noise_spectrum[k] > FLT_MIN ? noise_spectrum[k] : SPECTRAL_EPSILON;
    magnitude_spectrum[k] = smoothed_snr[k] * denom;
  }
}
