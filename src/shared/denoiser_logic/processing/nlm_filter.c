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
#include "shared/utils/thread_pool.h"

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
  } else if (self->config.h_parameter > NLM_MAX_H_PARAMETER) {
    self->config.h_parameter = NLM_MAX_H_PARAMETER;
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

  self->num_threads = self->config.num_threads > 0U ? self->config.num_threads
                                                    : NLM_NUM_THREADS_DEFAULT;
  if (self->num_threads > NLM_MAX_THREADS) {
    self->num_threads = NLM_MAX_THREADS;
  }
  if (self->num_threads > 1U) {
    // The calling thread participates in every dispatch, so the pool only
    // needs num_threads - 1 dedicated workers.
    self->pool = sb_thread_pool_create(self->num_threads - 1U);
    if (!self->pool) {
      self->num_threads = 1U; // Graceful fallback to sequential processing
    }
  } else {
    self->num_threads = 1U;
  }

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
                          self->config.search_range_time_future + 1 + 16;
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
    free((void*)filter->frame_ptrs);
  }

  sb_thread_pool_free(filter->pool);
  filter->pool = NULL;

  free(filter);
}

void nlm_filter_set_h_parameter(NlmFilter* filter, float h) {
  if (!filter) {
    return;
  }

  if (h <= 0.0F) {
    filter->config.h_parameter = 0.0F;
    filter->h_squared = 0.0F;
    filter->inv_h_squared = 0.0F;
    filter->distance_threshold_actual = 0.0F;
    return;
  }

  float target_h = h;
  if (target_h > NLM_MAX_H_PARAMETER) {
    target_h = NLM_MAX_H_PARAMETER;
  }
  filter->config.h_parameter = target_h;
  filter->h_squared = target_h * target_h;
  filter->inv_h_squared = 1.0F / filter->h_squared;

  if (filter->config.distance_threshold <= 0.0F) {
    filter->distance_threshold_actual = 4.0F * filter->h_squared;
  } else {
    filter->distance_threshold_actual = filter->config.distance_threshold;
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

bool nlm_filter_process_generic(NlmFilter* filter, float* smoothed_snr) {
  return nlm_filter_process_core(filter, smoothed_snr);
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

  sb_vec8_t noise_floor_min = sb_set8(NLM_SNR_NOISE_FLOOR_MIN);

  for (; k + 7 < spectrum_size; k += 8) {
    sb_vec8_t noise = sb_load8(noise_spectrum + k);
    sb_vec8_t mask = sb_gt8(noise, noise_floor_min);
    sb_vec8_t denom = sb_sel8(mask, noise, noise_floor_min);
    sb_vec8_t power_snr = sb_div8(sb_load8(reference_spectrum + k), denom);
    sb_vec8_t snr = sb_sqrt8(power_snr);
    sb_store8(snr_frame + k, snr);
  }

  for (; k < spectrum_size; k++) {
    float denom = noise_spectrum[k] > NLM_SNR_NOISE_FLOOR_MIN
                      ? noise_spectrum[k]
                      : NLM_SNR_NOISE_FLOOR_MIN;
    snr_frame[k] = sqrtf(reference_spectrum[k] / denom);
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

  sb_vec8_t noise_floor_min = sb_set8(NLM_SNR_NOISE_FLOOR_MIN);

  for (; k + 7 < spectrum_size; k += 8) {
    sb_vec8_t noise = sb_load8(noise_spectrum + k);
    sb_vec8_t mask = sb_gt8(noise, noise_floor_min);
    sb_vec8_t denom = sb_sel8(mask, noise, noise_floor_min);
    sb_vec8_t smoothed = sb_load8(smoothed_snr + k);
    sb_vec8_t smoothed_sq = sb_mul8(smoothed, smoothed);
    sb_vec8_t mag = sb_mul8(smoothed_sq, denom);
    sb_store8(magnitude_spectrum + k, mag);
  }

  for (; k < spectrum_size; k++) {
    float denom = noise_spectrum[k] > NLM_SNR_NOISE_FLOOR_MIN
                      ? noise_spectrum[k]
                      : NLM_SNR_NOISE_FLOOR_MIN;
    float smoothed = smoothed_snr[k];
    magnitude_spectrum[k] = smoothed * smoothed * denom;
  }
}
