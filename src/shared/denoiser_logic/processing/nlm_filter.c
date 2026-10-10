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

NlmFilter* nlm_filter_initialize(NlmFilterConfig config) {
  if (config.spectrum_size == 0) {
    return NULL;
  }

  NlmFilter* self = (NlmFilter*)calloc(1U, sizeof(NlmFilter));
  if (!self) {
    return NULL;
  }

  self->config = config;
  sb_resolve_patch_geometry_defaults(
      &self->config.patch_size, &self->config.paste_block_size,
      &self->config.search_range_freq, &self->config.search_range_time_past,
      &self->config.search_range_time_future);
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

  if (!patch_filter_context_initialize(
          &self->context, self->config.spectrum_size,
          self->config.time_buffer_size, self->config.search_range_time_past,
          self->config.search_range_time_future, self->config.num_threads)) {
    nlm_filter_free(self);
    return NULL;
  }

  self->weight_accum =
      (float*)calloc(self->config.spectrum_size, sizeof(float));
  if (!self->weight_accum) {
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

  patch_filter_context_free(&filter->context);
  if (filter->weight_accum) {
    free(filter->weight_accum);
  }

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

  patch_filter_context_push_frame(&filter->context, snr_frame);
}

bool nlm_filter_is_ready(const NlmFilter* filter) {
  if (!filter) {
    return false;
  }
  return patch_filter_context_is_ready(&filter->context);
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

  patch_filter_context_reset(&filter->context);
}

uint32_t nlm_filter_get_latency_frames(const NlmFilter* filter) {
  if (!filter) {
    return 0;
  }
  return patch_filter_context_get_latency_frames(&filter->context);
}

void nlm_filter_calculate_snr(const NlmFilter* filter,
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

void nlm_filter_reconstruct_magnitude(const NlmFilter* filter,
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
