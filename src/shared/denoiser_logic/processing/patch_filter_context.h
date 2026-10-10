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

#ifndef PATCH_FILTER_CONTEXT_H
#define PATCH_FILTER_CONTEXT_H

#include "shared/configurations.h"
#include "shared/utils/general_utils.h"
#include "shared/utils/thread_pool.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct PatchFilterContext {
  float** frame_buffer;
  uint32_t buffer_head;
  uint32_t frames_filled;
  float** frame_ptrs;
  uint32_t time_buffer_size;
  uint32_t spectrum_size;
  uint32_t search_range_time_past;
  uint32_t search_range_time_future;
  SbThreadPool* pool;
  uint32_t num_threads;
} PatchFilterContext;

bool patch_filter_context_initialize(PatchFilterContext* context,
                                     uint32_t spectrum_size,
                                     uint32_t time_buffer_size,
                                     uint32_t search_range_time_past,
                                     uint32_t search_range_time_future,
                                     uint32_t requested_threads);
void patch_filter_context_free(PatchFilterContext* context);
void patch_filter_context_push_frame(PatchFilterContext* context,
                                     const float* frame);
bool patch_filter_context_is_ready(const PatchFilterContext* context);
void patch_filter_context_reset(PatchFilterContext* context);
uint32_t patch_filter_context_get_latency_frames(
    const PatchFilterContext* context);

static inline SB_UNUSED uint32_t patch_filter_clamp_index(int32_t idx,
                                                          uint32_t max_val) {
  if (idx < 0) {
    return 0;
  }
  if ((uint32_t)idx >= max_val) {
    return max_val - 1;
  }
  return (uint32_t)idx;
}

// Keep ring/cache lookup with each filter; only the identical row sum is
// shared.
static inline SB_UNUSED float patch_filter_accumulate_patch_row_ssd(
    float distance, const float* target_row, const float* candidate_row,
    uint32_t target_freq, uint32_t candidate_freq, uint32_t patch_size,
    uint32_t half_patch_size, uint32_t spectrum_size) {
  for (uint32_t df = 0; df < patch_size; df++) {
    const uint32_t target_bin = patch_filter_clamp_index(
        (int32_t)target_freq + (int32_t)df - (int32_t)half_patch_size,
        spectrum_size);
    const uint32_t candidate_bin = patch_filter_clamp_index(
        (int32_t)candidate_freq + (int32_t)df - (int32_t)half_patch_size,
        spectrum_size);
    const float diff = target_row[target_bin] - candidate_row[candidate_bin];
    distance += diff * diff;
  }
  return distance;
}

static inline SB_UNUSED float* patch_filter_context_get_frame(
    PatchFilterContext* context, int32_t relative_offset) {
  const int32_t size = (int32_t)context->time_buffer_size;
  int32_t idx = (int32_t)context->buffer_head -
                (int32_t)context->search_range_time_future - 1 +
                relative_offset;
  idx = ((idx % size) + size) % size;
  return context->frame_buffer[idx];
}

static inline SB_UNUSED void patch_filter_context_populate_frame_ptrs(
    PatchFilterContext* context) {
  const int32_t past = (int32_t)context->search_range_time_past;
  const int32_t future = (int32_t)context->search_range_time_future;
  for (int32_t dt = -past - (int32_t)NLM_HALO_FRAMES;
       dt <= future + (int32_t)NLM_HALO_FRAMES; dt++) {
    context->frame_ptrs[past + (int32_t)NLM_HALO_FRAMES + dt] =
        patch_filter_context_get_frame(context, dt);
  }
}

static inline SB_UNUSED float* patch_filter_context_cached_get_frame(
    PatchFilterContext* context, int32_t dt) {
  return context->frame_ptrs[(int32_t)context->search_range_time_past +
                             (int32_t)NLM_HALO_FRAMES + dt];
}

#endif
