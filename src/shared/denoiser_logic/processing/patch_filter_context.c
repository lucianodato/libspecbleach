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

#include "patch_filter_context.h"
#include <stdlib.h>
#include <string.h>

bool patch_filter_context_initialize(PatchFilterContext* context,
                                     uint32_t spectrum_size,
                                     uint32_t time_buffer_size,
                                     uint32_t search_range_time_past,
                                     uint32_t search_range_time_future,
                                     uint32_t requested_threads) {
  if (!context) {
    return false;
  }

  memset(context, 0, sizeof(*context));
  context->spectrum_size = spectrum_size;
  context->time_buffer_size = time_buffer_size;
  context->search_range_time_past = search_range_time_past;
  context->search_range_time_future = search_range_time_future;

  context->num_threads =
      requested_threads > 0U ? requested_threads : NLM_NUM_THREADS_DEFAULT;
  if (context->num_threads > NLM_MAX_THREADS) {
    context->num_threads = NLM_MAX_THREADS;
  }
  if (context->num_threads > 1U) {
    context->pool = sb_thread_pool_create(context->num_threads - 1U);
    if (!context->pool) {
      context->num_threads = 1U;
    }
  } else {
    context->num_threads = 1U;
  }

  context->frame_buffer =
      (float**)calloc(context->time_buffer_size, sizeof(float*));
  if (!context->frame_buffer) {
    return false;
  }

  for (uint32_t i = 0; i < context->time_buffer_size; i++) {
    context->frame_buffer[i] =
        (float*)calloc(context->spectrum_size, sizeof(float));
    if (!context->frame_buffer[i]) {
      return false;
    }
  }

  const uint32_t total_time_span = context->search_range_time_past +
                                   context->search_range_time_future + 1U +
                                   (2U * NLM_HALO_FRAMES);
  context->frame_ptrs = (float**)calloc(total_time_span, sizeof(float*));
  return context->frame_ptrs != NULL;
}

void patch_filter_context_free(PatchFilterContext* context) {
  if (!context) {
    return;
  }
  if (context->frame_buffer) {
    for (uint32_t i = 0; i < context->time_buffer_size; i++) {
      free(context->frame_buffer[i]);
    }
    free((void*)context->frame_buffer);
  }
  free((void*)context->frame_ptrs);
  sb_thread_pool_free(context->pool);
  context->pool = NULL;
}

void patch_filter_context_push_frame(PatchFilterContext* context,
                                     const float* frame) {
  memcpy(context->frame_buffer[context->buffer_head], frame,
         context->spectrum_size * sizeof(float));
  context->buffer_head =
      (context->buffer_head + 1U) % context->time_buffer_size;
  if (context->frames_filled < context->time_buffer_size) {
    context->frames_filled++;
  }
}

bool patch_filter_context_is_ready(const PatchFilterContext* context) {
  return context && context->frames_filled >= context->time_buffer_size;
}

void patch_filter_context_reset(PatchFilterContext* context) {
  if (!context) {
    return;
  }
  for (uint32_t i = 0; i < context->time_buffer_size; i++) {
    memset(context->frame_buffer[i], 0, context->spectrum_size * sizeof(float));
  }
  context->buffer_head = 0;
  context->frames_filled = 0;
}

uint32_t patch_filter_context_get_latency_frames(
    const PatchFilterContext* context) {
  return context ? context->search_range_time_future : 0U;
}
