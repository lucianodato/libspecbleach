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
static inline __attribute__((unused)) float* cached_get_frame(
    NlmFilter* self, int32_t dt) {
  return self->frame_ptrs[(int32_t)self->config.search_range_time_past + 4 +
                          dt];
}

// Generic implementation (SSE/NEON/Scalar)
bool nlm_filter_process_generic(NlmFilter* filter, float* smoothed_snr);

// AVX implementation
bool nlm_filter_process_avx(NlmFilter* filter, float* smoothed_snr);

#endif /* NLM_FILTER_INTERNAL_H */
