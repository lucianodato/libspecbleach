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

#include "nlm_filter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Default NLM parameters - optimized for real-time performance
// (Reduced from Lukin's offline parameters for real-time use)
#define DEFAULT_PATCH_SIZE 4U
#define DEFAULT_PASTE_BLOCK_SIZE 2U
#define DEFAULT_SEARCH_RANGE_FREQ 4U
#define DEFAULT_SEARCH_RANGE_TIME_PAST 4U
#define DEFAULT_SEARCH_RANGE_TIME_FUTURE 1U
#define DEFAULT_H_PARAMETER 1.0F
#define MIN_WEIGHT 1e-10F

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
  float distance_threshold_actual;

  // Scratch buffer for processing (avoid realloc)
  float* weight_accum;
};

// Actually, let's use a simpler polynomial approximation for 0..10 range which
// is most relevant exp(-x) approx 1 / (1 + x + x^2/2)
static inline float fast_exp_neg(float x) {
  // x is distance/h^2, always positive
  // if x > 10, weight is negligible (< 4.5e-5)
  if (x > 10.0f)
    return 0.0f;
  return expf(-x); // Fallback to std lib for now, correctness first
}

// Helper: clamp index to valid range
static inline uint32_t clamp_index(int32_t idx, uint32_t max_val) {
  if (idx < 0) {
    return 0;
  }
  if ((uint32_t)idx >= max_val) {
    return max_val - 1;
  }
  return (uint32_t)idx;
}

// Helper: get frame from ring buffer (handles wrap-around)
static inline float* get_frame(NlmFilter* self, int32_t relative_offset) {
  int32_t idx = (int32_t)self->buffer_head -
                (int32_t)self->config.search_range_time_past + relative_offset;

  // Wrap around in ring buffer
  while (idx < 0) {
    idx += (int32_t)self->config.time_buffer_size;
  }
  idx %= (int32_t)self->config.time_buffer_size;

  return self->frame_buffer[idx];
}

// Helper: compute squared Euclidean distance between two patches
static float compute_patch_distance(NlmFilter* self, int32_t target_time,
                                    uint32_t target_freq,
                                    int32_t candidate_time,
                                    uint32_t candidate_freq) {
  float distance = 0.0F;
  const uint32_t patch_size = self->config.patch_size;
  const uint32_t half_patch = patch_size / 2;
  const uint32_t spectrum_size = self->config.spectrum_size;

  for (uint32_t dt = 0; dt < patch_size; dt++) {
    int32_t t_target = target_time + (int32_t)dt - (int32_t)half_patch;
    int32_t t_cand = candidate_time + (int32_t)dt - (int32_t)half_patch;

    // Clamp time indices to buffer bounds
    t_target = (int32_t)clamp_index(
        t_target + (int32_t)self->config.search_range_time_past,
        self->config.time_buffer_size);
    t_cand = (int32_t)clamp_index(
        t_cand + (int32_t)self->config.search_range_time_past,
        self->config.time_buffer_size);

    float* target_frame = get_frame(self, t_target);
    float* cand_frame = get_frame(self, t_cand);

    for (uint32_t df = 0; df < patch_size; df++) {
      uint32_t f_target =
          clamp_index((int32_t)target_freq + (int32_t)df - (int32_t)half_patch,
                      spectrum_size);
      uint32_t f_cand = clamp_index(
          (int32_t)candidate_freq + (int32_t)df - (int32_t)half_patch,
          spectrum_size);

      float diff = target_frame[f_target] - cand_frame[f_cand];
      distance += diff * diff;
    }
  }

  return distance;
}

NlmFilter* nlm_filter_initialize(NlmFilterConfig config) {
  // Validate configuration
  if (config.spectrum_size == 0) {
    return NULL;
  }

  NlmFilter* self = (NlmFilter*)calloc(1U, sizeof(NlmFilter));
  if (!self) {
    return NULL;
  }

  // Apply defaults for zero values
  self->config = config;
  if (self->config.patch_size == 0) {
    self->config.patch_size = DEFAULT_PATCH_SIZE;
  }
  if (self->config.paste_block_size == 0) {
    self->config.paste_block_size = DEFAULT_PASTE_BLOCK_SIZE;
  }
  if (self->config.search_range_freq == 0) {
    self->config.search_range_freq = DEFAULT_SEARCH_RANGE_FREQ;
  }
  if (self->config.search_range_time_past == 0) {
    self->config.search_range_time_past = DEFAULT_SEARCH_RANGE_TIME_PAST;
  }
  if (self->config.search_range_time_future == 0) {
    self->config.search_range_time_future = DEFAULT_SEARCH_RANGE_TIME_FUTURE;
  }
  if (self->config.h_parameter <= 0.0F) {
    self->config.h_parameter = DEFAULT_H_PARAMETER;
  }

  // Compute time buffer size from search ranges
  if (self->config.time_buffer_size == 0) {
    self->config.time_buffer_size = self->config.search_range_time_past +
                                    self->config.search_range_time_future + 1;
  }

  // Precompute values
  self->h_squared = self->config.h_parameter * self->config.h_parameter;
  if (self->config.distance_threshold <= 0.0F) {
    // Default: 4 * h² (patches 2× noise variance away contribute < 2%)
    self->distance_threshold_actual = 4.0F * self->h_squared;
  } else {
    self->distance_threshold_actual = self->config.distance_threshold;
  }

  // Target frame is at the "past" offset to allow look-ahead
  self->target_frame_offset = self->config.search_range_time_past;

  // Allocate ring buffer for frames
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

  // Allocate scratch buffer
  self->weight_accum =
      (float*)calloc(self->config.spectrum_size, sizeof(float));
  if (!self->weight_accum) {
    nlm_filter_free(self);
    return NULL;
  }

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

  free(filter);
}

void nlm_filter_set_h_parameter(NlmFilter* filter, float h) {
  if (!filter || h <= 0.0F) {
    return;
  }

  filter->config.h_parameter = h;
  filter->h_squared = h * h;

  // Update distance threshold if using default
  if (filter->config.distance_threshold <= 0.0F) {
    filter->distance_threshold_actual = 4.0F * filter->h_squared;
  }
}

void nlm_filter_push_frame(NlmFilter* filter, const float* snr_frame) {
  if (!filter || !snr_frame) {
    return;
  }

  // Copy frame to ring buffer
  memcpy(filter->frame_buffer[filter->buffer_head], snr_frame,
         filter->config.spectrum_size * sizeof(float));

  // Advance head
  filter->buffer_head =
      (filter->buffer_head + 1) % filter->config.time_buffer_size;

  // Track fill level
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

  // Need full buffer for processing
  if (!nlm_filter_is_ready(filter)) {
    return false;
  }

  const uint32_t spectrum_size = filter->config.spectrum_size;
  const uint32_t paste_size = filter->config.paste_block_size;
  const uint32_t search_freq = filter->config.search_range_freq;
  const int32_t search_time_past =
      (int32_t)filter->config.search_range_time_past;
  const int32_t search_time_future =
      (int32_t)filter->config.search_range_time_future;

  // Get target frame (center of buffer)
  float* target_frame = get_frame(filter, 0);

  // Initialize output with zeros
  memset(smoothed_snr, 0, spectrum_size * sizeof(float));

  // Use preallocated buffer
  float* weight_sum = filter->weight_accum;
  // Clear weight buffer for new frame
  memset(weight_sum, 0, spectrum_size * sizeof(float));

  // Process using pasted block optimization
  // For each paste block position in frequency
  for (uint32_t block_start = 0; block_start < spectrum_size;
       block_start += paste_size) {

    uint32_t block_center = block_start + (paste_size / 2);
    if (block_center >= spectrum_size) {
      block_center = spectrum_size - 1;
    }

    // Search over time and frequency window
    for (int32_t dt = -search_time_past; dt <= search_time_future; dt++) {
      // Optimization: inner loop over frequency often has more locality
      for (int32_t df = -(int32_t)search_freq; df <= (int32_t)search_freq;
           df++) {

        uint32_t cand_center =
            clamp_index((int32_t)block_center + df, spectrum_size);

        // Compute patch distance
        float distance =
            compute_patch_distance(filter, 0, block_center, dt, cand_center);

        // Distance thresholding (early termination optimization)
        if (distance > filter->distance_threshold_actual) {
          continue;
        }

        // Compute weight
        float weight = fast_exp_neg(distance / filter->h_squared);
        if (weight < MIN_WEIGHT) {
          continue;
        }

        // Get candidate frame
        float* cand_frame = get_frame(filter, dt);

        // Apply weight to all bins in the paste block
        // Unroll loop for small fixed paste size (2 or 4)
        for (uint32_t i = 0;
             i < paste_size && (block_start + i) < spectrum_size; i++) {
          uint32_t target_bin = block_start + i;
          uint32_t cand_bin =
              clamp_index((int32_t)target_bin + df, spectrum_size);

          smoothed_snr[target_bin] += weight * cand_frame[cand_bin];
          weight_sum[target_bin] += weight;
        }
      }
    }
  }

  // Normalize by weight sum
  for (uint32_t k = 0; k < spectrum_size; k++) {
    if (weight_sum[k] > MIN_WEIGHT) {
      smoothed_snr[k] /= weight_sum[k];
    } else {
      // Fallback to original value if no valid weights
      smoothed_snr[k] = target_frame[k];
    }
  }

  // Normalize by weight sum - no free(weight_sum) needed

  return true;
}

void nlm_filter_reset(NlmFilter* filter) {
  if (!filter) {
    return;
  }

  // Clear all frames
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
  // Latency is the number of look-ahead frames
  return filter->config.search_range_time_future;
}
