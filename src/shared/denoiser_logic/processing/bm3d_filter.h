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

#ifndef BM3D_FILTER_H
#define BM3D_FILTER_H

#include <stdbool.h>
#include <stdint.h>

typedef struct Bm3dFilter Bm3dFilter;

/**
 * Configuration for the 1-pass BM3D-lite collaborative filter.
 * Mirrors NlmFilterConfig geometry so both filters share window and latency.
 */
typedef struct Bm3dFilterConfig {
  uint32_t spectrum_size;          /**< Number of frequency bins */
  uint32_t time_buffer_size;       /**< Frames to buffer (past+future+1) */
  uint32_t patch_size;             /**< Patch size for matching */
  uint32_t paste_block_size;       /**< Bins per collaborative block */
  uint32_t search_range_freq;      /**< Search range in frequency bins */
  uint32_t search_range_time_past; /**< Search range in time - past frames */
  uint32_t
      search_range_time_future; /**< Search range in time - future frames */
  float h_parameter;            /**< Smoothing parameter (shared with NLM) */
  uint32_t num_threads;         /**< Total threads (0 = default) */
} Bm3dFilterConfig;

Bm3dFilter* bm3d_filter_initialize(Bm3dFilterConfig config);
void bm3d_filter_free(Bm3dFilter* filter);
void bm3d_filter_set_h_parameter(Bm3dFilter* filter, float h);
void bm3d_filter_push_frame(Bm3dFilter* filter, const float* snr_frame);
bool bm3d_filter_process(Bm3dFilter* filter, float* smoothed_snr);
bool bm3d_filter_is_ready(Bm3dFilter* filter);
void bm3d_filter_reset(Bm3dFilter* filter);
uint32_t bm3d_filter_get_latency_frames(Bm3dFilter* filter);
void bm3d_filter_calculate_snr(Bm3dFilter* filter,
                               const float* reference_spectrum,
                               const float* noise_spectrum, float* snr_frame);
void bm3d_filter_reconstruct_magnitude(Bm3dFilter* filter,
                                       const float* smoothed_snr,
                                       const float* noise_spectrum,
                                       float* magnitude_spectrum);

#endif /* BM3D_FILTER_H */
