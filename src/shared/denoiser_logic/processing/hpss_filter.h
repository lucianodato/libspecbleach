/*
libspecbleach - A spectral processing library

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

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

#ifndef SHARED_DENOISER_LOGIC_HPSS_FILTER_H
#define SHARED_DENOISER_LOGIC_HPSS_FILTER_H

#include "shared/configurations.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct HpssFilter HpssFilter;

typedef struct HpssConfig {
  uint32_t real_spectrum_size;
} HpssConfig;

HpssFilter* hpss_filter_initialize(HpssConfig config);
void hpss_filter_free(HpssFilter* self);

// Dynamic parameter update (allocator-free)
void hpss_filter_set_enabled(HpssFilter* self, bool enabled);

// Sliding HPSS is causal with 0 lookahead delay frames
uint32_t hpss_filter_get_latency_frames(const HpssFilter* self);
float hpss_filter_get_onset_ratio(const HpssFilter* self);
bool hpss_filter_is_transient_detected(const HpssFilter* self);
/**
 * Process a spectral magnitude frame with Sliding HPSS.
 *
 * @param self HPSS filter instance
 * @param current_magnitude Current frame magnitude spectrum
 * [real_spectrum_size]
 * @param noise_floor Optional noise floor spectrum (or NULL to skip sub-noise
 * gate)
 * @param current_magnitude_out Optional destination for current frame magnitude
 * copy (or NULL)
 * @param mask_harmonic_out Optional output array for harmonic soft mask (or
 * NULL)
 * @param mask_percussive_out Optional output array for percussive soft mask (or
 * NULL)
 * @return true on success, false on invalid arguments
 */
bool hpss_filter_process(HpssFilter* self, const float* current_magnitude,
                         const float* noise_floor, float* current_magnitude_out,
                         float* mask_harmonic_out, float* mask_percussive_out);

#endif // SHARED_DENOISER_LOGIC_HPSS_FILTER_H
