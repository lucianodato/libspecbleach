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

#ifndef HPSS_QUALITY_MODE_DEFINED
#define HPSS_QUALITY_MODE_DEFINED
typedef enum HpssQualityMode {
  HPSS_QUALITY_OFF = 0,
  HPSS_QUALITY_LOW = 1,
  HPSS_QUALITY_MEDIUM = 2,
  HPSS_QUALITY_HIGH = 3
} HpssQualityMode;
#endif

typedef struct HpssFilter HpssFilter;

typedef struct HpssConfig {
  uint32_t real_spectrum_size;
} HpssConfig;

HpssFilter* hpss_filter_initialize(HpssConfig config);
void hpss_filter_free(HpssFilter* self);

// Dynamic parameter update (allocator-free)
void hpss_filter_set_quality_mode(HpssFilter* self, HpssQualityMode mode);
void hpss_filter_set_sensitivity(HpssFilter* self, float sensitivity);

// Sliding HPSS is causal with 0 lookahead delay frames
uint32_t hpss_filter_get_latency_frames(const HpssFilter* self);
float hpss_filter_get_onset_ratio(const HpssFilter* self);
bool hpss_filter_is_transient_detected(const HpssFilter* self);

bool hpss_filter_process(HpssFilter* self, const float* current_magnitude,
                         const float* noise_floor, float* delayed_magnitude_out,
                         float* mask_harmonic_out, float* mask_percussive_out);

#endif // SHARED_DENOISER_LOGIC_HPSS_FILTER_H
