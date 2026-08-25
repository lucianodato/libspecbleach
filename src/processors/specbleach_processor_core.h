/*
libspecbleach - A spectral processing library

Copyright 2022-2026 Luciano Dato <lucianodato@gmail.com>

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

#ifndef SPECBLEACH_PROCESSOR_CORE_H
#define SPECBLEACH_PROCESSOR_CORE_H

#include "processors/denoiser/spectral_denoiser.h"
#include "processors/denoiser2d/spectral_2d_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/stft/stft_processor.h"
#include "specbleach_2d_denoiser.h"
#include "specbleach_denoiser.h"
#include <stdbool.h>
#include <stdint.h>

DenoiserParameters sb_denoiser_params_sanitize(
    const SpecbleachDenoiserParameters* parameters);

Denoiser2DParameters sb_denoiser_2d_params_sanitize(
    const Specbleach2DDenoiserParameters* parameters);

/**
 * Copies a caller-provided reduction curve bias into an internally owned
 * buffer with capacity reuse. Returns NULL when the curve is disabled, the
 * source is NULL, or required_size is zero; otherwise returns the owned
 * copy, which stays valid until the next load or until the instance is
 * freed.
 *
 * May allocate when required_size grows beyond current capacity; reuses the
 * existing buffer otherwise. Callers own *buffer and must free it during
 * instance teardown.
 */
const float* sb_curve_bias_copy(float** buffer, uint32_t* capacity,
                                uint32_t required_size, bool enabled,
                                const float* source);

typedef struct SbProcessorCore {
  uint32_t sample_rate;
  StftProcessor* stft_processor;
  NoiseProfile* noise_profile;
} SbProcessorCore;

SbProcessorCore* sb_processor_core_initialize(
    uint32_t sample_rate, float frame_size, uint32_t overlap_factor,
    ZeroPaddingType padding_type, uint32_t zeropadding_amount,
    WindowTypes input_window, WindowTypes output_window,
    uint32_t profile_spectrum_size);

void sb_processor_core_free(SbProcessorCore* core);

uint32_t sb_processor_core_get_noise_profile_size(SbProcessorCore* core);

bool sb_processor_core_load_noise_profile_for_mode(
    SbProcessorCore* core, const float* restored_profile, uint32_t profile_size,
    uint32_t block_count, int mode);

bool sb_processor_core_reset_noise_profile(SbProcessorCore* core);

uint32_t sb_processor_core_get_noise_profile_block_count_for_mode(
    SbProcessorCore* core, int mode);

float* sb_processor_core_get_noise_profile_for_mode(SbProcessorCore* core,
                                                    int mode);

bool sb_processor_core_noise_profile_available_for_mode(SbProcessorCore* core,
                                                        int mode);

#endif /* SPECBLEACH_PROCESSOR_CORE_H */
