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

#ifndef MASKING_VETO_H
#define MASKING_VETO_H

#include "../utils/spectral_features.h"
#include <stdbool.h>
#include <stdint.h>

/**
 * @file masking_veto.h
 * @brief Psychoacoustic veto module for spectral denoising.
 *
 * This module implements a masking veto mechanism that moderates noise
 * reduction based on psychoacoustic masking thresholds. It ensures that
 * noise reduction is only applied where it is perceptually necessary.
 */

typedef struct MaskingVeto MaskingVeto;

/**
 * Initialize a MaskingVeto instance.
 * @param fft_size Size of the FFT
 * @param sample_rate Audio sample rate
 * @param spectrum_type Scaling type (POWER or MAGNITUDE)
 * @return Initialized MaskingVeto instance or NULL on failure
 */
MaskingVeto* masking_veto_initialize(uint32_t fft_size, uint32_t sample_rate,
                                     SpectrumType spectrum_type);

/**
 * Free a MaskingVeto instance.
 * @param self Instance to free
 */
void masking_veto_free(MaskingVeto* self);

/**
 * Moderate noise reduction alpha values based on psychoacoustic masking.
 *
 * @param self MaskingVeto instance
 * @param spectrum Current spectral magnitude or power
 * @param noise_spectrum Estimated noise profile
 * @param alpha In/Out alpha map to be moderated
 * @param depth Veto strength (0.0: No veto, 1.0: Full psychoacoustic
 * protection)
 */
void masking_veto_apply(MaskingVeto* self, const float* smoothed_spectrum,
                        const float* noisy_spectrum,
                        const float* noise_spectrum, float* alpha,
                        float floor_alpha, float depth, float elasticity);

#endif
