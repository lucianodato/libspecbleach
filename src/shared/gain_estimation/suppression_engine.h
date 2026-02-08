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

#ifndef SUPPRESSION_ENGINE_H
#define SUPPRESSION_ENGINE_H

#include <stdint.h>

/**
 * @file suppression_engine.h
 * @brief SNR-dependent oversubtraction (Alpha) calculation engine.
 *
 * This module implements a Berouti-style oversubtraction factor calculation
 * where the suppression aggressiveness (Alpha) is determined per-bin based on
 * the local signal-to-noise ratio. This allows for aggressive cleaning of
 * noisy components while preserving high-SNR signal details.
 */

typedef struct SuppressionEngine SuppressionEngine;

/**
 * Initialize a SuppressionEngine instance.
 * @param real_spectrum_size Size of the real spectrum ((fft_size/2)+1)
 * @return Initialized instance or NULL on failure
 */
SuppressionEngine* suppression_engine_initialize(uint32_t real_spectrum_size);

/**
 * Free a SuppressionEngine instance.
 * @param self Instance to free
 */
void suppression_engine_free(SuppressionEngine* self);

/**
 * Calculate per-bin oversubtraction (alpha) and undersubtraction (beta)
 * factors.
 *
 * @param self SuppressionEngine instance
 * @param reference_spectrum Current spectral magnitude/power
 * @param noise_spectrum Estimated noise profile
 * @param suppression_strength Global aggressiveness factor (0.0 to 1.0)
 * @param alpha [Out] Calculated oversubtraction factors
 * @param beta [Out] Calculated undersubtraction factors
 */
void suppression_engine_calculate(SuppressionEngine* self,
                                  const float* reference_spectrum,
                                  const float* noise_spectrum,
                                  float suppression_strength, float* alpha,
                                  float* beta);

#endif
