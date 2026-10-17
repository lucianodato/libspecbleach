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

#ifndef SHARED_DENOISER_LOGIC_MASKING_VETO_H
#define SHARED_DENOISER_LOGIC_MASKING_VETO_H

#include "shared/utils/critical_bands.h"
#include "shared/utils/spectral_features.h"
#include <stdbool.h>
#include <stdint.h>

/**
 * @file masking_veto.h
 * @brief Psychoacoustic audibility floor for spectral denoising.
 *
 * This module computes a per-bin gain floor from a psychoacoustic masking
 * analysis: suppression is relieved where the denoised residual would fall
 * below the masking threshold (frequency spreading + temporal integration +
 * absolute-hearing floor). Pushing reduction deeper there is imperceptible
 * and only breeds musical noise.
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
                                     CriticalBandType critical_band_type,
                                     SpectrumType spectrum_type,
                                     bool use_absolute_threshold,
                                     bool use_temporal_masking);

/**
 * Free a MaskingVeto instance.
 * @param self Instance to free
 */
void masking_veto_free(MaskingVeto* self);

/**
 * Lift the suppression gain toward the psychoacoustic audibility floor.
 *
 * The masking threshold field is computed from a gated clean-signal
 * estimate (future spectrum adds backward masking), band NMR is converted
 * into an audibility gain floor (the gain at which the residual lands
 * exactly on the masking threshold) and the gain is lerped toward it:
 *
 *   gain = max(gain, gain + depth * (floor - gain))
 *
 * Where the original noise is already masked the relief floor saturates at
 * the fixed relief ceiling, so reduction is reduced but never fully
 * disabled ("soft suppression where noise is subjectively inaudible").
 * Relief only ever raises the gain, so the
 * maximum residual per band is bounded by the masked-audibility level.
 *
 * @param self MaskingVeto instance
 * @param spectrum Smoothed spectral power/magnitude of the masked frame
 * @param noise_spectrum Estimated noise profile (broadband)
 * @param future_spectrum Future frame for backward masking (or NULL)
 * @param gain In/Out suppression gains to be relieved
 * @param depth Veto strength (0.0: no relief, 1.0: full audibility floor)
 */
void masking_veto_apply_floor(MaskingVeto* self, const float* spectrum,
                              const float* noise_spectrum,
                              const float* future_spectrum, float* gain,
                              float depth);

/* Init-time only: retune smoothing + true-hop decays. */
void masking_veto_set_hop_sec(MaskingVeto* self, float hop_sec);

#endif // SHARED_DENOISER_LOGIC_MASKING_VETO_H
