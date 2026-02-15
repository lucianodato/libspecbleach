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

#ifndef MASKING_ESTIMATOR_H
#define MASKING_ESTIMATOR_H

#include "shared/utils/critical_bands.h"
#include "shared/utils/spectral_features.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct MaskingEstimator MaskingEstimator;

MaskingEstimator* masking_estimation_initialize(
    uint32_t fft_size, uint32_t sample_rate,
    CriticalBandType critical_band_type, SpectrumType spectrum_type);
void masking_estimation_free(MaskingEstimator* self);
bool compute_masking_thresholds(MaskingEstimator* self, const float* spectrum,
                                const float* future_spectrum,
                                float* masking_thresholds);
void masking_estimation_set_use_absolute_threshold(MaskingEstimator* self,
                                                   bool use_absolute_threshold);

/**
 * masking_estimation_set_temporal_masking: Toggle forward/backward masking
 * logic. Default: true (Enhanced Mode). Set to false for Pure Johnston (1988)
 * compliance.
 */
void masking_estimation_set_temporal_masking(MaskingEstimator* self,
                                             bool enabled);

/**
 * masking_estimation_set_spectral_additivity_exponent: Set the power-law
 * exponent for simultaneous masking (addition of spreading). Standard: 1.0
 * (Johnston 1988 / Linear Addition) Advanced: 0.4 (PEAQ / ITU-R BS.1387
 * Precision)
 */
void masking_estimation_set_spectral_additivity_exponent(MaskingEstimator* self,
                                                         float exponent);

#endif
