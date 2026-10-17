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
    CriticalBandType critical_band_type, SpectrumType spectrum_type,
    bool use_absolute_threshold, bool use_temporal_masking);
void masking_estimation_free(MaskingEstimator* self);
/* Compute per-bin masking thresholds from the estimated masker spectrum.
 *
 * @param future_spectrum Optional lookahead masker spectrum enabling
 *   backward masking (NULL for latency-free chains).
 * @param band_masker_mask Optional per-band 0/1 mask (critical-bands length)
 *   selecting bands whose masker evidence is trusted by the caller. Bands
 *   masked 0 contribute nothing to the per-frame frequency threshold, so
 *   the model's temporal memory is fed only by validated maskers (a biased
 *   clean estimate can otherwise accumulate phantom persistence in
 *   noise-only sections). NULL trusts every band.
 * @param masking_thresholds [Out] per-bin thresholds. */
bool compute_masking_thresholds(MaskingEstimator* self, const float* spectrum,
                                const float* future_spectrum,
                                const float* band_masker_mask,
                                float* masking_thresholds);

/* Read-only absolute-hearing-threshold spectrum (real_spectrum_size bins);
 * NULL when the absolute threshold is disabled. Valid while the instance
 * lives. */
const float* masking_estimation_get_absolute_thresholds(
    const MaskingEstimator* self);

/* Init-time only: recompute decays from the true hop (frame/overlap/sr). */
void masking_estimation_set_hop_sec(MaskingEstimator* self, float hop_sec);

/* Read-only per-band temporal masking claim (in band-power scale; length =
 * number of critical bands). Holds the integrated forward/backward masking
 * components of THIS frame without the fresh frequency estimate: nonzero
 * exactly where the model says a recent or incoming masker still masks.
 * Zeros when temporal masking is disabled. Valid while the instance
 * lives. */
const float* masking_estimation_get_temporal_persistence(
    const MaskingEstimator* self);

#endif
