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

#ifndef TRANSIENT_DETECTOR_H
#define TRANSIENT_DETECTOR_H

#include <stdbool.h>
#include <stdint.h>

typedef struct TransientDetector TransientDetector;

TransientDetector* transient_detector_initialize(uint32_t num_items);
void transient_detector_free(TransientDetector* self);

/**
 * Global transient detection for a full spectrum/buffer.
 * Used by SpectralSmoother to skip smoothing during onsets.
 */
bool transient_detector_run(TransientDetector* self, const float* spectrum);

/**
 * Process band energies and update transient weights.
 * Used by MaskingVeto for per-band psychoacoustic protection.
 * @param self TransientDetector instance
 * @param band_energies Current energy per critical band
 * @param onset_weights Output buffer for per-band onset weights (0.0:
 * steady, 1.0: transient)
 */
void transient_detector_process(TransientDetector* self,
                                const float* band_energies,
                                float* onset_weights);

#endif
