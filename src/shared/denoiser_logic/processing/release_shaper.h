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

#ifndef RELEASE_SHAPER_H
#define RELEASE_SHAPER_H

#include <stdbool.h>
#include <stdint.h>

typedef struct ReleaseShaper ReleaseShaper;

/**
 * Initialize the release shaper.
 *
 * Computes a per-bin release scale (0.0 = fast close, 1.0 = full slider
 * release) from band-level energy collapse evidence: when a band's
 * recent-energy envelope far exceeds its current energy, the signal in that
 * band has ended and its bins must close fast to avoid ghosting the signal's
 * spectral shape onto the residual noise. Random per-bin fluctuations never
 * move a band aggregate, so the full release is preserved in noise-only
 * stretches (anti-chirp protection).
 *
 * @param sample_rate Audio sample rate in Hz
 * @param fft_size FFT size used for analysis
 */
ReleaseShaper* release_shaper_initialize(uint32_t sample_rate,
                                         uint32_t fft_size);
void release_shaper_free(ReleaseShaper* self);
void release_shaper_reset(ReleaseShaper* self);

/* Init-time only: true STFT hop in seconds (0 = derive from fft/overlap). */
void release_shaper_set_hop_sec(ReleaseShaper* self, float hop_sec);

/**
 * Compute the per-bin release scale for the current frame.
 *
 * @param self ReleaseShaper instance
 * @param magnitude Current (pre-subtraction stabilized) magnitude spectrum
 * @param release_scale [Out] Per-bin scale in 0.0..1.0
 *                       (real_spectrum_size entries)
 * @return true on success
 */
bool release_shaper_compute(ReleaseShaper* self, const float* magnitude,
                            float* release_scale);

#endif
