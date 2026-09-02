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

#ifndef MINIMUM_PHASE_FIR_H
#define MINIMUM_PHASE_FIR_H

#include <stdbool.h>
#include <stdint.h>

#include "shared/stft/fft_transform.h"

/**
 * Synthesizes an M-tap minimum-phase FIR impulse response from a target
 * magnitude spectrum via homomorphic (cepstral) factorization:
 * log-magnitude -> IFFT -> cepstral fold -> FFT -> exp -> IFFT.
 *
 * Stateless and off-the-real-time-thread: run it from the analysis worker.
 * `fft` must be an FftTransform of exactly `fft_size` bins (owned by the
 * caller, reused across calls). `gain_bins` holds real_spectrum_size
 * magnitude gains. The precomputed `taper_window` (length `taps`,
 * typically flat with a decaying tail) smooths truncation to `taps`.
 *
 * @return true on success, false on NULL/invalid arguments.
 */
bool minimum_phase_fir_synthesize(uint32_t fft_size, uint32_t taps,
                                  const float* gain_bins, FftTransform* fft,
                                  const float* taper_window, float* out_ir);

/**
 * Builds the asymmetric truncation window: flat for the first
 * DSAF_ASYMMETRIC_TAPER_START fraction of taps, then a raised-cosine
 * decay to zero. Precompute once at initialization.
 */
bool minimum_phase_fir_taper_window(uint32_t taps, float* window);

#endif
