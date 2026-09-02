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

#ifndef FIR_HISTORY_CONVOLVER_H
#define FIR_HISTORY_CONVOLVER_H

#include <stdint.h>

typedef struct SbFirHistoryConvolver SbFirHistoryConvolver;

SbFirHistoryConvolver* fir_history_convolver_initialize(uint32_t taps);
void fir_history_convolver_free(SbFirHistoryConvolver* self);

/**
 * Zero-latency block convolution against the rolling sample history.
 *
 * RT-safe: no allocations, locks, or I/O. The history buffer is advanced
 * by `n` samples and the output is y[n] = sum_m coeffs[m] * x[n - m],
 * so an impulse passed with a delta-coefficient filter at index 0
 * emerges with zero delay.
 *
 * The coefficients buffer is owned by the caller (may be swapped
 * lock-free between calls) and must hold `taps` floats.
 */
void fir_history_convolver_process(SbFirHistoryConvolver* self, uint32_t n,
                                   const float* input, float* output,
                                   const float* coeffs);

#endif
