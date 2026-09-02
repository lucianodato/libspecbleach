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

#include "fir_history_convolver.h"
#include "../utils/simd_utils.h"
#include <stdlib.h>
#include <string.h>

struct SbFirHistoryConvolver {
  // Linear history: index 0 is the newest sample, index taps-1 the oldest.
  // Contiguous layout keeps the dot product SIMD-friendly.
  float* history;
  uint32_t taps;
};

SbFirHistoryConvolver* fir_history_convolver_initialize(const uint32_t taps) {
  if (taps == 0U) {
    return NULL;
  }

  SbFirHistoryConvolver* self =
      (SbFirHistoryConvolver*)calloc(1U, sizeof(SbFirHistoryConvolver));
  if (!self) {
    return NULL;
  }

  self->taps = taps;
  self->history = (float*)calloc(taps, sizeof(float));
  if (!self->history) {
    free(self);
    return NULL;
  }

  return self;
}

void fir_history_convolver_free(SbFirHistoryConvolver* self) {
  if (!self) {
    return;
  }
  free(self->history);
  free(self);
}

void fir_history_convolver_process(SbFirHistoryConvolver* self,
                                   const uint32_t n, const float* input,
                                   float* output, const float* coeffs) {
  if (!self || n == 0U || !input || !output || !coeffs) {
    return;
  }

  const uint32_t taps = self->taps;
  float* history = self->history;

  sb_simd_state_t simd_state = sb_simd_enable_ftz_daz();

  for (uint32_t i = 0U; i < n; i++) {
    // Advance history by one sample: newest at index 0
    memmove(history + 1U, history, (taps - 1U) * sizeof(float));
    history[0] = input[i];

    // y[i] = sum_m coeffs[m] * history[m]
    sb_vec8_t acc = sb_set8(0.0F);
    float tail_acc = 0.F;
    uint32_t m = 0U;
    for (; m + 7U < taps; m += 8U) {
      acc = sb_add8(acc, sb_mul8(sb_load8(coeffs + m), sb_load8(history + m)));
    }
    for (; m < taps; m++) {
      tail_acc += coeffs[m] * history[m];
    }

    float lanes[8];
    sb_store8(lanes, acc);
    output[i] = lanes[0] + lanes[1] + lanes[2] + lanes[3] + lanes[4] +
                lanes[5] + lanes[6] + lanes[7] + tail_acc;
  }

  sb_simd_restore_state(simd_state);
}
