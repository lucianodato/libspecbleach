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

#include "../src/shared/time_domain/fir_history_convolver.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TAPS 256U
#define BLOCK 64U
#define EPSILON 1e-3F

static void test_delta_coefficient_zero_latency(void) {
  SbFirHistoryConvolver* conv = fir_history_convolver_initialize(TAPS);
  assert(conv != NULL);

  float coeffs[TAPS];
  memset(coeffs, 0, sizeof(coeffs));
  coeffs[0] = 1.0F; // delta filter: y[n] = x[n]

  float input[BLOCK];
  float output[BLOCK];
  for (uint32_t i = 0U; i < BLOCK; i++) {
    input[i] = (float)(i + 1);
  }

  fir_history_convolver_process(conv, BLOCK, input, output, coeffs);
  for (uint32_t i = 0U; i < BLOCK; i++) {
    assert(fabsf(output[i] - input[i]) < EPSILON);
  }

  fir_history_convolver_free(conv);
}

static void test_matches_direct_convolution(void) {
  SbFirHistoryConvolver* conv = fir_history_convolver_initialize(TAPS);
  assert(conv != NULL);

  float coeffs[TAPS];
  for (uint32_t m = 0U; m < TAPS; m++) {
    coeffs[m] = 0.5F / (float)(m + 1U);
  }

  // Reference: maintain a scalar history over the whole signal
  const uint32_t total = 4U * BLOCK;
  float* input = (float*)calloc(total, sizeof(float));
  float* reference = (float*)calloc(total, sizeof(float));
  float* history = (float*)calloc(TAPS, sizeof(float)); // history[0] = newest
  for (uint32_t i = 0U; i < total; i++) {
    input[i] = sinf(0.05F * (float)i) + 0.25F * cosf(0.31F * (float)i);
    memmove(history + 1, history, (TAPS - 1) * sizeof(float));
    history[0] = input[i];
    for (uint32_t m = 0U; m < TAPS; m++) {
      reference[i] += coeffs[m] * history[m];
    }
  }

  // Run in blocks and compare
  for (uint32_t start = 0U; start < total; start += BLOCK) {
    float output[BLOCK];
    fir_history_convolver_process(conv, BLOCK, input + start, output, coeffs);
    for (uint32_t i = 0U; i < BLOCK; i++) {
      assert(fabsf(output[i] - reference[start + i]) < EPSILON);
    }
  }

  free(input);
  free(reference);
  free(history);
  fir_history_convolver_free(conv);
}

static void test_history_persistence_across_blocks(void) {
  // A single impulse in an earlier block must appear delayed in a later one
  float coeffs[TAPS];
  memset(coeffs, 0, sizeof(coeffs));
  coeffs[TAPS - 1] = 1.0F; // pure delay of TAPS-1 samples

  SbFirHistoryConvolver* conv = fir_history_convolver_initialize(TAPS);
  assert(conv != NULL);

  float zeros[BLOCK] = {0};
  float output[BLOCK];

  float impulse[BLOCK] = {0};
  impulse[0] = 1.0F;
  fir_history_convolver_process(conv, BLOCK, impulse, output, coeffs);
  for (uint32_t b = 1U; b < (TAPS / BLOCK); b++) {
    fir_history_convolver_process(conv, BLOCK, zeros, output, coeffs);
    if (b == (TAPS / BLOCK) - 1U) {
      // Impulse delay = TAPS - 1 samples = inside this block at index
      // (TAPS - 1) - b*BLOCK
      const uint32_t expected_pos = (TAPS - 1U) - b * BLOCK;
      assert(fabsf(output[expected_pos] - 1.0F) < EPSILON);
    }
  }

  fir_history_convolver_free(conv);
}

static void test_invalid_arguments(void) {
  assert(fir_history_convolver_initialize(0U) == NULL);

  SbFirHistoryConvolver* conv = fir_history_convolver_initialize(TAPS);
  assert(conv != NULL);

  float input[BLOCK] = {0};
  float output[BLOCK] = {0};
  float coeffs[TAPS] = {0};

  fir_history_convolver_process(NULL, BLOCK, input, output, coeffs);
  fir_history_convolver_process(conv, 0U, input, output, coeffs);
  fir_history_convolver_process(conv, BLOCK, NULL, output, coeffs);
  fir_history_convolver_process(conv, BLOCK, input, NULL, coeffs);
  fir_history_convolver_process(conv, BLOCK, input, output, NULL);

  fir_history_convolver_free(conv);
  fir_history_convolver_free(NULL);
}

int main(void) {
  test_delta_coefficient_zero_latency();
  test_matches_direct_convolution();
  test_history_persistence_across_blocks();
  test_invalid_arguments();

  printf("✓ FIR history convolver tests passed\n");
  return 0;
}
