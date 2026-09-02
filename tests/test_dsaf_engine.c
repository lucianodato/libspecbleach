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

#include "../src/processors/time_denoiser/dsaf_engine.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SAMPLE_RATE 48000U
#define BLOCK 64U
#define TOTAL_BLOCKS 400U

static void test_null_safety(void) {
  assert(dsaf_engine_initialize(0U) == NULL);
  assert(dsaf_engine_get_latency(NULL) == 0U);

  DsafParameters params = {0};
  dsaf_engine_load_parameters(NULL, &params);
  dsaf_engine_load_parameters(NULL, NULL);

  float input[BLOCK] = {0};
  float output[BLOCK] = {0};
  assert(!dsaf_engine_process(NULL, BLOCK, input, output));

  dsaf_engine_free(NULL);
}

static void test_silence_and_finite_output(void) {
  SbDsafEngine* engine = dsaf_engine_initialize(SAMPLE_RATE);
  assert(engine != NULL);

  assert(dsaf_engine_get_latency(engine) == 0U);

  float input[BLOCK];
  float output[BLOCK];

  // Feed silence: output must remain silent and finite while the worker
  // updates coefficients concurrently
  for (uint32_t b = 0U; b < TOTAL_BLOCKS; b++) {
    for (uint32_t i = 0U; i < BLOCK; i++) {
      input[i] = 0.F;
    }
    assert(dsaf_engine_process(engine, BLOCK, input, output));
    for (uint32_t i = 0U; i < BLOCK; i++) {
      assert(isfinite(output[i]));
      assert(fabsf(output[i]) < 1e-6F);
    }
  }

  // Feed noise: output must remain finite with the worker running
  srand(42U);
  for (uint32_t b = 0U; b < TOTAL_BLOCKS; b++) {
    for (uint32_t i = 0U; i < BLOCK; i++) {
      input[i] = ((float)rand() / (float)RAND_MAX) * 2.0F - 1.0F;
    }
    assert(dsaf_engine_process(engine, BLOCK, input, output));
    for (uint32_t i = 0U; i < BLOCK; i++) {
      assert(isfinite(output[i]));
    }
  }

  dsaf_engine_free(engine);
}

static void test_parameter_reactivity(void) {
  SbDsafEngine* engine = dsaf_engine_initialize(SAMPLE_RATE);
  assert(engine != NULL);

  DsafParameters params;
  params.reduction_gain = 0.063F; // ~-24 dB
  params.smoothing_factor = 0.5F;
  params.adaptive_noise = true;
  params.noise_estimation_method = SPP_MMSE_METHOD;
  params.suppression_strength = 0.5F;

  dsaf_engine_load_parameters(engine, &params);

  // Out-of-range values are clamped without breaking processing
  params.reduction_gain = 42.0F;
  params.smoothing_factor = -3.0F;
  params.suppression_strength = 100.0F;
  dsaf_engine_load_parameters(engine, &params);

  float input[BLOCK];
  float output[BLOCK];
  for (uint32_t i = 0U; i < BLOCK; i++) {
    input[i] = ((float)(i % 7) - 3.0F) * 0.1F;
  }
  for (uint32_t b = 0U; b < 100U; b++) {
    assert(dsaf_engine_process(engine, BLOCK, input, output));
    for (uint32_t i = 0U; i < BLOCK; i++) {
      assert(isfinite(output[i]));
    }
  }

  dsaf_engine_free(engine);
}

int main(void) {
  test_null_safety();
  test_silence_and_finite_output();
  test_parameter_reactivity();

  printf("✓ DSAF engine tests passed\n");
  return 0;
}
