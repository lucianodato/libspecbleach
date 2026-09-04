/*
libspecbleach - A spectral processing library

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

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

/* Unit tests for the DFTT-lite post-NLM filter (#155 phase 2). Fast and
 * exact: lifecycle, NULL guards, readiness gating, silence transparency,
 * dual-input transparency, determinism, and reset. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "shared/denoiser_logic/processing/dftt_filter.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define SPEC_SIZE 64U
#define TIME_SPAN 4U
#define BLOCK_FREQ 16U

static void test_lifecycle(void) {
  printf("Testing lifecycle...\n");
  TEST_ASSERT(dftt_filter_initialize(0U, TIME_SPAN, BLOCK_FREQ) == NULL,
              "zero spectrum must fail");
  TEST_ASSERT(dftt_filter_initialize(SPEC_SIZE, 0U, BLOCK_FREQ) == NULL,
              "zero time span must fail");
  TEST_ASSERT(dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, 0U) == NULL,
              "zero block must fail");
  DfttFilter* f = dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, BLOCK_FREQ);
  TEST_ASSERT(f != NULL, "valid config must init");
  TEST_ASSERT(!dftt_filter_is_ready(f), "fresh filter must not be ready");
  dftt_filter_free(f);
  dftt_filter_free(NULL);
  dftt_filter_reset(NULL);
  printf("✓ Lifecycle passed\n");
}

static void test_null_guards(void) {
  printf("Testing NULL guards...\n");
  DfttFilter* f = dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, BLOCK_FREQ);
  float frame[SPEC_SIZE];
  for (uint32_t k = 0U; k < SPEC_SIZE; k++) {
    frame[k] = 1.0F;
  }
  dftt_filter_push(NULL, frame, frame);
  dftt_filter_push(f, NULL, frame);
  dftt_filter_push(f, frame, NULL);
  TEST_ASSERT(!dftt_filter_process(NULL, frame), "NULL filter must fail");
  TEST_ASSERT(!dftt_filter_process(f, NULL), "NULL output must fail");
  TEST_ASSERT(!dftt_filter_process(f, frame), "must fail before ready");
  TEST_ASSERT(!dftt_filter_is_ready(NULL), "NULL must not be ready");
  dftt_filter_free(f);
  printf("✓ NULL guards passed\n");
}

static void push_pair(DfttFilter* f, const float* noisy,
                      const float* smoothed) {
  dftt_filter_push(f, noisy, smoothed);
}

static void test_uniform_response(void) {
  printf("Testing uniform response...\n");
  DfttFilter* f = dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, BLOCK_FREQ);
  float noisy[SPEC_SIZE];
  float smoothed[SPEC_SIZE];
  float out[SPEC_SIZE];
  for (uint32_t k = 0U; k < SPEC_SIZE; k++) {
    noisy[k] = 1.0F;
    smoothed[k] = 1.0F;
  }
  for (uint32_t t = 0U; t < TIME_SPAN; t++) {
    push_pair(f, noisy, smoothed);
  }
  TEST_ASSERT(dftt_filter_is_ready(f), "must be ready after span pushes");
  TEST_ASSERT(dftt_filter_process(f, out), "must process when ready");
  double mean = 0.0;
  for (uint32_t k = 0U; k < SPEC_SIZE; k++) {
    TEST_ASSERT(isfinite(out[k]), "output must be finite");
    mean += out[k];
  }
  mean /= SPEC_SIZE;
  TEST_ASSERT(mean > 0.3 && mean < 1.7, "uniform input must pass through");
  TEST_ASSERT(!dftt_filter_process(f, out), "stale frame must not reprocess");
  dftt_filter_free(f);
  printf("✓ Uniform response passed\n");
}

static void test_silence_transparency(void) {
  printf("Testing silence transparency...\n");
  DfttFilter* f = dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, BLOCK_FREQ);
  float zero[SPEC_SIZE] = {0};
  float out[SPEC_SIZE];
  for (uint32_t k = 0U; k < SPEC_SIZE; k++) {
    out[k] = 5.0F;
  }
  for (uint32_t t = 0U; t < TIME_SPAN; t++) {
    push_pair(f, zero, zero);
  }
  TEST_ASSERT(dftt_filter_process(f, out), "must process silence");
  for (uint32_t k = 0U; k < SPEC_SIZE; k++) {
    TEST_ASSERT(out[k] == 0.0F, "silence must stay silent");
  }
  dftt_filter_free(f);
  printf("✓ Silence transparency passed\n");
}

static void test_empty_reference_veto(void) {
  printf("Testing empty-reference veto...\n");
  /* Veto contract: the kill-map can only remove energy, never add. A zero
   * NLM verdict means zero output, even with a loud noisy tile. */
  DfttFilter* f = dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, BLOCK_FREQ);
  float noisy[SPEC_SIZE] = {0};
  float smoothed[SPEC_SIZE] = {0};
  float out[SPEC_SIZE] = {0};
  noisy[20] = 10.0F;
  for (uint32_t t = 0U; t < TIME_SPAN; t++) {
    push_pair(f, noisy, smoothed);
  }
  TEST_ASSERT(dftt_filter_process(f, out), "must process");
  TEST_ASSERT(fabsf(out[20]) < 1e-2F,
              "veto must never add energy above the NLM verdict");
  dftt_filter_free(f);
  printf("✓ Empty-reference veto passed\n");
}

static void test_veto_never_exceeds_reference(void) {
  printf("Testing veto never exceeds reference...\n");
  /* Hum-creep regression: a suppressed NLM bin (tonal noise during speech)
   * must stay suppressed no matter how loud the noisy tile is there. */
  DfttFilter* f = dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, BLOCK_FREQ);
  float noisy[SPEC_SIZE];
  float smoothed[SPEC_SIZE];
  float out[SPEC_SIZE] = {0};
  for (uint32_t k = 0U; k < SPEC_SIZE; k++) {
    noisy[k] = 2.0F + 0.5F * sinf(0.3F * (float)k);
    smoothed[k] = 1.0F + 0.1F * sinf(0.3F * (float)k);
  }
  smoothed[40] = 0.05F; /* NLM-suppressed hum line */
  noisy[40] = 3.0F;     /* ...that rings loud in the raw map */
  for (uint32_t t = 0U; t < TIME_SPAN; t++) {
    push_pair(f, noisy, smoothed);
  }
  TEST_ASSERT(dftt_filter_process(f, out), "must process");
  TEST_ASSERT(out[40] <= smoothed[40] + 1e-2F,
              "suppressed bin must stay suppressed");
  dftt_filter_free(f);
  printf("✓ Veto never exceeds reference passed\n");
}

static void test_determinism_and_reset(void) {
  printf("Testing determinism and reset...\n");
  DfttFilter* a = dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, BLOCK_FREQ);
  DfttFilter* b = dftt_filter_initialize(SPEC_SIZE, TIME_SPAN, BLOCK_FREQ);
  float noisy[SPEC_SIZE];
  float smoothed[SPEC_SIZE];
  float out_a[SPEC_SIZE];
  float out_b[SPEC_SIZE];
  for (uint32_t k = 0U; k < SPEC_SIZE; k++) {
    noisy[k] = 0.5F + 0.01F * (float)((k * 37U) % 11U);
    smoothed[k] = 0.4F + 0.01F * (float)((k * 53U) % 7U);
  }
  for (uint32_t t = 0U; t < TIME_SPAN + 2U; t++) {
    push_pair(a, noisy, smoothed);
    push_pair(b, noisy, smoothed);
  }
  TEST_ASSERT(dftt_filter_process(a, out_a), "a must process");
  TEST_ASSERT(dftt_filter_process(b, out_b), "b must process");
  for (uint32_t k = 0U; k < SPEC_SIZE; k++) {
    TEST_ASSERT(out_a[k] == out_b[k], "identical feeds must match exactly");
  }
  dftt_filter_reset(a);
  TEST_ASSERT(!dftt_filter_is_ready(a), "reset must clear readiness");
  TEST_ASSERT(!dftt_filter_process(a, out_a), "reset must block process");
  dftt_filter_free(a);
  dftt_filter_free(b);
  printf("✓ Determinism and reset passed\n");
}

int main(void) {
  printf("Running DFTT filter tests...\n\n");
  test_lifecycle();
  test_null_guards();
  test_uniform_response();
  test_silence_transparency();
  test_empty_reference_veto();
  test_veto_never_exceeds_reference();
  test_determinism_and_reset();
  printf("\n✅ All DFTT filter tests passed!\n");
  return 0;
}
