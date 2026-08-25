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

#include "specbleach_transition.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define SAMPLE_RATE 48000
#define MAX_BLOCK 512
#define BLOCK_SIZE 512
#define CHANNELS 2

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      return;                                                                  \
    }                                                                          \
  } while (0)

static void test_init_and_validation(void) {
  printf("Testing transition init and validation...\n");

  TEST_ASSERT(specbleach_transition_initialize(0, MAX_BLOCK, CHANNELS) == NULL,
              "sample_rate == 0 rejected");
  TEST_ASSERT(
      specbleach_transition_initialize(SAMPLE_RATE, 0, CHANNELS) == NULL,
      "max_block_size == 0 rejected");
  TEST_ASSERT(
      specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, 0) == NULL,
      "channels == 0 rejected");

  specbleach_transition* t =
      specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, CHANNELS);
  TEST_ASSERT(t != NULL, "valid init succeeds");
  TEST_ASSERT(!specbleach_transition_active(t), "fresh instance is idle");
  TEST_ASSERT(specbleach_transition_get_latency(t) == 0,
              "latency zero before begin");

  /* Invalid process calls are rejected */
  float src[CHANNELS][BLOCK_SIZE] = {{0}};
  float dst[CHANNELS][BLOCK_SIZE] = {{0}};
  const float* from_ptrs[CHANNELS] = {src[0], src[1]};
  const float* to_ptrs[CHANNELS] = {dst[0], dst[1]};
  float* blended_ptrs[CHANNELS] = {dst[0], dst[1]};
  TEST_ASSERT(!specbleach_transition_process(NULL, BLOCK_SIZE, from_ptrs,
                                             to_ptrs, blended_ptrs),
              "NULL instance rejected");
  TEST_ASSERT(
      !specbleach_transition_process(t, 0, from_ptrs, to_ptrs, blended_ptrs),
      "zero samples rejected");
  TEST_ASSERT(!specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, NULL,
                                             blended_ptrs),
              "NULL to_out rejected");
  TEST_ASSERT(
      !specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs, NULL),
      "NULL blended rejected");

  specbleach_transition_free(t);
  specbleach_transition_free(NULL); /* must be a no-op */
  printf("✓ Transition init and validation tests passed\n");
}

static void test_idle_passthrough(void) {
  printf("Testing idle passthrough...\n");

  specbleach_transition* t =
      specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, CHANNELS);

  float source[CHANNELS][BLOCK_SIZE];
  float target[CHANNELS][BLOCK_SIZE];
  float blended[CHANNELS][BLOCK_SIZE];
  const float* from_ptrs[CHANNELS] = {source[0], source[1]};
  const float* to_ptrs[CHANNELS] = {target[0], target[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    source[0][s] = 0.25f;
    source[1][s] = 0.5f;
    target[0][s] = -0.75f;
    target[1][s] = -0.5f;
  }

  TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                            blended_ptrs),
              "idle passthrough processed");
  TEST_ASSERT(blended[0][10] == target[0][10],
              "idle passes through to_out ch0");
  TEST_ASSERT(blended[1][20] == target[1][20],
              "idle passes through to_out ch1");

  /* Aliased destination must not memcpy over itself (guard against UB). */
  const float* alias_to[CHANNELS] = {blended[0], blended[1]};
  TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, alias_to,
                                            blended_ptrs),
              "aliased passthrough processed");
  TEST_ASSERT(blended[1][0] == target[1][0], "alias passthrough intact");

  /* NULL from_out is allowed while idle. */
  TEST_ASSERT(
      specbleach_transition_process(t, BLOCK_SIZE, NULL, to_ptrs, blended_ptrs),
      "NULL from_out allowed while idle");

  specbleach_transition_free(t);
  printf("✓ Idle passthrough tests passed\n");
}

static void test_equal_power_fade(void) {
  printf("Testing equal-power fade...\n");

  specbleach_transition* t =
      specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, CHANNELS);

  const uint32_t latency_from = 1024;
  const uint32_t latency_to = 2048;

  float source[CHANNELS][BLOCK_SIZE];
  float target[CHANNELS][BLOCK_SIZE];
  float blended[CHANNELS][BLOCK_SIZE];
  const float* from_ptrs[CHANNELS] = {source[0], source[1]};
  const float* to_ptrs[CHANNELS] = {target[0], target[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    source[0][s] = 1.0f;
    source[1][s] = 1.0f;
    target[0][s] = 0.5f;
    target[1][s] = 0.5f;
  }

  TEST_ASSERT(specbleach_transition_begin(t, latency_from, latency_to),
              "begin accepted");
  TEST_ASSERT(specbleach_transition_active(t), "active during fade");
  TEST_ASSERT(specbleach_transition_get_latency(t) == latency_to,
              "reported latency equals target immediately");

  // NULL source during fade must be rejected
  TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, NULL, to_ptrs,
                                            blended_ptrs) == false,
              "NULL from_out during fade rejected");

  /* Continuity: the very first faded sample must equal the source exactly
   * (w_from == cos(0) == 1), proving no rewind or offset is applied. */
  TEST_ASSERT(
      specbleach_transition_process(t, 1, from_ptrs, to_ptrs, blended_ptrs),
      "single sample processed");
  TEST_ASSERT(blended[0][0] == 1.0f, "fade starts at source sample exactly");
  TEST_ASSERT(blended[1][0] == 1.0f, "ch1 fade starts at source exactly");

  uint32_t blocks = 0;
  while (specbleach_transition_active(t)) {
    TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                              blended_ptrs),
                "fade block processed");
    if (++blocks > 100U) {
      TEST_ASSERT(0, "fade did not finish in reasonable time");
    }
  }

  /* After the fade the output converges to the pure target signal. */
  TEST_ASSERT(fabsf(blended[0][BLOCK_SIZE - 1] - 0.5f) < 1e-4f,
              "converged to target after fade");
  TEST_ASSERT(!specbleach_transition_active(t), "inactive after fade");

  // Reverse direction exercises the same path toward the other engine and
  // in-place blending (blended aliases target).
  TEST_ASSERT(specbleach_transition_begin(t, latency_to, latency_from),
              "reverse begin accepted");
  blocks = 0;
  while (specbleach_transition_active(t)) {
    TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, to_ptrs, to_ptrs,
                                              blended_ptrs),
                "reverse fade block processed (in-place)");
    if (++blocks > 200U) {
      TEST_ASSERT(0, "reverse transition did not finish");
    }
  }
  TEST_ASSERT(fabsf(blended[0][BLOCK_SIZE - 1] - 0.5f) < 1e-3f,
              "converged after reverse fade");

  // Post-transition: from may be NULL now, passthrough applies.
  memset(blended, 0, sizeof(blended));
  TEST_ASSERT(
      specbleach_transition_process(t, BLOCK_SIZE, NULL, to_ptrs, blended_ptrs),
      "post-transition passthrough");
  TEST_ASSERT(blended[0][0] == 0.5f, "passthrough after settle");

  specbleach_transition_free(t);
  printf("✓ Equal-power fade tests passed\n");
}

static void test_equal_latency_noop(void) {
  printf("Testing equal-latency transition...\n");

  specbleach_transition* t =
      specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, CHANNELS);

  float source[CHANNELS][BLOCK_SIZE];
  float target[CHANNELS][BLOCK_SIZE];
  float blended[CHANNELS][BLOCK_SIZE];
  const float* from_ptrs[CHANNELS] = {source[0], source[1]};
  const float* to_ptrs[CHANNELS] = {target[0], target[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    source[0][s] = 0.25f;
    source[1][s] = 0.25f;
    target[0][s] = 0.75f;
    target[1][s] = 0.75f;
  }

  TEST_ASSERT(specbleach_transition_begin(t, 1024, 1024), "begin ok");
  TEST_ASSERT(!specbleach_transition_active(t),
              "equal latencies need no transition");
  TEST_ASSERT(specbleach_transition_get_latency(t) == 1024,
              "equal-latency begin reports target");

  TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                            blended_ptrs),
              "process ok");
  TEST_ASSERT(blended[0][10] == 0.75f, "equal latency passes through target");
  TEST_ASSERT(blended[1][10] == 0.75f, "ch1 equal latency passes through");

  specbleach_transition_free(t);
  printf("✓ Equal-latency tests passed\n");
}

int main(void) {
  test_init_and_validation();
  test_idle_passthrough();
  test_equal_power_fade();
  test_equal_latency_noop();

  printf("✅ All specbleach transition tests passed!\n");
  return 0;
}
