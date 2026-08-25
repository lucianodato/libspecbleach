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
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define SAMPLE_RATE 48000
#define MAX_BLOCK 512
#define BLOCK_SIZE 256
#define CHANNELS 2

static void test_init_and_validation(void) {
  printf("Testing transition init and validation...\n");

  TEST_ASSERT(specbleach_transition_initialize(0, MAX_BLOCK, CHANNELS) == NULL,
              "sample_rate == 0 rejected");
  TEST_ASSERT(
      specbleach_transition_initialize(SAMPLE_RATE, 0, CHANNELS) == NULL,
      "max_block == 0 rejected");
  TEST_ASSERT(
      specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, 0) == NULL,
      "channels == 0 rejected");
  specbleach_transition_free(NULL);

  specbleach_transition* t =
      specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, CHANNELS);
  TEST_ASSERT(t != NULL, "valid init");
  TEST_ASSERT(!specbleach_transition_active(t), "inactive after init");
  TEST_ASSERT(specbleach_transition_get_latency(t) == 0,
              "zero latency before begin");

  TEST_ASSERT(specbleach_transition_process(NULL, BLOCK_SIZE, NULL, NULL,
                                            NULL) == false,
              "NULL instance process rejected");
  TEST_ASSERT(specbleach_transition_begin(NULL, 0, 0) == false,
              "NULL instance begin rejected");

  specbleach_transition_free(t);
  printf("✓ Transition init and validation tests passed\n");
}

static void test_idle_passthrough(void) {
  printf("Testing idle passthrough...\n");

  specbleach_transition* t =
      specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, CHANNELS);

  float to[CHANNELS][BLOCK_SIZE];
  float blended[CHANNELS][BLOCK_SIZE];
  const float* from_ptrs[CHANNELS] = {NULL};
  const float* to_ptrs[CHANNELS] = {to[0], to[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    to[0][s] = (float)s;
    to[1][s] = -(float)s;
  }

  // Before begin(): to_out is passed through untouched (from_out may be NULL)
  TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                            blended_ptrs),
              "idle passthrough succeeds");
  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    TEST_ASSERT(blended[0][s] == to[0][s], "ch0 passthrough identical");
    TEST_ASSERT(blended[1][s] == to[1][s], "ch1 passthrough identical");
  }

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
  float blended_alias[CHANNELS][BLOCK_SIZE];
  const float* from_ptrs[CHANNELS] = {source[0], source[1]};
  const float* alias_from_ptrs[CHANNELS] = {blended_alias[0], blended_alias[1]};
  const float* to_ptrs[CHANNELS] = {target[0], target[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};
  float* alias_blended_ptrs[CHANNELS] = {blended_alias[0], blended_alias[1]};

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

  uint32_t blocks = 0;
  while (specbleach_transition_active(t)) {
    TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                              blended_ptrs),
                "fade block processed");
    if (++blocks > 100U) {
      TEST_ASSERT(0, "fade did not finish in reasonable time");
    }
  }

  /* After a fade that lands on the HIGHER latency engine there is no slew
   * phase; output must converge to the pure target signal. */
  TEST_ASSERT(fabsf(blended[0][BLOCK_SIZE - 1] - 0.5f) < 1e-4f,
              "converged to target after fade");

  // Restart toward the LOWER latency engine: exercises the slew phase and
  // in-place blending (blended aliases target).
  TEST_ASSERT(specbleach_transition_begin(t, latency_to, latency_from),
              "reverse begin accepted");
  blocks = 0;
  while (specbleach_transition_active(t)) {
    TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, alias_from_ptrs,
                                              to_ptrs, alias_blended_ptrs),
                "reverse fade block processed (in-place)");
    if (++blocks > 200U) {
      TEST_ASSERT(0, "reverse transition did not finish");
    }
  }
  TEST_ASSERT(fabsf(blended_alias[0][BLOCK_SIZE - 1] - 0.5f) < 1e-3f,
              "converged after slew phase");

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

  float source[BLOCK_SIZE];
  float target[BLOCK_SIZE];
  float blended[BLOCK_SIZE];
  const float* from_ptrs[1] = {source};
  const float* to_ptrs[1] = {target};
  float* blended_ptrs[1] = {blended};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    source[s] = 0.25f;
    target[s] = 0.75f;
  }

  TEST_ASSERT(specbleach_transition_begin(t, 1024, 1024), "begin ok");
  TEST_ASSERT(!specbleach_transition_active(t),
              "equal latencies need no transition");

  TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                            blended_ptrs),
              "process ok");
  TEST_ASSERT(blended[10] == 0.75f, "equal latency passes through target");

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
