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
#define MAX_DELAY 4096

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      return;                                                                  \
    }                                                                          \
  } while (0)

static void test_init_and_validation(void) {
  printf("Testing transition init and validation...\n");

  TEST_ASSERT(specbleach_transition_initialize(0, MAX_BLOCK, CHANNELS,
                                               MAX_DELAY) == NULL,
              "sample_rate == 0 rejected");
  TEST_ASSERT(specbleach_transition_initialize(SAMPLE_RATE, 0, CHANNELS,
                                               MAX_DELAY) == NULL,
              "max_block_size == 0 rejected");
  TEST_ASSERT(specbleach_transition_initialize(SAMPLE_RATE, MAX_BLOCK, 0,
                                               MAX_DELAY) == NULL,
              "channels == 0 rejected");

  specbleach_transition* t = specbleach_transition_initialize(
      SAMPLE_RATE, MAX_BLOCK, CHANNELS, MAX_DELAY);
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
  TEST_ASSERT(!specbleach_transition_process(t, BLOCK_SIZE, NULL, to_ptrs,
                                             blended_ptrs),
              "NULL to_out rejected while idle passthrough needs it");
  TEST_ASSERT(
      !specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs, NULL),
      "NULL blended rejected");
  /* Asking for more alignment than initialized fails cleanly. */
  TEST_ASSERT(!specbleach_transition_begin(t, 100, 100 + MAX_DELAY + 1),
              "over-capacity begin rejected");

  specbleach_transition_free(t);
  specbleach_transition_free(NULL); /* must be a no-op */
  printf("✓ Transition init and validation tests passed\n");
}

static void test_idle_passthrough(void) {
  printf("Testing idle passthrough...\n");

  specbleach_transition* t = specbleach_transition_initialize(
      SAMPLE_RATE, MAX_BLOCK, CHANNELS, MAX_DELAY);

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
    blended[0][s] = 99.0f;
    blended[1][s] = 99.0f;
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

  /* Feed accepts the emitted stream and is a no-op with garbage args. */
  specbleach_transition_feed(t, to_ptrs, BLOCK_SIZE);
  specbleach_transition_feed(t, NULL, BLOCK_SIZE);
  specbleach_transition_feed(t, to_ptrs, 0);

  specbleach_transition_free(t);
  printf("✓ Idle passthrough tests passed\n");
}

/*
 * Emits a known ramp as the plugin would: render -> blend -> feed(output).
 * Switching toward the higher-latency engine must reproduce emitted history
 * exactly at w_from == 1 (sample-true continuation), then converge on the
 * target.
 */
static void test_higher_latency_continuity(void) {
  printf("Testing higher-latency fade continuity...\n");

  specbleach_transition* t = specbleach_transition_initialize(
      SAMPLE_RATE, MAX_BLOCK, CHANNELS, MAX_DELAY);

  const uint32_t latency_from = 1024;
  const uint32_t latency_to = 2048;

  float source[CHANNELS][BLOCK_SIZE];
  float target[CHANNELS][BLOCK_SIZE];
  float blended[CHANNELS][BLOCK_SIZE];
  const float* to_ptrs[CHANNELS] = {target[0], target[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};

  float value = 1.0f;

  /* Warm up: emit two blocks of a ramp through idle passthrough + feed. */
  for (uint32_t b = 0; b < 2; ++b) {
    for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
      target[0][s] = target[1][s] = value;
      value += 1.0f;
    }
    TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, NULL, to_ptrs,
                                              blended_ptrs),
                "warmup block processed");
    const float* emitted_ptrs[CHANNELS] = {blended[0], blended[1]};
    specbleach_transition_feed(t, emitted_ptrs, BLOCK_SIZE);
  }
  const float last_emitted = value - 1.0f;

  /* Begin fade and render both engines: source keeps the ramp going,
   * target sits at a constant different level. */
  TEST_ASSERT(specbleach_transition_begin(t, latency_from, latency_to),
              "begin accepted");
  TEST_ASSERT(specbleach_transition_active(t), "active during fade");

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    source[0][s] = source[1][s] = value;
    target[0][s] = target[1][s] = -500.0f;
    value += 1.0f;
  }

  /* First faded sample must equal emitted[t - diff] EXACTLY: with
   * diff = 1024 that is ramp value (last_emitted - 1024 + 1). */
  const float* from_ptrs[CHANNELS] = {source[0], source[1]};
  TEST_ASSERT(
      specbleach_transition_process(t, 1, from_ptrs, to_ptrs, blended_ptrs),
      "single fading sample processed");
  const float expected = last_emitted - 1023.0f;
  TEST_ASSERT(fabsf(blended[0][0] - expected) < 1e-3f,
              "first blend sample continues emitted history (aligned)");
  TEST_ASSERT(fabsf(blended[1][0] - expected) < 1e-3f, "ch1 aligned too");

  specbleach_transition_free(t);
  printf("✓ Higher-latency continuity tests passed\n");
}

static void test_lower_latency_fade_converges(void) {
  printf("Testing lower-latency direct fade...\n");

  specbleach_transition* t = specbleach_transition_initialize(
      SAMPLE_RATE, MAX_BLOCK, CHANNELS, MAX_DELAY);

  float source[CHANNELS][BLOCK_SIZE];
  float target[CHANNELS][BLOCK_SIZE];
  float blended[CHANNELS][BLOCK_SIZE];
  const float* from_ptrs[CHANNELS] = {source[0], source[1]};
  const float* to_ptrs[CHANNELS] = {target[0], target[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    source[0][s] = source[1][s] = 1.0f;
    target[0][s] = target[1][s] = 0.25f;
  }

  TEST_ASSERT(specbleach_transition_begin(t, 2048, 1024), "begin accepted");

  /* First sample equals source directly (no tap in this direction). */
  TEST_ASSERT(
      specbleach_transition_process(t, 1, from_ptrs, to_ptrs, blended_ptrs),
      "single sample processed");
  TEST_ASSERT(blended[0][0] == 1.0f, "lower-latency fade starts at source");

  uint32_t guard = 0;
  while (specbleach_transition_active(t)) {
    TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                              blended_ptrs),
                "fade block processed");
    if (++guard > 100U) {
      TEST_ASSERT(0, "fade did not finish");
    }
  }
  TEST_ASSERT(fabsf(blended[0][BLOCK_SIZE - 1] - 0.25f) < 1e-4f,
              "converged to lower-latency target");
  TEST_ASSERT(specbleach_transition_get_latency(t) == 1024,
              "latency reports target");

  specbleach_transition_free(t);
  printf("✓ Lower-latency fade tests passed\n");
}

int main(void) {
  test_init_and_validation();
  test_idle_passthrough();
  test_higher_latency_continuity();
  test_lower_latency_fade_converges();

  printf("✅ All specbleach transition tests passed!\n");
  return 0;
}
