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
  /* While idle, from_out may be NULL. */
  TEST_ASSERT(
      specbleach_transition_process(t, BLOCK_SIZE, NULL, to_ptrs, blended_ptrs),
      "NULL from_out allowed while idle");
  TEST_ASSERT(!specbleach_transition_begin(t, 100, 100 + MAX_DELAY + 1),
              "over-capacity begin rejected");

  specbleach_transition_free(t);
  specbleach_transition_free(NULL);
  printf("✓ Transition init and validation tests passed\n");
}

/*
 * Fading toward the HIGHER-latency engine aligns the SOURCE side through
 * the internal history: after `diff` samples, every output must equal
 * w_from * from[s - diff] + w_to * to[s], computed independently here.
 */
static void test_higher_latency_alignment(void) {
  printf("Testing higher-latency aligned fade...\n");

  specbleach_transition* t = specbleach_transition_initialize(
      SAMPLE_RATE, MAX_BLOCK, CHANNELS, MAX_DELAY);

  const uint32_t latency_from = 1024;
  const uint32_t latency_to = 1024 + 8; /* small tap for exact math */
  const uint32_t diff = latency_to - latency_from;

  float source[CHANNELS][BLOCK_SIZE];
  float target[CHANNELS][BLOCK_SIZE];
  float blended[CHANNELS][BLOCK_SIZE];
  const float* from_ptrs[CHANNELS] = {source[0], source[1]};
  const float* to_ptrs[CHANNELS] = {target[0], target[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    source[0][s] = source[1][s] = (float)s * 0.001f;
    target[0][s] = target[1][s] = -0.5f;
  }

  TEST_ASSERT(specbleach_transition_begin(t, latency_from, latency_to),
              "begin accepted");
  TEST_ASSERT(specbleach_transition_active(t), "active during fade");

  TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                            blended_ptrs),
              "fade block processed");

  /* Independently recompute the expected blend for s >= diff. */
  const float fade_samples =
      (float)SAMPLE_RATE * (SPECBLEACH_TRANSITION_FADE_TIME_MS / 1000.0f);
/* Exact check: track the history model locally. */
#define TAP_RING (8 + 1)
  float model[TAP_RING];
  memset(model, 0, sizeof(model));
  uint32_t model_pos = 0;
  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    model[model_pos] = source[0][s];
    const size_t read = (model_pos + TAP_RING - diff) % (TAP_RING);
    if (s >= diff) {
      const float p = (float)s / fade_samples;
      const float rad = (p > 1.0f ? 1.0f : p) * 1.5707963267948966f;
      const float expected = cosf(rad) * model[read] + sinf(rad) * target[0][s];
      TEST_ASSERT(fabsf(blended[0][s] - expected) < 1e-3f,
                  "aligned blend matches independent model (ch0)");
    }
    model_pos = (model_pos + 1U) % (TAP_RING);
  }

  specbleach_transition_free(t);
  printf("✓ Higher-latency alignment tests passed\n");
}

/*
 * Fading toward the LOWER-latency engine blends onto an aligned copy of it,
 * then the slew phase slides the delay back out until the direct target is
 * emitted verbatim.
 */
static void test_lower_latency_slew(void) {
  printf("Testing lower-latency fade and slew...\n");

  specbleach_transition* t = specbleach_transition_initialize(
      SAMPLE_RATE, MAX_BLOCK, CHANNELS, MAX_DELAY);

  const uint32_t latency_from = 1024 + 8;
  const uint32_t latency_to = 1024;

  float source[CHANNELS][BLOCK_SIZE];
  float target[CHANNELS][BLOCK_SIZE];
  float blended[CHANNELS][BLOCK_SIZE];
  const float* from_ptrs[CHANNELS] = {source[0], source[1]};
  const float* to_ptrs[CHANNELS] = {target[0], target[1]};
  float* blended_ptrs[CHANNELS] = {blended[0], blended[1]};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    source[0][s] = source[1][s] = 1.0f;
    target[0][s] = target[1][s] = 0.25f + 0.001f * (float)s;
  }

  TEST_ASSERT(specbleach_transition_begin(t, latency_from, latency_to),
              "begin accepted");
  TEST_ASSERT(specbleach_transition_get_latency(t) == latency_to,
              "target latency reported immediately");

  uint32_t guard = 0;
  while (specbleach_transition_active(t)) {
    TEST_ASSERT(specbleach_transition_process(t, BLOCK_SIZE, from_ptrs, to_ptrs,
                                              blended_ptrs),
                "transition block processed");
    if (++guard > 400U) {
      TEST_ASSERT(0, "fade+slew did not finish");
    }
  }

  /* After everything settles, output must be the pure direct target. */
  TEST_ASSERT(
      fabsf(blended[0][BLOCK_SIZE - 1] - target[0][BLOCK_SIZE - 1]) < 1e-4f,
      "slew converges to direct lower-latency target");
  TEST_ASSERT(
      fabsf(blended[1][BLOCK_SIZE - 1] - target[1][BLOCK_SIZE - 1]) < 1e-4f,
      "ch1 slew converges too");

  /* Post-transition passthrough. */
  memset(blended, 0, sizeof(blended));
  TEST_ASSERT(
      specbleach_transition_process(t, BLOCK_SIZE, NULL, to_ptrs, blended_ptrs),
      "idle passthrough after slew");
  TEST_ASSERT(blended[0][5] == target[0][5], "passthrough intact");

  specbleach_transition_free(t);
  printf("✓ Lower-latency fade/slew tests passed\n");
}

int main(void) {
  test_init_and_validation();
  test_higher_latency_alignment();
  test_lower_latency_slew();

  printf("✅ All specbleach transition tests passed!\n");
  return 0;
}
