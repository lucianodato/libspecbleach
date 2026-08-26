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

#include "specbleach_delay_line.h"
#include <stdint.h>
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

#define CHANNELS 2
#define MAX_DELAY 960

/**
 * Verifies that initialization rejects the maximum unsigned 32-bit delay.
 */
static void test_max_delay_rejected(void) {
  printf("Testing maximum delay rejected...\n");
  TEST_ASSERT(specbleach_delay_line_initialize(UINT32_MAX, 1) == NULL,
              "UINT32_MAX max_delay rejected");
  TEST_ASSERT(specbleach_delay_line_initialize(UINT32_MAX, CHANNELS) == NULL,
              "UINT32_MAX max_delay rejected (multi-channel)");
}

static void test_init_and_free(void) {
  printf("Testing delay line init and free...\n");
  TEST_ASSERT(specbleach_delay_line_initialize(MAX_DELAY, 0) == NULL,
              "channels == 0 rejected");
  specbleach_delay_line_free(NULL);

  specbleach_delay_line* dl =
      specbleach_delay_line_initialize(MAX_DELAY, CHANNELS);
  TEST_ASSERT(dl != NULL, "init");
  TEST_ASSERT(specbleach_delay_line_get_max_delay(dl) == MAX_DELAY,
              "max delay roundtrip");
  TEST_ASSERT(specbleach_delay_line_get_channels(dl) == CHANNELS,
              "channel count roundtrip");
  TEST_ASSERT(specbleach_delay_line_get_delay(dl) == 0, "starts at 0");
  specbleach_delay_line_free(dl);
}

static void test_set_delay_validation(void) {
  printf("Testing set_delay validation...\n");
  specbleach_delay_line* dl =
      specbleach_delay_line_initialize(MAX_DELAY, CHANNELS);
  TEST_ASSERT(dl != NULL, "init");
  TEST_ASSERT(!specbleach_delay_line_set_delay(NULL, 10),
              "NULL instance rejected");
  TEST_ASSERT(!specbleach_delay_line_set_delay(dl, MAX_DELAY + 1U),
              "delay above maximum rejected");
  TEST_ASSERT(specbleach_delay_line_get_delay(dl) == 0,
              "rejected set leaves state untouched");
  TEST_ASSERT(specbleach_delay_line_set_delay(dl, MAX_DELAY),
              "maximum accepted");
  TEST_ASSERT(specbleach_delay_line_get_delay(dl) == MAX_DELAY,
              "set roundtrip");
  specbleach_delay_line_free(dl);
}

static void test_passthrough_at_zero(void) {
  printf("Testing passthrough at zero delay...\n");
  specbleach_delay_line* dl = specbleach_delay_line_initialize(0, 1);
  TEST_ASSERT(dl != NULL, "zero-max init");
  float in[4] = {0.1f, -0.2f, 0.3f, -0.4f};
  float out[4] = {0};
  const float* in_p = in;
  float* out_p = out;
  TEST_ASSERT(specbleach_delay_line_process(dl, 4, &in_p, &out_p),
              "process ok");
  for (int i = 0; i < 4; ++i) {
    TEST_ASSERT(out[i] == in[i], "passthrough sample");
  }
  specbleach_delay_line_free(dl);
}

/**
 * Verifies that an impulse appears at the configured delay tap and that all other samples are zero.
 */
static void test_impulse_lands_at_tap(void) {
  printf("Testing impulse lands exactly at tap...\n");
  specbleach_delay_line* dl = specbleach_delay_line_initialize(MAX_DELAY, 1);
  TEST_ASSERT(dl != NULL, "init");
  const uint32_t tap = 7;
  TEST_ASSERT(specbleach_delay_line_set_delay(dl, tap), "set tap");

  enum { N = 32 };
  float out[N] = {0};
  float in[N] = {0};
  in[3] = 1.0f;
  const float* in_p = in;
  float* out_p = out;
  TEST_ASSERT(specbleach_delay_line_process(dl, N, &in_p, &out_p),
              "process ok");
  for (uint32_t i = 0; i < N; ++i) {
    const float expected = (i == 3 + tap) ? 1.0f : 0.0f;
    TEST_ASSERT(out[i] == expected, "impulse position/exclusivity");
  }
  specbleach_delay_line_free(dl);
}

/**
 * Verifies delay-line continuity across irregular processing block sizes.
 */
static void test_continuity_across_odd_blocks(void) {
  printf("Testing continuity across odd block sizes...\n");
  specbleach_delay_line* dl =
      specbleach_delay_line_initialize(MAX_DELAY, CHANNELS);
  TEST_ASSERT(dl != NULL, "init");
  const uint32_t tap = 13;

  enum { TOTAL = 257 }; /* prime on purpose */
  float reference[TOTAL] = {0};
  for (uint32_t i = 0; i < TOTAL; ++i) {
    reference[i] = ((i * 2654435761U) % 1000) / 1000.0f - 0.5f;
  }

  TEST_ASSERT(specbleach_delay_line_set_delay(dl, tap), "set tap");
  float streamed[TOTAL];
  memset(streamed, 0, sizeof(streamed));
  uint32_t done = 0;
  uint32_t chunk_sizes[] = {1, 5, 64, 3, 128, 56};
  uint32_t chunk = 0;
  while (done < TOTAL) {
    uint32_t n = chunk_sizes[chunk++ % 6];
    if (done + n > TOTAL) {
      n = TOTAL - done;
    }
    float out_buf[256] = {0};
    float* out_ptrs[CHANNELS] = {out_buf, out_buf};
    const float* in_ptrs[CHANNELS] = {&reference[done], &reference[done]};
    /* both channels fed the same stream; independence checked below */
    TEST_ASSERT(
        specbleach_delay_line_process(dl, n, (const float* const*)in_ptrs,
                                      (float* const*)out_ptrs),
        "process ok");
    memcpy(&streamed[done], out_buf, n * sizeof(float));
    done += n;
  }

  for (uint32_t i = 0; i < TOTAL; ++i) {
    const float expected = i >= tap ? reference[i - tap] : 0.0f;
    TEST_ASSERT(streamed[i] == expected,
                "stream continuity matches ideal delay");
  }
  specbleach_delay_line_free(dl);
}

/**
 * Verifies argument validation, channel independence, and silence before the configured delay.
 */
static void test_channel_independence_and_null_args(void) {
  printf("Testing channel independence and argument checks...\n");
  specbleach_delay_line* dl =
      specbleach_delay_line_initialize(MAX_DELAY, CHANNELS);
  TEST_ASSERT(dl != NULL, "init");
  TEST_ASSERT(specbleach_delay_line_set_delay(dl, 2), "set tap");

  float a[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  float b[4] = {2.0f, 2.0f, 2.0f, 2.0f};
  float oa[4] = {0};
  float ob[4] = {0};
  const float* in_p[CHANNELS] = {a, b};
  float* out_p[CHANNELS] = {oa, ob};
  TEST_ASSERT(!specbleach_delay_line_process(NULL, 4, in_p, out_p),
              "NULL instance rejected");
  TEST_ASSERT(!specbleach_delay_line_process(dl, 4, NULL, out_p),
              "NULL input array rejected");
  TEST_ASSERT(!specbleach_delay_line_process(dl, 4, in_p, NULL),
              "NULL output array rejected");
  TEST_ASSERT(specbleach_delay_line_process(dl, 4, in_p, out_p), "process ok");
  TEST_ASSERT(oa[3] == 1.0f && ob[3] == 2.0f, "channels stay independent");
  TEST_ASSERT(oa[0] == 0.0f && ob[0] == 0.0f, "pre-tap silence");

  specbleach_delay_line_free(dl);
}

int main(void) {
  test_max_delay_rejected();
  test_init_and_free();
  test_set_delay_validation();
  test_passthrough_at_zero();
  test_impulse_lands_at_tap();
  test_continuity_across_odd_blocks();
  test_channel_independence_and_null_args();

  printf("✅ All specbleach delay line tests passed!\n");
  return 0;
}
