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

#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/release_shaper.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s (%s:%d)\n", message, __FILE__,          \
              __LINE__);                                                       \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

#define FFT_SIZE 1024
#define SAMPLE_RATE 44100
#define NUM_BINS ((FFT_SIZE / 2U) + 1U)
#define MID_BIN 256U

static void fill_magnitude(float* magnitude, float value) {
  for (uint32_t k = 0; k < NUM_BINS; k++) {
    magnitude[k] = value;
  }
}

static float run_frame(ReleaseShaper* shaper, float* magnitude, float* scale) {
  TEST_ASSERT(release_shaper_compute(shaper, magnitude, scale),
              "release_shaper_compute should succeed");
  float mid_scale = scale[MID_BIN];
  for (uint32_t k = 0; k < NUM_BINS; k++) {
    TEST_ASSERT(scale[k] >= 0.0F && scale[k] <= 1.0F,
                "Release scale must be within [0,1]");
  }
  return mid_scale;
}

static void test_release_shaper(void) {
  printf("Testing Release Shaper...\n");

  // NULL / invalid safety
  TEST_ASSERT(release_shaper_initialize(0, FFT_SIZE) == NULL,
              "Initialization with sample_rate=0 should fail");
  TEST_ASSERT(release_shaper_initialize(SAMPLE_RATE, 0) == NULL,
              "Initialization with fft_size=0 should fail");
  release_shaper_free(NULL);
  release_shaper_reset(NULL);
  release_shaper_set_hop_sec(NULL, 0.01F);
  TEST_ASSERT(!release_shaper_compute(NULL, NULL, NULL),
              "Compute with NULL self should fail");

  ReleaseShaper* shaper = release_shaper_initialize(SAMPLE_RATE, FFT_SIZE);
  TEST_ASSERT(shaper != NULL, "Initialization should succeed");
  release_shaper_set_hop_sec(
      shaper, (float)(FFT_SIZE / OVERLAP_FACTOR) / (float)SAMPLE_RATE);

  float magnitude[NUM_BINS];
  float scale[NUM_BINS];

  // Frame 1 (seed): full release everywhere, no crash on empty spectrum
  fill_magnitude(magnitude, 0.0F);
  TEST_FLOAT_CLOSE(run_frame(shaper, magnitude, scale), 1.0F, 0.0001F);

  // Steady noise-only floor: envelope tracks current -> no collapse evidence
  for (int frame = 0; frame < 10; frame++) {
    TEST_FLOAT_CLOSE(run_frame(shaper, magnitude, scale), 1.0F, 0.0001F);
  }

  // Small noise dip (~1.5 dB in power terms): below the collapse onset,
  // release stays long
  fill_magnitude(magnitude, 0.7F);
  TEST_ASSERT(run_frame(shaper, magnitude, scale) > 0.8F,
              "Small band dip must keep the full release");

  // Fill back to steady state before the real signal test
  fill_magnitude(magnitude, 1.0F);
  for (int frame = 0; frame < 10; frame++) {
    run_frame(shaper, magnitude, scale);
  }

  // Signal present for a while, then ends: band collapses -> fast release
  fill_magnitude(magnitude, 100.0F); // signal 40 dB above the floor
  for (int frame = 0; frame < 5; frame++) {
    run_frame(shaper, magnitude, scale);
  }
  // Hard cut to exact zero: band_energy == 0 must still count as collapse
  // evidence (floored for the ratio), not fall back to the full release
  fill_magnitude(magnitude, 0.0F);
  TEST_ASSERT(run_frame(shaper, magnitude, scale) < 0.2F,
              "Zero-energy band after loud signal must close fast");
  fill_magnitude(magnitude, 1.0F); // signal stops (back to noise floor)
  TEST_ASSERT(run_frame(shaper, magnitude, scale) < 0.2F,
              "Band collapse after signal end must shorten the release");
  for (int frame = 0; frame < 3; frame++) {
    TEST_ASSERT(run_frame(shaper, magnitude, scale) < 0.2F,
                "Fast release persists while the collapse gap holds");
  }

  // Envelope decays back to the floor: full release is restored (anti-chirp
  // protection active again in noise-only stretches)
  bool restored = false;
  for (int frame = 0; frame < 200; frame++) {
    if (run_frame(shaper, magnitude, scale) > 0.95F) {
      restored = true;
      break;
    }
  }
  TEST_ASSERT(restored,
              "Full release must be restored after the envelope decays");

  // Signal returns: fast envelope attack keeps collapse evidence away
  fill_magnitude(magnitude, 10.0F);
  TEST_FLOAT_CLOSE(run_frame(shaper, magnitude, scale), 1.0F, 0.05F);

  release_shaper_reset(shaper);
  fill_magnitude(magnitude, 0.0F);
  TEST_FLOAT_CLOSE(run_frame(shaper, magnitude, scale), 1.0F, 0.0001F);

  release_shaper_free(shaper);

  printf("✓ Release Shaper full coverage tests passed\n");
}

int main(void) {
  test_release_shaper();
  return 0;
}
