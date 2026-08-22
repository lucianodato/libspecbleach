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

#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/hpss_filter.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

void test_hpss_filter(void) {
  printf("Testing HPSS Filter...\n");

  // 1. Invalid Initialization
  HpssConfig invalid_cfg = {
      .real_spectrum_size = 0U,
  };
  HpssFilter* invalid_f = hpss_filter_initialize(invalid_cfg);
  TEST_ASSERT(invalid_f == NULL, "Initialization with size 0 should fail");

  // 2. NULL free check
  hpss_filter_free(NULL);

  // 3. Valid Initialization
  HpssConfig config = {
      .real_spectrum_size = 257U,
  };
  HpssFilter* filter = hpss_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "HPSS filter initialization should succeed");

  // Check latency is 0 (sliding HPSS is causal)
  uint32_t latency = hpss_filter_get_latency_frames(filter);
  TEST_ASSERT(latency == 0U, "Latency should be 0 frames for sliding HPSS");

  float current_mag[257];
  float delayed_mag[257];
  float mask_h[257];
  float mask_p[257];

  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = 1.0f;
  }

  // 4. Process multiple frames with constant spectrum
  for (uint32_t frame = 0; frame < 20; ++frame) {
    bool ok = hpss_filter_process(filter, current_mag, NULL, delayed_mag,
                                  mask_h, mask_p);
    TEST_ASSERT(ok, "hpss_filter_process should succeed");
  }

  // Verify masks sum to ~1.0
  for (uint32_t k = 0; k < 257; ++k) {
    float sum_mask = mask_h[k] + mask_p[k];
    TEST_FLOAT_CLOSE(sum_mask, 1.0f, 0.001f);
    TEST_ASSERT(mask_h[k] >= 0.0f && mask_h[k] <= 1.0f, "mask_h in [0, 1]");
    TEST_ASSERT(mask_p[k] >= 0.0f && mask_p[k] <= 1.0f, "mask_p in [0, 1]");
  }

  hpss_filter_free(filter);

  // 5. Test transient detection immediately on transient frame (0 latency)
  filter = hpss_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "HPSS filter initialization should succeed");

  // Feed baseline quiet frames
  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = 0.01f;
  }
  for (uint32_t frame = 0; frame < 20; ++frame) {
    (void)hpss_filter_process(filter, current_mag, NULL, delayed_mag, mask_h,
                              mask_p);
  }
  TEST_ASSERT(!hpss_filter_is_transient_detected(filter),
              "Baseline quiet should be harmonic (no transient)");

  // Feed sharp transient frame (huge broadband jump)
  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = 5.0f;
  }
  (void)hpss_filter_process(filter, current_mag, NULL, delayed_mag, mask_h,
                            mask_p);

  // Transient detected immediately on the exact arrival frame!
  TEST_ASSERT(hpss_filter_is_transient_detected(filter),
              "Transient frame should report transient detected immediately");

  // Feed quiet frames to ensure it returns to harmonic state
  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = 0.01f;
  }
  for (uint32_t frame = 0; frame < 20; ++frame) {
    (void)hpss_filter_process(filter, current_mag, NULL, delayed_mag, mask_h,
                              mask_p);
  }
  TEST_ASSERT(!hpss_filter_is_transient_detected(filter),
              "Quiet frames after transient should return to harmonic state");

  hpss_filter_free(filter);

  // 6. Test Enable / Disable On/Off
  filter = hpss_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "HPSS filter initialization should succeed");
  TEST_ASSERT(hpss_filter_get_latency_frames(filter) == 0U,
              "Sliding HPSS latency is 0");

  // Disable HPSS
  hpss_filter_set_enabled(filter, false);
  for (uint32_t frame = 0; frame < 10; ++frame) {
    bool ok = hpss_filter_process(filter, current_mag, NULL, delayed_mag,
                                  mask_h, mask_p);
    TEST_ASSERT(ok, "Process on disabled filter should succeed");
    TEST_ASSERT(delayed_mag[10] == current_mag[10],
                "Delayed magnitude should match current magnitude when disabled");
    TEST_ASSERT(mask_h[10] == 1.0f, "Harmonic mask should be 1.0 when disabled");
    TEST_ASSERT(mask_p[10] == 0.0f, "Percussive mask should be 0.0 when disabled");
  }
  TEST_ASSERT(!hpss_filter_is_transient_detected(filter),
              "Disabled HPSS should not report transients");

  // Enable HPSS back
  hpss_filter_set_enabled(filter, true);
  for (uint32_t frame = 0; frame < 10; ++frame) {
    bool ok = hpss_filter_process(filter, current_mag, NULL, delayed_mag,
                                  mask_h, mask_p);
    TEST_ASSERT(ok, "Process on enabled filter should succeed");
  }

  hpss_filter_free(filter);

  // 7. Test NULL safety
  hpss_filter_set_enabled(NULL, false);
  TEST_ASSERT(hpss_filter_get_latency_frames(NULL) == 0U, "NULL latency is 0");
  TEST_ASSERT(hpss_filter_get_onset_ratio(NULL) == 0.0f, "NULL onset ratio 0");
  TEST_ASSERT(hpss_filter_is_transient_detected(NULL) == false,
              "NULL transient detected is false");
  TEST_ASSERT(hpss_filter_process(NULL, current_mag, NULL, delayed_mag, mask_h,
                                  mask_p) == false,
              "NULL filter process fails");
  filter = hpss_filter_initialize(config);
  TEST_ASSERT(hpss_filter_is_transient_detected(filter) == false,
              "Initial transient detected is false");
  TEST_ASSERT(hpss_filter_process(filter, NULL, NULL, delayed_mag, mask_h,
                                  mask_p) == false,
              "NULL magnitude process fails");
  hpss_filter_free(filter);

  printf("✓ HPSS Filter tests passed\n");
}

void test_hpss_noise_only(void) {
  printf("Testing HPSS noise-only classification...\n");

  HpssConfig config = {
      .real_spectrum_size = 257U,
  };
  HpssFilter* filter = hpss_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "HPSS filter initialization should succeed");

  float current_mag[257];
  float noise_floor[257];
  float delayed_mag[257];
  float mask_h[257];
  float mask_p[257];

  for (uint32_t i = 0; i < 257; ++i) {
    noise_floor[i] = 0.15f;
  }

  uint32_t rng_state = 12345U;
  uint32_t false_transient_count = 0;
  for (uint32_t frame = 0; frame < 60; ++frame) {
    for (uint32_t i = 0; i < 257; ++i) {
      rng_state = rng_state * 1103515245U + 12345U;
      float r = 0.05f + 0.20f * ((float)(rng_state >> 16) / 65535.0f);
      current_mag[i] = r;
    }
    (void)hpss_filter_process(filter, current_mag, noise_floor, delayed_mag,
                              mask_h, mask_p);

    if (frame >= 5U) {
      if (hpss_filter_is_transient_detected(filter)) {
        false_transient_count++;
      }
    }
  }

  TEST_ASSERT(false_transient_count == 0,
              "Stationary noise should never trigger transient detection");

  hpss_filter_free(filter);

  printf("✓ HPSS noise-only classification tests passed\n");
}

void test_hpss_subband_transient_in_noise(void) {
  printf("Testing HPSS subband transient detection in stationary noise...\n");

  HpssConfig config = {
      .real_spectrum_size = 257U,
  };
  HpssFilter* filter = hpss_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "HPSS filter initialization should succeed");

  float current_mag[257];
  float noise_floor[257];
  float delayed_mag[257];
  float mask_h[257];
  float mask_p[257];

  for (uint32_t i = 0; i < 257; ++i) {
    noise_floor[i] = (i < 40) ? 2.0f : 0.2f;
  }

  // 1. Feed baseline heavy stationary noise (e.g. low-frequency hum/rumble +
  // noise)
  for (uint32_t frame = 0; frame < 30; ++frame) {
    for (uint32_t i = 0; i < 257; ++i) {
      current_mag[i] = (i < 40) ? 2.0f : 0.2f; // Heavy bass noise
    }
    (void)hpss_filter_process(filter, current_mag, noise_floor, delayed_mag,
                              mask_h, mask_p);
  }
  TEST_ASSERT(!hpss_filter_is_transient_detected(filter),
              "Baseline noise with heavy rumble should be harmonic");

  // 2. Inject high-frequency localized transient (e.g. snare/hat/pluck attack
  // in bins 80-160)
  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = (i < 40) ? 2.0f : 0.2f;
    if (i >= 80 && i <= 160) {
      current_mag[i] += 1.5f; // Localized attack
    }
  }
  (void)hpss_filter_process(filter, current_mag, noise_floor, delayed_mag,
                            mask_h, mask_p);

  // Transient should be detected immediately on the attack frame!
  TEST_ASSERT(hpss_filter_is_transient_detected(filter),
              "Subband transient should be detected amidst heavy noise");

  hpss_filter_free(filter);
  printf("✓ HPSS subband transient in noise tests passed\n");
}

int main(void) {
  test_hpss_filter();
  test_hpss_noise_only();
  test_hpss_subband_transient_in_noise();
  return 0;
}
