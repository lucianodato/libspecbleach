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
  HpssConfig invalid_cfg = {0, 17, 17, 1.1f};
  HpssFilter* invalid_f = hpss_filter_initialize(invalid_cfg);
  TEST_ASSERT(invalid_f == NULL, "Initialization with size 0 should fail");

  // 2. NULL free check
  hpss_filter_free(NULL);

  // 3. Valid Initialization
  HpssConfig config = {
      .real_spectrum_size = 257U,
      .time_window_size = HPSS_TIME_WINDOW_MEDIUM,
      .freq_window_size = HPSS_FREQ_WINDOW_MEDIUM,
      .noise_oversubtraction = HPSS_NOISE_OVERSUBTRACTION_DEFAULT,
  };
  HpssFilter* filter = hpss_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "HPSS filter initialization should succeed");

  // Check latency
  uint32_t latency = hpss_filter_get_latency_frames(filter);
  TEST_ASSERT(latency == 8U, "Latency should be 8 frames for 17-frame window");

  float current_mag[257];
  float noise_profile[257];
  float delayed_mag[257];
  float mask_h[257];
  float mask_p[257];

  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = 1.0f;
    noise_profile[i] = 0.2f;
  }

  // 4. Process multiple frames with constant spectrum
  for (uint32_t frame = 0; frame < 20; ++frame) {
    bool ok = hpss_filter_process(filter, current_mag, noise_profile,
                                  delayed_mag, mask_h, mask_p);
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

  // 5. Test transient onset detection
  filter = hpss_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "HPSS filter initialization should succeed");

  // Feed baseline quiet frames
  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = 0.01f;
    noise_profile[i] = 0.005f;
  }
  for (uint32_t frame = 0; frame < 20; ++frame) {
    (void)hpss_filter_process(filter, current_mag, noise_profile, delayed_mag,
                              mask_h, mask_p);
  }

  // Feed sharp transient frame (huge broadband jump)
  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = 5.0f;
  }
  (void)hpss_filter_process(filter, current_mag, noise_profile, delayed_mag,
                            mask_h, mask_p);

  // Feed quiet frames until the transient reaches center frame (latency frames
  // later)
  for (uint32_t i = 0; i < 257; ++i) {
    current_mag[i] = 0.01f;
  }
  for (uint32_t frame = 0; frame < latency; ++frame) {
    (void)hpss_filter_process(filter, current_mag, noise_profile, delayed_mag,
                              mask_h, mask_p);
  }

  // When the transient is in the center frame, mask_p should dominate (high
  // percussive weight)
  float avg_p = 0.0f;
  for (uint32_t k = 0; k < 257; ++k) {
    avg_p += mask_p[k];
  }
  avg_p /= 257.0f;
  TEST_ASSERT(avg_p > 0.5f,
              "Percussive mask should dominate on transient frame");

  hpss_filter_free(filter);

  // 6. Test Dynamic Quality Mode Reconfiguration
  filter = hpss_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "HPSS filter initialization should succeed");
  TEST_ASSERT(hpss_filter_get_latency_frames(filter) == 8U,
              "Default latency is 8");

  // Switch to Low
  hpss_filter_set_quality_mode(filter, HPSS_QUALITY_LOW);
  TEST_ASSERT(hpss_filter_get_latency_frames(filter) == 4U,
              "Low quality latency should be 4 frames");
  for (uint32_t frame = 0; frame < 10; ++frame) {
    bool ok = hpss_filter_process(filter, current_mag, noise_profile,
                                  delayed_mag, mask_h, mask_p);
    TEST_ASSERT(ok, "Process on Low quality mode should succeed");
  }

  // Switch to High
  hpss_filter_set_quality_mode(filter, HPSS_QUALITY_HIGH);
  TEST_ASSERT(hpss_filter_get_latency_frames(filter) == 16U,
              "High quality latency should be 16 frames");
  for (uint32_t frame = 0; frame < 10; ++frame) {
    bool ok = hpss_filter_process(filter, current_mag, noise_profile,
                                  delayed_mag, mask_h, mask_p);
    TEST_ASSERT(ok, "Process on High quality mode should succeed");
  }

  // Switch to Off
  hpss_filter_set_quality_mode(filter, HPSS_QUALITY_OFF);
  TEST_ASSERT(hpss_filter_get_latency_frames(filter) == 0U,
              "Off quality latency should be 0 frames");
  for (uint32_t frame = 0; frame < 10; ++frame) {
    bool ok = hpss_filter_process(filter, current_mag, noise_profile,
                                  delayed_mag, mask_h, mask_p);
    TEST_ASSERT(ok, "Process on Off quality mode should succeed");
    TEST_ASSERT(delayed_mag[10] == current_mag[10],
                "Delayed magnitude should match current magnitude in Off mode");
    TEST_ASSERT(mask_h[10] == 1.0f, "Harmonic mask should be 1.0 in Off mode");
  }

  // Switch back to Medium
  hpss_filter_set_quality_mode(filter, HPSS_QUALITY_MEDIUM);
  TEST_ASSERT(hpss_filter_get_latency_frames(filter) == 8U,
              "Medium quality latency should be 8 frames");
  // Test redundant mode set
  hpss_filter_set_quality_mode(filter, HPSS_QUALITY_MEDIUM);
  TEST_ASSERT(hpss_filter_get_latency_frames(filter) == 8U,
              "Medium quality latency should remain 8 frames");

  for (uint32_t frame = 0; frame < 10; ++frame) {
    bool ok = hpss_filter_process(filter, current_mag, noise_profile,
                                  delayed_mag, mask_h, mask_p);
    TEST_ASSERT(ok, "Process on Medium quality mode should succeed");
  }

  // Switch to invalid/default mode
  hpss_filter_set_quality_mode(filter, (HpssQualityMode)999);
  TEST_ASSERT(hpss_filter_get_latency_frames(filter) == 8U,
              "Default quality latency should be 8 frames");

  hpss_filter_free(filter);

  // 7. Test even window configuration and default noise oversubtraction
  HpssConfig even_cfg = {
      .real_spectrum_size = 257U,
      .time_window_size = 16U,
      .freq_window_size = 16U,
      .noise_oversubtraction = -1.0f,
  };
  HpssFilter* even_filter = hpss_filter_initialize(even_cfg);
  TEST_ASSERT(even_filter != NULL, "Even config should be normalized to odd");
  TEST_ASSERT(hpss_filter_get_latency_frames(even_filter) == 8U,
              "16 window size rounds to 17, latency 8");
  hpss_filter_free(even_filter);

  // 8. Test NULL safety
  hpss_filter_set_quality_mode(NULL, HPSS_QUALITY_LOW);
  TEST_ASSERT(hpss_filter_get_latency_frames(NULL) == 0U, "NULL latency is 0");
  TEST_ASSERT(hpss_filter_get_onset_ratio(NULL) == 0.0f, "NULL onset ratio 0");
  TEST_ASSERT(hpss_filter_process(NULL, current_mag, noise_profile, delayed_mag,
                                  mask_h, mask_p) == false,
              "NULL filter process fails");
  filter = hpss_filter_initialize(config);
  TEST_ASSERT(hpss_filter_process(filter, NULL, noise_profile, delayed_mag,
                                  mask_h, mask_p) == false,
              "NULL magnitude process fails");
  hpss_filter_free(filter);

  printf("✓ HPSS Filter tests passed\n");
}

int main(void) {
  test_hpss_filter();
  return 0;
}
