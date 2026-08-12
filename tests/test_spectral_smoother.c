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
#include "shared/utils/spectral_smoother.h"
#include "shared/utils/spectral_utils.h"
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

void test_spectral_smoother(void) {
  printf("Testing Spectral Smoother...\n");

  // 1. Invalid Initialization (fft_size == 0)
  SpectralSmoother* invalid_ss =
      spectral_smoothing_initialize(0, 44100, OVERLAP_FACTOR_1D, FIXED);
  TEST_ASSERT(invalid_ss == NULL, "Initialization with fft_size=0 should fail");

  // 2. Free NULL instance (safety check)
  spectral_smoothing_free(NULL);

  uint32_t fft_size = 1024;
  uint32_t num_bins = (fft_size / 2U) + 1U;
  uint32_t custom_overlap = 8U;

  // 3. Test initialization and run across smoothing types
  for (int type = FIXED; type <= TRANSIENT_AWARE; type++) {
    SpectralSmoother* ss = spectral_smoothing_initialize(
        fft_size, 44100, custom_overlap, (TimeSmoothingType)type);
    TEST_ASSERT(ss != NULL, "Spectral smoother initialization should succeed");

    float gains[1024] = {0.0f};
    for (uint32_t i = 0; i < num_bins; i++) {
      gains[i] = 1.0f + (0.5f * sinf((float)i * 0.1f));
    }
    sb_apply_spectral_symmetry(gains, num_bins, fft_size);

    // Test NULL safety in spectral_smoothing_run
    TimeSmoothingParameters params = {.smoothing = 0.8f};
    TEST_ASSERT(!spectral_smoothing_run(NULL, params, gains),
                "Run with NULL self should return false");
    TEST_ASSERT(!spectral_smoothing_run(ss, params, NULL),
                "Run with NULL gains should return false");

    // First run (uninitialized -> initializes history with gains)
    TEST_ASSERT(spectral_smoothing_run(ss, params, gains),
                "First run should succeed and initialize state");

    // Second run with active smoothing
    float new_gains[1024] = {0.0f};
    TEST_ASSERT(spectral_smoothing_run(ss, params, new_gains),
                "Second run should succeed");

    // Verify release time smoothing calculation (0.8 factor with custom
    // overlap)
    float test_tau = GAIN_SMOOTHING_MIN_RELEASE_SEC +
                     (0.8f * (GAIN_SMOOTHING_MAX_RELEASE_SEC -
                              GAIN_SMOOTHING_MIN_RELEASE_SEC));
    float test_dt = ((float)fft_size / (float)custom_overlap) / 44100.0f;
    float test_alpha = expf(-test_dt / test_tau);
    for (uint32_t i = 0; i < num_bins; i++) {
      float expected = test_alpha * (1.0f + (0.5f * sinf((float)i * 0.1f)));
      TEST_FLOAT_CLOSE(new_gains[i], expected, 0.001f);
    }

    // Verify mirrored bins produced by sb_apply_spectral_symmetry
    sb_apply_spectral_symmetry(new_gains, num_bins, fft_size);
    for (uint32_t i = 1; i < num_bins - 1; i++) {
      TEST_FLOAT_CLOSE(new_gains[fft_size - i], new_gains[i], 0.001f);
    }

    // Test 100% smoothing factor does not lock memory (alpha < 1.0)
    TimeSmoothingParameters max_params = {.smoothing = 1.0f};
    float max_gains[1024] = {0.0f};
    TEST_ASSERT(spectral_smoothing_run(ss, max_params, max_gains),
                "Run with smoothing=1.0 should succeed");
    // Ensure gains updated (decayed towards 0.0, not frozen)
    TEST_ASSERT(max_gains[0] < new_gains[0],
                "100% smoothing must update gains rather than freeze");

    // Third run with bypass smoothing <= 0.0f
    TimeSmoothingParameters zero_params = {.smoothing = 0.0f};
    float bypass_gains[1024] = {0.0f};
    for (uint32_t i = 0; i < num_bins; i++) {
      bypass_gains[i] = 0.5f;
    }
    TEST_ASSERT(spectral_smoothing_run(ss, zero_params, bypass_gains),
                "Run with smoothing=0.0 should succeed");
    for (uint32_t i = 0; i < num_bins; i++) {
      TEST_FLOAT_CLOSE(bypass_gains[i], 0.5f, 0.001f);
    }

    spectral_smoothing_free(ss);
  }

  // 4. Test spectral_smoothing_apply_spatial
  // Edge cases: NULL data, size < 2
  spectral_smoothing_apply_spatial(NULL, 10);
  float single_data[1] = {1.0f};
  spectral_smoothing_apply_spatial(single_data, 1);
  TEST_FLOAT_CLOSE(single_data[0], 1.0f, 0.001f);

  // Normal case: size = 5
  float spatial_data[5] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
  spectral_smoothing_apply_spatial(spatial_data, 5);
  TEST_FLOAT_CLOSE(spatial_data[1], 2.0f, 0.001f);
  TEST_FLOAT_CLOSE(spatial_data[4], 4.75f, 0.001f);

  // 5. Test spectral_smoothing_apply_simple_temporal
  // Edge cases: NULL parameters, size 0
  float current_buf[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  float memory_buf[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  spectral_smoothing_apply_simple_temporal(NULL, memory_buf, 4, 0.5f);
  spectral_smoothing_apply_simple_temporal(current_buf, NULL, 4, 0.5f);
  spectral_smoothing_apply_simple_temporal(current_buf, memory_buf, 0, 0.5f);

  // Normal case
  spectral_smoothing_apply_simple_temporal(current_buf, memory_buf, 4, 0.5f);
  for (int i = 0; i < 4; i++) {
    TEST_FLOAT_CLOSE(current_buf[i], 0.5f, 0.001f);
    TEST_FLOAT_CLOSE(memory_buf[i], 0.5f, 0.001f);
  }

  printf("✓ Spectral Smoother full coverage tests passed\n");
}

int main(void) {
  test_spectral_smoother();
  return 0;
}
