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

/*
 * Unit tests for Brandt Noise Estimator
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"
#include "shared/denoiser_logic/estimators/brandt_noise_estimator.h"
#include "shared/utils/spectral_utils.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_brandt_initialization(void) {
  printf("Testing Brandt initialization...\n");

  uint32_t real_size = 257;
  float duration = 1000.0f;
  uint32_t sr = 44100;
  uint32_t fft_size = 512;

  BrandtNoiseEstimator* est =
      brandt_noise_estimator_initialize(real_size, duration, sr, fft_size);
  TEST_ASSERT(est != NULL, "Initialization should succeed");

  brandt_noise_estimator_set_hop_sec(NULL, 0.0115F); // must not crash
  brandt_noise_estimator_set_hop_sec(est, 0.0115F);
  brandt_noise_estimator_set_hop_sec(est, 0.0F); // edge: ignored

  brandt_noise_estimator_free(est);
  brandt_noise_estimator_free(NULL); // Should be safe
  printf("✓ Brandt initialization tests passed\n");
}

void test_brandt_run_logic(void) {
  printf("Testing Brandt run logic...\n");

  uint32_t real_size = 64;
  uint32_t sr = 44100;
  uint32_t fft_size = 128;

  BrandtNoiseEstimator* est =
      brandt_noise_estimator_initialize(real_size, 500.0f, sr, fft_size);

  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));

  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1.0f;
  }

  // Use set_state to bypass the learning period/confidence rejected start
  brandt_noise_estimator_set_state(est, spectrum);

  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed");
  // Check output matches the set state
  for (uint32_t i = 0; i < real_size; i++) {
    TEST_FLOAT_CLOSE(noise_spectrum[i], 1.0f, 1e-4f);
  }

  // Silence check logic
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 0.0f;
  }
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed in silence");
  // Should keep output from previous state (1.0) due to silence threshold skip
  for (uint32_t i = 0; i < real_size; i++) {
    TEST_FLOAT_CLOSE(noise_spectrum[i], 1.0f, 1e-4f);
  }

  // Exercise sorted-history replacement across multiple ring wraps, including
  // state and floor updates that must keep the sorted copy in sync.
  brandt_noise_estimator_set_hop_sec(est, 0.01F);
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1.0F;
  }
  brandt_noise_estimator_set_state(est, spectrum);
  float floor_profile[64];
  float state_profile[64];
  for (uint32_t frame = 0; frame < 120U; frame++) {
    for (uint32_t i = 0; i < real_size; i++) {
      spectrum[i] = 0.05F + ((float)((frame * 7U + i * 3U) % 31U) * 0.01F);
      floor_profile[i] = 0.02F + ((float)(i % 5U) * 0.005F);
      state_profile[i] = 0.1F + ((float)(i % 9U) * 0.01F);
    }
    if (frame == 37U) {
      brandt_noise_estimator_apply_floor(est, floor_profile);
    } else if (frame == 73U) {
      brandt_noise_estimator_set_state(est, state_profile);
    }
    TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
                "Run should succeed while replacing sorted history");
    for (uint32_t i = 0; i < real_size; i++) {
      TEST_ASSERT(isfinite(noise_spectrum[i]),
                  "Sorted-history updates should remain finite");
    }
  }

  // NULL checks
  TEST_ASSERT(!brandt_noise_estimator_run(NULL, spectrum, noise_spectrum),
              "Should fail with NULL estimator");
  TEST_ASSERT(!brandt_noise_estimator_run(est, NULL, noise_spectrum),
              "Should fail with NULL input");
  TEST_ASSERT(!brandt_noise_estimator_run(est, spectrum, NULL),
              "Should fail with NULL output");

  brandt_noise_estimator_free(est);
  free(spectrum);
  free(noise_spectrum);
  printf("✓ Brandt run logic tests passed\n");
}

void test_brandt_state_management(void) {
  printf("Testing Brandt state management...\n");

  uint32_t real_size = 64;
  BrandtNoiseEstimator* est =
      brandt_noise_estimator_initialize(real_size, 500.0f, 44100, 128);

  float* profile = (float*)calloc(real_size, sizeof(float));
  for (uint32_t i = 0; i < real_size; i++) {
    profile[i] = 0.5f;
  }

  // Test set state
  brandt_noise_estimator_set_state(est, profile);
  brandt_noise_estimator_set_state(NULL, profile);
  brandt_noise_estimator_set_state(est, NULL);

  // Test update seed
  brandt_noise_estimator_update_seed(est, profile);
  brandt_noise_estimator_update_seed(NULL, profile);
  brandt_noise_estimator_update_seed(est, NULL);

  // Test apply floor
  float* floor = (float*)calloc(real_size, sizeof(float));
  for (uint32_t i = 0; i < real_size; i++) {
    floor[i] = 0.8f;
  }
  brandt_noise_estimator_apply_floor(est, floor);
  brandt_noise_estimator_apply_floor(NULL, floor);
  brandt_noise_estimator_apply_floor(est, NULL);

  brandt_noise_estimator_free(est);
  free(profile);
  free(floor);
  printf("✓ Brandt state management tests passed\n");
}

void test_brandt_sorted_history_updates(void) {
  printf("Testing Brandt sorted-history updates...\n");

  BrandtNoiseEstimator* est =
      brandt_noise_estimator_initialize(1U, 500.0F, 48000U, 4096U);
  TEST_ASSERT(est != NULL, "Sorted-history estimator should initialize");
  brandt_noise_estimator_set_hop_sec(est, 0.01F);

  float spectrum = 0.0F;
  float noise = 0.0F;
  for (uint32_t frame = 0U; frame < 60U; frame++) {
    uint32_t rank = (frame * 17U) % 50U;
    float u = ((float)rank + 0.5F) / 50.0F;
    spectrum = -logf(1.0F - u);
    TEST_ASSERT(brandt_noise_estimator_run(est, &spectrum, &noise),
                "Run should succeed for exponential-history input");

    if (frame == 45U) {
      TEST_FLOAT_CLOSE(noise, 0.8500484F, 1e-4F);
    } else if (frame == 50U) {
      TEST_FLOAT_CLOSE(noise, 0.9930851F, 1e-4F);
    }
  }

  brandt_noise_estimator_free(est);
  printf("✓ Brandt sorted-history tests passed\n");
}

int main(void) {
  printf("Running Brandt Noise Estimator tests...\n\n");

  test_brandt_initialization();
  test_brandt_run_logic();
  test_brandt_state_management();
  test_brandt_sorted_history_updates();

  printf("\n✅ All Brandt tests passed!\n");
  return 0;
}
