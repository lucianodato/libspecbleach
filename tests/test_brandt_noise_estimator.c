/*
 * Unit tests for Brandt Noise Estimator
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/configurations.h"
#include "../src/shared/noise_estimation/brandt_noise_estimator.h"
#include "../src/shared/utils/spectral_utils.h"

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

  // First frame: Should snap history to spectrum
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed");
  // Due to correction factor applied to trimmed mean of identical values,
  // it should equal exactly 1.0f.
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

int main(void) {
  printf("Running Brandt Noise Estimator tests...\n\n");

  test_brandt_initialization();
  test_brandt_run_logic();
  test_brandt_state_management();

  printf("\n✅ All Brandt tests passed!\n");
  return 0;
}
