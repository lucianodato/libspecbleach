/*
 * Unit tests for Louizou Noise Estimator
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/configurations.h"
#include "../src/shared/noise_estimation/louizou_noise_estimator.h"
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

void test_louizou_initialization(void) {
  printf("Testing Louizou initialization...\n");

  uint32_t real_size = 257;
  uint32_t sr = 44100;
  uint32_t fft_size = 512;

  LouizouNoiseEstimator* est =
      louizou_noise_estimator_initialize(real_size, sr, fft_size);
  TEST_ASSERT(est != NULL, "Initialization should succeed");

  louizou_noise_estimator_free(est);
  louizou_noise_estimator_free(NULL); // Should be safe
  printf("✓ Louizou initialization tests passed\n");
}

void test_louizou_run_logic(void) {
  printf("Testing Louizou run logic...\n");

  uint32_t real_size = 64;
  uint32_t sr = 44100;
  uint32_t fft_size = 127; // Odd size to check crossover bins

  LouizouNoiseEstimator* est =
      louizou_noise_estimator_initialize(real_size, sr, fft_size);

  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));

  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1.0f;
  }

  // First frame
  TEST_ASSERT(louizou_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed");
  for (uint32_t i = 0; i < real_size; i++) {
    TEST_FLOAT_CLOSE(noise_spectrum[i], 1.0f, 1e-6f);
  }

  // Subsequent frames (increasing spectrum to trigger speech detection logic)
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 10.0f;
  }
  TEST_ASSERT(louizou_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed");

  // NULL checks
  TEST_ASSERT(!louizou_noise_estimator_run(NULL, spectrum, noise_spectrum),
              "Should fail with NULL estimator");
  TEST_ASSERT(!louizou_noise_estimator_run(est, NULL, noise_spectrum),
              "Should fail with NULL input");
  TEST_ASSERT(!louizou_noise_estimator_run(est, spectrum, NULL),
              "Should fail with NULL output");

  louizou_noise_estimator_free(est);
  free(spectrum);
  free(noise_spectrum);
  printf("✓ Louizou run logic tests passed\n");
}

void test_louizou_state_management(void) {
  printf("Testing Louizou state management...\n");

  uint32_t real_size = 64;
  LouizouNoiseEstimator* est =
      louizou_noise_estimator_initialize(real_size, 44100, 128);

  float* profile = (float*)calloc(real_size, sizeof(float));
  for (uint32_t i = 0; i < real_size; i++) {
    profile[i] = 0.5f;
  }

  // Test set state
  louizou_noise_estimator_set_state(est, profile);
  louizou_noise_estimator_set_state(NULL, profile);
  louizou_noise_estimator_set_state(est, NULL);

  // Test update seed
  louizou_noise_estimator_update_seed(est, profile);
  louizou_noise_estimator_update_seed(NULL, profile);
  louizou_noise_estimator_update_seed(est, NULL);

  // Test apply floor
  float* floor = (float*)calloc(real_size, sizeof(float));
  for (uint32_t i = 0; i < real_size; i++) {
    floor[i] = 0.8f;
  }
  louizou_noise_estimator_apply_floor(est, floor);
  louizou_noise_estimator_apply_floor(NULL, floor);
  louizou_noise_estimator_apply_floor(est, NULL);

  louizou_noise_estimator_free(est);
  free(profile);
  free(floor);
  printf("✓ Louizou state management tests passed\n");
}

int main(void) {
  printf("Running Louizou Noise Estimator tests...\n\n");

  test_louizou_initialization();
  test_louizou_run_logic();
  test_louizou_state_management();

  printf("\n✅ All Louizou tests passed!\n");
  return 0;
}
