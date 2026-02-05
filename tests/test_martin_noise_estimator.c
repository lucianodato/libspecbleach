/*
 * Unit tests for Martin (2001) Minimum Statistics Noise Estimator
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/configurations.h"
#include "../src/shared/noise_estimation/martin_noise_estimator.h"
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

void test_martin_initialization(void) {
  printf("Testing Martin initialization...\n");

  uint32_t real_size = 257;
  uint32_t sr = 44100;
  uint32_t fft_size = 512;

  MartinNoiseEstimator* est =
      martin_noise_estimator_initialize(real_size, sr, fft_size);
  TEST_ASSERT(est != NULL, "Initialization should succeed");

  martin_noise_estimator_free(est);
  martin_noise_estimator_free(NULL); // Should be safe
  printf("✓ Martin initialization tests passed\n");
}

void test_martin_run_logic(void) {
  printf("Testing Martin run logic...\n");

  uint32_t real_size = 64;
  uint32_t sr = 44100;
  uint32_t fft_size = 128;

  MartinNoiseEstimator* est =
      martin_noise_estimator_initialize(real_size, sr, fft_size);

  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));

  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1.0f;
  }

  // First frame: Should snap to spectrum
  TEST_ASSERT(martin_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed");
  for (uint32_t i = 0; i < real_size; i++) {
    TEST_FLOAT_CLOSE(noise_spectrum[i], 1.0f, 1e-6f);
  }

  // Silence check logic
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 0.0f;
  }
  TEST_ASSERT(martin_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed in silence");
  // Should keep previous output (snapped 1.0) because of silence threshold
  for (uint32_t i = 0; i < real_size; i++) {
    TEST_FLOAT_CLOSE(noise_spectrum[i], 1.0f, 1e-6f);
  }

  // NULL checks
  TEST_ASSERT(!martin_noise_estimator_run(NULL, spectrum, noise_spectrum),
              "Should fail with NULL estimator");
  TEST_ASSERT(!martin_noise_estimator_run(est, NULL, noise_spectrum),
              "Should fail with NULL input");
  TEST_ASSERT(!martin_noise_estimator_run(est, spectrum, NULL),
              "Should fail with NULL output");

  martin_noise_estimator_free(est);
  free(spectrum);
  free(noise_spectrum);
  printf("✓ Martin run logic tests passed\n");
}

void test_martin_state_management(void) {
  printf("Testing Martin state management...\n");

  uint32_t real_size = 64;
  MartinNoiseEstimator* est =
      martin_noise_estimator_initialize(real_size, 44100, 128);

  float* profile = (float*)calloc(real_size, sizeof(float));
  for (uint32_t i = 0; i < real_size; i++) {
    profile[i] = 0.5f;
  }

  // Test set state
  martin_noise_estimator_set_state(est, profile);
  martin_noise_estimator_set_state(NULL, profile);
  martin_noise_estimator_set_state(est, NULL);

  // Test update seed
  martin_noise_estimator_update_seed(est, profile);
  martin_noise_estimator_update_seed(NULL, profile);
  martin_noise_estimator_update_seed(est, NULL);

  // Test apply floor
  float* floor = (float*)calloc(real_size, sizeof(float));
  for (uint32_t i = 0; i < real_size; i++) {
    floor[i] = 0.8f;
  }
  martin_noise_estimator_apply_floor(est, floor);
  martin_noise_estimator_apply_floor(NULL, floor);
  martin_noise_estimator_apply_floor(est, NULL);

  martin_noise_estimator_free(est);
  free(profile);
  free(floor);
  printf("✓ Martin state management tests passed\n");
}

int main(void) {
  printf("Running Martin Noise Estimator tests...\n\n");

  test_martin_initialization();
  test_martin_run_logic();
  test_martin_state_management();

  printf("\n✅ All Martin tests passed!\n");
  return 0;
}
