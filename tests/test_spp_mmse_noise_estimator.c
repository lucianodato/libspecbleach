/*
 * Unit tests for SPP-MMSE Noise Estimator
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/configurations.h"
#include "../src/shared/noise_estimation/spp_mmse_noise_estimator.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_spp_initialization(void) {
  printf("Testing SPP-MMSE initialization...\n");

  uint32_t real_size = 257;
  SppMmseNoiseEstimator* est =
      spp_mmse_noise_estimator_initialize(real_size, 44100, 512);
  TEST_ASSERT(est != NULL, "Initialization should succeed");

  spp_mmse_noise_estimator_free(est);
  spp_mmse_noise_estimator_free(NULL);
  printf("✓ SPP-MMSE initialization tests passed\n");
}

void test_spp_run_logic(void) {
  printf("Testing SPP-MMSE run logic...\n");

  uint32_t real_size = 64;
  SppMmseNoiseEstimator* est =
      spp_mmse_noise_estimator_initialize(real_size, 44100, 128);

  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));

  // 1. Silence threshold on first frame
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = ESTIMATOR_SILENCE_THRESHOLD * 0.5f;
  }
  TEST_ASSERT(spp_mmse_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed");
  TEST_FLOAT_CLOSE(noise_spectrum[0], 0.0f, 1e-6f);

  // 2. Normal first frame
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1.0f;
  }
  TEST_ASSERT(spp_mmse_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed");
  TEST_FLOAT_CLOSE(noise_spectrum[0], 1.0f, 1e-6f);

  // 3. Silence threshold on subsequent frame
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = ESTIMATOR_SILENCE_THRESHOLD * 0.1f;
  }
  TEST_ASSERT(spp_mmse_noise_estimator_run(est, spectrum, noise_spectrum),
              "Run should succeed");
  // Should keep previous noise
  TEST_FLOAT_CLOSE(noise_spectrum[0], 1.0f, 1e-6f);

  // 4. Stagnation cap branch
  // Process many frames with very high energy to trigger stagnation
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 100.0f;
  }
  for (int i = 0; i < 50; i++) {
    spp_mmse_noise_estimator_run(est, spectrum, noise_spectrum);
  }

  // NULL checks
  TEST_ASSERT(!spp_mmse_noise_estimator_run(NULL, spectrum, noise_spectrum),
              "NULL estimator");
  TEST_ASSERT(!spp_mmse_noise_estimator_run(est, NULL, noise_spectrum),
              "NULL input");
  TEST_ASSERT(!spp_mmse_noise_estimator_run(est, spectrum, NULL),
              "NULL output");

  spp_mmse_noise_estimator_free(est);
  free(spectrum);
  free(noise_spectrum);
  printf("✓ SPP-MMSE run logic tests passed\n");
}

void test_spp_state_management(void) {
  printf("Testing SPP-MMSE state management...\n");

  uint32_t real_size = 64;
  SppMmseNoiseEstimator* est =
      spp_mmse_noise_estimator_initialize(real_size, 44100, 128);

  float* profile = (float*)calloc(real_size, sizeof(float));
  for (uint32_t i = 0; i < real_size; i++) {
    profile[i] = 0.5f;
  }

  spp_mmse_noise_estimator_set_state(est, profile);
  spp_mmse_noise_estimator_set_state(NULL, profile);
  spp_mmse_noise_estimator_set_state(est, NULL);

  spp_mmse_noise_estimator_update_seed(est, profile);
  spp_mmse_noise_estimator_update_seed(NULL, profile);
  spp_mmse_noise_estimator_update_seed(est, NULL);

  spp_mmse_noise_estimator_apply_floor(est, profile);
  spp_mmse_noise_estimator_apply_floor(NULL, profile);
  spp_mmse_noise_estimator_apply_floor(est, NULL);

  spp_mmse_noise_estimator_free(est);
  free(profile);
  printf("✓ SPP-MMSE state management tests passed\n");
}

int main(void) {
  printf("Running SPP-MMSE Noise Estimator tests...\n\n");

  test_spp_initialization();
  test_spp_run_logic();
  test_spp_state_management();

  printf("\n✅ All SPP-MMSE tests passed!\n");
  return 0;
}
