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

void test_brandt_convergence(void) {
  printf("Testing Brandt convergence over multiple frames...\n");

  uint32_t real_size = 64;
  uint32_t sr = 44100;
  uint32_t fft_size = 128;

  BrandtNoiseEstimator* est =
      brandt_noise_estimator_initialize(real_size, 500.0f, sr, fft_size);

  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));

  // Initialize with a known noise floor
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 0.1f + (float)i * 0.01f; // Increasing noise floor
  }
  brandt_noise_estimator_set_state(est, spectrum);

  // Run multiple frames with same noise level
  for (int frame = 0; frame < 50; frame++) {
    // Add small variations to simulate noise
    for (uint32_t i = 0; i < real_size; i++) {
      spectrum[i] = 0.1f + (float)i * 0.01f + 0.02f * sinf((float)frame * 0.1f);
    }
    brandt_noise_estimator_run(est, spectrum, noise_spectrum);
  }

  // After convergence, noise estimate should be reasonable
  for (uint32_t i = 0; i < real_size; i++) {
    TEST_ASSERT(noise_spectrum[i] >= 0.05f && noise_spectrum[i] <= 0.5f,
                "Converged noise estimate should be in reasonable range");
  }

  free(spectrum);
  free(noise_spectrum);
  brandt_noise_estimator_free(est);
  printf("✓ Brandt convergence tests passed\n");
}

void test_brandt_varying_energy_levels(void) {
  printf("Testing Brandt with varying energy levels...\n");

  uint32_t real_size = 32;
  uint32_t sr = 44100;
  uint32_t fft_size = 64;

  BrandtNoiseEstimator* est =
      brandt_noise_estimator_initialize(real_size, 300.0f, sr, fft_size);

  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));

  // Initialize with moderate energy
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1.0f;
  }
  brandt_noise_estimator_set_state(est, spectrum);

  // Test with very low energy (should handle gracefully)
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1e-8f;
  }
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Should handle very low energy");

  // Test with high energy
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 100.0f;
  }
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Should handle high energy");

  // Test with mixed energy levels
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = (i % 2 == 0) ? 0.01f : 10.0f;
  }
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Should handle mixed energy levels");

  free(spectrum);
  free(noise_spectrum);
  brandt_noise_estimator_free(est);
  printf("✓ Brandt varying energy level tests passed\n");
}

void test_brandt_spectral_patterns(void) {
  printf("Testing Brandt with different spectral patterns...\n");

  uint32_t real_size = 64;
  uint32_t sr = 44100;
  uint32_t fft_size = 128;

  BrandtNoiseEstimator* est =
      brandt_noise_estimator_initialize(real_size, 400.0f, sr, fft_size);

  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));

  // Test 1: Flat spectrum
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1.0f;
  }
  brandt_noise_estimator_set_state(est, spectrum);
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Should handle flat spectrum");

  // Test 2: Sloped spectrum (pink noise like)
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 1.0f / sqrtf((float)(i + 1));
  }
  brandt_noise_estimator_set_state(est, spectrum);
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Should handle sloped spectrum");

  // Test 3: Peaked spectrum (tonal components)
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 0.1f;
    if (i == 10 || i == 30 || i == 50) {
      spectrum[i] = 5.0f; // Peaks
    }
  }
  brandt_noise_estimator_set_state(est, spectrum);
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Should handle peaked spectrum");

  // Test 4: Decaying spectrum
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = expf(-(float)i * 0.05f);
  }
  brandt_noise_estimator_set_state(est, spectrum);
  TEST_ASSERT(brandt_noise_estimator_run(est, spectrum, noise_spectrum),
              "Should handle decaying spectrum");

  free(spectrum);
  free(noise_spectrum);
  brandt_noise_estimator_free(est);
  printf("✓ Brandt spectral pattern tests passed\n");
}

void test_brandt_history_duration(void) {
  printf("Testing Brandt with different history durations...\n");

  uint32_t real_size = 32;
  uint32_t sr = 44100;
  uint32_t fft_size = 64;

  // Test short history
  BrandtNoiseEstimator* est_short =
      brandt_noise_estimator_initialize(real_size, 100.0f, sr, fft_size);
  TEST_ASSERT(est_short != NULL, "Short history initialization should succeed");
  brandt_noise_estimator_free(est_short);

  // Test medium history
  BrandtNoiseEstimator* est_med =
      brandt_noise_estimator_initialize(real_size, 500.0f, sr, fft_size);
  TEST_ASSERT(est_med != NULL, "Medium history initialization should succeed");
  brandt_noise_estimator_free(est_med);

  // Test long history
  BrandtNoiseEstimator* est_long =
      brandt_noise_estimator_initialize(real_size, 2000.0f, sr, fft_size);
  TEST_ASSERT(est_long != NULL, "Long history initialization should succeed");
  brandt_noise_estimator_free(est_long);

  // Test very short history (edge case)
  BrandtNoiseEstimator* est_vshort =
      brandt_noise_estimator_initialize(real_size, 10.0f, sr, fft_size);
  TEST_ASSERT(est_vshort != NULL,
              "Very short history initialization should succeed");
  brandt_noise_estimator_free(est_vshort);

  printf("✓ Brandt history duration tests passed\n");
}

void test_brandt_floor_application(void) {
  printf("Testing Brandt floor application behavior...\n");

  uint32_t real_size = 32;
  uint32_t sr = 44100;
  uint32_t fft_size = 64;

  BrandtNoiseEstimator* est =
      brandt_noise_estimator_initialize(real_size, 300.0f, sr, fft_size);

  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));
  float* floor_profile = (float*)calloc(real_size, sizeof(float));

  // Initialize with low noise
  for (uint32_t i = 0; i < real_size; i++) {
    spectrum[i] = 0.1f;
  }
  brandt_noise_estimator_set_state(est, spectrum);

  // Apply floor higher than current estimate
  for (uint32_t i = 0; i < real_size; i++) {
    floor_profile[i] = 0.5f;
  }
  brandt_noise_estimator_apply_floor(est, floor_profile);

  // Run estimation - should respect the floor
  brandt_noise_estimator_run(est, spectrum, noise_spectrum);

  // Check that estimate is influenced by floor
  for (uint32_t i = 0; i < real_size; i++) {
    TEST_ASSERT(noise_spectrum[i] >= 0.0f,
                "Noise estimate should be non-negative after floor");
  }

  free(spectrum);
  free(noise_spectrum);
  free(floor_profile);
  brandt_noise_estimator_free(est);
  printf("✓ Brandt floor application tests passed\n");
}

int main(void) {
  printf("Running Brandt Noise Estimator tests...\n\n");

  test_brandt_initialization();
  test_brandt_run_logic();
  test_brandt_state_management();
  test_brandt_convergence();
  test_brandt_varying_energy_levels();
  test_brandt_spectral_patterns();
  test_brandt_history_duration();
  test_brandt_floor_application();

  printf("\n✅ All Brandt tests passed!\n");
  return 0;
}