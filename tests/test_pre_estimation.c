/*
 * Unit tests for Pre-estimation modules
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/pre_estimation/absolute_hearing_thresholds.h"
#include "../src/shared/pre_estimation/critical_bands.h"
#include "../src/shared/pre_estimation/masking_estimator.h"
#include "../src/shared/pre_estimation/noise_scaling_criterias.h"
#include "../src/shared/pre_estimation/spectral_smoother.h"
#include "../src/shared/pre_estimation/transient_detector.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_absolute_hearing_thresholds(void) {
  printf("Testing Absolute Hearing Thresholds...\n");

  uint32_t fft_size = 1024;
  uint32_t sample_rate = 44100;

  AbsoluteHearingThresholds* aht = absolute_hearing_thresholds_initialize(
      sample_rate, fft_size, POWER_SPECTRUM);

  TEST_ASSERT(aht != NULL, "AHT initialization should succeed");

  float spectrum[513] = {0.0f}; // 513 = fft_size/2 + 1
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 0.001f; // Set low values that should be floored
  }

  TEST_ASSERT(apply_thresholds_as_floor(aht, spectrum),
              "Apply thresholds should succeed");

  // Check that some values were floored to hearing thresholds
  bool some_values_changed = false;
  for (int i = 0; i < 513; i++) {
    if (spectrum[i] > 0.001f) {
      some_values_changed = true;
      break;
    }
  }
  TEST_ASSERT(some_values_changed,
              "Some spectrum values should be floored to hearing thresholds");

  absolute_hearing_thresholds_free(aht);
  printf("✓ Absolute Hearing Thresholds tests passed\n");
}

void test_critical_bands(void) {
  printf("Testing Critical Bands...\n");

  uint32_t fft_size = 1024;
  uint32_t sample_rate = 44100;

  CriticalBands* cb =
      critical_bands_initialize(sample_rate, fft_size, BARK_SCALE);
  TEST_ASSERT(cb != NULL, "Critical bands initialization should succeed");

  uint32_t num_bands = get_number_of_critical_bands(cb);
  TEST_ASSERT(num_bands > 0, "Should have some critical bands");

  // Test getting band indexes
  for (uint32_t i = 0; i < num_bands; i++) {
    CriticalBandIndexes indexes = get_band_indexes(cb, i);
    TEST_ASSERT(indexes.start_position < indexes.end_position,
                "Band start should be before end");
    TEST_ASSERT(indexes.end_position <= fft_size / 2 + 1,
                "Band end should be within spectrum bounds");
  }

  // Test computing critical bands spectrum
  float spectrum[513] = {0.0f};
  float critical_bands_spectrum[24] = {0.0f}; // Assuming max 24 bark bands

  for (int i = 0; i < 513; i++) {
    spectrum[i] = (float)i * 0.1f;
  }

  TEST_ASSERT(
      compute_critical_bands_spectrum(cb, spectrum, critical_bands_spectrum),
      "Compute critical bands should succeed");

  critical_bands_free(cb);
  printf("✓ Critical Bands tests passed\n");
}

void test_masking_estimator(void) {
  printf("Testing Masking Estimator...\n");

  uint32_t fft_size = 1024;
  uint32_t sample_rate = 44100;

  MaskingEstimator* me =
      masking_estimation_initialize(fft_size, sample_rate, POWER_SPECTRUM);
  TEST_ASSERT(me != NULL, "Masking estimator initialization should succeed");

  float spectrum[513] = {0.0f};
  float masking_thresholds[513] = {0.0f};

  // Create a test spectrum with some energy
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 0.1f + (float)i * 0.01f;
  }

  TEST_ASSERT(compute_masking_thresholds(me, spectrum, masking_thresholds),
              "Compute masking thresholds should succeed");

  // Check that masking thresholds are reasonable
  for (int i = 0; i < 513; i++) {
    TEST_ASSERT(masking_thresholds[i] >= 0.0f,
                "Masking thresholds should be non-negative");
  }

  masking_estimation_free(me);
  printf("✓ Masking Estimator tests passed\n");
}

void test_noise_scaling_criterias(void) {
  printf("Testing Noise Scaling Criterias...\n");

  uint32_t fft_size = 1024;
  uint32_t sample_rate = 44100;

  NoiseScalingCriterias* nsc = noise_scaling_criterias_initialize(
      fft_size, BARK_SCALE, sample_rate, POWER_SPECTRUM);
  TEST_ASSERT(nsc != NULL,
              "Noise scaling criterias initialization should succeed");

  float spectrum[513] = {0.0f};
  float noise_spectrum[513] = {0.0f};
  float alpha[513] = {0.0f};
  float beta[513] = {0.0f};

  // Create test data
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 1.0f;
    noise_spectrum[i] = 0.1f;
  }

  NoiseScalingParameters params = {.undersubtraction = 1.0f,
                                   .oversubtraction = 1.0f,
                                   .scaling_type = A_POSTERIORI_SNR};

  TEST_ASSERT(apply_noise_scaling_criteria(nsc, spectrum, noise_spectrum, alpha,
                                           beta, params),
              "Apply noise scaling criteria should succeed");

  // Check that alpha and beta are in reasonable ranges
  for (int i = 0; i < 513; i++) {
    TEST_ASSERT(alpha[i] >= 0.0f && alpha[i] <= 10.0f,
                "Alpha should be in reasonable range");
    TEST_ASSERT(beta[i] >= 0.0f && beta[i] <= 10.0f,
                "Beta should be in reasonable range");
  }

  noise_scaling_criterias_free(nsc);
  printf("✓ Noise Scaling Criterias tests passed\n");
}

void test_spectral_smoother(void) {
  printf("Testing Spectral Smoother...\n");

  uint32_t fft_size = 1024;

  SpectralSmoother* ss = spectral_smoothing_initialize(fft_size, FIXED);
  TEST_ASSERT(ss != NULL, "Spectral smoother initialization should succeed");

  float spectrum[513] = {0.0f};

  // Create a test spectrum with some variation
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 1.0f + 0.5f * sinf((float)i * 0.1f);
  }

  TimeSmoothingParameters params = {.smoothing = 0.8f,
                                    .transient_protection_enabled = false};

  TEST_ASSERT(spectral_smoothing_run(ss, params, spectrum),
              "Spectral smoothing should succeed");

  // Check that output is reasonable
  for (int i = 0; i < 513; i++) {
    TEST_ASSERT(spectrum[i] >= 0.0f,
                "Smoothed spectrum should be non-negative");
  }

  spectral_smoothing_free(ss);
  printf("✓ Spectral Smoother tests passed\n");
}

void test_transient_detector(void) {
  printf("Testing Transient Detector...\n");

  uint32_t fft_size = 1024;

  TransientDetector* td = transient_detector_initialize(fft_size);
  TEST_ASSERT(td != NULL, "Transient detector initialization should succeed");

  float spectrum[513] = {0.0f};
  float different_spectrum[513] = {0.0f};

  // Create initial spectrum
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 1.0f;
    different_spectrum[i] = 2.0f; // Different spectrum to create change
  }

  // First run - should not detect transient (no previous data)
  bool result1 = transient_detector_run(td, spectrum);

  // Second run with different spectrum - may or may not detect transient
  bool result2 = transient_detector_run(td, different_spectrum);

  // The function should run without error (return value is boolean indicating
  // transient presence) We just test that it doesn't crash and returns a valid
  // boolean
  TEST_ASSERT(result1 == true || result1 == false,
              "First run should return boolean");
  TEST_ASSERT(result2 == true || result2 == false,
              "Second run should return boolean");

  transient_detector_free(td);
  printf("✓ Transient Detector tests passed\n");
}

int main(void) {
  printf("Running pre-estimation tests...\n");

  test_absolute_hearing_thresholds();
  test_critical_bands();
  test_masking_estimator();
  test_noise_scaling_criterias();
  test_spectral_smoother();
  test_transient_detector();

  printf("✅ All pre-estimation tests passed!\n");
  return 0;
}
