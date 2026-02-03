/*
 * Unit tests for Gain Estimation
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/gain_estimation/gain_estimators.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_gain_estimation_wiener(void) {
  printf("Testing Wiener gain estimation...\n");

  uint32_t fft_size = 8;
  uint32_t real_spectrum_size = fft_size / 2 + 1;

  float spectrum[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 0.0f, 0.0f, 0.0f};
  float noise_spectrum[8] = {0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 0.0f, 0.0f, 0.0f};
  float gain_spectrum[8] = {0.0f};
  float alpha[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  float beta[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

  estimate_gains(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                 gain_spectrum, alpha, beta, WIENER);

  // Wiener in this lib: gain = (spectrum - noise) / spectrum
  // For each bin: gain = (spectrum - noise) / spectrum
  // gain[0] = (1.0 - 0.5) / 1.0 = 0.5
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.5f, 0.01f);

  printf("✓ Wiener gain estimation tests passed\n");
}

void test_gain_estimation_gates(void) {
  printf("Testing Gates gain estimation...\n");

  uint32_t fft_size = 8;
  uint32_t real_spectrum_size = fft_size / 2 + 1;

  float spectrum[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 0.0f, 0.0f, 0.0f};
  float noise_spectrum[8] = {0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 0.0f, 0.0f, 0.0f};
  float gain_spectrum[8] = {0.0f};
  float alpha[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  float beta[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

  estimate_gains(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                 gain_spectrum, alpha, beta, GATES);

  // Gates: gain = 1 if spectrum > noise, 0 otherwise
  // gain[0] = 1.0 > 0.5 = 1.0
  // gain[1] = 2.0 > 1.0 = 1.0
  // gain[2] = 3.0 > 1.5 = 1.0
  // gain[3] = 4.0 > 2.0 = 1.0
  // gain[4] = 5.0 > 2.5 = 1.0
  for (int i = 0; i < 5; i++) {
    TEST_FLOAT_CLOSE(gain_spectrum[i], 1.0f, 0.001f);
  }

  printf("✓ Gates gain estimation tests passed\n");
}

void test_gain_estimation_spectral_subtraction(void) {
  printf("Testing Generalized Spectral Subtraction gain estimation...\n");

  uint32_t fft_size = 8;
  uint32_t real_spectrum_size = fft_size / 2 + 1;

  float spectrum[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 0.0f, 0.0f, 0.0f};
  float noise_spectrum[8] = {0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 0.0f, 0.0f, 0.0f};
  float gain_spectrum[8] = {0.0f};
  float alpha[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  float beta[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

  estimate_gains(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                 gain_spectrum, alpha, beta, GENERALIZED_SPECTRALSUBTRACTION);

  // Generalized spectral subtraction: gain = sqrt(max(0, (spectrum^2 -
  // alpha*noise^2) / spectrum^2)) For each bin: gain = sqrt(max(0, (spectrum^2
  // - alpha*noise^2) / spectrum^2)) gain[0] = sqrt(max(0, (1.0^2 - 1.0*0.5^2)
  // / 1.0^2)) = sqrt(max(0, (1.0 - 0.25) / 1.0)) = sqrt(0.75) ≈ 0.866
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.866f, 0.01f);

  printf("✓ Generalized Spectral Subtraction gain estimation tests passed\n");
}

void test_gain_estimation_edge_cases(void) {
  printf("Testing gain estimation edge cases...\n");

  uint32_t fft_size = 8;
  uint32_t real_spectrum_size = fft_size / 2 + 1;

  float spectrum[8] = {0.0f, 0.0f, 1.0f, 2.0f, 3.0f, 0.0f, 0.0f, 0.0f};
  float noise_spectrum[8] = {0.0f, 1.0f, 1.0f, 2.0f, 2.0f, 0.0f, 0.0f, 0.0f};
  float gain_spectrum[8] = {0.0f};
  float alpha[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  float beta[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

  // Test with zero spectrum values
  estimate_gains(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                 gain_spectrum, alpha, beta, WIENER);

  // Should handle division by zero gracefully
  TEST_ASSERT(gain_spectrum[0] >= 0.0f && gain_spectrum[0] <= 1.0f,
              "Gain should be in valid range");
  TEST_ASSERT(gain_spectrum[1] >= 0.0f && gain_spectrum[1] <= 1.0f,
              "Gain should be in valid range");

  printf("✓ Gain estimation edge cases tests passed\n");
}

int main(void) {
  printf("Running gain estimation tests...\n");

  test_gain_estimation_wiener();
  test_gain_estimation_gates();
  test_gain_estimation_spectral_subtraction();
  test_gain_estimation_edge_cases();

  printf("✅ All gain estimation tests passed!\n");
  return 0;
}
