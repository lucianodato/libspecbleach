/*
 * Unit tests for Spectral Whitening
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/processing/spectral_whitening.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_whitening_lifecycle(void) {
  printf("Testing Spectral Whitening lifecycle...\n");

  SpectralWhitening* sw = spectral_whitening_initialize(1024);
  TEST_ASSERT(sw != NULL, "Initialization should succeed");

  spectral_whitening_free(NULL); // Coverage: free NULL
  spectral_whitening_free(sw);
  printf("✓ Spectral Whitening lifecycle tests passed\n");
}

void test_whitening_get_ideal_reduction_db(void) {
  printf("Testing Spectral Whitening get_ideal_reduction_db...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  SpectralWhitening* sw = spectral_whitening_initialize(fft_size);

  float* weights = (float*)calloc(real_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  // Initialize noise profile with a ramp
  for (uint32_t k = 0; k < real_size; k++) {
    noise_profile[k] = 1.0f + (float)k;
  }

  // Coverage: NULL inputs
  spectral_whitening_get_ideal_reduction_db(NULL, 0.1778f, noise_profile, weights);
  spectral_whitening_get_ideal_reduction_db(sw, 0.1778f, NULL, weights);
  spectral_whitening_get_ideal_reduction_db(sw, 0.1778f, noise_profile, NULL);

  // Test with standard reduction limit (15 dB / 0.1778f)
  spectral_whitening_get_ideal_reduction_db(sw, 0.1778f, noise_profile, weights);

  // Values above median (e.g. k=512, val=513 > median 257) get reduction weight > 0.0.
  // Values below median (e.g. k=0, val=1 < median 257) get reduction weight == 0.0 (left alone).
  TEST_ASSERT(weights[real_size - 1] > 0.0f,
              "Bins above median should get reduction weight > 0.0");
  TEST_ASSERT(weights[0] == 0.0f,
              "Bins below median should get reduction weight == 0.0 (left alone)");

  // Test with zero noise profile (1e-12 floor handling)
  memset(noise_profile, 0, real_size * sizeof(float));
  spectral_whitening_get_ideal_reduction_db(sw, 0.1778f, noise_profile, weights);

  free(weights);
  free(noise_profile);
  spectral_whitening_free(sw);
  printf("✓ Spectral Whitening get_ideal_reduction_db tests passed\n");
}

void test_whitening_tonal_peak(void) {
  printf("Testing Spectral Whitening with tonal peak...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  SpectralWhitening* sw = spectral_whitening_initialize(fft_size);

  float* weights = (float*)calloc(real_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  // Create a profile: Mostly quiet (0.1), one BIG peak (10.0) at index 100
  for (uint32_t i = 0; i < real_size; i++) {
    noise_profile[i] = 0.1f;
  }
  noise_profile[100] = 10.0f; // The hum/tonal noise

  // Standard Reduction (15 dB / 0.1778f)
  spectral_whitening_get_ideal_reduction_db(sw, 0.1778f, noise_profile, weights);

  // Quiet bins at or under median anchor level (e.g. 200, val=0.1) should have 0.0 reduction depth weight at whitening=1.0 (left alone)
  TEST_ASSERT(weights[200] == 0.0f,
              "Quiet bins at anchor level should have 0 reduction depth weight");

  // The Hum spike (index 100, val=10.0) is far above median anchor level (0.1), so h = 1 - 0.1/10 = 0.9
  // Weight should be close to 0.99
  TEST_ASSERT(weights[100] > 0.9f,
              "Hum spike should receive high reduction depth weight");

  free(weights);
  free(noise_profile);
  spectral_whitening_free(sw);
  printf("✓ Spectral Whitening tonal peak tests passed\n");
}

void test_whitening_band_limited(void) {
  printf("Testing Spectral Whitening with band-limited audio...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  SpectralWhitening* sw = spectral_whitening_initialize(fft_size);

  float* weights = (float*)calloc(real_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  // Band-limited audio: Only bins 0..100 have noise energy (0.1f).
  // Upper sub-bands (bins 101..512, >80% of spectrum) are zero.
  for (uint32_t i = 0; i <= 100; i++) {
    noise_profile[i] = 0.1f;
  }
  // Add a hum spike at bin 30 inside the passband
  noise_profile[30] = 1.0f;

  spectral_whitening_get_ideal_reduction_db(sw, 0.1778f, noise_profile, weights);

  // Bins at baseline passband level (0.1f) are at/under median, so reduction depth weight = 0.0f (left alone)
  TEST_ASSERT(weights[50] == 0.0f,
              "Active passband bins at median should have 0 reduction weight");

  // Hum spike in active passband (1.0f vs 0.1f) should receive high reduction depth weight (h = 0.9)
  TEST_ASSERT(weights[30] > 0.8f,
              "Tonal peak in band-limited profile should receive high reduction weight");

  // Inactive stopband bins should have 0.0 reduction depth weight
  TEST_ASSERT(weights[300] == 0.0f,
              "Stopband bins should have 0 reduction weight");

  free(weights);
  free(noise_profile);
  spectral_whitening_free(sw);
  printf("✓ Spectral Whitening band-limited tests passed\n");
}

int main(void) {
  printf("Running Spectral Whitening tests...\n\n");

  test_whitening_lifecycle();
  test_whitening_get_ideal_reduction_db();
  test_whitening_tonal_peak();
  test_whitening_band_limited();

  printf("\n✅ All Spectral Whitening tests passed!\n");
  return 0;
}
