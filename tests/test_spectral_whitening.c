/*
 * Unit tests for Spectral Whitening
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/post_estimation/spectral_whitening.h"

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

void test_whitening_get_weights(void) {
  printf("Testing Spectral Whitening get_weights...\n");

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
  spectral_whitening_get_weights(NULL, 0.5f, noise_profile, weights);
  spectral_whitening_get_weights(sw, 0.5f, NULL, weights);
  spectral_whitening_get_weights(sw, 0.5f, noise_profile, NULL);

  // Test with whitening factor 0 (should be just tapering window)
  spectral_whitening_get_weights(sw, 0.0f, noise_profile, weights);
  TEST_ASSERT(weights[0] > 0.0f, "Weight should be positive");

  // Test with whitening factor 1.0 (full whitening)
  spectral_whitening_get_weights(sw, 1.0f, noise_profile, weights);

  // With Median anchoring, values above median get weights < 1.0, and vice
  // versa. Profile is a ramp [1..513]. Median is ~257. Weight at k=0 (val=1)
  // should be 257/1 = 257. Weight at k=512 (val=513) should be 257/513 = 0.5.
  // crucially, the results are NOT tapered by Hamming.
  TEST_ASSERT(weights[real_size - 1] < 1.0f,
              "HF should not be boosted if above median");
  TEST_ASSERT(weights[0] > 1.0f, "LF should be boosted if below median");

  // Test with very low whitening factor
  spectral_whitening_get_weights(sw, 0.01f, noise_profile, weights);

  // Test with zero noise profile (1e-12 floor handling)
  memset(noise_profile, 0, real_size * sizeof(float));
  spectral_whitening_get_weights(sw, 1.0f, noise_profile, weights);

  free(weights);
  free(noise_profile);
  spectral_whitening_free(sw);
  printf("✓ Spectral Whitening get_weights tests passed\n");
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

  // Whitening Factor 1.0 (100% - function uses linear scale)
  spectral_whitening_get_weights(sw, 1.0f, noise_profile, weights);

  // Check floor level at a quiet bin (e.g., 200)
  // Floor = Profile * Weight
  // With Median anchoring (0.1), the weight at index 200 (val=0.1) should
  // be 1.0. With Tapering REMOVED, the weight should be exactly 1.0 even at HF.
  TEST_ASSERT(weights[200] >= 0.99f && weights[200] <= 1.01f,
              "Quiet bins at anchor level should have unit weight");

  // Check the Hum spike (index 100, val=10.0)
  // Weight should be Median/Profile = 0.1 / 10.0 = 0.01 (approx)
  // This 'shaves' the peak in the floor calculation.
  float floor_at_hum = noise_profile[100] * weights[100];
  TEST_ASSERT(floor_at_hum <= 0.2f,
              "Hum spike should be flattened in the floor");

  // Check the overall flat floor level (should be approx the median = 0.1)
  TEST_ASSERT(floor_at_hum >= 0.05f,
              "Hum spike should not be over-suppressed in floor");
  TEST_FLOAT_CLOSE(noise_profile[200] * weights[200], 0.1f, 0.02f);

  free(weights);
  free(noise_profile);
  spectral_whitening_free(sw);
  printf("✓ Spectral Whitening tonal peak tests passed\n");
}

int main(void) {
  printf("Running Spectral Whitening tests...\n\n");

  test_whitening_lifecycle();
  test_whitening_get_weights();
  test_whitening_tonal_peak();

  printf("\n✅ All Spectral Whitening tests passed!\n");
  return 0;
}
