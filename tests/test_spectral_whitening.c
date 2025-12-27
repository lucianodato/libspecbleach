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
  uint32_t real_size = fft_size / 2 + 1;
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

int main(void) {
  printf("Running Spectral Whitening tests...\n\n");

  test_whitening_lifecycle();
  test_whitening_get_weights();

  printf("\n✅ All Spectral Whitening tests passed!\n");
  return 0;
}
