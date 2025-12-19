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

  SpectralWhitening* sw = spectral_whitening_initialize(1024, 44100, 256);
  TEST_ASSERT(sw != NULL, "Initialization should succeed");

  spectral_whitening_free(sw);
  printf("✓ Spectral Whitening lifecycle tests passed\n");
}

void test_whitening_run(void) {
  printf("Testing Spectral Whitening run...\n");

  uint32_t fft_size = 8;
  SpectralWhitening* sw = spectral_whitening_initialize(fft_size, 44100, 256);

  float spectrum[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};

  // Test with whitening factor 0 (no effect expected on gain calculation logic,
  // but it should run)
  TEST_ASSERT(spectral_whitening_run(sw, 0.0f, spectrum), "Run should succeed");

  // Test with whitening factor 1.0
  TEST_ASSERT(spectral_whitening_run(sw, 1.0f, spectrum), "Run should succeed");

  spectral_whitening_free(sw);
  printf("✓ Spectral Whitening run tests passed\n");
}

int main(void) {
  printf("Running Spectral Whitening tests...\n\n");

  test_whitening_lifecycle();
  test_whitening_run();

  printf("\n✅ All Spectral Whitening tests passed!\n");
  return 0;
}
