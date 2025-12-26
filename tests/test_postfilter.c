/*
 * Unit tests for Post-filter
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/post_estimation/postfilter.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_postfilter_lifecycle(void) {
  printf("Testing Post-filter lifecycle...\n");

  uint32_t fft_size = 1024;

  PostFilter* pf = postfilter_initialize(fft_size);
  TEST_ASSERT(pf != NULL, "Post-filter initialization should succeed");

  postfilter_free(pf);
  printf("✓ Post-filter lifecycle tests passed\n");
}

void test_postfilter_apply(void) {
  printf("Testing Post-filter apply...\n");

  uint32_t fft_size = 1024;

  PostFilter* pf = postfilter_initialize(fft_size);
  TEST_ASSERT(pf != NULL, "Post-filter initialization should succeed");

  float spectrum[1024] = {0.0f};
  float gain_spectrum[1024] = {0.0f};

  // Create test data (only fill real spectrum part)
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 1.0f + (float)i * 0.01f;
    gain_spectrum[i] = 0.8f; // 80% gain
  }

  PostFiltersParameters params = {
      .snr_threshold = -10.0f // -10 dB SNR threshold
  };

  TEST_ASSERT(postfilter_apply(pf, spectrum, gain_spectrum, params),
              "Post-filter apply should succeed");

  // Check that gains are in valid range [0, 1]
  for (int i = 0; i < 1024; i++) {
    TEST_ASSERT(gain_spectrum[i] >= 0.0f && gain_spectrum[i] <= 1.0f,
                "Post-filter gain should be in [0, 1] range");
  }

  postfilter_free(pf);
  printf("✓ Post-filter apply tests passed\n");
}

void test_postfilter_parameters(void) {
  printf("Testing Post-filter with different parameters...\n");

  uint32_t fft_size = 1024;

  PostFilter* pf = postfilter_initialize(fft_size);
  TEST_ASSERT(pf != NULL, "Post-filter initialization should succeed");

  float spectrum[1024] = {0.0f};
  float gain_spectrum[1024] = {0.0f};

  // Create test spectrum with high SNR (only fill real spectrum part)
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 10.0f; // High signal
    gain_spectrum[i] = 1.0f;
  }

  // Test with very low SNR threshold (should apply more filtering)
  PostFiltersParameters params_strict = {.snr_threshold = -20.0f};

  TEST_ASSERT(postfilter_apply(pf, spectrum, gain_spectrum, params_strict),
              "Post-filter apply with strict threshold should succeed");

  // Test with high SNR threshold (should apply less filtering)
  PostFiltersParameters params_lenient = {.snr_threshold = 10.0f};

  // Reset gain spectrum
  for (int i = 0; i < 1024; i++) {
    gain_spectrum[i] = 1.0f;
  }

  TEST_ASSERT(postfilter_apply(pf, spectrum, gain_spectrum, params_lenient),
              "Post-filter apply with lenient threshold should succeed");

  postfilter_free(pf);
  printf("✓ Post-filter parameter tests passed\n");
}

int main(void) {
  printf("Running post-filter tests...\n");

  test_postfilter_lifecycle();
  test_postfilter_apply();
  test_postfilter_parameters();

  // Test edge cases for coverage
  printf("Testing Post-filter edge cases...\n");
  uint32_t fft_size = 1024;
  PostFilter* pf = postfilter_initialize(fft_size);
  float spectrum[513] = {0.0f};
  float gain_spectrum[513] = {1.0f};
  PostFiltersParameters params = {.snr_threshold = 0.5f};

  // 1. Zero energy spectrum (hits line 81 in postfilter.c)
  TEST_ASSERT(postfilter_apply(pf, spectrum, gain_spectrum, params),
              "Zero energy apply should succeed");

  // 2. High SNR / Unit gain (hits line 88 in postfilter.c)
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 10.0f;
    gain_spectrum[i] = 1.0f;
  }
  TEST_ASSERT(postfilter_apply(pf, spectrum, gain_spectrum, params),
              "Unit gain apply should succeed");

  // 3. NULL arguments (hits line 141 in postfilter.c)
  TEST_ASSERT(postfilter_apply(NULL, spectrum, gain_spectrum, params) == false,
              "NULL filter should fail");
  TEST_ASSERT(postfilter_apply(pf, NULL, gain_spectrum, params) == false,
              "NULL spectrum should fail");
  // 4. Large window size / Small spectrum (hits moving_average boundaries)
  // We can't directly call moving_average (static), so we trigger it via
  // postfilter_apply
  PostFilter* pf_small = postfilter_initialize(16); // small real_spectrum_size
  float spectrum_small[16] = {1.0f};
  float gain_small[16] = {0.5f};
  PostFiltersParameters params_strict = {.snr_threshold = 1.0f};
  TEST_ASSERT(
      postfilter_apply(pf_small, spectrum_small, gain_small, params_strict),
      "Large window apply should succeed");
  postfilter_free(pf_small);

  // 5. Gain floor test
  for (int i = 0; i < 513; i++) {
    gain_spectrum[i] = 0.0f; // Force min gain floor
  }
  TEST_ASSERT(postfilter_apply(pf, spectrum, gain_spectrum, params),
              "Gain floor apply should succeed");
  for (int i = 0; i < 513; i++) {
    TEST_ASSERT(gain_spectrum[i] > 0.0f, "Gain should be floored");
  }

  postfilter_free(pf);
  printf("✓ Post-filter edge case tests passed\n");

  printf("✅ All post-filter tests passed!\n");
  return 0;
}
