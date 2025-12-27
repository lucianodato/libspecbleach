/*
 * Unit tests for Noise Floor Manager
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/post_estimation/noise_floor_manager.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_noise_floor_manager_lifecycle(void) {
  printf("Testing Noise Floor Manager lifecycle...\n");

  NoiseFloorManager* nfm = noise_floor_manager_initialize(1024, 44100, 256);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager lifecycle tests passed\n");
}

void test_noise_floor_manager_apply(void) {
  printf("Testing Noise Floor Manager apply...\n");

  uint32_t fft_size = 8;
  uint32_t real_size = fft_size / 2 + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size, 44100, 256);

  float gain_spectrum[8] = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
  float noise_profile[5] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

  // Test with reduction amount 0 and whitening factor 0
  // Gains should remain unchanged (except for tapering)
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.0f, 0.0f);

  // Since reduction_amount is 0 and whitening is 0, floor is 0.
  // gain = 0 + (1-0)*gain = gain.
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.5f, 0.01f);

  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager apply tests passed\n");
}

int main(void) {
  printf("Running Noise Floor Manager tests...\n\n");

  test_noise_floor_manager_lifecycle();
  test_noise_floor_manager_apply();

  printf("\n✅ All Noise Floor Manager tests passed!\n");
  return 0;
}
