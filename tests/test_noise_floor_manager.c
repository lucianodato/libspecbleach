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

  noise_floor_manager_free(NULL); // Coverage: free NULL
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager lifecycle tests passed\n");
}

void test_noise_floor_manager_apply(void) {
  printf("Testing Noise Floor Manager apply...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = fft_size / 2 + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size, 44100, 256);

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  for (uint32_t k = 0; k < fft_size; k++)
    gain_spectrum[k] = 0.5f;
  for (uint32_t k = 0; k < real_size; k++)
    noise_profile[k] = 1.0f;

  // Coverage: NULL inputs
  noise_floor_manager_apply(NULL, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.5f);
  noise_floor_manager_apply(nfm, real_size, fft_size, NULL, noise_profile, 0.1f,
                            0.5f);
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum, NULL, 0.1f,
                            0.5f);

  // Test with large reduction amount
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 2.0f, 0.5f);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 1.0f, 0.001f); // Should be capped at 1.0

  // Test standard application
  for (uint32_t k = 0; k < fft_size; k++)
    gain_spectrum[k] = 0.1f;
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.5f);

  free(gain_spectrum);
  free(noise_profile);
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
