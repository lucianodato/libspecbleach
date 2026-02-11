/*
 * Unit tests for Noise Floor Manager
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/denoiser_logic/core/noise_floor_manager.h"

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
  uint32_t real_size = (fft_size / 2) + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size, 44100, 256);

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.5f;
  }
  for (uint32_t k = 0; k < real_size; k++) {
    noise_profile[k] = 1.0f;
  }

  // Coverage: NULL inputs
  noise_floor_manager_apply(NULL, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 0.5f);
  noise_floor_manager_apply(nfm, real_size, fft_size, NULL, noise_profile, 0.1f,
                            0.1f, NULL, 0.5f);
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum, NULL, 0.1f,
                            0.1f, NULL, 0.5f);

  // Test with 1.0 linear reduction (should result in 1.0 floor/gain)
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.5f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 1.0f, 1.0f, NULL, 0.0f);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 1.0f, 0.001f);

  // Test with 0.1 linear reduction (equivalent to 20dB)
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f; // Original gain is 0
  }
  // No whitening (0.0) -> floor should be 0.1
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 0.0f);
  // floor + (1-floor)*gain = 0.1 + 0.9*0 = 0.1
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);

  // Test with whitening (Max anchoring)
  // Profile is all 1.0, so Max=1.0. Whitney factor 100% -> weight=1.0 for all.
  // Result should still be 0.1 because profile is already flat.
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 1.0f);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);

  // Test Clamping: Gain is 0.05, Floor is 0.1 -> Result should be 0.1
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.05f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 0.0f);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);

  // Test Transparency Guard: Even if original gain is low, 0dB reduction
  // FORCES 1.0
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.5f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 1.0f, 1.0f, NULL, 1.0f);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 1.0f, 0.001f);

  // Test RX Alignment: 100% Whitening overrides tonal path
  // Broadband reduction = 20dB (0.1), Tonal reduction = 0dB (1.0)
  // Mask is 1.0 (pure tone). At 100% whitening, it should follow broadband
  // (0.1)
  float* tonal_mask_test = (float*)calloc(real_size, sizeof(float));
  for (uint32_t k = 0; k < real_size; k++) {
    tonal_mask_test[k] = 1.0f;
  }
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 1.0f, tonal_mask_test, 1.0f);
  // target_reduction should be 0.1 because whitening is 1.0
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);
  free(tonal_mask_test);

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
