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

  NoiseFloorManager* nfm = noise_floor_manager_initialize(1024);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

  noise_floor_manager_free(NULL); // Coverage: free NULL
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager lifecycle tests passed\n");
}

void test_noise_floor_manager_apply(void) {
  printf("Testing Noise Floor Manager apply...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

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

  // Coverage: Mismatched sizes
  // Force 0.0 gain, apply 0.1 reduction, with wrong sizes given.
  // It should still process using the correct internal sizes and output 0.1
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }
  noise_floor_manager_apply(nfm, 999, 999, gain_spectrum, noise_profile, 0.1f,
                            0.1f, NULL, 0.0f);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);
  TEST_FLOAT_CLOSE(gain_spectrum[fft_size - 1], 0.1f,
                   0.001f); // symmetric copy checks out

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

void test_noise_floor_manager_spectral_shaping(void) {
  printf("Testing Noise Floor Manager spectral shaping...\n");

  uint32_t fft_size = 512;
  uint32_t real_size = (fft_size / 2) + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  // Create a non-flat noise profile (shaped spectrum)
  for (uint32_t k = 0; k < real_size; k++) {
    // Decreasing noise floor with frequency (like pink noise)
    noise_profile[k] = 1.0f / sqrtf((float)(k + 1));
  }

  // Zero initial gain
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }

  // Apply with whitening to see spectral shaping
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 1.0f);

  // Verify floor was applied
  TEST_ASSERT(gain_spectrum[0] >= 0.09f && gain_spectrum[0] <= 0.11f,
              "Floor should be applied at low frequency");

  // Verify symmetry in FFT
  for (uint32_t k = 1; k < fft_size - k; k++) {
    TEST_FLOAT_CLOSE(gain_spectrum[k], gain_spectrum[fft_size - k], 1e-5f);
  }

  free(gain_spectrum);
  free(noise_profile);
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager spectral shaping tests passed\n");
}

void test_noise_floor_manager_tonal_handoff(void) {
  printf("Testing Noise Floor Manager tonal handoff logic...\n");

  uint32_t fft_size = 256;
  uint32_t real_size = (fft_size / 2) + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));
  float* tonal_mask = (float*)calloc(real_size, sizeof(float));

  for (uint32_t k = 0; k < real_size; k++) {
    noise_profile[k] = 1.0f;
    gain_spectrum[k] = 0.0f;
    // Create mask with some tonal components
    tonal_mask[k] = (k % 10 == 0) ? 1.0f : 0.0f;
  }

  // Test with 0% whitening (tonal reduction should be preserved)
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.2f, 0.5f, tonal_mask, 0.0f);

  // Bins with mask=1.0 should follow tonal reduction (0.5)
  // Bins with mask=0.0 should follow broadband reduction (0.2)
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.5f, 0.01f);  // Tonal bin
  TEST_FLOAT_CLOSE(gain_spectrum[1], 0.2f, 0.01f);  // Non-tonal bin

  // Test with 100% whitening (should flatten to broadband)
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.2f, 0.5f, tonal_mask, 1.0f);

  // All bins should converge to broadband reduction
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.2f, 0.01f);
  TEST_FLOAT_CLOSE(gain_spectrum[1], 0.2f, 0.01f);

  free(gain_spectrum);
  free(noise_profile);
  free(tonal_mask);
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager tonal handoff tests passed\n");
}

void test_noise_floor_manager_dynamic_whitening(void) {
  printf("Testing Noise Floor Manager dynamic whitening depth...\n");

  uint32_t fft_size = 128;
  uint32_t real_size = (fft_size / 2) + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  for (uint32_t k = 0; k < real_size; k++) {
    noise_profile[k] = 1.0f;
  }

  // Test: At 0dB reduction (1.0), whitening should have no effect
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.8f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 1.0f, 1.0f, NULL, 1.0f);

  // Transparency guard should force unity gain
  TEST_FLOAT_CLOSE(gain_spectrum[0], 1.0f, 0.001f);

  // Test: At maximum reduction, whitening should be at full strength
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.0f, 0.0f, NULL, 1.0f);

  // Should apply full whitening effect
  TEST_ASSERT(gain_spectrum[0] >= 0.0f && gain_spectrum[0] <= 0.1f,
              "Full reduction should apply strong floor");

  free(gain_spectrum);
  free(noise_profile);
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager dynamic whitening tests passed\n");
}

void test_noise_floor_manager_edge_cases(void) {
  printf("Testing Noise Floor Manager edge cases...\n");

  uint32_t fft_size = 64;
  uint32_t real_size = (fft_size / 2) + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  // Test with very low noise profile values
  for (uint32_t k = 0; k < real_size; k++) {
    noise_profile[k] = 1e-10f;
    gain_spectrum[k] = 0.5f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.5f, 0.5f, NULL, 0.5f);
  TEST_ASSERT(gain_spectrum[0] >= 0.4f && gain_spectrum[0] <= 0.6f,
              "Should handle very low noise profile");

  // Test with very high noise profile values
  for (uint32_t k = 0; k < real_size; k++) {
    noise_profile[k] = 1e10f;
    gain_spectrum[k] = 0.5f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.5f, 0.5f, NULL, 0.5f);
  TEST_ASSERT(gain_spectrum[0] >= 0.0f && gain_spectrum[0] <= 1.5f,
              "Should handle very high noise profile");

  // Test with gain already above floor
  for (uint32_t k = 0; k < real_size; k++) {
    noise_profile[k] = 1.0f;
    gain_spectrum[k] = 0.9f; // Already high
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 0.0f);
  // Gain should remain unchanged (no clamping up needed)
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.9f, 0.001f);

  free(gain_spectrum);
  free(noise_profile);
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager edge case tests passed\n");
}

int main(void) {
  printf("Running Noise Floor Manager tests...\n\n");

  test_noise_floor_manager_lifecycle();
  test_noise_floor_manager_apply();
  test_noise_floor_manager_spectral_shaping();
  test_noise_floor_manager_tonal_handoff();
  test_noise_floor_manager_dynamic_whitening();
  test_noise_floor_manager_edge_cases();

  printf("\n✅ All Noise Floor Manager tests passed!\n");
  return 0;
}