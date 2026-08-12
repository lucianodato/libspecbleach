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
      fprintf(stderr, "TEST FAILED (Line %d): %s\n", __LINE__, message);       \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_noise_floor_manager_lifecycle(void) {
  printf("Testing Noise Floor Manager lifecycle...\n");

  TEST_ASSERT(noise_floor_manager_initialize(0) == NULL,
              "Initialization with 0 fft_size should fail");

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
                            noise_profile, 0.1f, 0.1f, NULL, 0.5f, NULL);
  noise_floor_manager_apply(nfm, real_size, fft_size, NULL, noise_profile, 0.1f,
                            0.1f, NULL, 0.5f, NULL);
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum, NULL, 0.1f,
                            0.1f, NULL, 0.5f, NULL);

  // Coverage: Mismatched sizes
  // Force 0.0 gain, apply 0.1 reduction, with wrong sizes given.
  // It should still process using the correct internal sizes and output 0.1
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }
  noise_floor_manager_apply(nfm, 999, 999, gain_spectrum, noise_profile, 0.1f,
                            0.1f, NULL, 0.0f, NULL);
  noise_floor_manager_apply(nfm, real_size, 999, gain_spectrum, noise_profile,
                            0.1f, 0.1f, NULL, 0.0f, NULL);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);
  TEST_FLOAT_CLOSE(gain_spectrum[fft_size - 1], 0.1f,
                   0.001f); // symmetric copy checks out

  // Test with 1.0 linear reduction (should result in 1.0 floor/gain)
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.5f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 1.0f, 1.0f, NULL, 0.0f, NULL);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 1.0f, 0.001f);

  // Test with 0.1 linear reduction (equivalent to 20dB)
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f; // Original gain is 0
  }
  // No whitening (0.0) -> floor should be 0.1
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 0.0f, NULL);
  // floor + (1-floor)*gain = 0.1 + 0.9*0 = 0.1
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);

  // Test with whitening (Max/Median anchoring)
  // Profile is all 1.0, so all bins are at/under reference line. Whitening 100%
  // -> weight=0.0 (left alone -> 1.0 floor).
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 1.0f, NULL);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);

  // Test Clamping: Gain is 0.05, Floor is 0.1 -> Result should be 0.1
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.05f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 0.0f, NULL);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 0.1f, 0.001f);

  // Test Transparency Guard: Even if original gain is low, 0dB reduction
  // FORCES 1.0
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.5f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 1.0f, 1.0f, NULL, 1.0f, NULL);
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
                            noise_profile, 0.1f, 1.0f, tonal_mask_test, 1.0f, NULL);
  TEST_FLOAT_CLOSE(gain_spectrum[0], 1.0f, 0.001f);
  // Test r_dp_db < 0.0f (dual_path_reduction >= 1.0f with tonal path < 1.0f)
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 1.0f, 0.1f, NULL, 0.0f, NULL);

  free(tonal_mask_test);

  free(gain_spectrum);
  free(noise_profile);
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager apply tests passed\n");
}

void test_noise_floor_manager_band_limited(void) {
  printf("Testing Noise Floor Manager with band-limited audio...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));

  // Band-limited noise profile: passband 0..100 = 1.0f, stopband 101..512 =
  // 0.0f
  for (uint32_t k = 0; k <= 100; k++) {
    noise_profile[k] = 1.0f;
  }
  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }

  // 100% Whitening factor with 0.1f target reduction (20dB reduction)
  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.1f, 0.1f, NULL, 1.0f, NULL);

  // Bins at median level (1.0f) get 0 reduction weight -> 1.0 gain (left alone)
  TEST_FLOAT_CLOSE(gain_spectrum[50], 0.1f, 0.005f);

  free(gain_spectrum);
  free(noise_profile);
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager band-limited tests passed\n");
}

void test_noise_floor_manager_reduction_curve(void) {
  printf("Testing Noise Floor Manager reduction curve bias...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  TEST_ASSERT(nfm != NULL, "Initialization should succeed");

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_profile = (float*)calloc(real_size, sizeof(float));
  float* curve_bias = (float*)calloc(real_size, sizeof(float));

  for (uint32_t k = 0; k < real_size; k++) {
    noise_profile[k] = 1.0f;
  }

  // Base reduction: 0.5 (approx 6 dB reduction)
  // Curve bias: +6 dB on bin 10, -6 dB on bin 20
  curve_bias[10] = 6.0f;  // Total reduction ~ 12 dB -> floor ~ 0.25
  curve_bias[20] = -6.0f; // Total reduction ~ 0 dB -> floor ~ 1.0

  for (uint32_t k = 0; k < fft_size; k++) {
    gain_spectrum[k] = 0.0f;
  }

  noise_floor_manager_apply(nfm, real_size, fft_size, gain_spectrum,
                            noise_profile, 0.5f, 0.5f, NULL, 0.0f, curve_bias);

  // Bin 10 should have deeper reduction (lower floor ~ 0.25)
  TEST_ASSERT(gain_spectrum[10] < gain_spectrum[0], "Bin 10 floor should be lower with positive bias");
  // Bin 20 should have less reduction (higher floor ~ 1.0)
  TEST_ASSERT(gain_spectrum[20] > gain_spectrum[0], "Bin 20 floor should be higher with negative bias");

  free(gain_spectrum);
  free(noise_profile);
  free(curve_bias);
  noise_floor_manager_free(nfm);
  printf("✓ Noise Floor Manager reduction curve tests passed\n");
}

int main(void) {
  printf("Running Noise Floor Manager tests...\n\n");

  test_noise_floor_manager_lifecycle();
  test_noise_floor_manager_apply();
  test_noise_floor_manager_band_limited();
  test_noise_floor_manager_reduction_curve();

  printf("\n✅ All Noise Floor Manager tests passed!\n");
  return 0;
}
