/*
 * Unit tests for Masking Veto
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/utils/spectral_features.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_masking_veto_logic(void) {
  printf("Testing MaskingVeto logic...\n");

  uint32_t fft_size = 1024;
  MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE,
                                            POWER_SPECTRUM, false, true);
  TEST_ASSERT(mv != NULL, "Init failed");

  float smoothed[513];
  float noisy[513];
  float noise[513];
  float alpha[513];
  for (uint32_t i = 0; i < 513; i++) {
    smoothed[i] = 1.0F;
    noisy[i] = 1.0F;
    noise[i] = 0.1F;
    alpha[i] = 2.0F;
  }

  // 1. Normal run
  masking_veto_apply(mv, smoothed, noisy, noise, NULL, alpha, 1.0F, 1.0F, 0.2F);

  // 2. NMR branches: Low noise vs threshold
  // Set noise very low so NMR <= 0
  for (uint32_t i = 0; i < 513; i++) {
    noise[i] = 1e-15F;
  }
  masking_veto_apply(mv, smoothed, noisy, noise, NULL, alpha, 1.0F, 1.0F, 0.2F);

  // 3. Edge/Invalid Depth
  masking_veto_apply(mv, smoothed, noisy, noise, NULL, alpha, 1.0F, -1.0F,
                     0.2F);

  // 4. NULL guards
  masking_veto_apply(NULL, smoothed, noisy, noise, NULL, alpha, 1.0F, 1.0F,
                     0.2F);
  masking_veto_apply(mv, NULL, noisy, noise, NULL, alpha, 1.0F, 1.0F, 0.2F);

  masking_veto_free(mv);
  printf("✓ Masking Veto tests passed\n");
}

void test_dot_artifact_removal(void) {
  printf("Testing Veto Interpolation and Protection...\n");

  uint32_t fft_size = 1024; // Real size = 513
  MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE,
                                            POWER_SPECTRUM, false, true);
  TEST_ASSERT(mv != NULL, "Init failed");

  uint32_t real_size = (fft_size / 2) + 1;
  float* smoothed = (float*)calloc(real_size, sizeof(float));
  float* noisy = (float*)calloc(real_size, sizeof(float));
  float* noise = (float*)calloc(real_size, sizeof(float));
  float* alpha = (float*)calloc(real_size, sizeof(float));

  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 1000.0F; // High signal -> High threshold
    noisy[i] = 1000.0F;
    alpha[i] = 1.0F; // Initial alpha (high reduction)

    // Default noise: Very low (Masked)
    noise[i] = 1.0F;
  }

  // A wide spike ensures we hit a few Bark band centers in the OPUS scale
  for (uint32_t i = 5; i <= 15; i++) {
    noise[i] = 10000.0F;
  }

  // Run Veto
  // alpha input is 1.0.
  // floor_alpha is 0.1.
  // depth = 1.0 (Full protection).
  // elasticity = 0.0 (Rigid).
  // Run Veto
  // alpha input is 1.0.
  // floor_alpha is 0.1.
  // depth = 1.0 (Full protection).
  // elasticity = 0.0 (Rigid).
  masking_veto_apply(mv, smoothed, noisy, noise, NULL, alpha, 0.1F, 1.0F, 0.0F);

  // Verify Bin 10 is rescued (high alpha) and its neighbors are at least
  // partially rescued (> floor) due to interpolation.
  TEST_ASSERT(alpha[10] > 0.5F, "Bin 10 should be rescued");
  TEST_ASSERT(alpha[9] > 0.1F, "Bin 9 should be above floor");
  TEST_ASSERT(alpha[11] > 0.1F, "Bin 11 should be above floor");

  free(smoothed);
  free(noisy);
  free(noise);
  free(alpha);
  masking_veto_free(mv);

  printf("✓ Interpolation protection tests passed\n");
}

void test_temporal_masking(void) {
  printf("Testing Temporal Masking (Forward and Backward)...\n");

  uint32_t fft_size = 1024;
  MaskingVeto* mv = masking_veto_initialize(fft_size, 48000, OPUS_SCALE,
                                            POWER_SPECTRUM, false, true);
  TEST_ASSERT(mv != NULL, "Init failed");

  uint32_t real_size = 513;
  float* smoothed = (float*)calloc(real_size, sizeof(float));
  float* noisy = (float*)calloc(real_size, sizeof(float));
  float* future = (float*)calloc(real_size, sizeof(float));
  float* noise = (float*)calloc(real_size, sizeof(float));
  float* alpha = (float*)calloc(real_size, sizeof(float));

  // Initialize: No signal
  for (uint32_t i = 0; i < real_size; i++) {
    alpha[i] = 4.0F; // High reduction
    noise[i] = 0.1F; // Medium noise
  }

  // 1. Forward Masking: Impulse in previous frame
  // We simulate this by first running with an impulse
  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 1000.0F;
  }
  masking_veto_apply(mv, smoothed, smoothed, noise, NULL, alpha, 0.1F, 1.0F,
                     0.0F);

  // Now run with a quiet signal (15dB SNR).
  // Noise is 0.1, so 15dB SNR is approx 3.16.
  // 3.16 is much lower than the previous 1000.0, so it would normally
  // be reduced, but forward masking should rescue it.
  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 3.2F;
    alpha[i] = 4.0F;
  }
  masking_veto_apply(mv, smoothed, smoothed, noise, NULL, alpha, 0.1F, 1.0F,
                     0.0F);

  // Verify bin 100 is rescued due to forward masking
  TEST_ASSERT(alpha[100] < 4.0F,
              "Forward masking should rescue alpha for quiet signals");

  // 2. Backward Masking: Impulse in future frame
  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 3.2F;
    future[i] = 1000.0F; // Impulse in future
    alpha[i] = 4.0F;
  }
  masking_veto_apply(mv, smoothed, smoothed, noise, future, alpha, 0.1F, 1.0F,
                     0.0F);

  // Verify bin 100 is rescued due to backward masking
  TEST_ASSERT(alpha[100] < 4.0F,
              "Backward masking should rescue alpha for quiet signals");

  free(smoothed);
  free(noisy);
  free(future);
  free(noise);
  free(alpha);
  masking_veto_free(mv);

  printf("✓ Temporal masking tests passed\n");
}

int main(void) {
  printf("Running Masking Veto tests...\n");
  test_masking_veto_logic();
  test_dot_artifact_removal();
  test_temporal_masking();
  printf("✅ All Masking Veto tests passed!\n");
  return 0;
}
