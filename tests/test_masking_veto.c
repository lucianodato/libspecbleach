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
  float noise[513];
  float alpha[513];
  for (uint32_t i = 0; i < 513; i++) {
    smoothed[i] = 1.0F;
    noise[i] = 0.1F;
    alpha[i] = 2.0F;
  }

  // 1. Normal run
  masking_veto_apply(mv, smoothed, noise, NULL, alpha, 1.0F);

  // 2. NMR branches: Low noise vs threshold
  // Set noise very low so NMR <= 0
  for (uint32_t i = 0; i < 513; i++) {
    noise[i] = 1e-15F;
  }
  masking_veto_apply(mv, smoothed, noise, NULL, alpha, 1.0F);

  // 3. Edge/Invalid Depth
  masking_veto_apply(mv, smoothed, noise, NULL, alpha, -1.0F);

  // 4. NULL guards
  masking_veto_apply(NULL, smoothed, noise, NULL, alpha, 1.0F);
  masking_veto_apply(mv, NULL, noise, NULL, alpha, 1.0F);

  masking_veto_free(mv);
  printf("✓ Masking Veto tests passed\n");
}

void test_nmr_protection(void) {
  printf("Testing NMR-based Protection...\n");

  uint32_t fft_size = 1024; // Real size = 513
  MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE,
                                            POWER_SPECTRUM, false, true);
  TEST_ASSERT(mv != NULL, "Init failed");

  uint32_t real_size = (fft_size / 2) + 1;
  float* smoothed = (float*)calloc(real_size, sizeof(float));
  float* noise = (float*)calloc(real_size, sizeof(float));
  float* alpha = (float*)calloc(real_size, sizeof(float));

  // Scenario: Strong signal (1000.0) with low noise (1.0)
  // The masking threshold should be high -> noise (1.0) is masked
  // -> NMR < 0dB -> full protection -> alpha lerped toward 1.0
  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 1000.0F; // Strong signal
    noise[i] = 1.0F;       // Low noise (masked by signal)
    alpha[i] = 4.0F;       // High oversubtraction
  }

  masking_veto_apply(mv, smoothed, noise, NULL, alpha, 1.0F);

  // With noise masked, alpha should be lerped toward 1.0
  // (less than original 4.0)
  TEST_ASSERT(alpha[100] < 4.0F, "Masked noise should reduce alpha toward 1.0");

  free(smoothed);
  free(noise);
  free(alpha);
  masking_veto_free(mv);

  // Scenario 2: Noise-only (noise >> signal -> NMR high -> no protection)
  // Use a fresh instance to avoid state from previous scenario
  mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE, POWER_SPECTRUM,
                               false, true);
  TEST_ASSERT(mv != NULL, "Init failed for noise-only scenario");

  smoothed = (float*)calloc(real_size, sizeof(float));
  noise = (float*)calloc(real_size, sizeof(float));
  alpha = (float*)calloc(real_size, sizeof(float));

  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 0.01F; // Near silence
    noise[i] = 100.0F;   // Very loud noise
    alpha[i] = 4.0F;     // High oversubtraction
  }

  // Run multiple frames to let the clean signal estimation stabilize at ~0
  for (int frame = 0; frame < 10; frame++) {
    for (uint32_t i = 0; i < real_size; i++) {
      alpha[i] = 4.0F;
    }
    masking_veto_apply(mv, smoothed, noise, NULL, alpha, 1.0F);
  }

  // With noise clearly audible (NMR >> 0dB), alpha should stay near 4.0
  TEST_ASSERT(alpha[100] > 3.0F,
              "Audible noise should leave alpha mostly unchanged");

  free(smoothed);
  free(noise);
  free(alpha);
  masking_veto_free(mv);

  printf("✓ NMR protection tests passed\n");
}

void test_temporal_masking(void) {
  printf("Testing Temporal Masking (Forward and Backward)...\n");

  uint32_t fft_size = 1024;
  MaskingVeto* mv = masking_veto_initialize(fft_size, 48000, OPUS_SCALE,
                                            POWER_SPECTRUM, false, true);
  TEST_ASSERT(mv != NULL, "Init failed");

  uint32_t real_size = 513;
  float* smoothed = (float*)calloc(real_size, sizeof(float));
  float* future = (float*)calloc(real_size, sizeof(float));
  float* noise = (float*)calloc(real_size, sizeof(float));
  float* alpha = (float*)calloc(real_size, sizeof(float));

  // Initialize
  for (uint32_t i = 0; i < real_size; i++) {
    alpha[i] = 4.0F; // High reduction
    noise[i] = 0.1F; // Medium noise
  }

  // 1. Forward Masking: Impulse in previous frame
  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 1000.0F;
  }
  masking_veto_apply(mv, smoothed, noise, NULL, alpha, 1.0F);

  // Now run with a quiet signal that follows the impulse.
  // Forward masking should still provide some protection.
  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 3.2F;
    alpha[i] = 4.0F;
  }
  masking_veto_apply(mv, smoothed, noise, NULL, alpha, 1.0F);

  // Forward masking should reduce alpha (provide some protection)
  TEST_ASSERT(alpha[100] < 4.0F,
              "Forward masking should provide some protection");

  // 2. Backward Masking: Impulse in future frame
  for (uint32_t i = 0; i < real_size; i++) {
    smoothed[i] = 3.2F;
    future[i] = 1000.0F;
    alpha[i] = 4.0F;
  }
  masking_veto_apply(mv, smoothed, noise, future, alpha, 1.0F);

  // Backward masking should also provide protection
  TEST_ASSERT(alpha[100] < 4.0F,
              "Backward masking should provide some protection");

  free(smoothed);
  free(future);
  free(noise);
  free(alpha);
  masking_veto_free(mv);

  printf("✓ Temporal masking tests passed\n");
}

int main(void) {
  printf("Running Masking Veto tests...\n");
  test_masking_veto_logic();
  test_nmr_protection();
  test_temporal_masking();
  printf("✅ All Masking Veto tests passed!\n");
  return 0;
}
