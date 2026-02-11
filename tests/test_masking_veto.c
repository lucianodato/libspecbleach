/*
 * Unit tests for Masking Veto
 */

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
  MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, POWER_SPECTRUM);
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
  masking_veto_apply(mv, smoothed, noisy, noise, alpha, 1.0F, 1.0F, 0.2F);

  // 2. NMR branches: Low noise vs threshold
  // Set noise very low so NMR <= 0
  for (uint32_t i = 0; i < 513; i++) {
    noise[i] = 1e-15F;
  }
  masking_veto_apply(mv, smoothed, noisy, noise, alpha, 1.0F, 1.0F, 0.2F);

  // 3. Edge/Invalid Depth
  masking_veto_apply(mv, smoothed, noisy, noise, alpha, 1.0F, -1.0F, 0.2F);

  // 4. NULL guards
  masking_veto_apply(NULL, smoothed, noisy, noise, alpha, 1.0F, 1.0F, 0.2F);
  masking_veto_apply(mv, NULL, noisy, noise, alpha, 1.0F, 1.0F, 0.2F);

  masking_veto_free(mv);
  printf("✓ Masking Veto tests passed\n");
}

int main(void) {
  printf("Running Masking Veto tests...\n");
  test_masking_veto_logic();
  printf("✅ All Masking Veto tests passed!\n");
  return 0;
}
