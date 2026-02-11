/*
 * Unit tests for Suppression Engine
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/suppression_engine.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_suppression_engine_branches(void) {
  printf("Testing SuppressionEngine branches...\n");

  uint32_t size = 10;
  uint32_t sample_rate = 44100;
  SuppressionEngine* se = suppression_engine_initialize(
      size, sample_rate, OPUS_SCALE, POWER_SPECTRUM);
  TEST_ASSERT(se != NULL, "Init failed");

  float ref[10];
  float noise[10];
  float alpha[10];
  float beta[10];

  SuppressionParameters params = {.type = SUPPRESSION_BEROUTI_PER_BIN,
                                  .strength = 0.5F,
                                  .undersubtraction = 0.0F};

  // 1. Berouti Per Bin - Low SNR (<= -5dB) -> alpha should be alpha_max_user
  for (uint32_t i = 0; i < size; i++) {
    ref[i] = 0.01F;
    noise[i] = 1.0F;
  } // SNR ~ -20dB
  suppression_engine_calculate(se, ref, noise, params, alpha, beta);
  float alpha_max_05 = ALPHA_MIN + (0.5F * (ALPHA_MAX - ALPHA_MIN));
  for (uint32_t i = 0; i < size; i++) {
    TEST_FLOAT_CLOSE(alpha[i], alpha_max_05, 0.001F);
  }

  // 2. Berouti Per Bin - High SNR (>= 20dB) -> alpha should be ALPHA_MIN
  for (uint32_t i = 0; i < size; i++) {
    ref[i] = 1000.0F;
    noise[i] = 1.0F;
  } // SNR ~ 30dB
  suppression_engine_calculate(se, ref, noise, params, alpha, beta);
  for (uint32_t i = 0; i < size; i++) {
    TEST_FLOAT_CLOSE(alpha[i], ALPHA_MIN, 0.001F);
  }

  // 3. Global SNR
  params.type = SUPPRESSION_GLOBAL_SNR;
  params.strength = 4.0F;
  // Low SNR -> Alpha = strength
  for (uint32_t i = 0; i < size; i++) {
    ref[i] = 0.01F;
    noise[i] = 1.0F;
  }
  suppression_engine_calculate(se, ref, noise, params, alpha, beta);
  TEST_FLOAT_CLOSE(alpha[0], 4.0F, 0.001F);

  // 4. Critical Bands SNR
  params.type = SUPPRESSION_CRITICAL_BANDS_SNR;
  suppression_engine_calculate(se, ref, noise, params, alpha, beta);
  TEST_FLOAT_CLOSE(alpha[0], 4.0F, 0.001F);

  // 5. Masking Thresholds
  params.type = SUPPRESSION_MASKING_THRESHOLDS;
  // High SNR/Masked -> High NMR -> Alpha = strength
  for (uint32_t i = 0; i < size; i++) {
    ref[i] = 10.0F;
    noise[i] = 1e-12F;
  }
  suppression_engine_calculate(se, ref, noise, params, alpha, beta);
  // Conservative masked alpha should be close to ALPHA_MIN + gap*0.1 = 1.0 +
  // 3.0*0.1 = 1.3
  TEST_FLOAT_CLOSE(alpha[0], 1.3F, 0.001F);

  // 6. NULL guards
  suppression_engine_calculate(NULL, ref, noise, params, alpha, beta);
  suppression_engine_calculate(se, NULL, noise, params, alpha, beta);

  suppression_engine_free(se);
  printf("✓ Suppression Engine tests passed\n");
}

int main(void) {
  printf("Running Suppression Engine tests...\n");
  test_suppression_engine_branches();
  printf("✅ All Suppression Engine tests passed!\n");
  return 0;
}
