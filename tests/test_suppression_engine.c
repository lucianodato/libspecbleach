/*
 * Unit tests for Suppression Engine
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/configurations.h"
#include "../src/shared/gain_estimation/suppression_engine.h"

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
  SuppressionEngine* se = suppression_engine_initialize(size);
  TEST_ASSERT(se != NULL, "Init failed");

  float ref[10];
  float noise[10];
  float alpha[10];
  float beta[10];

  // 1. Low SNR (<= -5dB) -> alpha should be alpha_max_user
  for (uint32_t i = 0; i < size; i++) {
    ref[i] = 0.01F;
    noise[i] = 1.0F;
  } // SNR ~ -20dB
  suppression_engine_calculate(se, ref, noise, 0.5F, alpha, beta);
  float alpha_max_05 = ALPHA_MIN + (0.5F * (ALPHA_MAX - ALPHA_MIN));
  for (uint32_t i = 0; i < size; i++) {
    TEST_FLOAT_CLOSE(alpha[i], alpha_max_05, 0.001F);
  }

  // 2. High SNR (>= 20dB) -> alpha should be ALPHA_MIN
  for (uint32_t i = 0; i < size; i++) {
    ref[i] = 1000.0F;
    noise[i] = 1.0F;
  } // SNR ~ 30dB
  suppression_engine_calculate(se, ref, noise, 0.5F, alpha, beta);
  for (uint32_t i = 0; i < size; i++) {
    TEST_FLOAT_CLOSE(alpha[i], ALPHA_MIN, 0.001F);
  }

  // 3. NULL guards
  suppression_engine_calculate(NULL, ref, noise, 0.5F, alpha, beta);
  suppression_engine_calculate(se, NULL, noise, 0.5F, alpha, beta);

  suppression_engine_free(se);
  printf("✓ Suppression Engine tests passed\n");
}

int main(void) {
  printf("Running Suppression Engine tests...\n");
  test_suppression_engine_branches();
  printf("✅ All Suppression Engine tests passed!\n");
  return 0;
}
