/*
 * Unit tests for Transient Detector
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/utils/transient_detector.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

void test_per_band_transient_detection(void) {
  printf("Testing per-band transient detection...\n");

  uint32_t num_bands = 25;
  TransientDetector* td = transient_detector_initialize(num_bands);
  TEST_ASSERT(td != NULL, "Init failed");

  float energies[25];
  float weights[25];

  // 1. Initial run: steady noise
  for (uint32_t i = 0; i < 25; i++) {
    energies[i] = 1.0F;
  }
  transient_detector_process(td, energies, weights);

  // Weights should be 0 because innovation is low (or we are initializing)
  // Actually on first run it might be 1.0 if smoothed is 0
  // Let's run twice to stabilize
  transient_detector_process(td, energies, weights);
  for (uint32_t i = 0; i < 25; i++) {
    TEST_ASSERT(weights[i] < 0.1F, "Weights should be low for steady signal");
  }

  // 2. Sudden burst on band 10
  energies[10] = 100.0F;
  transient_detector_process(td, energies, weights);

  TEST_ASSERT(weights[10] > 0.9F, "Band 10 should be detected as transient");
  TEST_ASSERT(weights[0] < 0.1F, "Band 0 should NOT be detected as transient");

  // 3. Steady state after burst (allow it to settle over a few frames)
  for (int i = 0; i < 5; i++) {
    transient_detector_process(td, energies, weights);
  }
  // Reference should have increased, so weight should drop
  TEST_ASSERT(weights[10] < 0.5F,
              "Weight should drop after signal becomes steady");

  transient_detector_free(td);
  printf("✓ Per-band transient detection passed\n");
}

void test_global_transient_detection(void) {
  printf("Testing global transient detection...\n");

  uint32_t size = 1024;
  TransientDetector* td = transient_detector_initialize(size);
  TEST_ASSERT(td != NULL, "Init failed");

  float spectrum[1024];

  // 1. Steady
  for (uint32_t i = 0; i < 1024; i++) {
    spectrum[i] = 1.0F;
  }
  transient_detector_run(td, spectrum);
  transient_detector_run(td, spectrum);

  TEST_ASSERT(transient_detector_run(td, spectrum) == false,
              "Should not detect transient in steady signal");

  // 2. Burst
  for (uint32_t i = 0; i < 1024; i++) {
    spectrum[i] = 100.0F;
  }
  TEST_ASSERT(transient_detector_run(td, spectrum) == true,
              "Should detect global transient");

  transient_detector_free(td);
  printf("✓ Global transient detection passed\n");
}

int main(void) {
  printf("Running Transient Detector tests...\n");
  test_per_band_transient_detection();
  test_global_transient_detection();
  printf("✅ All Transient Detector tests passed!\n");
  return 0;
}
