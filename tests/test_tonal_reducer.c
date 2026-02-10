/*
libspecbleach - A spectral processing library
Test suite for the Tonal Reducer module.
*/

#include "shared/noise_estimation/tonal_reducer.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_SAMPLE_RATE 48000
#define TEST_FFT_SIZE 2048
#define TEST_SPECTRUM_SIZE ((TEST_FFT_SIZE / 2) + 1)
#define BIN_HZ ((float)TEST_SAMPLE_RATE / (float)TEST_FFT_SIZE)

// Helper: find the bin closest to a target frequency
static int freq_to_bin(float freq_hz) {
  return (int)roundf(freq_hz / BIN_HZ);
}

void test_initialization(void) {
  printf("Testing initialization...\n");
  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);
  if (!reducer) {
    fprintf(stderr, "FAIL: Failed to initialize tonal_reducer\n");
    exit(1);
  }

  const float* mask = tonal_reducer_get_mask(reducer);
  if (!mask) {
    fprintf(stderr, "FAIL: get_mask returned NULL\n");
    exit(1);
  }

  tonal_reducer_free(reducer);
  printf("✓ Initialization passed\n");
}

void test_flat_noise_no_boost(void) {
  printf("Testing flat noise (no boost expected)...\n");
  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);

  float alpha[TEST_SPECTRUM_SIZE];
  float noise_spectrum[TEST_SPECTRUM_SIZE];
  float max_profile[TEST_SPECTRUM_SIZE];
  float median_profile[TEST_SPECTRUM_SIZE];

  // Initialize with flat noise
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    alpha[i] = 1.0f;           // Default alpha (no oversubtraction yet)
    noise_spectrum[i] = 0.01f; // Flat noise
    max_profile[i] = 0.01f;
    median_profile[i] = 0.01f;
  }

  // Run with heavy reduction requested (-60dB -> alpha target ~1.0)
  // But wait, alpha is oversubtraction factor. Usually alpha > 1.0 means more
  // reduction. Wait, let's check tonal_reducer.c logic. "alpha[k] =
  // fmaxf(alpha[k], mask[k] * target_alpha)" If reduction is 0.0 (0dB),
  // target_alpha = 1.0. If reduction is 0.9 (20dB), target_alpha = ? Let's
  // assume standard Berouti alpha usage or similar. The code says: float
  // target_alpha = alpha_needed * tonal_mask[k]; Actually, let's trust the
  // module logic: it boosts alpha. If no tonal components, mask should be 0. So
  // alpha should remain 1.0.

  // Run with no reduction requested (gain 1.0)
  tonal_reducer_run(reducer, noise_spectrum, max_profile, median_profile, alpha,
                    1.0f);

  const float* mask = tonal_reducer_get_mask(reducer);
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    if (mask[i] > 0.0f) {
      fprintf(stderr, "FAIL: False detection at bin %d\n", i);
      exit(1);
    }
    if (alpha[i] != 1.0f) {
      fprintf(stderr, "FAIL: Alpha modified at bin %d (expected 1.0, got %f)\n",
              i, alpha[i]);
      exit(1);
    }
  }

  tonal_reducer_free(reducer);
  printf("✓ Flat noise passed (no changes)\n");
}

void test_tonal_boost(void) {
  printf("Testing tonal boost (alpha increase)...\n");
  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);

  float alpha[TEST_SPECTRUM_SIZE];
  float noise_spectrum[TEST_SPECTRUM_SIZE];
  float max_profile[TEST_SPECTRUM_SIZE];
  float median_profile[TEST_SPECTRUM_SIZE];

  // Initialize flat with realistic alpha=1.0 (ALPHA_MIN)
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    alpha[i] = 1.0f;
    noise_spectrum[i] = 0.01f;
    max_profile[i] = 0.01f;
    median_profile[i] = 0.01f;
  }

  // Add a tone at 1 kHz (stable -> max ~ median)
  int bin = freq_to_bin(1000.0f);
  noise_spectrum[bin] = 0.1f;
  max_profile[bin] = 0.1f;
  median_profile[bin] = 0.1f;
  // Neighbors to simulate peak shape for parabolic interp
  max_profile[bin - 1] = 0.03f;
  max_profile[bin + 1] = 0.03f;
  median_profile[bin - 1] = 0.03f;
  median_profile[bin + 1] = 0.03f;

  float reduction_gain = 0.00398f; // ~48dB reduction (gain close to 0)

  tonal_reducer_run(reducer, noise_spectrum, max_profile, median_profile, alpha,
                    reduction_gain);

  const float* mask = tonal_reducer_get_mask(reducer);

  if (mask[bin] <= 0.0f) {
    fprintf(stderr, "FAIL: Tone not detected at bin %d\n", bin);
    exit(1);
  }
  printf("  Tone detected at bin %d (mask=%.3f)\n", bin, mask[bin]);

  // Alpha should be boosted.
  // With ALPHA_MAX_TONAL=10.0, we expect it to be close to 10.0.
  if (alpha[bin] <= 9.0f) {
    fprintf(stderr,
            "FAIL: Alpha not boosted aggressively (<=9.0) at tonal bin %d (got "
            "%f)\n",
            bin, alpha[bin]);
    exit(1);
  }
  printf("  Alpha at bin %d: %.3f (boosted aggressively) ✓\n", bin, alpha[bin]);

  tonal_reducer_free(reducer);
  printf("✓ Tonal boost passed\n");
}

int main(void) {
  printf("=== Tonal Reducer Tests ===\n\n");
  test_initialization();
  test_flat_noise_no_boost();
  test_tonal_boost();
  printf("\n=== All tonal reducer tests passed ===\n");
  return 0;
}
