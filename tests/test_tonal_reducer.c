/*
libspecbleach - A spectral processing library
Test suite for the Tonal Reducer module.
*/

#include "shared/denoiser_logic/processing/tonal_reducer.h"
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
  float cv_mask[TEST_SPECTRUM_SIZE];

  // Initialize with flat noise
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    alpha[i] = 1.0f;           // Default alpha (no oversubtraction yet)
    noise_spectrum[i] = 0.01f; // Flat noise
    cv_mask[i] = 0.0f;
  }

  // Run with no reduction requested (gain 1.0)
  tonal_reducer_run(reducer, noise_spectrum, cv_mask, alpha, 1.0f);

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
  float cv_mask[TEST_SPECTRUM_SIZE];

  // Initialize flat with realistic alpha=1.0 (ALPHA_MIN)
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    alpha[i] = 1.0f;
    noise_spectrum[i] = 0.01f;
    cv_mask[i] = 0.0f;
  }

  // Add a tone at 1 kHz in CV mask
  int bin = freq_to_bin(1000.0f);
  noise_spectrum[bin] = 0.1f;
  cv_mask[bin] = 1.0f;
  cv_mask[bin - 1] = 0.5f;
  cv_mask[bin + 1] = 0.5f;

  float reduction_gain = 0.00398f; // ~48dB reduction (gain close to 0)

  tonal_reducer_run(reducer, noise_spectrum, cv_mask, alpha, reduction_gain);

  const float* mask = tonal_reducer_get_mask(reducer);

  if (mask[bin] <= 0.0f) {
    fprintf(stderr, "FAIL: Tone not detected at bin %d\n", bin);
    exit(1);
  }
  printf("  Tone detected at bin %d (mask=%.3f)\n", bin, mask[bin]);

  // Alpha should be boosted aggressively.
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

void test_caching_and_adaptive_support(void) {
  printf("Testing mask caching and adaptive support in tonal reducer...\n");
  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);

  float alpha[TEST_SPECTRUM_SIZE];
  float noise_spectrum[TEST_SPECTRUM_SIZE];

  // 1. Initialize empty CV mask (adaptive mode)
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    alpha[i] = 1.0f;
    noise_spectrum[i] = 0.01f;
  }

  // Add a tone at bin 100
  int bin1 = 100;
  noise_spectrum[bin1] = 0.1f;
  noise_spectrum[bin1 - 1] = 0.03f;
  noise_spectrum[bin1 + 1] = 0.03f;

  float reduction_gain = 0.0f; // max reduction strength

  // Run 1: Should detect the tone at bin1 in adaptive mode
  tonal_reducer_run(reducer, noise_spectrum, NULL, alpha, reduction_gain);

  const float* mask = tonal_reducer_get_mask(reducer);
  if (mask[bin1] <= 0.0f) {
    fprintf(stderr,
            "FAIL: Reducer in adaptive mode did not detect tone at bin %d\n",
            bin1);
    exit(1);
  }
  printf("  Run 1: Tone at bin %d detected in adaptive mode (mask=%.3f) ✓\n",
         bin1, mask[bin1]);

  // Reset alpha
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    alpha[i] = 1.0f;
  }

  // Move tone to bin 200 but keep the sum identical
  int bin2 = 200;
  // Restore bin1 to noise floor
  noise_spectrum[bin1] = 0.01f;
  noise_spectrum[bin1 - 1] = 0.01f;
  noise_spectrum[bin1 + 1] = 0.01f;
  // Put tone at bin2
  noise_spectrum[bin2] = 0.1f;
  noise_spectrum[bin2 - 1] = 0.03f;
  noise_spectrum[bin2 + 1] = 0.03f;

  // Run 2: Since sum is identical, caching should skip detection.
  // The cached mask should still have tone at bin1 and NOT at bin2.
  tonal_reducer_run(reducer, noise_spectrum, NULL, alpha, reduction_gain);
  mask = tonal_reducer_get_mask(reducer);

  if (mask[bin1] <= 0.0f) {
    fprintf(stderr, "FAIL: Caching failed, cached mask at bin %d cleared\n",
            bin1);
    exit(1);
  }
  if (mask[bin2] > 0.0f) {
    fprintf(stderr,
            "FAIL: Redetected tone at bin %d despite identical sum (caching "
            "failed)\n",
            bin2);
    exit(1);
  }
  printf(
      "  Run 2: Tonal mask successfully cached (old tone remains, new tone "
      "skipped) ✓\n");

  // Reset alpha
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    alpha[i] = 1.0f;
  }

  // Run 3: Modify sum slightly (change bin 5 to 0.02)
  noise_spectrum[5] = 0.02f;

  // Should trigger recalculation. Old tone at bin1 should be cleared, new tone
  // at bin2 detected.
  tonal_reducer_run(reducer, noise_spectrum, NULL, alpha, reduction_gain);
  mask = tonal_reducer_get_mask(reducer);

  if (mask[bin1] > 0.0f) {
    fprintf(stderr,
            "FAIL: Cache did not invalidate after sum changed; old tone at bin "
            "%d still present\n",
            bin1);
    exit(1);
  }
  if (mask[bin2] <= 0.0f) {
    fprintf(stderr,
            "FAIL: Cache invalidation did not detect new tone at bin %d\n",
            bin2);
    exit(1);
  }
  printf(
      "  Run 3: Cache successfully invalidated and updated (new tone detected, "
      "old tone cleared) ✓\n");

  // Run 4: Manual CV mask mode
  float manual_cv_mask[TEST_SPECTRUM_SIZE];
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    manual_cv_mask[i] = 0.0f;
    alpha[i] = 1.0f;
  }
  int cv_bin = 50;
  manual_cv_mask[cv_bin] = 1.0f;

  tonal_reducer_run(reducer, noise_spectrum, manual_cv_mask, alpha,
                    reduction_gain);
  mask = tonal_reducer_get_mask(reducer);

  if (mask[cv_bin] <= 0.0f) {
    fprintf(stderr,
            "FAIL: Reducer in manual CV mask mode did not set mask at bin %d\n",
            cv_bin);
    exit(1);
  }
  printf("  Run 4: CV mask directly applied at bin %d ✓\n", cv_bin);

  tonal_reducer_free(reducer);
  printf("✓ Cache and adaptive support tests passed\n");
}

void test_tonal_reducer_peaks(void) {
  printf("Testing tonal_reducer_get_peaks and edge cases...\n");
  float freqs[10];
  if (tonal_reducer_get_peaks(NULL, freqs, 10) != 0 ||
      tonal_reducer_get_mask(NULL) != NULL) {
    fprintf(stderr, "FAIL: Null reducer handling failed\n");
    exit(1);
  }

  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);

  if (tonal_reducer_get_peaks(reducer, freqs, 10) != 0 ||
      tonal_reducer_get_peaks(reducer, freqs, 0) != 0) {
    fprintf(stderr,
            "FAIL: Uninitialized mask or max_peaks=0 should return 0 peaks\n");
    exit(1);
  }

  float alpha[TEST_SPECTRUM_SIZE];
  float noise[TEST_SPECTRUM_SIZE] = {0.01f};
  float cv_mask[TEST_SPECTRUM_SIZE] = {0.0f};
  // Gain 1.0f (no reduction requested, short-circuits early)
  tonal_reducer_run(reducer, noise, cv_mask, alpha, 1.0f);

  tonal_reducer_free(reducer);
  printf("✓ Tonal reducer peaks test passed\n");
}

int main(void) {
  printf("=== Tonal Reducer Tests ===\n\n");
  test_initialization();
  test_flat_noise_no_boost();
  test_tonal_boost();
  test_caching_and_adaptive_support();
  test_tonal_reducer_peaks();
  printf("\n=== All tonal reducer tests passed ===\n");
  return 0;
}
