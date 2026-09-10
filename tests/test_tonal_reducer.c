/*
libspecbleach - A spectral processing library

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
libspecbleach - A spectral processing library
Test suite for the Tonal Reducer module (dual-path split + parallel tonal
gain path).
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

void test_flat_noise_identity_split(void) {
  printf("Testing flat noise (identity split, no tonal path)...\n");
  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);

  float noise_spectrum[TEST_SPECTRUM_SIZE];
  float noise_bb[TEST_SPECTRUM_SIZE];
  float noise_tonal[TEST_SPECTRUM_SIZE];
  float gain_tonal[TEST_SPECTRUM_SIZE];
  float smoothed[TEST_SPECTRUM_SIZE];

  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    noise_spectrum[i] = 0.01f; // Flat noise
    smoothed[i] = 0.2f;        // Somewhat above the noise
  }

  // Max tonal reduction requested, adaptive mode (no CV profile)
  tonal_reducer_compute_split(reducer, noise_spectrum, NULL, false, 0.0f,
                              noise_bb, noise_tonal);

  const float* mask = tonal_reducer_get_mask(reducer);
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    if (mask[i] > 0.0f) {
      fprintf(stderr, "FAIL: False tonal detection at bin %d\n", i);
      exit(1);
    }
    if (noise_bb[i] != noise_spectrum[i]) {
      fprintf(stderr, "FAIL: bb profile modified at bin %d (%f vs %f)\n", i,
              noise_bb[i], noise_spectrum[i]);
      exit(1);
    }
    if (noise_tonal[i] != 0.0f) {
      fprintf(stderr, "FAIL: Non-zero tonal residual at bin %d\n", i);
      exit(1);
    }
  }

  // Tonal gains must be unity everywhere (no residual)
  tonal_reducer_compute_tonal_gains(reducer, 0U, smoothed, noise_tonal, mask,
                                    0.0f, gain_tonal);
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    if (gain_tonal[i] != 1.0f) {
      fprintf(stderr, "FAIL: Tonal gain not unity at bin %d\n", i);
      exit(1);
    }
  }

  tonal_reducer_free(reducer);
  printf("✓ Flat noise passed (identity split)\n");
}

void test_tonal_split_and_gains(void) {
  printf("Testing tonal split and parallel gains (CV mask)...\n");
  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);

  float noise_spectrum[TEST_SPECTRUM_SIZE];
  float cv_mask[TEST_SPECTRUM_SIZE];
  float noise_bb[TEST_SPECTRUM_SIZE];
  float noise_tonal[TEST_SPECTRUM_SIZE];
  float gain_tonal[TEST_SPECTRUM_SIZE];
  float smoothed[TEST_SPECTRUM_SIZE];

  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    noise_spectrum[i] = 0.01f; // Flat noise floor
    smoothed[i] = 0.2f;        // Signal somewhat above the noise
    cv_mask[i] = 0.0f;
  }

  // Add a tone at 1 kHz in CV mask (with energy shoulders on the masked
  // side bins so the partial-mask softening is actually exercised)
  int bin = freq_to_bin(1000.0f);
  noise_spectrum[bin] = 0.1f;
  noise_spectrum[bin - 1] = 0.05f;
  noise_spectrum[bin + 1] = 0.05f;
  cv_mask[bin] = 1.0f;
  cv_mask[bin - 1] = 0.5f;
  cv_mask[bin + 1] = 0.5f;

  float reduction_gain = 0.0f; // max reduction strength

  tonal_reducer_compute_split(reducer, noise_spectrum, cv_mask, true,
                              reduction_gain, noise_bb, noise_tonal);

  const float* mask = tonal_reducer_get_mask(reducer);
  if (mask[bin] <= 0.0f) {
    fprintf(stderr, "FAIL: Tone not detected at bin %d\n", bin);
    exit(1);
  }
  printf("  Tone detected at bin %d (mask=%.3f)\n", bin, mask[bin]);

  // At the mask bin the residual must be positive and the bb profile must be
  // pushed toward the envelope (below the raw profile).
  if (noise_tonal[bin] <= 0.0f) {
    fprintf(stderr, "FAIL: No tonal residual at bin %d (%f)\n", bin,
            noise_tonal[bin]);
    exit(1);
  }
  if (noise_bb[bin] >= noise_spectrum[bin]) {
    fprintf(stderr,
            "FAIL: bb profile not flattened at tonal bin %d (%f >= %f)\n", bin,
            noise_bb[bin], noise_spectrum[bin]);
    exit(1);
  }
  // A mask-free bin keeps the profile exactly and carries no residual.
  if (noise_bb[bin + 10] != noise_spectrum[bin + 10] ||
      noise_tonal[bin + 10] != 0.0f) {
    fprintf(stderr, "FAIL: Identity split violated at mask-free bin %d\n",
            bin + 10);
    exit(1);
  }

  // Tonal gains: deep notch at the full-strength bin, softer at partial-mask
  // neighbors, unity far away.
  tonal_reducer_compute_tonal_gains(reducer, 0U, smoothed, noise_tonal, mask,
                                    reduction_gain, gain_tonal);
  if (gain_tonal[bin] > 0.5f) {
    fprintf(stderr, "FAIL: Tonal gain not notched at bin %d (%f)\n", bin,
            gain_tonal[bin]);
    exit(1);
  }
  if (gain_tonal[bin + 10] != 1.0f) {
    fprintf(stderr, "FAIL: Tonal gain not unity at non-tonal bin %d (%f)\n",
            bin + 10, gain_tonal[bin + 10]);
    exit(1);
  }
  if (gain_tonal[bin - 1] >= 1.0f || gain_tonal[bin - 1] <= gain_tonal[bin]) {
    fprintf(stderr,
            "FAIL: Partial-mask bin %d not softer than full bin %d (got %f vs "
            "%f)\n",
            bin - 1, bin, gain_tonal[bin - 1], gain_tonal[bin]);
    exit(1);
  }
  printf("  Gains: full=%f partial=%f far=1.0 ✓\n", gain_tonal[bin],
         gain_tonal[bin - 1]);

  // Second call: feed a changed magnitude so the raw target moves up. The
  // one-pole must land strictly between the seeded gain and the new raw target
  // (a fresh reducer's first call exposes that unstabilized target).
  const int probe = bin - 1; // partial-mask bin: unsaturated and movable
  float smoothed_second[TEST_SPECTRUM_SIZE];
  memcpy(smoothed_second, smoothed, sizeof(smoothed));
  smoothed_second[probe] = 0.5f; // more signal -> higher raw gain
  float gain_second[TEST_SPECTRUM_SIZE];
  tonal_reducer_compute_tonal_gains(reducer, 0U, smoothed_second, noise_tonal,
                                    mask, reduction_gain, gain_second);
  if (gain_second[bin] < gain_tonal[bin] - 1e-6f || gain_second[bin] > 1.0f) {
    fprintf(stderr, "FAIL: One-pole trajectory invalid at bin %d (%f -> %f)\n",
            bin, gain_tonal[bin], gain_second[bin]);
    exit(1);
  }

  TonalReducer* fresh = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);
  float fresh_bb[TEST_SPECTRUM_SIZE];
  float fresh_tonal[TEST_SPECTRUM_SIZE];
  tonal_reducer_compute_split(fresh, noise_spectrum, cv_mask, true,
                              reduction_gain, fresh_bb, fresh_tonal);
  float target[TEST_SPECTRUM_SIZE];
  tonal_reducer_compute_tonal_gains(fresh, 0U, smoothed_second, fresh_tonal,
                                    tonal_reducer_get_mask(fresh),
                                    reduction_gain, target);
  tonal_reducer_free(fresh);
  if (!(gain_second[probe] > gain_tonal[probe] + 1e-6f &&
        gain_second[probe] < target[probe] - 1e-6f)) {
    fprintf(stderr,
            "FAIL: One-pole not between seeded and raw target at bin %d "
            "(%f -> %f, target %f)\n",
            probe, gain_tonal[probe], gain_second[probe], target[probe]);
    exit(1);
  }

  tonal_reducer_free(reducer);
  printf("✓ Tonal split and gains passed\n");
}

void test_disabled_reduction_is_legacy(void) {
  printf("Testing disabled tonal reduction (legacy single path)...\n");
  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);

  float noise_spectrum[TEST_SPECTRUM_SIZE];
  float cv_mask[TEST_SPECTRUM_SIZE];
  float noise_bb[TEST_SPECTRUM_SIZE];
  float noise_tonal[TEST_SPECTRUM_SIZE];

  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    noise_spectrum[i] = 0.01f;
    cv_mask[i] = 0.0f;
  }
  int bin = freq_to_bin(1000.0f);
  noise_spectrum[bin] = 0.1f;
  cv_mask[bin] = 1.0f;

  // Reduction disabled (gain 1.0): mask still published but the split must
  // be a pure copy (legacy behavior).
  tonal_reducer_compute_split(reducer, noise_spectrum, cv_mask, true, 1.0f,
                              noise_bb, noise_tonal);
  const float* mask = tonal_reducer_get_mask(reducer);
  if (mask[bin] <= 0.0f) {
    fprintf(stderr, "FAIL: Mask not published when reduction disabled\n");
    exit(1);
  }
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    if (noise_bb[i] != noise_spectrum[i]) {
      fprintf(stderr,
              "FAIL: Disabled reduction modified bb profile at bin %d\n", i);
      exit(1);
    }
    if (noise_tonal[i] != 0.0f) {
      fprintf(stderr, "FAIL: Disabled reduction produced residual at bin %d\n",
              i);
      exit(1);
    }
  }

  tonal_reducer_free(reducer);
  printf("✓ Disabled reduction passed (legacy path)\n");
}

void test_mask_refresh_and_adaptive_support(void) {
  printf("Testing adaptive mask refresh and manual mask support...\n");
  TonalReducer* reducer = tonal_reducer_initialize(
      TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE);

  float noise_spectrum[TEST_SPECTRUM_SIZE];
  float noise_bb[TEST_SPECTRUM_SIZE];
  float noise_tonal[TEST_SPECTRUM_SIZE];

  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    noise_spectrum[i] = 0.01f;
  }

  // Add a tone at bin 100
  int bin1 = 100;
  noise_spectrum[bin1] = 0.1f;
  noise_spectrum[bin1 - 1] = 0.03f;
  noise_spectrum[bin1 + 1] = 0.03f;

  // Run 1: Should detect the tone at bin1 in adaptive mode
  tonal_reducer_compute_split(reducer, noise_spectrum, NULL, false, 0.0f,
                              noise_bb, noise_tonal);
  const float* mask = tonal_reducer_get_mask(reducer);
  if (mask[bin1] <= 0.0f) {
    fprintf(stderr,
            "FAIL: Reducer in adaptive mode did not detect tone at bin %d\n",
            bin1);
    exit(1);
  }
  if (noise_tonal[bin1] <= 0.0f) {
    fprintf(stderr, "FAIL: No residual extracted for adaptive tone at bin %d\n",
            bin1);
    exit(1);
  }
  printf("  Run 1: Tone at bin %d detected in adaptive mode ✓\n", bin1);

  // Move tone to bin 200 but keep the structure similar
  int bin2 = 200;
  noise_spectrum[bin1] = 0.01f;
  noise_spectrum[bin1 - 1] = 0.01f;
  noise_spectrum[bin1 + 1] = 0.01f;
  noise_spectrum[bin2] = 0.1f;
  noise_spectrum[bin2 - 1] = 0.03f;
  noise_spectrum[bin2 + 1] = 0.03f;

  tonal_reducer_compute_split(reducer, noise_spectrum, NULL, false, 0.0f,
                              noise_bb, noise_tonal);
  mask = tonal_reducer_get_mask(reducer);
  if (mask[bin1] > 0.0f) {
    fprintf(stderr,
            "FAIL: Old tone at bin %d was not cleared on spectrum update\n",
            bin1);
    exit(1);
  }
  if (mask[bin2] <= 0.0f) {
    fprintf(stderr,
            "FAIL: Moved tone at bin %d was not detected on spectrum update\n",
            bin2);
    exit(1);
  }
  printf("  Run 2: Tonal mask refreshed on spectrum update ✓\n", bin2);

  // Add a third tone at bin 50
  noise_spectrum[50] = 0.1f;
  noise_spectrum[49] = 0.03f;
  noise_spectrum[51] = 0.03f;

  tonal_reducer_compute_split(reducer, noise_spectrum, NULL, false, 0.0f,
                              noise_bb, noise_tonal);
  mask = tonal_reducer_get_mask(reducer);
  if (mask[50] <= 0.0f) {
    fprintf(stderr, "FAIL: New tone at bin 50 was not detected\n");
    exit(1);
  }
  if (mask[bin2] <= 0.0f) {
    fprintf(stderr, "FAIL: Existing tone at bin %d was not detected\n", bin2);
    exit(1);
  }
  printf("  Run 3: Multiple tones detected ✓\n");

  // Run 4: Manual CV mask mode with high-frequency bins > 99
  float manual_cv_mask[TEST_SPECTRUM_SIZE];
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    manual_cv_mask[i] = 0.0f;
  }
  int cv_bin = 150; // Above index 99 to ensure full spectrum scanning
  manual_cv_mask[cv_bin] = 1.0f;
  noise_spectrum[cv_bin] = 0.1f; // Tonality needs profile energy to notch
  int cv_bin2 = 180;             // Second high-frequency bin above index 99
  manual_cv_mask[cv_bin2] = 0.5f;

  tonal_reducer_compute_split(reducer, noise_spectrum, manual_cv_mask, true,
                              0.0f, noise_bb, noise_tonal);
  mask = tonal_reducer_get_mask(reducer);

  if (mask[cv_bin] <= 0.0f) {
    fprintf(stderr,
            "FAIL: Reducer in manual CV mask mode did not set mask at bin %d\n",
            cv_bin);
    exit(1);
  }
  if (fabsf(mask[cv_bin2] - 0.5f) > 1e-6f) {
    fprintf(stderr,
            "FAIL: Reducer in manual CV mask mode did not preserve fractional "
            "mask at bin %d (got %f)\n",
            cv_bin2, mask[cv_bin2]);
    exit(1);
  }
  if (noise_tonal[cv_bin] <= 0.0f) {
    fprintf(stderr, "FAIL: No residual at manual CV bin %d\n", cv_bin);
    exit(1);
  }
  printf("  Run 4: CV mask directly applied at bins %d and %d ✓\n", cv_bin,
         cv_bin2);

  // Run 5: Available all-zero CV mask honored without adaptive detection
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    manual_cv_mask[i] = 0.0f;
  }
  tonal_reducer_compute_split(reducer, noise_spectrum, manual_cv_mask, true,
                              0.0f, noise_bb, noise_tonal);
  mask = tonal_reducer_get_mask(reducer);
  for (int i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    if (mask[i] > 0.0f) {
      fprintf(stderr,
              "FAIL: Available all-zero CV mask produced non-zero mask at bin "
              "%d\n",
              i);
      exit(1);
    }
    if (noise_bb[i] != noise_spectrum[i]) {
      fprintf(
          stderr,
          "FAIL: Available all-zero CV mask modified bb profile at bin %d\n",
          i);
      exit(1);
    }
  }
  printf("  Run 5: Available all-zero CV mask honored ✓\n");

  tonal_reducer_free(reducer);
  printf("✓ Adaptive support and manual mask tests passed\n");
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

  float noise_spectrum[TEST_SPECTRUM_SIZE] = {0.01f};
  float noise_bb[TEST_SPECTRUM_SIZE];
  float noise_tonal[TEST_SPECTRUM_SIZE];
  // Gain 1.0f (no reduction requested, split short-circuits)
  tonal_reducer_compute_split(reducer, noise_spectrum, NULL, false, 1.0f,
                              noise_bb, noise_tonal);

  tonal_reducer_free(reducer);
  printf("✓ Tonal reducer peaks test passed\n");
}

int main(void) {
  printf("=== Tonal Reducer Tests ===\n\n");
  test_initialization();
  test_flat_noise_identity_split();
  test_tonal_split_and_gains();
  test_disabled_reduction_is_legacy();
  test_mask_refresh_and_adaptive_support();
  test_tonal_reducer_peaks();
  printf("\n=== All tonal reducer tests passed ===\n");
  return 0;
}
