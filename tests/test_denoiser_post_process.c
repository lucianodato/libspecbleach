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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/denoiser_logic/core/denoiser_post_process.h"
#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/processing/tonal_reducer.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_post_process_null_guards(void) {
  printf("Testing denoiser post process NULL guards...\n");

  uint32_t fft_size = 8;
  uint32_t real_size = (fft_size / 2) + 1;

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));
  float* fft_spectrum = (float*)calloc(fft_size, sizeof(float));

  for (uint32_t i = 0; i < fft_size; i++) {
    gain_spectrum[i] = 0.5f;
    fft_spectrum[i] = 1.0f;
  }
  for (uint32_t i = 0; i < real_size; i++) {
    noise_spectrum[i] = 0.1f;
  }

  // Case 1: Both components are NULL
  DenoiserPostProcessParams params1 = {
      .fft_size = fft_size,
      .real_spectrum_size = real_size,
      .reduction_amount = 0.1f,
      .tonal_reduction = 0.1f,
      .whitening_factor = 0.0f,
      .residual_listen = false,
      .noise_floor_manager = NULL,
      .tonal_reducer = NULL,
      .gain_spectrum = gain_spectrum,
      .noise_spectrum = noise_spectrum,
      .fft_spectrum = fft_spectrum,
  };

  denoiser_post_process_apply(params1);

  // Since both are NULL, it should bypass noise floor manager and just do
  // mixing fft_spectrum should be original * gain = 1.0 * 0.5 = 0.5
  for (uint32_t i = 0; i < fft_size; i++) {
    TEST_FLOAT_CLOSE(fft_spectrum[i], 0.5f, 0.0001f);
  }

  // Reset fft_spectrum
  for (uint32_t i = 0; i < fft_size; i++) {
    fft_spectrum[i] = 1.0f;
  }

  // Case 2: noise_spectrum is NULL
  DenoiserPostProcessParams params2 = params1;
  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  params2.noise_floor_manager = nfm;
  params2.noise_spectrum = NULL;

  denoiser_post_process_apply(params2);

  // Should also bypass noise floor manager (due to NULL noise_spectrum guard)
  // and do mixing
  for (uint32_t i = 0; i < fft_size; i++) {
    TEST_FLOAT_CLOSE(fft_spectrum[i], 0.5f, 0.0001f);
  }

  noise_floor_manager_free(nfm);
  free(gain_spectrum);
  free(noise_spectrum);
  free(fft_spectrum);
  printf("✓ NULL guards tests passed\n");
}

void test_post_process_mixing(void) {
  printf("Testing denoiser post process mixing modes...\n");

  uint32_t fft_size = 4;
  uint32_t real_size = (fft_size / 2) + 1;

  float gain_spectrum[4] = {0.8f, 0.5f, 0.5f, 0.8f};
  float noise_spectrum[3] = {0.1f, 0.1f, 0.1f};
  float fft_spectrum[4] = {1.0f, 2.0f, 3.0f, 4.0f};

  // Normal Mixing Mode
  DenoiserPostProcessParams params = {
      .fft_size = fft_size,
      .real_spectrum_size = real_size,
      .reduction_amount = 0.1f,
      .tonal_reduction = 0.1f,
      .whitening_factor = 0.0f,
      .residual_listen = false,
      .noise_floor_manager = NULL,
      .tonal_reducer = NULL,
      .gain_spectrum = gain_spectrum,
      .noise_spectrum = noise_spectrum,
      .fft_spectrum = fft_spectrum,
  };

  denoiser_post_process_apply(params);

  TEST_FLOAT_CLOSE(fft_spectrum[0], 0.8f, 0.0001f);
  TEST_FLOAT_CLOSE(fft_spectrum[1], 1.0f, 0.0001f);
  TEST_FLOAT_CLOSE(fft_spectrum[2], 1.5f, 0.0001f);
  TEST_FLOAT_CLOSE(fft_spectrum[3], 3.2f, 0.0001f);

  // Reset fft_spectrum
  fft_spectrum[0] = 1.0f;
  fft_spectrum[1] = 2.0f;
  fft_spectrum[2] = 3.0f;
  fft_spectrum[3] = 4.0f;

  // Residual Listen Mode
  params.residual_listen = true;
  denoiser_post_process_apply(params);

  // original - original * gain = original * (1 - gain)
  // 1.0 * (1 - 0.8) = 0.2
  // 2.0 * (1 - 0.5) = 1.0
  // 3.0 * (1 - 0.5) = 1.5
  // 4.0 * (1 - 0.8) = 0.8
  TEST_FLOAT_CLOSE(fft_spectrum[0], 0.2f, 0.0001f);
  TEST_FLOAT_CLOSE(fft_spectrum[1], 1.0f, 0.0001f);
  TEST_FLOAT_CLOSE(fft_spectrum[2], 1.5f, 0.0001f);
  TEST_FLOAT_CLOSE(fft_spectrum[3], 0.8f, 0.0001f);

  printf("✓ Mixing modes tests passed\n");
}

void test_post_process_gain_patching(void) {
  printf("Testing denoiser post process spectral gain patching...\n");

  // We need a larger spectrum size to test neighbors up to 6 bins
  uint32_t fft_size = 32;
  uint32_t real_size = (fft_size / 2) + 1;

  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* noise_spectrum = (float*)calloc(real_size, sizeof(float));
  float* fft_spectrum = (float*)calloc(fft_size, sizeof(float));

  for (uint32_t i = 0; i < fft_size; i++) {
    gain_spectrum[i] = 0.5f;
    fft_spectrum[i] = 1.0f;
  }
  for (uint32_t i = 0; i < real_size; i++) {
    noise_spectrum[i] = 0.1f;
  }

  NoiseFloorManager* nfm = noise_floor_manager_initialize(fft_size);
  TonalReducer* tr = tonal_reducer_initialize(real_size, 44100, fft_size);

  // We will run tonal reducer to produce a tonal mask.
  // We want a tone at bin 8.
  // In adaptive mode (all profiles 0):
  float reducer_noise[17];
  float reducer_max[17];
  float reducer_median[17];
  float reducer_alpha[17];
  for (int i = 0; i < 17; i++) {
    reducer_noise[i] = 0.01f;
    reducer_max[i] = 0.0f;
    reducer_median[i] = 0.0f;
    reducer_alpha[i] = 1.0f;
  }
  reducer_noise[8] = 0.1f;

  tonal_reducer_run(tr, reducer_noise, reducer_max, reducer_median,
                    reducer_alpha, 0.0f);

  // Verify that we got a mask at bin 8
  const float* mask = tonal_reducer_get_mask(tr);
  TEST_ASSERT(mask[8] > 0.0f, "Should have tonal mask at bin 8");
  TEST_ASSERT(mask[7] == 0.0f, "Should not have tonal mask at bin 7");

  // Case 1: Interpolation with both left and right neighbors
  // Let's set neighbor gains at 7 and 9
  gain_spectrum[7] = 0.4f;
  gain_spectrum[9] = 0.6f;
  gain_spectrum[8] = 0.1f; // original gain at tonal bin

  DenoiserPostProcessParams params = {
      .fft_size = fft_size,
      .real_spectrum_size = real_size,
      .reduction_amount = 0.1f,
      .tonal_reduction = 0.0f,
      .whitening_factor = 0.0f,
      .residual_listen = false,
      .noise_floor_manager = nfm,
      .tonal_reducer = tr,
      .gain_spectrum = gain_spectrum,
      .noise_spectrum = noise_spectrum,
      .fft_spectrum = fft_spectrum,
  };

  denoiser_post_process_apply(params);

  // Bin 8 should be patched using interpolation of 7 and 9.
  // t = (8-7)/(9-7) = 0.5.
  // interpolated g_floor = 0.4 + 0.5 * (0.6 - 0.4) = 0.5.
  // Since s_noise/s_in and w_sig will affect target_factor:
  // With some signal presence or default, patched_gain will be blended with
  // original. Let's verify that the patched gain at bin 8 is strictly greater
  // than the original 0.1.
  TEST_ASSERT(gain_spectrum[8] > 0.15f,
              "Gain patching should increase gain at tonal bin 8");
  // Symmetry check: fft_size - k = 32 - 8 = 24.
  TEST_FLOAT_CLOSE(gain_spectrum[24], gain_spectrum[8], 0.0001f);

  // Case 2: Only left neighbor available (tonal bins on the right)
  // Let's reinitialize
  for (uint32_t i = 0; i < fft_size; i++) {
    gain_spectrum[i] = 0.5f;
  }
  gain_spectrum[8] = 0.1f;
  gain_spectrum[7] = 0.4f;

  // Let's force a mask on bins 9, 10, 11, 12, 13, 14 (so right neighbor is at
  // least 6 bins away i.e. not found) We can do this by setting them as tonal.
  // But since we use real tonal reducer, let's just make sure they are tonal.
  // Wait! We can initialize a mock mask by running tonal_reducer_run with peaks
  // at 8, 9, 10, 11, 12, 13, 14.
  for (int i = 8; i <= 14; i++) {
    reducer_noise[i] = 0.1f;
  }
  tonal_reducer_run(tr, reducer_noise, reducer_max, reducer_median,
                    reducer_alpha, 0.0f);

  denoiser_post_process_apply(params);

  // For bin 8:
  // Left neighbor is 7 (gain_spectrum[7] = 0.4).
  // Right neighbor is unavailable (>= 6 bins away).
  // g_floor should be gain_spectrum[7] = 0.4.
  // Let's check that patching is active and gain_spectrum[8] is updated.
  TEST_ASSERT(gain_spectrum[8] > 0.15f,
              "Gain patching with only left neighbor should be active");

  // Case 2b: Only right neighbor available (tonal bins on the left)
  // Let's reinitialize
  for (uint32_t i = 0; i < fft_size; i++) {
    gain_spectrum[i] = 0.5f;
  }
  gain_spectrum[8] = 0.1f;
  gain_spectrum[9] = 0.6f;

  // Let's force a mask on bins 2..8 (so left neighbor is 7 bins away)
  for (int i = 0; i < 17; i++) {
    reducer_noise[i] = 0.01f;
  }
  for (int i = 2; i <= 8; i++) {
    reducer_noise[i] = 0.1f;
  }
  tonal_reducer_run(tr, reducer_noise, reducer_max, reducer_median,
                    reducer_alpha, 0.0f);

  denoiser_post_process_apply(params);

  // For bin 8:
  // Left neighbor is unavailable (>= 6 bins away).
  // Right neighbor is 9 (gain_spectrum[9] = 0.6).
  // g_floor should be gain_spectrum[9] = 0.6.
  // Let's check that patching is active and gain_spectrum[8] is updated.
  TEST_ASSERT(gain_spectrum[8] > 0.15f,
              "Gain patching with only right neighbor should be active");

  // Case 3: No neighbor available fallback (both left and right unavailable)
  // Let's clear neighbors.
  // Run tonal reducer with tones covering 1 to 14, so for bin 7, there are no
  // neighbors within 6 bins.
  for (int i = 1; i <= 14; i++) {
    reducer_noise[i] = 0.1f;
  }
  tonal_reducer_run(tr, reducer_noise, reducer_max, reducer_median,
                    reducer_alpha, 0.0f);

  for (uint32_t i = 0; i < fft_size; i++) {
    gain_spectrum[i] = 0.35f;
  }
  // Bin 7 gain
  gain_spectrum[7] = 0.35f;

  denoiser_post_process_apply(params);

  // Since both left and right neighbors are unavailable, patched_gain should
  // remain g_floor (original gain, 0.35f). And the final gain should remain
  // exactly 0.35f (it shouldn't drop to 0.0f!).
  TEST_FLOAT_CLOSE(gain_spectrum[7], 0.35f, 0.0001f);

  noise_floor_manager_free(nfm);
  tonal_reducer_free(tr);
  free(gain_spectrum);
  free(noise_spectrum);
  free(fft_spectrum);
  printf("✓ Spectral gain patching tests passed\n");
}

int main(void) {
  printf("=== Denoiser Post Process Tests ===\n\n");

  test_post_process_null_guards();
  test_post_process_mixing();
  test_post_process_gain_patching();

  printf("\n=== All denoiser post process tests passed ===\n");
  return 0;
}
