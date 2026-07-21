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

int main(void) {
  printf("=== Denoiser Post Process Tests ===\n\n");

  test_post_process_null_guards();
  test_post_process_mixing();

  printf("\n=== All denoiser post process tests passed ===\n");
  return 0;
}
