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

#include "shared/utils/spectral_smoother.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

void test_spectral_smoother(void) {
  printf("Testing Spectral Smoother...\n");

  uint32_t fft_size = 1024;

  for (int type = FIXED; type <= TRANSIENT_AWARE; type++) {
    SpectralSmoother* ss =
        spectral_smoothing_initialize(fft_size, 44100, (TimeSmoothingType)type);
    TEST_ASSERT(ss != NULL, "Spectral smoother initialization should succeed");

    float spectrum[513] = {0.0f};
    for (int i = 0; i < 513; i++) {
      spectrum[i] = 1.0f + (0.5f * sinf((float)i * 0.1f));
    }

    TimeSmoothingParameters params = {.smoothing = 0.8f};
    TEST_ASSERT(spectral_smoothing_run(ss, params, spectrum),
                "Spectral smoothing should succeed");

    // Run again to test previous spectrum logic
    TEST_ASSERT(spectral_smoothing_run(ss, params, spectrum),
                "Spectral smoothing should succeed on second run");

    // Check that output is reasonable
    for (int i = 0; i < 513; i++) {
      TEST_ASSERT(spectrum[i] >= 0.0f,
                  "Smoothed spectrum should be non-negative");
    }

    spectral_smoothing_free(ss);
  }

  printf("✓ Spectral Smoother tests passed\n");
}

int main(void) {
  test_spectral_smoother();
  return 0;
}
