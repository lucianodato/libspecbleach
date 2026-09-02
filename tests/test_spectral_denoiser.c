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

#include "processors/denoiser/spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "Assertion failed: %s (%s:%d)\n", message, __FILE__,     \
              __LINE__);                                                       \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

int main(void) {
  printf("Testing spectral denoiser internal API validation...\n");

  // Valid noise profile for a 2048-point FFT (real spectrum size is 1025)
  NoiseProfile* profile = noise_profile_initialize(1025U);
  TEST_ASSERT(profile != NULL, "Noise profile allocation should succeed");

  // Invalid initialization arguments: each invalid field isolated with a
  // valid profile
  TEST_ASSERT(spectral_denoiser_initialize(0U, 2048U, 4U, profile) == NULL,
              "Zero sample rate must fail");
  TEST_ASSERT(spectral_denoiser_initialize(44100U, 0U, 4U, profile) == NULL,
              "Zero fft size must fail");
  TEST_ASSERT(spectral_denoiser_initialize(44100U, 2048U, 0U, profile) == NULL,
              "Zero overlap must fail");
  TEST_ASSERT(spectral_denoiser_initialize(44100U, 2048U, 4U, NULL) == NULL,
              "NULL noise profile must fail");

  // NULL-handle guards
  float spectrum[1026] = {0};
  TEST_ASSERT(load_reduction_parameters(NULL, (DenoiserParameters){0}) == false,
              "NULL load parameters must fail");
  TEST_ASSERT(spectral_denoiser_run(NULL, spectrum) == false,
              "NULL instance run must fail");
  TEST_ASSERT(spectral_denoiser_get_tonal_mask(NULL) == NULL,
              "NULL tonal mask must be NULL");
  TEST_ASSERT(spectral_denoiser_get_active_noise_profile(NULL) == NULL,
              "NULL active profile must be NULL");
  TEST_ASSERT(spectral_denoiser_get_latency_frames(NULL) == 0,
              "NULL latency must be 0");
  TEST_ASSERT(spectral_denoiser_is_transient_detected(NULL) == false,
              "NULL transient detected must be false");
  TEST_ASSERT(spectral_denoiser_get_transient_intensity(NULL) == 0.0f,
              "NULL transient intensity must be 0");
  float peaks[8];
  TEST_ASSERT(spectral_denoiser_get_peaks(NULL, peaks, 8) == 0,
              "NULL peaks must return 0");
  spectral_denoiser_reset_noise_profile(NULL); // must not crash

  // Valid handle: NULL spectrum must fail and reset must be safe
  SpectralProcessorHandle handle =
      spectral_denoiser_initialize(44100U, 2048U, OVERLAP_FACTOR_1D, profile);
  TEST_ASSERT(handle != NULL, "Initialization should succeed");

  TEST_ASSERT(spectral_denoiser_run(handle, NULL) == false,
              "NULL spectrum with valid handle must fail");
  TEST_ASSERT(spectral_denoiser_get_latency_frames(handle) > 0,
              "Latency frames must be positive");
  spectral_denoiser_reset_noise_profile(handle);

  spectral_denoiser_free(handle);
  noise_profile_free(profile);

  printf("✅ All spectral denoiser internal tests passed!\n");
  return 0;
}
