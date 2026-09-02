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

#include "../src/shared/time_domain/minimum_phase_fir.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FFT_SIZE 512U
#define TAPS 256U
#define REAL_SIZE (FFT_SIZE / 2U + 1U)

static void test_unity_gains_passthrough(void) {
  // Unity magnitude spectrum -> minimum-phase IR ~ delta at n=0
  FftTransform* fft = fft_transform_initialize_bins(FFT_SIZE);
  assert(fft != NULL);

  float gains[REAL_SIZE];
  for (uint32_t k = 0U; k < REAL_SIZE; k++) {
    gains[k] = 1.0F;
  }

  float taper[TAPS];
  assert(minimum_phase_fir_taper_window(TAPS, taper));

  float ir[TAPS];
  assert(minimum_phase_fir_synthesize(FFT_SIZE, TAPS, gains, fft, taper, ir));

  // Energy must concentrate at t=0 (delta passthrough)
  const float energy_first = fabsf(ir[0]);
  float energy_rest = 0.F;
  for (uint32_t i = 1U; i < TAPS; i++) {
    energy_rest += fabsf(ir[i]);
  }
  assert(energy_first > 0.9F);
  assert(energy_rest < 0.1F);
  assert(energy_first > energy_rest);
  printf("  delta energy first=%.6f rest=%.6f\n", energy_first, energy_rest);

  fft_transform_free(fft);
}

static void test_lowpass_minimum_phase(void) {
  FftTransform* fft = fft_transform_initialize_bins(FFT_SIZE);
  assert(fft != NULL);

  // Simple lowpass magnitude: 1.0 below bin 32, rolling to 0.1 above
  float gains[REAL_SIZE];
  for (uint32_t k = 0U; k < REAL_SIZE; k++) {
    gains[k] = (k < 32U) ? 1.0F : 0.1F;
  }

  float taper[TAPS];
  assert(minimum_phase_fir_taper_window(TAPS, taper));

  float ir[TAPS];
  assert(minimum_phase_fir_synthesize(FFT_SIZE, TAPS, gains, fft, taper, ir));

  // Finite values only
  for (uint32_t i = 0U; i < TAPS; i++) {
    assert(isfinite(ir[i]));
  }

  // Minimum phase: energy concentrated in the first quarter of taps
  float energy_early = 0.F;
  float energy_late = 0.F;
  for (uint32_t i = 0U; i < TAPS / 4U; i++) {
    energy_early += ir[i] * ir[i];
  }
  for (uint32_t i = TAPS / 4U; i < TAPS; i++) {
    energy_late += ir[i] * ir[i];
  }
  assert(energy_early > energy_late);
  assert(energy_early > 0.F);
  printf("  lowpass energy early=%.6f late=%.6f\n", energy_early, energy_late);

  fft_transform_free(fft);
}

static void test_constant_gain_delta(void) {
  // Constant gain c -> minimum-phase IR ~ c * delta at n=0
  FftTransform* fft = fft_transform_initialize_bins(FFT_SIZE);
  assert(fft != NULL);

  float gains[REAL_SIZE];
  const float gain = 0.1F;
  for (uint32_t k = 0U; k < REAL_SIZE; k++) {
    gains[k] = gain;
  }

  float taper[TAPS];
  assert(minimum_phase_fir_taper_window(TAPS, taper));

  float ir[TAPS];
  assert(minimum_phase_fir_synthesize(FFT_SIZE, TAPS, gains, fft, taper, ir));

  assert(fabsf(ir[0] - gain) < 1e-3F);
  for (uint32_t i = 1U; i < TAPS; i++) {
    assert(fabsf(ir[i]) < 1e-3F);
  }

  fft_transform_free(fft);
}

static void test_invalid_arguments(void) {
  FftTransform* fft = fft_transform_initialize_bins(FFT_SIZE);
  assert(fft != NULL);

  float gains[REAL_SIZE] = {0};
  float taper[TAPS] = {0};
  float ir[TAPS] = {0};

  assert(!minimum_phase_fir_synthesize(0U, TAPS, gains, fft, taper, ir));
  assert(!minimum_phase_fir_synthesize(FFT_SIZE, 0U, gains, fft, taper, ir));
  assert(!minimum_phase_fir_synthesize(FFT_SIZE, TAPS, NULL, fft, taper, ir));
  assert(!minimum_phase_fir_synthesize(FFT_SIZE, TAPS, gains, NULL, taper, ir));
  assert(!minimum_phase_fir_synthesize(FFT_SIZE, TAPS, gains, fft, NULL, ir));
  assert(
      !minimum_phase_fir_synthesize(FFT_SIZE, TAPS, gains, fft, taper, NULL));
  assert(!minimum_phase_fir_taper_window(0U, taper));
  assert(!minimum_phase_fir_taper_window(TAPS, NULL));

  fft_transform_free(fft);
}

int main(void) {
  test_unity_gains_passthrough();
  test_constant_gain_delta();
  test_lowpass_minimum_phase();
  test_invalid_arguments();

  printf("✓ Minimum phase FIR tests passed\n");
  return 0;
}
