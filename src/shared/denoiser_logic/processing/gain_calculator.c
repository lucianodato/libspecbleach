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

#include "shared/denoiser_logic/processing/gain_calculator.h"
#include "shared/configurations.h"
#include "shared/utils/simd_utils.h"
#include "shared/utils/spectral_utils.h"
#include <float.h>
#include <math.h>

static void wiener_subtraction(const uint32_t real_spectrum_size,
                               const uint32_t fft_size, const float* spectrum,
                               const float* noise_spectrum, const float* alpha,
                               float* gain_spectrum) {
  uint32_t k = 0;
  sb_vec8_t flt_min = sb_set8(FLT_MIN);
  sb_vec8_t zero = sb_set8(0.0f);
  sb_vec8_t one = sb_set8(1.0f);

  for (; k + 7 < real_spectrum_size; k += 8) {
    sb_vec8_t s = sb_load8(spectrum + k);
    sb_vec8_t n = sb_load8(noise_spectrum + k);
    sb_vec8_t a = sb_load8(alpha + k);
    sb_vec8_t num = sb_sub8(s, sb_mul8(a, n));
    num = sb_sel8(sb_gt8(num, zero), num, zero);
    sb_vec8_t denom = sb_sel8(sb_gt8(s, flt_min), s, flt_min);
    sb_vec8_t gain = sb_div8(num, denom);
    gain = sb_sel8(sb_gt8(gain, one), one, gain);
    sb_store8(gain_spectrum + k, gain);
  }

  for (; k < real_spectrum_size; k++) {
    float s = spectrum[k];
    float n = noise_spectrum[k];
    float a = alpha[k];
    float num = s - (a * n);
    if (num < 0.0f) num = 0.0f;
    float denom = s > FLT_MIN ? s : FLT_MIN;
    float gain = num / denom;
    if (gain > 1.0f) gain = 1.0f;
    gain_spectrum[k] = gain;
  }

  for (uint32_t i = 1; i < real_spectrum_size - 1; i++) {
    gain_spectrum[i] = 0.25f * gain_spectrum[i - 1] + 0.5f * gain_spectrum[i] + 0.25f * gain_spectrum[i + 1];
  }

  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}

static void gates_subtraction(const uint32_t real_spectrum_size,
                              const uint32_t fft_size, const float* spectrum,
                              const float* noise_spectrum, const float* alpha,
                              float* gain_spectrum) {
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float s = spectrum[k];
    float n = noise_spectrum[k];
    float a = alpha[k];
    if (s >= a * n) {
      gain_spectrum[k] = 1.0f;
    } else {
      gain_spectrum[k] = 0.0f;
    }
  }
  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}

static void generalized_spectral_subtraction(const uint32_t real_spectrum_size,
                                             const uint32_t fft_size,
                                             const float* spectrum,
                                             const float* noise_spectrum,
                                             const float* alpha,
                                             const float* beta,
                                             float* gain_spectrum) {
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float s = spectrum[k];
    float n = noise_spectrum[k];
    float a = alpha[k];
    float b = beta[k];

    float num = s - (a * powf(n, GSS_EXPONENT));
    if (num < 0.0f) num = b * powf(n, GSS_EXPONENT);
    float denom = s > FLT_MIN ? s : FLT_MIN;
    float gain = powf(num / denom, 1.0f / GSS_EXPONENT);
    if (gain > 1.0f) gain = 1.0f;
    if (gain < 0.0f) gain = 0.0f;
    gain_spectrum[k] = gain;
  }
  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}

void calculate_gains(uint32_t real_spectrum_size, uint32_t fft_size,
                     const float* spectrum, const float* noise_spectrum,
                     float* gain_spectrum, const float* alpha,
                     const float* beta, GainCalculationType type) {
  if (!spectrum || !noise_spectrum || !gain_spectrum || !alpha) {
    return;
  }

  switch (type) {
    case WIENER:
      wiener_subtraction(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                           alpha, gain_spectrum);
      break;
    case GATES:
      gates_subtraction(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                        alpha, gain_spectrum);
      break;
    case GENERALIZED_SPECTRALSUBTRACTION:
      generalized_spectral_subtraction(real_spectrum_size, fft_size, spectrum,
                                         noise_spectrum, alpha, beta,
                                         gain_spectrum);
      break;
    default:
      wiener_subtraction(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                           alpha, gain_spectrum);
      break;
    }
}
