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
    sb_vec8_t scaled_noise = sb_mul8(n, a);

    sb_vec8_t mask_noise = sb_gt8(scaled_noise, flt_min);
    sb_vec8_t mask_gain = sb_gt8(s, scaled_noise);

    sb_vec8_t gain = sb_div8(sb_sub8(s, scaled_noise), s);
    gain = sb_sel8(mask_gain, gain, zero);
    gain = sb_sel8(mask_noise, gain, one);

    sb_store8(gain_spectrum + k, gain);
  }

  for (; k < real_spectrum_size; k++) {
    float scaled_noise = noise_spectrum[k] * alpha[k];
    if (scaled_noise > FLT_MIN) {
      if (spectrum[k] > scaled_noise) {
        gain_spectrum[k] = (spectrum[k] - (scaled_noise)) / spectrum[k];
      } else {
        gain_spectrum[k] = 0.F;
      }
    } else {
      gain_spectrum[k] = 1.F;
    }
  }

  // Symmetry
  for (uint32_t i = 1U; i < real_spectrum_size && i < (fft_size - i); i++) {
    gain_spectrum[fft_size - i] = gain_spectrum[i];
  }
}

static void spectral_gating(const uint32_t real_spectrum_size,
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
    sb_vec8_t scaled_noise = sb_mul8(n, a);

    sb_vec8_t mask_noise = sb_gt8(scaled_noise, flt_min);
    sb_vec8_t mask_gate_fail = sb_gt8(scaled_noise, s);

    sb_vec8_t gain = sb_sel8(mask_gate_fail, zero, one);
    gain = sb_sel8(mask_noise, gain, one);

    sb_store8(gain_spectrum + k, gain);
  }

  for (; k < real_spectrum_size; k++) {
    float scaled_noise = noise_spectrum[k] * alpha[k];
    if (scaled_noise > FLT_MIN) {
      if (spectrum[k] >= scaled_noise) {
        gain_spectrum[k] = 1.F;
      } else {
        gain_spectrum[k] = 0.F;
      }
    } else {
      gain_spectrum[k] = 1.F;
    }
  }

  // Symmetry
  for (uint32_t i = 1U; i < real_spectrum_size && i < (fft_size - i); i++) {
    gain_spectrum[fft_size - i] = gain_spectrum[i];
  }
}

static void generalized_spectral_subtraction(
    const uint32_t real_spectrum_size, const uint32_t fft_size,
    const float* spectrum, const float* noise_spectrum, float* gain_spectrum,
    const float* alpha, const float* beta) {
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    if (spectrum[k] > FLT_MIN) {
      // Use multiplications instead of powf for exponent 2.0
      float ratio = noise_spectrum[k] / spectrum[k];
      float ratio_sq = ratio * ratio;
      if (ratio_sq < (1.F / (alpha[k] + beta[k]))) {
        float val = fmaxf(1.F - (alpha[k] * ratio_sq), 0.0f);
        gain_spectrum[k] = fmaxf(sqrtf(val), 0.F);
      } else {
        float val = fmaxf(beta[k] * ratio_sq, 0.0f);
        gain_spectrum[k] = fmaxf(sqrtf(val), 0.F);
      }
    } else {
      gain_spectrum[k] = 1.F;
    }
  }

  // Symmetry
  for (uint32_t i = 1U; i < real_spectrum_size && i < (fft_size - i); i++) {
    gain_spectrum[fft_size - i] = gain_spectrum[i];
  }
}

void calculate_gains(uint32_t real_spectrum_size, uint32_t fft_size,
                     const float* spectrum, const float* noise_spectrum,
                     float* gain_spectrum, const float* alpha,
                     const float* beta, GainCalculationType type) {
  switch (type) {
    case GATES:
      spectral_gating(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                      alpha, gain_spectrum);
      break;
    case WIENER:
      wiener_subtraction(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                         alpha, gain_spectrum);
      break;
    case GENERALIZED_SPECTRALSUBTRACTION:
      generalized_spectral_subtraction(real_spectrum_size, fft_size, spectrum,
                                       noise_spectrum, gain_spectrum, alpha,
                                       beta);
      break;

    default:
      break;
  }
}
