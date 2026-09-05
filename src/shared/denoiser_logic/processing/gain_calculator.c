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
                               float* gain_spectrum,
                               const float* knee_spectrum) {
  uint32_t k = 0;
  sb_vec8_t flt_min = sb_set8(FLT_MIN);
  sb_vec8_t zero = sb_set8(0.0f);
  sb_vec8_t one = sb_set8(1.0f);

  for (; k + 7 < real_spectrum_size; k += 8) {
    sb_vec8_t s = sb_load8(spectrum + k);
    sb_vec8_t n = sb_load8(noise_spectrum + k);
    sb_vec8_t a = sb_load8(alpha + k);
    sb_vec8_t scaled_noise = sb_mul8(n, a);

    // Soft knee: lift the numerator by knee*noise so bins slightly below the
    // oversubtraction threshold keep a small positive gain (finite slope at
    // the floor crossing instead of an abrupt cut). knee = 0 reproduces the
    // plain subtraction curve.
    sb_vec8_t knee_vec = knee_spectrum ? sb_load8(knee_spectrum + k) : zero;
    sb_vec8_t knee_noise = sb_mul8(n, knee_vec);
    sb_vec8_t denom = sb_add8(s, knee_noise);

    sb_vec8_t mask_noise = sb_gt8(scaled_noise, flt_min);
    sb_vec8_t mask_gain = sb_gt8(denom, scaled_noise);

    sb_vec8_t safe_denom = sb_sel8(mask_gain, denom, one);
    sb_vec8_t diff = sb_sel8(mask_gain, sb_sub8(denom, scaled_noise), zero);
    sb_vec8_t power_gain = sb_div8(diff, safe_denom);
    sb_vec8_t gain = sb_sqrt8(power_gain);
    gain = sb_sel8(mask_gain, gain, zero);
    gain = sb_sel8(mask_noise, gain, one);

    sb_store8(gain_spectrum + k, gain);
  }

  for (; k < real_spectrum_size; k++) {
    float knee = knee_spectrum ? knee_spectrum[k] : 0.0F;
    float scaled_noise = noise_spectrum[k] * alpha[k];
    if (scaled_noise > FLT_MIN) {
      float denom = spectrum[k] + (noise_spectrum[k] * knee);
      float diff = denom - scaled_noise;
      if (diff > 0.0F) {
        gain_spectrum[k] = sqrtf(diff / denom);
      } else {
        gain_spectrum[k] = 0.F;
      }
    } else {
      gain_spectrum[k] = 1.F;
    }
  }

  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
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

  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}

static void generalized_spectral_subtraction(
    const uint32_t real_spectrum_size, const uint32_t fft_size,
    const float* spectrum, const float* noise_spectrum, float* gain_spectrum,
    const float* alpha, const float* beta) {
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    if (spectrum[k] > FLT_MIN) {
      float ratio_sq = noise_spectrum[k] / spectrum[k];
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

  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}

void calculate_gains(uint32_t real_spectrum_size, uint32_t fft_size,
                     const float* spectrum, const float* noise_spectrum,
                     float* gain_spectrum, const float* alpha,
                     const float* beta, GainCalculationType type,
                     const float* knee) {
  switch (type) {
    case GATES:
      spectral_gating(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                      alpha, gain_spectrum);
      break;
    case WIENER:
      wiener_subtraction(real_spectrum_size, fft_size, spectrum, noise_spectrum,
                         alpha, gain_spectrum, knee);
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
