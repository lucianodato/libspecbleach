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
#include <stdbool.h>
#include <stdint.h>

static void wiener_subtraction(const uint32_t real_spectrum_size,
                               const float* spectrum,
                               const float* noise_spectrum,
                               const float* alpha,
                               const float* beta,
                               float* gain_spectrum) {
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float s = spectrum[k];
    float n = noise_spectrum[k];
    float a = alpha ? alpha[k] : 1.0f;
    float b = beta ? beta[k] : 0.0f;

    float sub = s - (a * n);
    float g = (s > SPECTRAL_EPSILON) ? (sub / s) : 0.0f;

    if (g < b) {
      g = b;
    }
    if (g > 1.0f) {
      g = 1.0f;
    }
    gain_spectrum[k] = g;
  }
}

static void gates_subtraction(const uint32_t real_spectrum_size,
                              const float* spectrum,
                              const float* noise_spectrum,
                              const float* alpha,
                              const float* beta,
                              float* gain_spectrum) {
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float s = spectrum[k];
    float n = noise_spectrum[k];
    float a = alpha ? alpha[k] : 1.0f;
    float b = beta ? beta[k] : 0.0f;

    float snr_lin = (n > SPECTRAL_EPSILON) ? (s / n) : 100.0f;
    float diff = snr_lin - a;
    float gate = 1.0f / (1.0f + expf(-2.0f * diff));
    float g = b + (1.0f - b) * gate;

    if (g < b) {
      g = b;
    }
    if (g > 1.0f) {
      g = 1.0f;
    }
    gain_spectrum[k] = g;
  }
}

static void generalized_spectral_subtraction(const uint32_t real_spectrum_size,
                                            const float* spectrum,
                                            const float* noise_spectrum,
                                            const float* alpha,
                                            const float* beta,
                                            float* gain_spectrum) {
  float exp_val = GSS_EXPONENT;
  float half_exp = exp_val * 0.5f;

  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float s = spectrum[k];
    float n = noise_spectrum[k];
    float a = alpha ? alpha[k] : 1.0f;
    float b = beta ? beta[k] : 0.0f;

    if (s <= SPECTRAL_EPSILON) {
      gain_spectrum[k] = b;
      continue;
    }

    float s_pow = powf(s, half_exp);
    float n_pow = powf(n, half_exp);
    float sub = s_pow - (a * n_pow);

    float g_pow = sub / s_pow;
    if (g_pow < 0.0f) {
      g_pow = 0.0f;
    }

    float g = powf(g_pow, 1.0f / half_exp);

    if (g < b) {
      g = b;
    }
    if (g > 1.0f) {
      g = 1.0f;
    }
    gain_spectrum[k] = g;
  }
}

static void apply_neighborhood_gain_coupling(const uint32_t real_spectrum_size,
                                              const float* spectrum,
                                              float* gain_spectrum) {
  if (real_spectrum_size < 3U) {
    return;
  }

  for (uint32_t k = 1U; k < real_spectrum_size - 1U; k++) {
    float g_prev = gain_spectrum[k - 1U];
    float g_curr = gain_spectrum[k];
    float g_next = gain_spectrum[k + 1U];

    float g_max_neighbor = fmaxf(g_prev, g_next);
    float g_min_neighbor = fminf(g_prev, g_next);

    float s_prev = spectrum[k - 1U];
    float s_curr = spectrum[k];
    float s_next = spectrum[k + 1U];
    float s_max_neighbor = fmaxf(s_prev, s_next) + SPECTRAL_EPSILON;
    float tonality_ratio = s_curr / s_max_neighbor;

    if (g_curr > g_max_neighbor) {
      if (tonality_ratio < 1.35f) {
        float blend = (tonality_ratio > 1.0f) ? ((tonality_ratio - 1.0f) / 0.35f) : 0.0f;
        gain_spectrum[k] = g_max_neighbor + (g_curr - g_max_neighbor) * blend;
      }
    } else if (g_curr < g_min_neighbor) {
      if (tonality_ratio > 0.85f) {
        float blend = (tonality_ratio < 1.0f) ? ((1.0f - tonality_ratio) / 0.15f) : 0.0f;
        gain_spectrum[k] = g_min_neighbor - (g_min_neighbor - g_curr) * blend;
      }
    }
  }
}

void calculate_gains(uint32_t real_spectrum_size, uint32_t fft_size,
                     const float* spectrum, const float* noise_spectrum,
                     float* gain_spectrum, const float* alpha,
                     const float* beta, GainCalculationType type) {
  if (!spectrum || !noise_spectrum || !gain_spectrum || real_spectrum_size == 0U) {
    return;
  }

  switch (type) {
    case WIENER:
      wiener_subtraction(real_spectrum_size, spectrum, noise_spectrum, alpha, beta, gain_spectrum);
      break;
    case GATES:
      gates_subtraction(real_spectrum_size, spectrum, noise_spectrum, alpha, beta, gain_spectrum);
      break;
    case GENERALIZED_SPECTRALSUBTRACTION:
      generalized_spectral_subtraction(real_spectrum_size, spectrum, noise_spectrum, alpha, beta, gain_spectrum);
      break;
    default:
      wiener_subtraction(real_spectrum_size, spectrum, noise_spectrum, alpha, beta, gain_spectrum);
      break;
  }

  apply_neighborhood_gain_coupling(real_spectrum_size, spectrum, gain_spectrum);

  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}
