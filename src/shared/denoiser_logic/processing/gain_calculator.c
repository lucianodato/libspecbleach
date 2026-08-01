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

    float noise_est = a * n;
    float diff = s - noise_est;
    float floor_val = b * n;

    float soft_diff = 0.5f * (diff + sqrtf(diff * diff + 1e-8f));
    float clean_est = fmaxf(soft_diff, floor_val);

    float gain = clean_est / (s + SPECTRAL_EPSILON);
    if (gain > 1.0f) gain = 1.0f;
    if (gain < 0.0f) gain = 0.0f;
    gain_spectrum[k] = gain;
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

    float threshold = a * n;
    float snr = s / (threshold + SPECTRAL_EPSILON);

    float log_snr = logf(snr + SPECTRAL_EPSILON);
    float g_exp = 1.0f / (1.0f + expf(-2.0f * log_snr));
    float gain = b + (1.0f - b) * g_exp;

    if (gain > 1.0f) gain = 1.0f;
    if (gain < 0.0f) gain = 0.0f;
    gain_spectrum[k] = gain;
  }
}

static void generalized_subtraction(const uint32_t real_spectrum_size,
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

    float p = GSS_EXPONENT;
    float s_p = powf(s + SPECTRAL_EPSILON, p);
    float n_p = powf(n + SPECTRAL_EPSILON, p);

    float diff_p = s_p - (a * n_p);
    float floor_p = b * n_p;

    float soft_diff_p = 0.5f * (diff_p + sqrtf(diff_p * diff_p + 1e-8f));
    float clean_p = fmaxf(soft_diff_p, floor_p);

    float gain = powf(clean_p / (s_p + SPECTRAL_EPSILON), 1.0f / p);
    if (gain > 1.0f) gain = 1.0f;
    if (gain < 0.0f) gain = 0.0f;
    gain_spectrum[k] = gain;
  }
}

static void apply_spatial_gain_smoothing(float* gain_spectrum, uint32_t size) {
  if (size < 3U) return;
  float prev = gain_spectrum[0];
  for (uint32_t k = 1U; k < size - 1U; k++) {
    float curr = gain_spectrum[k];
    float next = gain_spectrum[k + 1];
    float smoothed = (0.25f * prev) + (0.50f * curr) + (0.25f * next);
    prev = curr;
    gain_spectrum[k] = smoothed;
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
      generalized_subtraction(real_spectrum_size, spectrum, noise_spectrum, alpha, beta, gain_spectrum);
      break;
    default:
      wiener_subtraction(real_spectrum_size, spectrum, noise_spectrum, alpha, beta, gain_spectrum);
      break;
  }

  apply_spatial_gain_smoothing(gain_spectrum, real_spectrum_size);
  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}
