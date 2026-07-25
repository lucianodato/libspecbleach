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

#include "general_utils.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

float sanitize_denormal(float value) {
  if (!isnormal(value)) {
    value = 0.F;
  }
  return value;
}

float from_db_to_coefficient(const float value_db) {
  return powf(10.0F, value_db / 20.0F);
}

float remap_percentage_log_like_unity(const float value) {
  return 1.F - expf(-3.F * (value));
}

int get_next_divisible_two(int number) {
  return number + (number & 1);
}

int get_next_power_two(int number) {
  if (number <= 0) {
    return 1;
  }
  number--;
  number |= number >> 1;
  number |= number >> 2;
  number |= number >> 4;
  number |= number >> 8;
  number |= number >> 16;
  number++;
  return number;
}

DenoiserParameters sb_denoiser_params_sanitize(
    SpectralBleachDenoiserParameters parameters) {
  return (DenoiserParameters){
      .learn_noise = parameters.learn_noise,
      .residual_listen = parameters.residual_listen,
      .reduction_amount =
          from_db_to_coefficient(parameters.reduction_amount * -1.F),
      .smoothing_factor =
          remap_percentage_log_like_unity(parameters.smoothing_factor / 100.F),
      .whitening_factor = parameters.whitening_factor / 100.F,
      .adaptive_noise = parameters.adaptive_noise,
      .noise_estimation_method = parameters.noise_estimation_method,
      .masking_depth = parameters.masking_depth,
      .suppression_strength = parameters.suppression_strength / 100.F,
      .aggressiveness = parameters.aggressiveness,
      .tonal_reduction =
          from_db_to_coefficient(parameters.tonal_reduction * -1.F),
  };
}

Denoiser2DParameters sb_denoiser_2d_params_sanitize(
    SpectralBleach2DDenoiserParameters parameters) {
  return (Denoiser2DParameters){
      .learn_noise = parameters.learn_noise,
      .residual_listen = parameters.residual_listen,
      .reduction_amount =
          from_db_to_coefficient(parameters.reduction_amount * -1.F),
      .smoothing_factor = parameters.smoothing_factor,
      .whitening_factor = parameters.whitening_factor / 100.F,
      .adaptive_noise = parameters.adaptive_noise,
      .noise_estimation_method = parameters.noise_estimation_method,
      .nlm_masking_protection = parameters.nlm_masking_protection,
      .suppression_strength = parameters.suppression_strength / 100.F,
      .aggressiveness = parameters.aggressiveness,
      .tonal_reduction =
          from_db_to_coefficient(parameters.tonal_reduction * -1.F),
  };
}
