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

#include "processors/time_denoiser/dsaf_engine.h"
#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"
#include "specbleach_time_denoiser.h"
#include <stdlib.h>
#include <string.h>

struct specbleach_time_denoiser {
  SbDsafEngine* engine;
};

specbleach_time_denoiser* specbleach_time_denoiser_initialize(
    const uint32_t sample_rate) {
  specbleach_time_denoiser* self =
      (specbleach_time_denoiser*)calloc(1U, sizeof(specbleach_time_denoiser));
  if (!self) {
    return NULL;
  }

  self->engine = dsaf_engine_initialize(sample_rate);
  if (!self->engine) {
    free(self);
    return NULL;
  }

  return self;
}

void specbleach_time_denoiser_free(specbleach_time_denoiser* instance) {
  if (!instance) {
    return;
  }

  dsaf_engine_free(instance->engine);
  free(instance);
}

bool specbleach_time_denoiser_load_parameters(
    specbleach_time_denoiser* instance,
    const SpecbleachTimeDenoiserParameters* parameters,
    const uint32_t parameters_size) {
  if (!instance || !parameters ||
      parameters_size != sizeof(SpecbleachTimeDenoiserParameters)) {
    return false;
  }

  DsafParameters engine_parameters;
  engine_parameters.reduction_gain = parameters->reduction_gain;
  engine_parameters.smoothing_factor = parameters->smoothing_factor;
  engine_parameters.adaptive_noise = parameters->adaptive_noise;
  engine_parameters.noise_estimation_method =
      (AdaptiveNoiseEstimationMethod)parameters->noise_estimation_method;
  engine_parameters.suppression_strength = parameters->suppression_strength;

  dsaf_engine_load_parameters(instance->engine, &engine_parameters);
  return true;
}

bool specbleach_time_denoiser_process(specbleach_time_denoiser* instance,
                                      const uint32_t number_of_samples,
                                      const float* input, float* output) {
  if (!instance) {
    return false;
  }

  return dsaf_engine_process(instance->engine, number_of_samples, input,
                             output);
}

uint32_t specbleach_time_denoiser_get_latency(
    specbleach_time_denoiser* instance) {
  if (!instance) {
    return 0U;
  }

  return dsaf_engine_get_latency(instance->engine);
}
