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

#include "adaptive_noise_estimator.h"
#include "louizou_noise_estimator.h"
#include "spp_mmse_noise_estimator.h"
#include <stdlib.h>

struct AdaptiveNoiseEstimator {
  AdaptiveNoiseEstimationMethod method;
  void* internal_estimator;
};

AdaptiveNoiseEstimator* louizou_estimator_initialize(
    uint32_t noise_spectrum_size, uint32_t sample_rate, uint32_t fft_size) {
  AdaptiveNoiseEstimator* self =
      (AdaptiveNoiseEstimator*)calloc(1U, sizeof(AdaptiveNoiseEstimator));
  if (!self) {
    return NULL;
  }

  self->method = LOUIZOU_METHOD;
  self->internal_estimator = louizou_noise_estimator_initialize(
      noise_spectrum_size, sample_rate, fft_size);

  if (!self->internal_estimator) {
    free(self);
    return NULL;
  }

  return self;
}

AdaptiveNoiseEstimator* spp_mmse_estimator_initialize(
    uint32_t noise_spectrum_size, uint32_t sample_rate, uint32_t fft_size) {
  AdaptiveNoiseEstimator* self =
      (AdaptiveNoiseEstimator*)calloc(1U, sizeof(AdaptiveNoiseEstimator));
  if (!self) {
    return NULL;
  }

  self->method = SPP_MMSE_METHOD;
  self->internal_estimator = spp_mmse_noise_estimator_initialize(
      noise_spectrum_size, sample_rate, fft_size);

  if (!self->internal_estimator) {
    free(self);
    return NULL;
  }

  return self;
}

void spp_mmse_estimator_free(AdaptiveNoiseEstimator* self) {
  louizou_estimator_free(self); // Dispatcher handles both
}

void louizou_estimator_free(AdaptiveNoiseEstimator* self) {
  if (!self) {
    return;
  }

  if (self->method == LOUIZOU_METHOD) {
    louizou_noise_estimator_free(
        (LouizouNoiseEstimator*)self->internal_estimator);
  } else {
    spp_mmse_noise_estimator_free(
        (SppMmseNoiseEstimator*)self->internal_estimator);
  }

  free(self);
}

bool louizou_estimator_run(AdaptiveNoiseEstimator* self, const float* spectrum,
                           float* noise_spectrum) {
  if (!self || self->method != LOUIZOU_METHOD) {
    return false;
  }
  return louizou_noise_estimator_run(
      (LouizouNoiseEstimator*)self->internal_estimator, spectrum,
      noise_spectrum);
}

bool spp_mmse_estimator_run(AdaptiveNoiseEstimator* self, const float* spectrum,
                            float* noise_spectrum) {
  if (!self || self->method != SPP_MMSE_METHOD) {
    return false;
  }
  return spp_mmse_noise_estimator_run(
      (SppMmseNoiseEstimator*)self->internal_estimator, spectrum,
      noise_spectrum);
}

void adaptive_estimator_set_state(AdaptiveNoiseEstimator* self,
                                  const float* initial_profile,
                                  int estimation_method) {
  (void)estimation_method;
  if (!self) {
    return;
  }

  if (self->method == LOUIZOU_METHOD) {
    louizou_noise_estimator_set_state(
        (LouizouNoiseEstimator*)self->internal_estimator, initial_profile);
  } else {
    spp_mmse_noise_estimator_set_state(
        (SppMmseNoiseEstimator*)self->internal_estimator, initial_profile);
  }
}

void adaptive_estimator_apply_floor(AdaptiveNoiseEstimator* self,
                                    const float* floor_profile) {
  if (!self) {
    return;
  }

  if (self->method == LOUIZOU_METHOD) {
    louizou_noise_estimator_apply_floor(
        (LouizouNoiseEstimator*)self->internal_estimator, floor_profile);
  } else {
    spp_mmse_noise_estimator_apply_floor(
        (SppMmseNoiseEstimator*)self->internal_estimator, floor_profile);
  }
}

void adaptive_estimator_update_seed(AdaptiveNoiseEstimator* self,
                                    const float* seed_profile) {
  if (!self) {
    return;
  }

  if (self->method == LOUIZOU_METHOD) {
    louizou_noise_estimator_update_seed(
        (LouizouNoiseEstimator*)self->internal_estimator, seed_profile);
  } else {
    spp_mmse_noise_estimator_update_seed(
        (SppMmseNoiseEstimator*)self->internal_estimator, seed_profile);
  }
}
