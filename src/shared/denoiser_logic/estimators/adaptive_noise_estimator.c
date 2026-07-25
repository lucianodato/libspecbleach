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

/* removed old path include */
#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/estimators/brandt_noise_estimator.h"
#include "shared/denoiser_logic/estimators/martin_noise_estimator.h"
#include "shared/denoiser_logic/estimators/spp_mmse_noise_estimator.h"
#include <stdlib.h>

typedef struct AdaptiveEstimatorOps {
  bool (*run)(void* inst, const float* spectrum, float* noise);
  void (*set_state)(void* inst, const float* initial);
  void (*apply_floor)(void* inst, const float* floor);
  void (*update_seed)(void* inst, const float* seed);
  void (*free)(void* inst);
} AdaptiveEstimatorOps;

struct AdaptiveNoiseEstimator {
  AdaptiveNoiseEstimationMethod method;
  uint32_t spectrum_size;
  void* internal_estimator;
  AdaptiveEstimatorOps ops;
};

static bool spp_mmse_run_wrapper(void* inst, const float* s, float* n) {
  return spp_mmse_noise_estimator_run((SppMmseNoiseEstimator*)inst, s, n);
}
static void spp_mmse_set_state_wrapper(void* inst, const float* p) {
  spp_mmse_noise_estimator_set_state((SppMmseNoiseEstimator*)inst, p);
}
static void spp_mmse_apply_floor_wrapper(void* inst, const float* f) {
  spp_mmse_noise_estimator_apply_floor((SppMmseNoiseEstimator*)inst, f);
}
static void spp_mmse_update_seed_wrapper(void* inst, const float* s) {
  spp_mmse_noise_estimator_update_seed((SppMmseNoiseEstimator*)inst, s);
}
static void spp_mmse_free_wrapper(void* inst) {
  spp_mmse_noise_estimator_free((SppMmseNoiseEstimator*)inst);
}

static bool brandt_run_wrapper(void* inst, const float* s, float* n) {
  return brandt_noise_estimator_run((BrandtNoiseEstimator*)inst, s, n);
}
static void brandt_set_state_wrapper(void* inst, const float* p) {
  brandt_noise_estimator_set_state((BrandtNoiseEstimator*)inst, p);
}
static void brandt_apply_floor_wrapper(void* inst, const float* f) {
  brandt_noise_estimator_apply_floor((BrandtNoiseEstimator*)inst, f);
}
static void brandt_update_seed_wrapper(void* inst, const float* s) {
  brandt_noise_estimator_update_seed((BrandtNoiseEstimator*)inst, s);
}
static void brandt_free_wrapper(void* inst) {
  brandt_noise_estimator_free((BrandtNoiseEstimator*)inst);
}

static bool martin_run_wrapper(void* inst, const float* s, float* n) {
  return martin_noise_estimator_run((MartinNoiseEstimator*)inst, s, n);
}
static void martin_set_state_wrapper(void* inst, const float* p) {
  martin_noise_estimator_set_state((MartinNoiseEstimator*)inst, p);
}
static void martin_apply_floor_wrapper(void* inst, const float* f) {
  martin_noise_estimator_apply_floor((MartinNoiseEstimator*)inst, f);
}
static void martin_update_seed_wrapper(void* inst, const float* s) {
  martin_noise_estimator_update_seed((MartinNoiseEstimator*)inst, s);
}
static void martin_free_wrapper(void* inst) {
  martin_noise_estimator_free((MartinNoiseEstimator*)inst);
}

AdaptiveNoiseEstimator* adaptive_estimator_initialize(
    uint32_t noise_spectrum_size, uint32_t sample_rate, uint32_t fft_size,
    AdaptiveNoiseEstimationMethod method) {
  AdaptiveNoiseEstimator* self =
      (AdaptiveNoiseEstimator*)calloc(1U, sizeof(AdaptiveNoiseEstimator));
  if (!self) {
    return NULL;
  }

  self->method = method;
  self->spectrum_size = noise_spectrum_size;

  if (method == SPP_MMSE_METHOD) {
    self->internal_estimator = spp_mmse_noise_estimator_initialize(
        noise_spectrum_size, sample_rate, fft_size);
    self->ops = (AdaptiveEstimatorOps){
        .run = spp_mmse_run_wrapper,
        .set_state = spp_mmse_set_state_wrapper,
        .apply_floor = spp_mmse_apply_floor_wrapper,
        .update_seed = spp_mmse_update_seed_wrapper,
        .free = spp_mmse_free_wrapper,
    };
  } else if (method == BRANDT_METHOD) {
    self->internal_estimator = brandt_noise_estimator_initialize(
        noise_spectrum_size, BRANDT_DEFAULT_HISTORY_MS, sample_rate, fft_size);
    self->ops = (AdaptiveEstimatorOps){
        .run = brandt_run_wrapper,
        .set_state = brandt_set_state_wrapper,
        .apply_floor = brandt_apply_floor_wrapper,
        .update_seed = brandt_update_seed_wrapper,
        .free = brandt_free_wrapper,
    };
  } else {
    self->method = MARTIN_METHOD;
    self->internal_estimator = martin_noise_estimator_initialize(
        noise_spectrum_size, sample_rate, fft_size);
    self->ops = (AdaptiveEstimatorOps){
        .run = martin_run_wrapper,
        .set_state = martin_set_state_wrapper,
        .apply_floor = martin_apply_floor_wrapper,
        .update_seed = martin_update_seed_wrapper,
        .free = martin_free_wrapper,
    };
  }

  if (!self->internal_estimator) {
    free(self);
    return NULL;
  }

  return self;
}

AdaptiveNoiseEstimationMethod adaptive_estimator_get_method(
    const AdaptiveNoiseEstimator* self) {
  return self ? self->method : MARTIN_METHOD;
}

void adaptive_estimator_free(AdaptiveNoiseEstimator* self) {
  if (!self) {
    return;
  }
  if (self->internal_estimator && self->ops.free) {
    self->ops.free(self->internal_estimator);
  }
  free(self);
}

bool adaptive_estimator_run(AdaptiveNoiseEstimator* self, const float* spectrum,
                            float* noise_spectrum, float* aggressiveness,
                            float param_aggressiveness) {
  if (!self || !self->internal_estimator || !self->ops.run) {
    return false;
  }

  bool success = self->ops.run(self->internal_estimator, spectrum, noise_spectrum);

  if (success) {
    if (aggressiveness) {
      *aggressiveness = param_aggressiveness;
    }
    interpolate_spectrum_gaps(noise_spectrum, self->spectrum_size,
                              NOISE_ESTIMATION_INTERPOLATION_THRESHOLD);
    smooth_spectrum(noise_spectrum, self->spectrum_size,
                    ADAPTIVE_NOISE_FLOOR_SMOOTHING);
  }

  return success;
}

void adaptive_estimator_set_state(AdaptiveNoiseEstimator* self,
                                  const float* initial_profile,
                                  int estimation_method) {
  (void)estimation_method;
  if (self && self->internal_estimator && self->ops.set_state) {
    self->ops.set_state(self->internal_estimator, initial_profile);
  }
}

void adaptive_estimator_apply_floor(AdaptiveNoiseEstimator* self,
                                    const float* floor_profile) {
  if (self && self->internal_estimator && self->ops.apply_floor) {
    self->ops.apply_floor(self->internal_estimator, floor_profile);
  }
}

void adaptive_estimator_update_seed(AdaptiveNoiseEstimator* self,
                                    const float* seed_profile) {
  if (self && self->internal_estimator && self->ops.update_seed) {
    self->ops.update_seed(self->internal_estimator, seed_profile);
  }
}
