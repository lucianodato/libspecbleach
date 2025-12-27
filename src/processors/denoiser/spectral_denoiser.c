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

#include "spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/gain_estimation/gain_estimators.h"
#include "shared/noise_estimation/noise_estimator.h"
#include "shared/post_estimation/noise_floor_manager.h"
#include "shared/post_estimation/postfilter.h"
#include "shared/pre_estimation/critical_bands.h"
#include "shared/pre_estimation/noise_scaling_criterias.h"
#include "shared/pre_estimation/spectral_smoother.h"
#include "shared/utils/denoise_mixer.h"
#include "shared/utils/spectral_features.h"
#include "shared/utils/spectral_utils.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct SbSpectralDenoiser {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t hop;
  float default_oversubtraction;
  float default_undersubtraction;

  float* gain_spectrum;
  float* alpha;
  float* beta;
  float* noise_spectrum;

  SpectrumType spectrum_type;
  CriticalBandType band_type;
  DenoiserParameters denoise_parameters;
  GainEstimationType gain_estimation_type;
  TimeSmoothingType time_smoothing_type;
  NoiseEstimatorType noise_estimator_type;

  NoiseEstimator* noise_estimator;
  PostFilter* postfiltering;
  NoiseProfile* noise_profile;
  SpectralFeatures* spectral_features;
  DenoiseMixer* mixer;
  NoiseScalingCriterias* noise_scaling_criteria;
  SpectralSmoother* spectrum_smoothing;
  NoiseFloorManager* noise_floor_manager;
  bool postfiltering_enabled;
  bool whitening_enabled;
} SbSpectralDenoiser;

SpectralProcessorHandle spectral_denoiser_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, NoiseProfile* noise_profile) {

  if (!noise_profile || sample_rate == 0 || fft_size == 0 ||
      overlap_factor == 0) {
    return NULL;
  }

  SbSpectralDenoiser* self =
      (SbSpectralDenoiser*)calloc(1U, sizeof(SbSpectralDenoiser));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->hop = self->fft_size / overlap_factor;
  self->sample_rate = sample_rate;
  self->spectrum_type = SPECTRAL_TYPE_GENERAL;
  self->band_type = CRITICAL_BANDS_TYPE;
  self->default_oversubtraction = DEFAULT_OVERSUBTRACTION;
  self->default_undersubtraction = DEFAULT_UNDERSUBTRACTION;
  self->gain_estimation_type = GAIN_ESTIMATION_TYPE;
  self->time_smoothing_type = TIME_SMOOTHING_TYPE;
  self->postfiltering_enabled = POSTFILTER_ENABLED_GENERAL;
  self->whitening_enabled = WHITENING_ENABLED_GENERAL;

  self->gain_spectrum = (float*)calloc(self->fft_size, sizeof(float));
  if (!self->gain_spectrum) {
    spectral_denoiser_free(self);
    return NULL;
  }
  (void)initialize_spectrum_with_value(self->gain_spectrum, self->fft_size,
                                       1.F);

  self->alpha = (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->alpha) {
    spectral_denoiser_free(self);
    return NULL;
  }
  (void)initialize_spectrum_with_value(self->alpha, self->real_spectrum_size,
                                       1.F);

  self->beta = (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->beta) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->noise_profile = noise_profile;
  self->noise_spectrum =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->noise_spectrum) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->noise_estimator =
      noise_estimation_initialize(self->fft_size, noise_profile);
  if (!self->noise_estimator) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->spectral_features =
      spectral_features_initialize(self->real_spectrum_size);
  if (!self->spectral_features) {
    spectral_denoiser_free(self);
    return NULL;
  }

  if (self->postfiltering_enabled) {
    self->postfiltering = postfilter_initialize(self->fft_size);
    if (!self->postfiltering) {
      spectral_denoiser_free(self);
      return NULL;
    }
  }

  self->spectrum_smoothing =
      spectral_smoothing_initialize(self->fft_size, self->time_smoothing_type);
  if (!self->spectrum_smoothing) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->noise_scaling_criteria = noise_scaling_criterias_initialize(
      self->fft_size, self->band_type, self->sample_rate, self->spectrum_type);
  if (!self->noise_scaling_criteria) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->mixer =
      denoise_mixer_initialize(self->fft_size, self->sample_rate, self->hop);
  if (!self->mixer) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->noise_floor_manager = noise_floor_manager_initialize(
      self->fft_size, self->sample_rate, self->hop);
  if (!self->noise_floor_manager) {
    spectral_denoiser_free(self);
    return NULL;
  }

  return self;
}

void spectral_denoiser_free(SpectralProcessorHandle instance) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  if (!self) {
    return;
  }

  // Don't free noise profile used as reference here

  if (self->noise_estimator) {
    noise_estimation_free(self->noise_estimator);
  }
  if (self->spectral_features) {
    spectral_features_free(self->spectral_features);
  }
  if (self->spectrum_smoothing) {
    spectral_smoothing_free(self->spectrum_smoothing);
  }
  if (self->noise_scaling_criteria) {
    noise_scaling_criterias_free(self->noise_scaling_criteria);
  }
  if (self->postfiltering) {
    postfilter_free(self->postfiltering);
  }
  if (self->mixer) {
    denoise_mixer_free(self->mixer);
  }
  if (self->gain_spectrum) {
    free(self->gain_spectrum);
  }

  if (self->noise_floor_manager) {
    noise_floor_manager_free(self->noise_floor_manager);
  }
  if (self->alpha) {
    free(self->alpha);
  }
  if (self->beta) {
    free(self->beta);
  }
  if (self->noise_spectrum) {
    free(self->noise_spectrum);
  }

  free(self);
}

bool load_reduction_parameters(SpectralProcessorHandle instance,
                               DenoiserParameters parameters) {
  if (!instance) {
    return false;
  }

  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;
  self->denoise_parameters = parameters;

  return true;
}

bool spectral_denoiser_run(SpectralProcessorHandle instance,
                           float* fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  float* reference_spectrum =
      get_spectral_feature(self->spectral_features, fft_spectrum,
                           self->fft_size, self->spectrum_type);

  if (self->denoise_parameters.learn_noise > 0) {
    // Learn all modes simultaneously
    for (int mode = ROLLING_MEAN; mode <= MAX; mode++) {
      noise_estimation_run(self->noise_estimator, (NoiseEstimatorType)mode,
                           reference_spectrum);
    }
  } else if (is_noise_estimation_available(
                 self->noise_profile,
                 self->denoise_parameters.noise_reduction_mode)) {
    memcpy(self->noise_spectrum,
           get_noise_profile(self->noise_profile,
                             self->denoise_parameters.noise_reduction_mode),
           self->real_spectrum_size * sizeof(float));

    NoiseScalingParameters oversubtraction_parameters =
        (NoiseScalingParameters){
            .oversubtraction = (self->default_oversubtraction +
                                self->denoise_parameters.noise_rescale),
            .undersubtraction = self->denoise_parameters.reduction_amount,
            .scaling_type = self->denoise_parameters.noise_scaling_type,
        };

    float whitening_factor = self->whitening_enabled
                                 ? self->denoise_parameters.whitening_factor
                                 : 0.0f;

    apply_noise_scaling_criteria(
        self->noise_scaling_criteria, reference_spectrum, self->noise_spectrum,
        self->alpha, self->beta, oversubtraction_parameters);

    TimeSmoothingParameters spectral_smoothing_parameters =
        (TimeSmoothingParameters){
            .smoothing = self->denoise_parameters.smoothing_factor,
        };
    spectral_smoothing_run(self->spectrum_smoothing,
                           spectral_smoothing_parameters, reference_spectrum);

    estimate_gains(self->real_spectrum_size, self->fft_size, reference_spectrum,
                   self->noise_spectrum, self->gain_spectrum, self->alpha,
                   self->beta, self->gain_estimation_type);

    noise_floor_manager_apply(
        self->noise_floor_manager, self->real_spectrum_size, self->fft_size,
        self->gain_spectrum, self->noise_spectrum,
        self->denoise_parameters.reduction_amount, whitening_factor);

    if (self->postfiltering_enabled) {
      PostFiltersParameters post_filter_parameters = (PostFiltersParameters){
          .snr_threshold = self->denoise_parameters.post_filter_threshold,
          .gain_floor = self->denoise_parameters.reduction_amount,
      };
      postfilter_apply(self->postfiltering, fft_spectrum, self->gain_spectrum,
                       post_filter_parameters);
    }

    DenoiseMixerParameters mixer_parameters = (DenoiseMixerParameters){
        .noise_level = self->denoise_parameters.reduction_amount,
        .residual_listen = self->denoise_parameters.residual_listen,
        .whitening_amount = whitening_factor,
    };

    denoise_mixer_run(self->mixer, fft_spectrum, self->gain_spectrum,
                      mixer_parameters);
  }

  return true;
}
