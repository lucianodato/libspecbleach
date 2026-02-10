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
#include "shared/gain_estimation/suppression_engine.h"
#include "shared/noise_estimation/adaptive_noise_estimator.h"
#include "shared/noise_estimation/noise_estimator.h"
#include "shared/noise_estimation/tonal_reducer.h"
#include "shared/post_estimation/masking_veto.h"
#include "shared/post_estimation/noise_floor_manager.h"
#include "shared/pre_estimation/critical_bands.h"
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
  float* manual_noise_floor;
  float* noisy_reference;

  TonalReducer* tonal_reducer;

  SpectrumType spectrum_type;
  CriticalBandType band_type;
  DenoiserParameters denoise_parameters;
  GainEstimationType gain_estimation_type;
  TimeSmoothingType time_smoothing_type;
  NoiseEstimatorType noise_estimator_type;

  NoiseEstimator* noise_estimator;
  AdaptiveNoiseEstimator* adaptive_estimator;
  NoiseProfile* noise_profile;
  SpectralFeatures* spectral_features;
  DenoiseMixer* mixer;
  SpectralSmoother* spectrum_smoothing;
  NoiseFloorManager* noise_floor_manager;
  MaskingVeto* masking_veto;
  SuppressionEngine* suppression_engine;

  int last_adaptive_state;
  int last_noise_estimation_method;
  float aggressiveness;
  bool was_learning;
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
  self->real_spectrum_size = (self->fft_size / 2U) + 1U;
  self->hop = self->fft_size / overlap_factor;
  self->sample_rate = sample_rate;
  self->spectrum_type = SPECTRAL_TYPE_GENERAL;
  self->band_type = CRITICAL_BANDS_TYPE;
  self->default_oversubtraction = DEFAULT_OVERSUBTRACTION;
  self->default_undersubtraction = DEFAULT_UNDERSUBTRACTION;
  self->gain_estimation_type = GAIN_ESTIMATION_TYPE;
  self->time_smoothing_type = TIME_SMOOTHING_TYPE;

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

  self->manual_noise_floor =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noisy_reference =
      (float*)calloc(self->real_spectrum_size, sizeof(float));

  if (!self->manual_noise_floor || !self->noisy_reference) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->tonal_reducer = tonal_reducer_initialize(
      self->real_spectrum_size, self->sample_rate, self->fft_size);
  if (!self->tonal_reducer) {
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

  self->spectrum_smoothing =
      spectral_smoothing_initialize(self->fft_size, self->time_smoothing_type);
  if (!self->spectrum_smoothing) {
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

  self->was_learning = false;
  self->aggressiveness = 0.0f;
  self->denoise_parameters.tonal_reduction = 0.0f;

  self->masking_veto = masking_veto_initialize(
      self->fft_size, self->sample_rate, self->spectrum_type);
  self->suppression_engine =
      suppression_engine_initialize(self->real_spectrum_size);

  if (!self->noise_floor_manager || !self->masking_veto ||
      !self->suppression_engine) {
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
  if (self->adaptive_estimator) {
    adaptive_estimator_free(self->adaptive_estimator);
  }
  if (self->spectral_features) {
    spectral_features_free(self->spectral_features);
  }
  if (self->spectrum_smoothing) {
    spectral_smoothing_free(self->spectrum_smoothing);
  }
  if (self->masking_veto) {
    masking_veto_free(self->masking_veto);
  }
  if (self->suppression_engine) {
    suppression_engine_free(self->suppression_engine);
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
  if (self->manual_noise_floor) {
    free(self->manual_noise_floor);
  }
  if (self->noisy_reference) {
    free(self->noisy_reference);
  }
  if (self->tonal_reducer) {
    tonal_reducer_free(self->tonal_reducer);
  }

  free(self);
}

bool load_reduction_parameters(SpectralProcessorHandle instance,
                               DenoiserParameters parameters) {
  if (!instance) {
    return false;
  }

  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  // Check if we need to initialize or re-initialize the adaptive estimator
  if (parameters.adaptive_noise) {
    AdaptiveNoiseEstimationMethod requested_method =
        (AdaptiveNoiseEstimationMethod)parameters.noise_estimation_method;

    bool needs_init = !self->adaptive_estimator ||
                      adaptive_estimator_get_method(self->adaptive_estimator) !=
                          requested_method;

    if (needs_init) {
      adaptive_estimator_free(self->adaptive_estimator);
      self->adaptive_estimator = adaptive_estimator_initialize(
          self->real_spectrum_size, self->sample_rate, self->fft_size,
          requested_method);
      self->last_adaptive_state = 0;
    }
  }

  self->denoise_parameters = parameters;
  self->aggressiveness = parameters.aggressiveness;

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
    for (int mode = ROLLING_MEAN; mode <= MINIMUM; mode++) {
      noise_estimation_run(self->noise_estimator, (NoiseEstimatorType)mode,
                           reference_spectrum);
    }
    self->was_learning = true;
    return true;
  }

  if (self->was_learning) {
    // User just stopped learning -> Finalize all captures
    for (int mode = ROLLING_MEAN; mode <= MINIMUM; mode++) {
      noise_estimation_finalize(self->noise_estimator,
                                (NoiseEstimatorType)mode);
    }
    self->was_learning = false;
  }

  // --- Denoising Path ---

  if (self->denoise_parameters.adaptive_noise && self->adaptive_estimator) {
    // ... (Adaptive logic remains similar but uses morphed profile as base)
    // Check for state transitions
    bool state_changed = !self->last_adaptive_state;
    bool mode_changed = fabsf(self->aggressiveness -
                              self->denoise_parameters.aggressiveness) > 0.01f;

    if (state_changed || mode_changed) {
      // Calculate morphed base profile
      get_morphed_profile(self->manual_noise_floor,
                          get_noise_profile(self->noise_profile, ROLLING_MEAN),
                          get_noise_profile(self->noise_profile, MEDIAN),
                          get_noise_profile(self->noise_profile, MAX),
                          get_noise_profile(self->noise_profile, MINIMUM),
                          self->real_spectrum_size,
                          self->denoise_parameters.aggressiveness);

      adaptive_estimator_update_seed(self->adaptive_estimator,
                                     self->manual_noise_floor);

      self->last_adaptive_state = 1;
      self->aggressiveness = self->denoise_parameters.aggressiveness;
    }

    // Run adaptive estimator
    adaptive_estimator_run(self->adaptive_estimator, reference_spectrum,
                           self->noise_spectrum);

    // Apply morphed profile as a floor
    adaptive_estimator_apply_floor(self->adaptive_estimator,
                                   self->manual_noise_floor);
    for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
      if (self->noise_spectrum[k] < self->manual_noise_floor[k]) {
        self->noise_spectrum[k] = self->manual_noise_floor[k];
      }
    }

    // Smooth the morphed/refined floor every frame to eliminate
    // residual musical noise in specific steering modes (e.g. Median)
    smooth_spectrum(self->noise_spectrum, self->real_spectrum_size, 0.5f);
  } else {
    // Manual Denoising Mode
    self->last_adaptive_state = 0;

    // Use morphed profile
    get_morphed_profile(self->noise_spectrum,
                        get_noise_profile(self->noise_profile, ROLLING_MEAN),
                        get_noise_profile(self->noise_profile, MEDIAN),
                        get_noise_profile(self->noise_profile, MAX),
                        get_noise_profile(self->noise_profile, MINIMUM),
                        self->real_spectrum_size,
                        self->denoise_parameters.aggressiveness);
  }

  // --- Common Processing Path ---

  // 2. Calculate SNR-dependent oversubtraction factors (Alpha/Beta)
  // The SuppressionEngine implement Berouti-style scaling.
  // The Veto engine will later moderate this based on audibility/transients.
  suppression_engine_calculate(
      self->suppression_engine, reference_spectrum, self->noise_spectrum,
      self->denoise_parameters.suppression_strength, self->alpha, self->beta);

  // 3. Detect tonal components and boost alpha at tonal bins.
  // Runs before the masking veto so the veto can protect signal harmonics.
  tonal_reducer_run(self->tonal_reducer, self->noise_spectrum,
                    get_noise_profile(self->noise_profile, MAX),
                    get_noise_profile(self->noise_profile, MEDIAN), self->alpha,
                    self->denoise_parameters.tonal_reduction);

  // Preserve 'noisy' reference before temporal smoothing for Veto comparison
  memcpy(self->noisy_reference, reference_spectrum,
         self->real_spectrum_size * sizeof(float));

  TimeSmoothingParameters spectral_smoothing_parameters =
      (TimeSmoothingParameters){
          .smoothing = self->denoise_parameters.smoothing_factor,
      };
  spectral_smoothing_run(self->spectrum_smoothing,
                         spectral_smoothing_parameters, reference_spectrum);

  // Apply Structural Veto to rescue transients and moderate artifacts
  masking_veto_apply(self->masking_veto, reference_spectrum,
                     self->noisy_reference, self->noise_spectrum, self->alpha,
                     1.0F, self->denoise_parameters.masking_depth,
                     self->denoise_parameters.masking_elasticity);

  estimate_gains(self->real_spectrum_size, self->fft_size, reference_spectrum,
                 self->noise_spectrum, self->gain_spectrum, self->alpha,
                 self->beta, self->gain_estimation_type);

  // Apply noise floor management
  noise_floor_manager_apply(self->noise_floor_manager, self->real_spectrum_size,
                            self->fft_size, self->gain_spectrum,
                            self->noise_spectrum,
                            self->denoise_parameters.reduction_amount,
                            self->denoise_parameters.tonal_reduction,
                            tonal_reducer_get_mask(self->tonal_reducer),
                            self->denoise_parameters.whitening_factor);

  DenoiseMixerParameters mixer_parameters = (DenoiseMixerParameters){
      .noise_level = self->denoise_parameters.reduction_amount,
      .residual_listen = self->denoise_parameters.residual_listen,
      .whitening_amount = self->denoise_parameters.whitening_factor,
  };

  denoise_mixer_run(self->mixer, fft_spectrum, self->gain_spectrum,
                    mixer_parameters);

  return true;
}
