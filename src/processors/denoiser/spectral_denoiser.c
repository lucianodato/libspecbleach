#include "spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/denoiser_post_process.h"
#include "shared/denoiser_logic/core/denoiser_profile_core.h"
#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/denoiser_logic/processing/gain_calculator.h"
#include "shared/denoiser_logic/processing/hpss_filter.h"
#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/denoiser_logic/processing/suppression_engine.h"
#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/stft/stft_processor.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/spectral_circular_buffer.h"
#include "shared/utils/spectral_features.h"
#include "shared/utils/spectral_smoother.h"
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

  float* gain_spectrum;
  float* alpha;
  float* beta;
  float* noise_spectrum;
  float* manual_noise_floor;

  TonalReducer* tonal_reducer;

  SpectrumType spectrum_type;
  CriticalBandType band_type;
  DenoiserParameters denoise_parameters;
  GainCalculationType gain_calculation_type;
  TimeSmoothingType time_smoothing_type;
  NoiseEstimatorType noise_estimator_type;

  NoiseEstimator* noise_estimator;
  AdaptiveNoiseEstimator* adaptive_estimator;
  NoiseProfile* noise_profile;
  SpectralFeatures* spectral_features;
  SpectralSmoother* spectrum_smoothing;
  NoiseFloorManager* noise_floor_manager;
  MaskingVeto* masking_veto;
  SuppressionEngine* suppression_engine;
  HpssFilter* hpss_filter;

  SbSpectralCircularBuffer* circular_buffer;
  uint32_t layer_fft;
  uint32_t layer_noise;

  float* delayed_magnitude;
  float* mask_harmonic;
  float* mask_percussive;
  float* g_h;

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
  self->spectrum_type = SPECTRAL_TYPE_1D;
  self->band_type = CRITICAL_BANDS_TYPE_1D;
  self->gain_calculation_type = GAIN_ESTIMATION_TYPE_1D;
  self->time_smoothing_type = FIXED;

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
  if (!self->manual_noise_floor) {
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
      spectral_smoothing_initialize(self->fft_size, self->sample_rate,
                                    overlap_factor, self->time_smoothing_type);
  if (!self->spectrum_smoothing) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->noise_floor_manager = noise_floor_manager_initialize(self->fft_size);

  self->was_learning = false;
  self->aggressiveness = 0.0f;
  self->denoise_parameters.tonal_reduction = 0.0f;

  self->masking_veto = masking_veto_initialize(
      self->fft_size, self->sample_rate, CRITICAL_BANDS_TYPE_1D,
      self->spectrum_type, false, USE_TEMPORAL_MASKING_1D_DEFAULT);
  self->suppression_engine = suppression_engine_initialize(
      self->real_spectrum_size, self->sample_rate, self->band_type,
      self->spectrum_type, true, USE_TEMPORAL_MASKING_1D_DEFAULT);

  self->circular_buffer = spectral_circular_buffer_create(DELAY_BUFFER_FRAMES);
  if (!self->circular_buffer) {
    spectral_denoiser_free(self);
    return NULL;
  }

  self->layer_fft =
      spectral_circular_buffer_add_layer(self->circular_buffer, self->fft_size);
  self->layer_noise = spectral_circular_buffer_add_layer(
      self->circular_buffer, self->real_spectrum_size);

  if (self->layer_fft == 0xFFFFFFFFU || self->layer_noise == 0xFFFFFFFFU) {
    spectral_denoiser_free(self);
    return NULL;
  }

  HpssConfig hpss_cfg = {
      .real_spectrum_size = self->real_spectrum_size,
  };
  self->hpss_filter = hpss_filter_initialize(hpss_cfg);
  self->delayed_magnitude =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->mask_harmonic = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->mask_percussive =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->g_h = (float*)calloc(self->fft_size, sizeof(float));

  if (!self->noise_floor_manager || !self->masking_veto ||
      !self->suppression_engine || !self->hpss_filter ||
      !self->delayed_magnitude || !self->mask_harmonic ||
      !self->mask_percussive || !self->g_h) {
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

  if (self->circular_buffer) {
    spectral_circular_buffer_free(self->circular_buffer);
  }
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

  if (self->tonal_reducer) {
    tonal_reducer_free(self->tonal_reducer);
  }

  if (self->hpss_filter) {
    hpss_filter_free(self->hpss_filter);
  }
  if (self->delayed_magnitude) {
    free(self->delayed_magnitude);
  }
  if (self->mask_harmonic) {
    free(self->mask_harmonic);
  }
  if (self->mask_percussive) {
    free(self->mask_percussive);
  }
  if (self->g_h) {
    free(self->g_h);
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

  if (self->hpss_filter) {
    hpss_filter_set_quality_mode(self->hpss_filter,
                                 (HpssQualityMode)parameters.hpss_quality_mode);
    hpss_filter_set_sensitivity(self->hpss_filter, parameters.hpss_sensitivity);
  }

  return true;
}

bool spectral_denoiser_run(SpectralProcessorHandle instance,
                           float* fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;

  // 1. Preparation: Get reference spectrum and handle learning mode
  float* reference_spectrum =
      get_spectral_feature(self->spectral_features, fft_spectrum,
                           self->fft_size, self->spectrum_type);

  if (denoiser_profile_core_handle_learning_mode(
          self->noise_estimator, reference_spectrum,
          self->denoise_parameters.learn_noise, &self->was_learning)) {
    return true;
  }

  // 2. Noise Estimation: Update noise profile (Adaptive or Manual)
  DenoiserProfileCoreParams profile_params = {
      .adaptive_enabled = self->denoise_parameters.adaptive_noise,
      .spectrum_size = self->real_spectrum_size,
      .aggressiveness = &self->aggressiveness,
      .param_aggressiveness = self->denoise_parameters.aggressiveness,
      .last_adaptive_state = &self->last_adaptive_state,
      .adaptive_estimator = self->adaptive_estimator,
      .noise_profile = self->noise_profile,
      .manual_noise_floor = self->manual_noise_floor,
      .noise_spectrum = self->noise_spectrum,
      .noise_estimator = self->noise_estimator,
      .noise_profile_offset_linear =
          self->denoise_parameters.noise_profile_offset_linear,
  };
  denoiser_profile_core_update(profile_params, reference_spectrum);

  // 2.1 Align internal state and output to the delayed frame (temporal plumbing)
  const uint32_t delay_frames =
      hpss_filter_get_latency_frames(self->hpss_filter);

  // 1. Push incoming frames
  spectral_circular_buffer_push(self->circular_buffer, self->layer_fft,
                                fft_spectrum);
  spectral_circular_buffer_push(self->circular_buffer, self->layer_noise,
                                self->noise_spectrum);

  // 2. Process HPSS (internally produces masks aligned at delay_frames)
  if (!hpss_filter_process(self->hpss_filter, reference_spectrum,
                           self->noise_spectrum, self->delayed_magnitude,
                           self->mask_harmonic, self->mask_percussive)) {
    memcpy(self->delayed_magnitude, reference_spectrum,
           self->real_spectrum_size * sizeof(float));
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      self->mask_harmonic[k] = 1.0f;
      self->mask_percussive[k] = 0.0f;
    }
  }

  // 3. Retrieve aligned FFT and noise at the exact same delay
  float* delayed_spectrum = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_fft, delay_frames);
  float* delayed_noise = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->layer_noise, delay_frames);

  if (!delayed_spectrum) {
    delayed_spectrum = fft_spectrum;
  }
  if (!delayed_noise) {
    delayed_noise = self->noise_spectrum;
  }

  // Align output to delayed frame for post-processing
  if (delayed_spectrum != fft_spectrum) {
    memcpy(fft_spectrum, delayed_spectrum, self->fft_size * sizeof(float));
  }

  float* input_magnitude = self->delayed_magnitude;

  // 3.1. Calculate SNR-dependent oversubtraction factors (Alpha/Beta)
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->denoise_parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->suppression_engine, input_magnitude,
                               delayed_noise, suppression_params, self->alpha,
                               self->beta);

  // 3.2. Detect tonal components and boost alpha at tonal bins
  tonal_reducer_run(self->tonal_reducer, delayed_noise,
                    get_noise_profile(self->noise_profile, CV_MASK),
                    is_noise_estimation_available(self->noise_profile, CV_MASK),
                    self->alpha, self->denoise_parameters.tonal_reduction);

  // 3.3. Apply Structural Veto to rescue transients and moderate artifacts
  masking_veto_apply(self->masking_veto, input_magnitude, delayed_noise, NULL,
                     self->alpha, self->denoise_parameters.masking_depth);

  // 3.4. Gain Calculation & HPSS Combination
  // 1. Calculate heavily smoothed Harmonic Gain (G_H) on stationary path
  calculate_gains(self->real_spectrum_size, self->fft_size, input_magnitude,
                  delayed_noise, self->g_h, self->alpha, self->beta,
                  self->gain_calculation_type);

  TimeSmoothingParameters spectral_smoothing_parameters =
      (TimeSmoothingParameters){
          .smoothing = self->denoise_parameters.smoothing_factor,
      };
  spectral_smoothing_run(self->spectrum_smoothing,
                         spectral_smoothing_parameters, self->g_h);
  if (self->denoise_parameters.smoothing_factor > 0.0f) {
    int passes = 1 + (int)(self->denoise_parameters.smoothing_factor * 2.0f);
    for (int p = 0; p < passes; ++p) {
      spectral_smoothing_apply_spatial(self->g_h, self->real_spectrum_size);
    }
  }

  // 2. Recombine: Harmonic gets G_H suppression, Percussive gets unity pass-through
  for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
    float w_h = self->mask_harmonic[k];
    float w_p = self->mask_percussive[k];

    // G_final blends Wiener suppression on stationary content with unity gain on transient bursts
    float combined_gain = (w_h * self->g_h[k]) + w_p;
    self->gain_spectrum[k] = fminf(fmaxf(combined_gain, 0.0f), 1.0f);
  }
  sb_apply_spectral_symmetry(self->gain_spectrum, self->real_spectrum_size,
                             self->fft_size);

  // 4. Post-Processing: Final gain management and mixing
  DenoiserPostProcessParams post_params = {
      .fft_size = self->fft_size,
      .real_spectrum_size = self->real_spectrum_size,
      .reduction_amount = self->denoise_parameters.reduction_amount,
      .tonal_reduction = self->denoise_parameters.tonal_reduction,
      .whitening_factor = self->denoise_parameters.whitening_factor,
      .residual_listen = self->denoise_parameters.residual_listen,
      .noise_floor_manager = self->noise_floor_manager,
      .tonal_reducer = self->tonal_reducer,
      .gain_spectrum = self->gain_spectrum,
      .noise_spectrum = delayed_noise,
      .fft_spectrum = fft_spectrum,
      .reduction_curve_bias = self->denoise_parameters.reduction_curve_bias,
  };
  denoiser_post_process_apply(post_params);

  // Advance circular buffer write index
  spectral_circular_buffer_advance(self->circular_buffer);

  return true;
}

const float* spectral_denoiser_get_tonal_mask(
    SpectralProcessorHandle instance) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;
  return (self && self->tonal_reducer)
             ? tonal_reducer_get_mask(self->tonal_reducer)
             : NULL;
}

uint32_t spectral_denoiser_get_peaks(SpectralProcessorHandle instance,
                                     float* peak_freqs_hz, uint32_t max_peaks) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;
  return self ? tonal_reducer_get_peaks(self->tonal_reducer, peak_freqs_hz,
                                        max_peaks)
              : 0;
}

const float* spectral_denoiser_get_active_noise_profile(
    SpectralProcessorHandle instance) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;
  return self ? self->noise_spectrum : NULL;
}

void spectral_denoiser_reset_noise_profile(SpectralProcessorHandle instance) {
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;
  if (!self) {
    return;
  }
  if (self->noise_estimator) {
    noise_estimation_reset(self->noise_estimator);
  }
  if (self->tonal_reducer) {
    tonal_reducer_reset(self->tonal_reducer);
  }
  if (self->manual_noise_floor) {
    memset(self->manual_noise_floor, 0,
           self->real_spectrum_size * sizeof(float));
  }
  self->was_learning = false;
  self->last_adaptive_state = 0;
}

uint32_t spectral_denoiser_get_latency_frames(
    SpectralProcessorHandle instance) {
  if (!instance) {
    return 0;
  }
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;
  return hpss_filter_get_latency_frames(self->hpss_filter);
}

bool spectral_denoiser_is_transient_detected(SpectralProcessorHandle instance) {
  if (!instance) {
    return false;
  }
  SbSpectralDenoiser* self = (SbSpectralDenoiser*)instance;
  return hpss_filter_is_transient_detected(self->hpss_filter);
}
