#include "spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/denoiser_post_process.h"
#include "shared/denoiser_logic/core/denoiser_profile_core.h"
#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/denoiser_logic/processing/gain_calculator.h"
#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/denoiser_logic/processing/suppression_engine.h"
#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/stft/stft_processor.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/spectral_features.h"
#include "shared/utils/spectral_smoother.h"
#include "shared/utils/spectral_utils.h"
#include "shared/utils/transient_detector.h"
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
  CriticalBands* critical_bands;
  TransientDetector* transient_detector;

  float* band_energies;
  float* onset_weights;
  float* transient_mask;
  float* smoothed_magnitude;
  float* clean_magnitude;

  bool is_transient_detected;
  float transient_intensity;

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
                                       1.0F);

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

  self->critical_bands = critical_bands_initialize(
      self->sample_rate, self->fft_size, self->band_type);
  uint32_t num_bands = self->critical_bands
                           ? get_number_of_critical_bands(self->critical_bands)
                           : 0U;
  self->transient_detector = transient_detector_initialize(num_bands);

  self->band_energies =
      (float*)calloc(num_bands > 0U ? num_bands : 1U, sizeof(float));
  self->onset_weights =
      (float*)calloc(num_bands > 0U ? num_bands : 1U, sizeof(float));
  self->transient_mask =
      (float*)calloc(self->real_spectrum_size, sizeof(float));

  self->smoothed_magnitude =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->clean_magnitude =
      (float*)calloc(self->real_spectrum_size, sizeof(float));

  if (!self->noise_floor_manager || !self->masking_veto ||
      !self->suppression_engine || !self->critical_bands ||
      !self->transient_detector || !self->band_energies ||
      !self->onset_weights || !self->transient_mask ||
      !self->smoothed_magnitude || !self->clean_magnitude) {
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
  if (self->critical_bands) {
    critical_bands_free(self->critical_bands);
  }
  if (self->transient_detector) {
    transient_detector_free(self->transient_detector);
  }
  if (self->band_energies) {
    free(self->band_energies);
  }
  if (self->onset_weights) {
    free(self->onset_weights);
  }
  if (self->transient_mask) {
    free(self->transient_mask);
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

  if (self->smoothed_magnitude) {
    free(self->smoothed_magnitude);
  }
  if (self->clean_magnitude) {
    free(self->clean_magnitude);
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
      .tonal_noise_profile_offset_linear =
          self->denoise_parameters.tonal_noise_profile_offset_linear,
      // Previous frame's tonal mask (one-frame latency) is used to apply the
      // tonal threshold offset at detected tonal bins
      .tonal_mask = tonal_reducer_get_mask(self->tonal_reducer),
  };
  denoiser_profile_core_update(profile_params, reference_spectrum);

  // 2.1 Transient Detection via Transient Detector across Critical Bands
  // Transient detection runs on a clean signal estimate with scaled-up noise
  // subtraction to avoid false triggering from residual musical noise.
  bool transient_enabled = (self->denoise_parameters.hpss_enable != 0);
  if (transient_enabled && self->critical_bands && self->transient_detector) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      // Scale noise up using TRANSIENT_CLEAN_NOISE_SCALE to eliminate spurious
      // noise peaks
      float clean = fmaxf(reference_spectrum[k] - (TRANSIENT_CLEAN_NOISE_SCALE *
                                                   self->noise_spectrum[k]),
                          0.0F);
      self->clean_magnitude[k] = clean;
    }

    compute_critical_bands_spectrum(self->critical_bands, self->clean_magnitude,
                                    self->band_energies);
    self->is_transient_detected = transient_detector_process(
        self->transient_detector, self->band_energies, self->onset_weights,
        &self->transient_intensity);

    // Expand critical band onset weights to per-bin transient mask
    uint32_t num_bands = get_number_of_critical_bands(self->critical_bands);
    memset(self->transient_mask, 0, self->real_spectrum_size * sizeof(float));
    for (uint32_t b = 0; b < num_bands; ++b) {
      float bw = self->onset_weights[b];
      if (bw > 0.0F) {
        CriticalBandIndexes idx = get_band_indexes(self->critical_bands, b);
        uint32_t end = (idx.end_position < self->real_spectrum_size)
                           ? idx.end_position
                           : self->real_spectrum_size;
        for (uint32_t k = idx.start_position; k < end; ++k) {
          self->transient_mask[k] = fmaxf(self->transient_mask[k], bw);
        }
      }
    }
  } else {
    self->is_transient_detected = false;
    self->transient_intensity = 0.0f;
    memset(self->transient_mask, 0, self->real_spectrum_size * sizeof(float));
  }

  // 3. Pre-Subtraction Smoothing:
  // Pre-subtraction temporal and spatial smoothing on input magnitude to
  // suppress musical noise. Transient bins are spared bin-by-bin during gain
  // calculation and time smoothing.
  const float* input_magnitude = reference_spectrum;
  const float* effective_magnitude;

  if (self->denoise_parameters.smoothing_factor > 0.0F) {
    float smooth_alpha =
        fminf(0.95F, self->denoise_parameters.smoothing_factor * 0.90F);

    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float raw = input_magnitude[k];
      float prev = self->smoothed_magnitude[k];
      // Bin-by-bin adaptive smoothing: open immediately on transient bins while
      // keeping full smoothing elsewhere
      float t_w = self->transient_mask[k];
      float adapt_alpha = (1.0F - t_w) * smooth_alpha;
      self->smoothed_magnitude[k] =
          (adapt_alpha * prev) + ((1.0F - adapt_alpha) * raw);
    }

    int spatial_passes =
        (int)(self->denoise_parameters.smoothing_factor * 2.0F);
    for (int p = 0; p < spatial_passes; ++p) {
      spectral_smoothing_apply_spatial(self->smoothed_magnitude,
                                       self->real_spectrum_size);
    }

    effective_magnitude = self->smoothed_magnitude;
  } else {
    memcpy(self->smoothed_magnitude, input_magnitude,
           self->real_spectrum_size * sizeof(float));
    effective_magnitude = input_magnitude;
  }

  // 3.1. Calculate SNR-dependent oversubtraction factors (Alpha/Beta)
  SuppressionParameters suppression_params = {
      .type = SUPPRESSION_BEROUTI_PER_BIN,
      .strength = self->denoise_parameters.suppression_strength,
      .undersubtraction = 0.0F};
  suppression_engine_calculate(self->suppression_engine, effective_magnitude,
                               self->noise_spectrum, suppression_params,
                               self->alpha, self->beta);

  // 3.2. Detect tonal components and boost alpha at tonal bins
  tonal_reducer_run(self->tonal_reducer, self->noise_spectrum,
                    get_noise_profile(self->noise_profile, CV_MASK),
                    is_noise_estimation_available(self->noise_profile, CV_MASK),
                    self->alpha, self->denoise_parameters.tonal_reduction);

  // 3.3. Apply Structural Veto to rescue psychoacoustically masked signal
  masking_veto_apply(self->masking_veto, effective_magnitude,
                     self->noise_spectrum, NULL, self->alpha,
                     self->denoise_parameters.masking_depth);

  // 3.4. When transients are detected and enabled, drop alphas firmly to
  // ALPHA_MIN (1.0) strictly on the specific frequencies where transient energy
  // was detected. Minimum alpha is never allowed below ALPHA_MIN (1.0).
  if (self->is_transient_detected) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_mask[k];
      if (t_weight > 0.0F) {
        float prot_factor = sqrtf(t_weight);
        self->alpha[k] =
            (self->alpha[k] * (1.0F - prot_factor)) + (ALPHA_MIN * prot_factor);
        self->alpha[k] = fmaxf(ALPHA_MIN, self->alpha[k]);
      }
    }
  }

  // 3.5. Gain Calculation & Smoothing
  calculate_gains(self->real_spectrum_size, self->fft_size, effective_magnitude,
                  self->noise_spectrum, self->gain_spectrum, self->alpha,
                  self->beta, self->gain_calculation_type);

  // Transient Protection: ensure transient bins have gain near 1.0 (preserved
  // in residual listen)
  if (self->is_transient_detected) {
    for (uint32_t k = 0U; k < self->real_spectrum_size; ++k) {
      float t_weight = self->transient_mask[k];
      if (t_weight > 0.0F) {
        // Protect transient bins from being attenuated: gain -> 1.0
        self->gain_spectrum[k] = fmaxf(self->gain_spectrum[k], t_weight);
      }
    }
  }

  // On detected transients, time smoother attacks instantly (0 attack time)
  // strictly on transient bins, while non-transient bins receive 100% of full
  // time smoothing to prevent musical noise.
  TimeSmoothingParameters spectral_smoothing_parameters =
      (TimeSmoothingParameters){
          .smoothing = self->denoise_parameters.smoothing_factor,
          .transient_mask = self->transient_mask,
      };
  spectral_smoothing_run(self->spectrum_smoothing,
                         spectral_smoothing_parameters, self->gain_spectrum);
  if (self->denoise_parameters.smoothing_factor > 0.0f) {
    int passes = 1 + (int)(self->denoise_parameters.smoothing_factor * 2.0f);
    for (int p = 0; p < passes; ++p) {
      spectral_smoothing_apply_spatial(self->gain_spectrum,
                                       self->real_spectrum_size);
    }
  }

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
      .noise_spectrum = self->noise_spectrum,
      .fft_spectrum = fft_spectrum,
      .reduction_curve_bias = self->denoise_parameters.reduction_curve_bias,
  };
  denoiser_post_process_apply(post_params);

  return true;
}

const float* spectral_denoiser_get_tonal_mask(
    SpectralProcessorHandle instance) {
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  return (self && self->tonal_reducer)
             ? tonal_reducer_get_mask(self->tonal_reducer)
             : NULL;
}

uint32_t spectral_denoiser_get_peaks(SpectralProcessorHandle instance,
                                     float* peak_freqs_hz, uint32_t max_peaks) {
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  return self ? tonal_reducer_get_peaks(self->tonal_reducer, peak_freqs_hz,
                                        max_peaks)
              : 0;
}

const float* spectral_denoiser_get_active_noise_profile(
    SpectralProcessorHandle instance) {
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
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
  (void)instance;
  return 0U;
}

bool spectral_denoiser_is_transient_detected(SpectralProcessorHandle instance) {
  if (!instance) {
    return false;
  }
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  return self->is_transient_detected;
}

float spectral_denoiser_get_transient_intensity(
    SpectralProcessorHandle instance) {
  if (!instance) {
    return 0.0f;
  }
  const SbSpectralDenoiser* self = (const SbSpectralDenoiser*)instance;
  return self->transient_intensity;
}
