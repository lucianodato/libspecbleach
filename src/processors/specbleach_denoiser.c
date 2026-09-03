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

#include "specbleach_denoiser.h"
#include "denoiser/spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/stft/stft_processor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct specbleach_denoiser { // NOLINT(readability-identifier-naming)
  uint32_t hop;
  uint32_t sample_rate;
  float frame_size_ms;
  uint32_t frame_size_samples;
  uint32_t fft_size;
  StftProcessor* stft_processor;
  NoiseProfile* noise_profile;
  SpectralProcessorHandle spectral_denoiser;
  float* reduction_curve_copy;
  uint32_t reduction_curve_capacity;
} SbDenoiserInstance;

static DenoiserParameters sanitize_denoiser_parameters(
    const SpecbleachDenoiserParameters* parameters) {
  return (DenoiserParameters){
      .learn_noise = (int)parameters->learn_noise,
      .residual_listen = parameters->residual_listen,
      .reduction_amount = fmaxf(0.0f, fminf(1.0f, parameters->reduction_gain)),
      .smoothing_factor =
          fmaxf(0.0f, fminf(1.0f, parameters->smoothing_factor)),
      .whitening_factor =
          fmaxf(0.0f, fminf(1.0f, parameters->whitening_factor)),
      .adaptive_noise = parameters->adaptive_noise ? 1 : 0,
      .noise_estimation_method = (int)parameters->noise_estimation_method,
      .smoothing_mode =
          (int)parameters->smoothing_mode == (int)SPECBLEACH_SMOOTHING_NLM_2D
              ? (int)SPECBLEACH_SMOOTHING_NLM_2D
              : (int)SPECBLEACH_SMOOTHING_TEMPORAL,
      .masking_depth = fmaxf(0.0f, fminf(1.0f, parameters->masking_depth)),
      .suppression_strength =
          fmaxf(0.0f, fminf(1.0f, parameters->suppression_strength)),
      .aggressiveness = fmaxf(-1.0f, fminf(1.0f, parameters->aggressiveness)),
      .tonal_reduction =
          fmaxf(0.0f, fminf(1.0f, parameters->tonal_reduction_gain)),
      .hpss_enable = parameters->hpss_enable ? 1 : 0,
      .noise_profile_offset_linear =
          fmaxf(0.01f, fminf(100.0f, parameters->noise_profile_scale > 0.0f
                                         ? parameters->noise_profile_scale
                                         : 1.0f)),
      .tonal_noise_profile_offset_linear = fmaxf(
          0.01f, fminf(100.0f, parameters->tonal_noise_profile_scale > 0.0f
                                   ? parameters->tonal_noise_profile_scale
                                   : 1.0f)),
      .reduction_curve_bias = NULL,
      .reduction_curve_enabled = parameters->reduction_curve_enabled,
  };
}

static const float* curve_bias_copy(float** buffer, uint32_t* capacity,
                                    const uint32_t required_size,
                                    const bool enabled, const float* source) {
  if (!buffer || !capacity || !enabled || !source || required_size == 0) {
    return NULL;
  }

  if (!*buffer || *capacity < required_size) {
    float* resized = (float*)realloc(*buffer, required_size * sizeof(float));
    if (!resized) {
      return NULL;
    }
    *buffer = resized;
    *capacity = required_size;
  }

  memcpy(*buffer, source, required_size * sizeof(float));
  return *buffer;
}

SpecbleachDenoiserParameters specbleach_denoiser_get_default_parameters(void) {
  SpecbleachDenoiserParameters p;
  memset(&p, 0, sizeof(p));
  p.learn_noise = SPECBLEACH_LEARN_OFF;
  p.residual_listen = false;
  p.reduction_gain = 0.1f;
  p.smoothing_factor = 0.5f;
  p.smoothing_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  p.whitening_factor = 0.0f;
  p.adaptive_noise = false;
  p.noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE;
  p.masking_depth = 0.5f;
  p.suppression_strength = 1.0f;
  p.aggressiveness = 0.0f;
  p.tonal_reduction_gain = 0.1f;
  p.hpss_enable = false;
  p.noise_profile_scale = 1.0f;
  p.reduction_curve_bias = NULL;
  p.reduction_curve_enabled = false;
  p.reduction_curve_size = 0;
  p.tonal_noise_profile_scale = 1.0f;
  return p;
}

static bool rebuild_engines(SbDenoiserInstance* self) {
  if (self->spectral_denoiser) {
    spectral_denoiser_free(self->spectral_denoiser);
    self->spectral_denoiser = NULL;
  }
  if (self->noise_profile) {
    noise_profile_free(self->noise_profile);
    self->noise_profile = NULL;
  }
  if (self->stft_processor) {
    stft_processor_free(self->stft_processor);
    self->stft_processor = NULL;
  }
  self->stft_processor = stft_processor_initialize(
      self->sample_rate, self->frame_size_ms, OVERLAP_FACTOR,
      PADDING_CONFIGURATION, ZEROPADDING_AMOUNT, INPUT_WINDOW_TYPE,
      OUTPUT_WINDOW_TYPE);
  if (!self->stft_processor) {
    return false;
  }
  const uint32_t fft_size = get_stft_fft_size(self->stft_processor);
  self->fft_size = fft_size;
  self->hop = get_stft_hop_size(self->stft_processor);
  self->noise_profile = noise_profile_initialize(
      get_stft_real_spectrum_size(self->stft_processor));
  if (!self->noise_profile) {
    return false;
  }
  self->spectral_denoiser = spectral_denoiser_initialize(
      self->sample_rate, fft_size, OVERLAP_FACTOR, self->noise_profile);
  return self->spectral_denoiser != NULL;
}

specbleach_denoiser* specbleach_denoiser_initialize(uint32_t sample_rate,
                                                    float frame_size_ms) {
  if (sample_rate == 0 || !isfinite(frame_size_ms) || frame_size_ms <= 0.0f) {
    return NULL;
  }

  const double frame_samples =
      ((double)frame_size_ms / 1000.0) * (double)sample_rate;
  if (!isfinite(frame_samples) || frame_samples < 1.0 ||
      frame_samples > (double)UINT32_MAX) {
    return NULL;
  }

  SbDenoiserInstance* self =
      (SbDenoiserInstance*)calloc(1U, sizeof(SbDenoiserInstance));
  if (!self) {
    return NULL;
  }

  self->sample_rate = sample_rate;
  self->frame_size_ms = frame_size_ms;
  self->frame_size_samples = (uint32_t)frame_samples;
  self->stft_processor = stft_processor_initialize(
      sample_rate, frame_size_ms, OVERLAP_FACTOR, PADDING_CONFIGURATION,
      ZEROPADDING_AMOUNT, INPUT_WINDOW_TYPE, OUTPUT_WINDOW_TYPE);

  if (!self->stft_processor) {
    specbleach_denoiser_free(self);
    return NULL;
  }

  self->fft_size = get_stft_fft_size(self->stft_processor);
  self->hop = get_stft_hop_size(self->stft_processor);

  self->noise_profile = noise_profile_initialize(
      get_stft_real_spectrum_size(self->stft_processor));
  if (!self->noise_profile) {
    specbleach_denoiser_free(self);
    return NULL;
  }

  self->spectral_denoiser = spectral_denoiser_initialize(
      self->sample_rate, self->fft_size, OVERLAP_FACTOR, self->noise_profile);

  if (!self->spectral_denoiser) {
    specbleach_denoiser_free(self);
    return NULL;
  }

  return self;
}

void specbleach_denoiser_free(specbleach_denoiser* instance) {
  SbDenoiserInstance* self = instance;

  if (!self) {
    return;
  }

  if (self->spectral_denoiser) {
    spectral_denoiser_free(self->spectral_denoiser);
  }
  if (self->noise_profile) {
    noise_profile_free(self->noise_profile);
  }
  if (self->stft_processor) {
    stft_processor_free(self->stft_processor);
  }
  if (self->reduction_curve_copy) {
    free(self->reduction_curve_copy);
  }

  free(self);
}

uint32_t specbleach_denoiser_get_latency(const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;

  if (!self || !self->stft_processor) {
    return 0;
  }

  uint32_t stft_latency = get_stft_latency(self->stft_processor);
  uint32_t denoiser_latency_frames =
      spectral_denoiser_get_latency_frames(self->spectral_denoiser);
  uint32_t denoiser_latency_samples = denoiser_latency_frames * self->hop;

  return stft_latency + denoiser_latency_samples;
}

bool specbleach_denoiser_process(specbleach_denoiser* instance,
                                 uint32_t number_of_samples, const float* input,
                                 float* output) {
  SbDenoiserInstance* self = instance;
  if (!instance || !input || !output) {
    return false;
  }
  if (number_of_samples == 0) {
    return false;
  }

  stft_processor_run(self->stft_processor, number_of_samples, input, output,
                     &spectral_denoiser_run, self->spectral_denoiser);

  return true;
}

uint32_t specbleach_denoiser_get_sample_rate(
    const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return self ? self->sample_rate : 0;
}

uint32_t specbleach_denoiser_get_frame_size(
    const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return self ? self->frame_size_samples : 0;
}

uint32_t specbleach_denoiser_get_fft_size(const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return self ? self->fft_size : 0;
}

uint32_t specbleach_denoiser_get_hop_size(const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return self ? self->hop : 0;
}

uint32_t specbleach_denoiser_get_noise_profile_size(
    const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return (self && self->noise_profile)
             ? get_noise_profile_size(self->noise_profile)
             : 0;
}

bool specbleach_denoiser_load_noise_profile_for_mode(
    specbleach_denoiser* instance, const float* restored_profile,
    const uint32_t profile_size, const uint32_t block_count,
    const SpecbleachProfileMode mode) {
  SbDenoiserInstance* self = instance;
  if (!self || !self->noise_profile || !restored_profile) {
    return false;
  }
  if (mode < SPECBLEACH_PROFILE_MODE_FIRST ||
      mode > SPECBLEACH_PROFILE_MODE_LAST) {
    return false;
  }

  if (profile_size != get_noise_profile_size(self->noise_profile)) {
    return false;
  }

  return set_noise_profile(self->noise_profile, (int)mode, restored_profile,
                           profile_size, block_count);
}

bool specbleach_denoiser_load_noise_profile_resampled(
    specbleach_denoiser* instance, const float* restored_profile,
    const uint32_t source_size, const uint32_t block_count,
    const SpecbleachProfileMode mode) {
  SbDenoiserInstance* self = instance;
  if (!self || !self->noise_profile || !restored_profile || source_size == 0) {
    return false;
  }
  if (mode < SPECBLEACH_PROFILE_MODE_FIRST ||
      mode > SPECBLEACH_PROFILE_MODE_LAST) {
    return false;
  }
  const uint32_t native_size = get_noise_profile_size(self->noise_profile);
  if (source_size == native_size) {
    return specbleach_denoiser_load_noise_profile_for_mode(
        instance, restored_profile, source_size, block_count, mode);
  }
  float* tmp = (float*)malloc(native_size * sizeof(float));
  if (!tmp) {
    return false;
  }
  // ponytail: linear interpolation, caller's curve if higher quality needed
  for (uint32_t i = 0; i < native_size; ++i) {
    const float pos =
        (native_size == 1)
            ? 0.0f
            : ((float)i * (float)(source_size - 1) / (float)(native_size - 1));
    const uint32_t lo = (uint32_t)pos;
    const uint32_t hi = lo + 1 < source_size ? lo + 1 : lo;
    const float frac = pos - (float)lo;
    tmp[i] = restored_profile[lo] * (1.0f - frac) + restored_profile[hi] * frac;
  }
  const bool ok = specbleach_denoiser_load_noise_profile_for_mode(
      instance, tmp, native_size, block_count, mode);
  free(tmp);
  return ok;
}

void specbleach_denoiser_reset_noise_profile(specbleach_denoiser* instance) {
  SbDenoiserInstance* self = instance;
  if (!self) {
    return;
  }
  if (self->spectral_denoiser) {
    spectral_denoiser_reset_noise_profile(self->spectral_denoiser);
  }
  if (self->noise_profile) {
    reset_noise_profile(self->noise_profile);
  }
}

bool specbleach_denoiser_reset_dsp_state(specbleach_denoiser* instance) {
  SbDenoiserInstance* self = instance;
  if (!self) {
    return false;
  }
  return rebuild_engines(self);
}

uint32_t specbleach_denoiser_get_noise_profile_block_count_for_mode(
    const specbleach_denoiser* instance, SpecbleachProfileMode mode) {
  const SbDenoiserInstance* self = instance;
  if (!self || !self->noise_profile || mode < SPECBLEACH_PROFILE_MODE_FIRST ||
      mode > SPECBLEACH_PROFILE_MODE_LAST) {
    return 0;
  }
  return get_noise_profile_block_count(self->noise_profile, (int)mode);
}

const float* specbleach_denoiser_get_noise_profile_for_mode(
    const specbleach_denoiser* instance, SpecbleachProfileMode mode) {
  const SbDenoiserInstance* self = instance;
  if (!self || !self->noise_profile || mode < SPECBLEACH_PROFILE_MODE_FIRST ||
      mode > SPECBLEACH_PROFILE_MODE_LAST) {
    return NULL;
  }
  return get_noise_profile(self->noise_profile, (int)mode);
}

bool specbleach_denoiser_noise_profile_available_for_mode(
    const specbleach_denoiser* instance, SpecbleachProfileMode mode) {
  const SbDenoiserInstance* self = instance;
  if (!self || !self->noise_profile || mode < SPECBLEACH_PROFILE_MODE_FIRST ||
      mode > SPECBLEACH_PROFILE_MODE_LAST) {
    return false;
  }
  return is_noise_estimation_available(self->noise_profile, (int)mode);
}

bool specbleach_denoiser_has_any_profile(const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  if (!self || !self->noise_profile) {
    return false;
  }
  for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
    if (is_noise_estimation_available(self->noise_profile, mode)) {
      return true;
    }
  }
  return false;
}

bool specbleach_denoiser_load_parameters(
    specbleach_denoiser* instance,
    const SpecbleachDenoiserParameters* parameters, uint32_t parameters_size) {
  SbDenoiserInstance* self = instance;
  if (!instance || !parameters) {
    return false;
  }
  if (parameters_size != sizeof(SpecbleachDenoiserParameters)) {
    return false;
  }

  const uint32_t profile_size = get_noise_profile_size(self->noise_profile);
  if (parameters->reduction_curve_enabled &&
      parameters->reduction_curve_size != profile_size) {
    return false;
  }

  if (parameters->reduction_curve_enabled &&
      parameters->reduction_curve_bias == NULL) {
    return false;
  }

  const float* owned_bias = curve_bias_copy(
      &self->reduction_curve_copy, &self->reduction_curve_capacity,
      profile_size, parameters->reduction_curve_enabled,
      parameters->reduction_curve_bias);
  if (parameters->reduction_curve_enabled && !owned_bias) {
    return false;
  }

  DenoiserParameters denoise_parameters =
      sanitize_denoiser_parameters(parameters);
  denoise_parameters.reduction_curve_bias = owned_bias;

  return load_reduction_parameters(self->spectral_denoiser,
                                   denoise_parameters);
}

const float* specbleach_denoiser_get_tonal_mask(
    const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return self ? spectral_denoiser_get_tonal_mask(self->spectral_denoiser)
              : NULL;
}

uint32_t specbleach_denoiser_get_tonal_peaks(
    const specbleach_denoiser* instance, float* peak_freqs_hz,
    uint32_t max_peaks) {
  const SbDenoiserInstance* self = instance;
  return self ? spectral_denoiser_get_peaks(self->spectral_denoiser,
                                            peak_freqs_hz, max_peaks)
              : 0;
}

const float* specbleach_denoiser_get_active_noise_profile(
    const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return self ? spectral_denoiser_get_active_noise_profile(
                    self->spectral_denoiser)
              : NULL;
}

bool specbleach_denoiser_is_transient_detected(
    const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return self ? spectral_denoiser_is_transient_detected(self->spectral_denoiser)
              : false;
}

float specbleach_denoiser_get_transient_intensity(
    const specbleach_denoiser* instance) {
  const SbDenoiserInstance* self = instance;
  return self ? spectral_denoiser_get_transient_intensity(
                    self->spectral_denoiser)
              : 0.0f;
}
