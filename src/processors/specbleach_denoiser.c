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
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/stft/stft_processor.h"
#include "shared/utils/general_utils.h"
#include "shared/utils/tonal_detector.h"
#include "specbleach_processor_core.h"
#include <stdlib.h>
#include <string.h>

typedef struct specbleach_denoiser { // NOLINT(readability-identifier-naming)
  uint32_t hop;
  SbProcessorCore* core;
  SpectralProcessorHandle spectral_denoiser;
  float* reduction_curve_copy;
  uint32_t reduction_curve_capacity;
} SbSpectralDenoiser;

specbleach_denoiser* specbleach_denoiser_initialize(const uint32_t sample_rate,
                                                    float frame_size) {
  SbSpectralDenoiser* self =
      (SbSpectralDenoiser*)calloc(1U, sizeof(SbSpectralDenoiser));
  if (!self) {
    return NULL;
  }

  self->core = sb_processor_core_initialize(
      sample_rate, frame_size, OVERLAP_FACTOR_1D, PADDING_CONFIGURATION_1D,
      ZEROPADDING_AMOUNT_1D, INPUT_WINDOW_TYPE_1D, OUTPUT_WINDOW_TYPE_1D, 0);

  if (!self->core) {
    specbleach_denoiser_free(self);
    return NULL;
  }

  const uint32_t fft_size = get_stft_fft_size(self->core->stft_processor);
  self->hop = get_stft_hop_size(self->core->stft_processor);

  self->spectral_denoiser = spectral_denoiser_initialize(
      self->core->sample_rate, fft_size, OVERLAP_FACTOR_1D,
      self->core->noise_profile);

  if (!self->spectral_denoiser) {
    specbleach_denoiser_free(self);
    return NULL;
  }

  return self;
}

void specbleach_denoiser_free(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;

  if (!self) {
    return;
  }

  if (self->spectral_denoiser) {
    spectral_denoiser_free(self->spectral_denoiser);
  }
  if (self->core) {
    sb_processor_core_free(self->core);
  }
  if (self->reduction_curve_copy) {
    free(self->reduction_curve_copy);
  }

  free(self);
}

uint32_t specbleach_denoiser_get_latency(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;

  if (!self || !self->core || !self->core->stft_processor) {
    return 0;
  }

  uint32_t stft_latency = get_stft_latency(self->core->stft_processor);
  uint32_t denoiser_latency_frames =
      spectral_denoiser_get_latency_frames(self->spectral_denoiser);
  uint32_t denoiser_latency_samples = denoiser_latency_frames * self->hop;

  return stft_latency + denoiser_latency_samples;
}

bool specbleach_denoiser_process(specbleach_denoiser* instance,
                                 const uint32_t number_of_samples,
                                 const float* input, float* output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  SbSpectralDenoiser* self = instance;

  stft_processor_run(self->core->stft_processor, number_of_samples, input,
                     output, &spectral_denoiser_run, self->spectral_denoiser);

  return true;
}

uint32_t specbleach_denoiser_get_noise_profile_size(
    specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_size(self->core) : 0;
}

bool specbleach_denoiser_load_noise_profile_for_mode(
    specbleach_denoiser* instance, const float* restored_profile,
    const uint32_t profile_size, const uint32_t block_count, const int mode) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_load_noise_profile_for_mode(
                    self->core, restored_profile, profile_size, block_count,
                    mode)
              : false;
}

bool specbleach_denoiser_reset_noise_profile(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  if (!self) {
    return false;
  }
  if (self->spectral_denoiser) {
    spectral_denoiser_reset_noise_profile(self->spectral_denoiser);
  }
  return sb_processor_core_reset_noise_profile(self->core);
}

uint32_t specbleach_denoiser_get_noise_profile_block_count_for_mode(
    specbleach_denoiser* instance, int mode) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_block_count_for_mode(
                    self->core, mode)
              : 0;
}

float* specbleach_denoiser_get_noise_profile_for_mode(
    specbleach_denoiser* instance, int mode) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_for_mode(self->core, mode)
              : NULL;
}

bool specbleach_denoiser_noise_profile_available_for_mode(
    specbleach_denoiser* instance, int mode) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_noise_profile_available_for_mode(self->core,
                                                                   mode)
              : false;
}

bool specbleach_denoiser_load_parameters(
    specbleach_denoiser* instance,
    const SpecbleachDenoiserParameters* parameters,
    const uint32_t parameters_size) {
  if (!instance || !parameters ||
      parameters_size != sizeof(SpecbleachDenoiserParameters)) {
    return false;
  }

  SbSpectralDenoiser* self = instance;

  const uint32_t profile_size =
      sb_processor_core_get_noise_profile_size(self->core);
  if (parameters->reduction_curve_enabled &&
      parameters->reduction_curve_size != profile_size) {
    return false;
  }

  const float* owned_bias = sb_curve_bias_copy(
      &self->reduction_curve_copy, &self->reduction_curve_capacity,
      profile_size, parameters->reduction_curve_enabled,
      parameters->reduction_curve_bias);
  if (parameters->reduction_curve_enabled && !owned_bias) {
    return false;
  }

  DenoiserParameters denoise_parameters =
      sb_denoiser_params_sanitize(parameters);
  denoise_parameters.reduction_curve_bias = owned_bias;

  return load_reduction_parameters(self->spectral_denoiser, denoise_parameters);
}

const float* specbleach_denoiser_get_tonal_mask(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_get_tonal_mask(self->spectral_denoiser)
              : NULL;
}

uint32_t specbleach_denoiser_get_tonal_peaks(specbleach_denoiser* instance,
                                             float* peak_freqs_hz,
                                             uint32_t max_peaks) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_get_peaks(self->spectral_denoiser,
                                            peak_freqs_hz, max_peaks)
              : 0;
}

uint32_t specbleach_denoiser_get_tonal_peaks_for_profile(
    specbleach_denoiser* instance, const float* profile, uint32_t profile_size,
    float* peak_freqs_hz, uint32_t max_peaks) {
  SbSpectralDenoiser* self = instance;
  if (!self || !self->core || !profile || profile_size == 0 ||
      profile_size != sb_processor_core_get_noise_profile_size(self->core)) {
    return 0;
  }
  uint32_t fft_size = get_stft_fft_size(self->core->stft_processor);

  float* median_profile =
      sb_processor_core_get_noise_profile_for_mode(self->core, MEDIAN);
  const float* eval_profile = median_profile ? median_profile : profile;

  return tonal_detector_get_peaks_from_profile(
      eval_profile, profile_size, self->core->sample_rate, fft_size,
      peak_freqs_hz, max_peaks);
}

const float* specbleach_denoiser_get_active_noise_profile(
    specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_get_active_noise_profile(
                    self->spectral_denoiser)
              : NULL;
}

bool specbleach_denoiser_is_transient_detected(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_is_transient_detected(self->spectral_denoiser)
              : false;
}

float specbleach_denoiser_get_transient_intensity(
    specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_get_transient_intensity(
                    self->spectral_denoiser)
              : 0.0f;
}
