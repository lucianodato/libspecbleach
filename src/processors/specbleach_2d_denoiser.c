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

#include "specbleach_2d_denoiser.h"
#include "denoiser2d/spectral_2d_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/stft/stft_processor.h"
#include "shared/utils/general_utils.h"
#include "shared/utils/tonal_detector.h"
#include "specbleach_processor_core.h"
#include <stdlib.h>
#include <string.h>

typedef struct Sb2DDenoiser {
  uint32_t hop;
  SbProcessorCore* core;
  SpectralProcessorHandle spectral_2d_denoiser;
} Sb2DDenoiser;

SpectralBleachHandle specbleach_2d_initialize(const uint32_t sample_rate,
                                              float frame_size) {
  if (sample_rate == 0 || frame_size <= 0.0f) {
    return NULL;
  }

  Sb2DDenoiser* self = (Sb2DDenoiser*)calloc(1U, sizeof(Sb2DDenoiser));
  if (!self) {
    return NULL;
  }

  self->core = sb_processor_core_initialize(
      sample_rate, frame_size, OVERLAP_FACTOR_2D, PADDING_CONFIGURATION_2D,
      ZEROPADDING_AMOUNT_2D, INPUT_WINDOW_TYPE_2D, OUTPUT_WINDOW_TYPE_2D,
      SB_PROCESSOR_CORE_DEFAULT_REAL_SPECTRUM);

  if (!self->core) {
    specbleach_2d_free(self);
    return NULL;
  }

  const uint32_t fft_size = get_stft_fft_size(self->core->stft_processor);
  self->hop = fft_size / OVERLAP_FACTOR_2D;

  self->spectral_2d_denoiser = spectral_2d_denoiser_initialize(
      sample_rate, fft_size, OVERLAP_FACTOR_2D, self->core->noise_profile);

  if (!self->spectral_2d_denoiser) {
    specbleach_2d_free(self);
    return NULL;
  }

  return self;
}

void specbleach_2d_free(SpectralBleachHandle instance) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;

  if (!self) {
    return;
  }

  if (self->spectral_2d_denoiser) {
    spectral_2d_denoiser_free(self->spectral_2d_denoiser);
  }
  if (self->core) {
    sb_processor_core_free(self->core);
  }

  free(self);
}

uint32_t specbleach_2d_get_latency(SpectralBleachHandle instance) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;

  if (!self || !self->core || !self->core->stft_processor) {
    return 0;
  }

  // Base STFT latency
  uint32_t stft_latency = get_stft_latency(self->core->stft_processor);

  // Additional NLM look-ahead latency
  uint32_t nlm_latency_frames =
      spectral_2d_denoiser_get_latency_frames(self->spectral_2d_denoiser);
  uint32_t nlm_latency_samples = nlm_latency_frames * self->hop;

  return stft_latency + nlm_latency_samples;
}

bool specbleach_2d_process(SpectralBleachHandle instance,
                           const uint32_t number_of_samples, const float* input,
                           float* output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;

  stft_processor_run(self->core->stft_processor, number_of_samples, input,
                     output, &spectral_2d_denoiser_run,
                     self->spectral_2d_denoiser);

  return true;
}

uint32_t specbleach_2d_get_noise_profile_size(SpectralBleachHandle instance) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? sb_processor_core_get_noise_profile_size(self->core) : 0;
}

bool specbleach_2d_load_noise_profile_for_mode(SpectralBleachHandle instance,
                                               const float* restored_profile,
                                               const uint32_t profile_size,
                                               const uint32_t block_count,
                                               const int mode) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? sb_processor_core_load_noise_profile_for_mode(
                    self->core, restored_profile, profile_size, block_count,
                    mode)
              : false;
}

bool specbleach_2d_reset_noise_profile(SpectralBleachHandle instance) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? sb_processor_core_reset_noise_profile(self->core) : false;
}

uint32_t specbleach_2d_get_noise_profile_block_count_for_mode(
    SpectralBleachHandle instance, int mode) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? sb_processor_core_get_noise_profile_block_count_for_mode(
                    self->core, mode)
              : 0;
}

float* specbleach_2d_get_noise_profile_for_mode(SpectralBleachHandle instance,
                                                int mode) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? sb_processor_core_get_noise_profile_for_mode(self->core, mode)
              : NULL;
}

bool specbleach_2d_noise_profile_available_for_mode(
    SpectralBleachHandle instance, int mode) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? sb_processor_core_noise_profile_available_for_mode(self->core,
                                                                   mode)
              : false;
}

bool specbleach_2d_load_parameters(
    SpectralBleachHandle instance,
    SpectralBleach2DDenoiserParameters parameters) {
  if (!instance) {
    return false;
  }

  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  Denoiser2DParameters denoise_parameters =
      sb_denoiser_2d_params_sanitize(parameters);

  load_2d_reduction_parameters(self->spectral_2d_denoiser, denoise_parameters);

  return true;
}

const float* specbleach_2d_get_tonal_mask(SpectralBleachHandle instance) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? spectral_2d_denoiser_get_tonal_mask(self->spectral_2d_denoiser)
              : NULL;
}

uint32_t specbleach_2d_get_tonal_peaks(SpectralBleachHandle instance,
                                       float* peak_freqs_hz,
                                       uint32_t max_peaks) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? spectral_2d_denoiser_get_peaks(self->spectral_2d_denoiser,
                                               peak_freqs_hz, max_peaks)
              : 0;
}

uint32_t specbleach_2d_get_tonal_peaks_for_profile(
    SpectralBleachHandle instance, const float* profile, uint32_t profile_size,
    float* peak_freqs_hz, uint32_t max_peaks) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  if (!self || !self->core || !profile || profile_size == 0) {
    return 0;
  }
  uint32_t fft_size = get_stft_fft_size(self->core->stft_processor);
  return tonal_detector_get_peaks_from_profile(
      profile, profile_size, self->core->sample_rate, fft_size, peak_freqs_hz,
      max_peaks);
}

float* specbleach_2d_get_active_noise_profile(SpectralBleachHandle instance) {
  Sb2DDenoiser* self = (Sb2DDenoiser*)instance;
  return self ? spectral_2d_denoiser_get_active_noise_profile(
                    self->spectral_2d_denoiser)
              : NULL;
}
