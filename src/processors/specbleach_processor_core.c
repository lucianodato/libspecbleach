/*
libspecbleach - A spectral processing library

Copyright 2022-2026 Luciano Dato <lucianodato@gmail.com>

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

#include "specbleach_processor_core.h"
#include "shared/configurations.h"
#include "shared/utils/general_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

DenoiserParameters sb_denoiser_params_sanitize(
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
          (int)parameters->smoothing_mode == (int)SB_SMOOTHING_NLM_2D
              ? (int)SB_SMOOTHING_NLM_2D
              : (int)SB_SMOOTHING_TEMPORAL,
      .masking_depth = fmaxf(0.0f, fminf(1.0f, parameters->masking_depth)),
      .suppression_strength =
          fmaxf(0.0f, fminf(1.0f, parameters->suppression_strength)),
      .aggressiveness = parameters->aggressiveness,
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

const float* sb_curve_bias_copy(float** buffer, uint32_t* capacity,
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

SbProcessorCore* sb_processor_core_initialize(
    const uint32_t sample_rate, const float frame_size,
    const uint32_t overlap_factor, const ZeroPaddingType padding_type,
    const uint32_t zeropadding_amount, const WindowTypes input_window,
    const WindowTypes output_window, const uint32_t profile_spectrum_size) {
  if (sample_rate == 0 || frame_size <= 0.0f) {
    return NULL;
  }

  SbProcessorCore* core = (SbProcessorCore*)calloc(1U, sizeof(SbProcessorCore));
  if (!core) {
    return NULL;
  }

  core->sample_rate = sample_rate;
  core->stft_processor = stft_processor_initialize(
      sample_rate, frame_size, overlap_factor, padding_type, zeropadding_amount,
      input_window, output_window);

  if (!core->stft_processor) {
    sb_processor_core_free(core);
    return NULL;
  }

  uint32_t spectrum_size = profile_spectrum_size;
  if (spectrum_size == SB_PROCESSOR_CORE_DEFAULT_REAL_SPECTRUM) {
    spectrum_size = get_stft_real_spectrum_size(core->stft_processor);
  } else if (spectrum_size == SB_PROCESSOR_CORE_FULL_FFT_SPECTRUM) {
    spectrum_size = get_stft_fft_size(core->stft_processor);
  }

  core->noise_profile = noise_profile_initialize(spectrum_size);
  if (!core->noise_profile) {
    sb_processor_core_free(core);
    return NULL;
  }

  return core;
}

void sb_processor_core_free(SbProcessorCore* core) {
  if (!core) {
    return;
  }

  if (core->noise_profile) {
    noise_profile_free(core->noise_profile);
  }
  if (core->stft_processor) {
    stft_processor_free(core->stft_processor);
  }

  free(core);
}

uint32_t sb_processor_core_get_noise_profile_size(SbProcessorCore* core) {
  if (!core || !core->noise_profile) {
    return 0;
  }
  return get_noise_profile_size(core->noise_profile);
}

bool sb_processor_core_load_noise_profile_for_mode(
    SbProcessorCore* core, const float* restored_profile,
    const uint32_t profile_size, const uint32_t block_count, const int mode) {
  if (!core || !core->noise_profile || !restored_profile || mode < 1 ||
      mode > 4) {
    return false;
  }

  if (profile_size != get_noise_profile_size(core->noise_profile)) {
    return false;
  }

  return set_noise_profile(core->noise_profile, mode, restored_profile,
                           profile_size, block_count);
}

bool sb_processor_core_reset_noise_profile(SbProcessorCore* core) {
  if (!core || !core->noise_profile) {
    return false;
  }
  return reset_noise_profile(core->noise_profile);
}

uint32_t sb_processor_core_get_noise_profile_block_count_for_mode(
    SbProcessorCore* core, int mode) {
  if (!core || !core->noise_profile || mode < 1 || mode > 4) {
    return 0;
  }
  return get_noise_profile_block_count(core->noise_profile, mode);
}

float* sb_processor_core_get_noise_profile_for_mode(SbProcessorCore* core,
                                                    int mode) {
  if (!core || !core->noise_profile || mode < 1 || mode > 4) {
    return NULL;
  }
  return get_noise_profile(core->noise_profile, mode);
}

bool sb_processor_core_noise_profile_available_for_mode(SbProcessorCore* core,
                                                        int mode) {
  if (!core || !core->noise_profile || mode < 1 || mode > 4) {
    return false;
  }
  return is_noise_estimation_available(core->noise_profile, mode);
}
