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

#include "specbleach_stereo.h"
#include <stdlib.h>

typedef struct specbleach_stereo { // NOLINT(readability-identifier-naming)
  uint32_t channels;
  void** instances; /* one opaque engine handle per channel */
  SpecbleachError last_error;
} SpecbleachStereoState;

specbleach_stereo* specbleach_stereo_initialize(uint32_t sample_rate,
                                                float frame_size_ms,
                                                uint32_t channels) {
  if (sample_rate == 0 || frame_size_ms <= 0.0f || channels == 0) {
    return NULL;
  }

  SpecbleachStereoState* self =
      (SpecbleachStereoState*)calloc(1U, sizeof(SpecbleachStereoState));
  if (!self) {
    return NULL;
  }
  self->last_error = SPECBLEACH_OK;

  self->channels = channels;
  self->instances = calloc((size_t)channels, sizeof(void*));
  if (!self->instances) {
    specbleach_stereo_free(self);
    return NULL;
  }

  for (uint32_t ch = 0; ch < channels; ++ch) {
    self->instances[ch] =
        specbleach_denoiser_initialize(sample_rate, frame_size_ms);
    if (!self->instances[ch]) {
      specbleach_stereo_free(self);
      return NULL;
    }
  }

  return self;
}

void specbleach_stereo_free(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return;
  }

  if (self->instances) {
    for (uint32_t ch = 0; ch < self->channels; ++ch) {
      if (self->instances[ch]) {
        specbleach_denoiser_free((specbleach_denoiser*)self->instances[ch]);
      }
    }
    free(self->instances);
  }

  free(self);
}

bool specbleach_stereo_load_parameters(
    specbleach_stereo* instance, const SpecbleachDenoiserParameters* parameters,
    uint32_t parameters_size) {
  SpecbleachStereoState* self = instance;

  if (!self || !parameters ||
      parameters_size != sizeof(SpecbleachDenoiserParameters)) {
    return false;
  }

  bool result = true;
  SpecbleachError first_error = SPECBLEACH_OK;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (!specbleach_denoiser_load_parameters(
            (specbleach_denoiser*)self->instances[ch], parameters,
            parameters_size)) {
      if (result) {
        first_error = specbleach_denoiser_get_last_error(
            (specbleach_denoiser*)self->instances[ch]);
      }
      result = false;
    }
  }
  self->last_error = result ? SPECBLEACH_OK : first_error;
  return result;
}

bool specbleach_stereo_process(specbleach_stereo* instance,
                               uint32_t number_of_samples,
                               const float** input, float** output) {
  SpecbleachStereoState* self = instance;

  if (!self || number_of_samples == 0 || !input || !output) {
    if (self) {
      self->last_error = (!self) ? SPECBLEACH_ERR_NULL_ARG
                         : number_of_samples == 0 ? SPECBLEACH_ERR_EMPTY
                                                  : SPECBLEACH_ERR_NULL_ARG;
    }
    return false;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (!input[ch] || !output[ch]) {
      self->last_error = SPECBLEACH_ERR_NULL_ARG;
      return false;
    }
    if (!specbleach_denoiser_process((specbleach_denoiser*)self->instances[ch],
                                     number_of_samples, input[ch],
                                     output[ch])) {
      self->last_error = specbleach_denoiser_get_last_error(
          (specbleach_denoiser*)self->instances[ch]);
      return false;
    }
  }

  self->last_error = SPECBLEACH_OK;
  return true;
}

uint32_t specbleach_stereo_get_channel_count(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  return self ? self->channels : 0;
}

uint32_t specbleach_stereo_get_sample_rate(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  if (!self || self->channels == 0) {
    return 0;
  }
  return specbleach_denoiser_get_sample_rate(
      (const specbleach_denoiser*)self->instances[0]);
}

uint32_t specbleach_stereo_get_frame_size(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  if (!self || self->channels == 0) {
    return 0;
  }
  return specbleach_denoiser_get_frame_size(
      (const specbleach_denoiser*)self->instances[0]);
}

uint32_t specbleach_stereo_get_fft_size(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  if (!self || self->channels == 0) {
    return 0;
  }
  return specbleach_denoiser_get_fft_size(
      (const specbleach_denoiser*)self->instances[0]);
}

uint32_t specbleach_stereo_get_hop_size(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  if (!self || self->channels == 0) {
    return 0;
  }
  return specbleach_denoiser_get_hop_size(
      (const specbleach_denoiser*)self->instances[0]);
}

uint32_t specbleach_stereo_get_latency(const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;

  if (!self || self->channels == 0) {
    return 0;
  }

  return specbleach_denoiser_get_latency(
      (const specbleach_denoiser*)self->instances[0]);
}

bool specbleach_stereo_sync_profiles(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  bool complete = true;

  for (SpecbleachProfileMode mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST;
       mode = (SpecbleachProfileMode)((int)mode + 1)) {
    int reference = -1;
    for (uint32_t ch = 0; ch < self->channels && reference < 0; ++ch) {
      if (specbleach_denoiser_noise_profile_available_for_mode(
              (specbleach_denoiser*)self->instances[ch], mode)) {
        reference = (int)ch;
      }
    }
    if (reference < 0) {
      complete = false;
      continue;
    }

    uint32_t size = specbleach_denoiser_get_noise_profile_size(
        (specbleach_denoiser*)self->instances[reference]);
    uint32_t blocks =
        specbleach_denoiser_get_noise_profile_block_count_for_mode(
            (specbleach_denoiser*)self->instances[reference], mode);
    const float* profile = specbleach_denoiser_get_noise_profile_for_mode(
        (const specbleach_denoiser*)self->instances[reference], mode);

    for (uint32_t ch = 0; ch < self->channels; ++ch) {
      if ((int)ch == reference ||
          specbleach_denoiser_noise_profile_available_for_mode(
              (specbleach_denoiser*)self->instances[ch], mode)) {
        continue;
      }
      if (!profile || !specbleach_denoiser_load_noise_profile_for_mode(
                          (specbleach_denoiser*)self->instances[ch], profile,
                          size, blocks, mode)) {
        complete = false;
      }
    }
  }

  return complete;
}

void specbleach_stereo_reset_profiles(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    specbleach_denoiser_reset_noise_profile(
        (specbleach_denoiser*)self->instances[ch]);
  }
  self->last_error = SPECBLEACH_OK;
}

bool specbleach_stereo_reset_dsp_state(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  bool result = true;
  SpecbleachError first_error = SPECBLEACH_OK;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (!specbleach_denoiser_reset_dsp_state(
            (specbleach_denoiser*)self->instances[ch])) {
      if (result) {
        first_error = specbleach_denoiser_get_last_error(
            (specbleach_denoiser*)self->instances[ch]);
      }
      result = false;
    }
  }
  self->last_error = result ? SPECBLEACH_OK : first_error;
  return result;
}

const float* specbleach_stereo_get_noise_profile_for_channel(
    const specbleach_stereo* instance, const uint32_t channel,
    const SpecbleachProfileMode mode) {
  const SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return specbleach_denoiser_get_noise_profile_for_mode(
      (const specbleach_denoiser*)self->instances[channel], mode);
}

bool specbleach_stereo_load_noise_profile_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const float* profile,
    const uint32_t profile_size, const uint32_t block_count,
    const SpecbleachProfileMode mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    if (self) {
      self->last_error = SPECBLEACH_ERR_INVALID_CHANNEL;
    }
    return false;
  }

  const bool ok = specbleach_denoiser_load_noise_profile_for_mode(
      (specbleach_denoiser*)self->instances[channel], profile, profile_size,
      block_count, mode);
  self->last_error = ok ? SPECBLEACH_OK
                        : specbleach_denoiser_get_last_error(
                              (specbleach_denoiser*)self->instances[channel]);
  return ok;
}

bool specbleach_stereo_load_noise_profile_resampled_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const float* profile,
    const uint32_t source_size, const uint32_t block_count,
    const SpecbleachProfileMode mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    if (self) {
      self->last_error = SPECBLEACH_ERR_INVALID_CHANNEL;
    }
    return false;
  }

  const bool ok = specbleach_denoiser_load_noise_profile_resampled(
      (specbleach_denoiser*)self->instances[channel], profile, source_size,
      block_count, mode);
  self->last_error = ok ? SPECBLEACH_OK
                        : specbleach_denoiser_get_last_error(
                              (specbleach_denoiser*)self->instances[channel]);
  return ok;
}

bool specbleach_stereo_profile_available_for_channel(
    const specbleach_stereo* instance, const uint32_t channel,
    const SpecbleachProfileMode mode) {
  const SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return false;
  }

  return specbleach_denoiser_noise_profile_available_for_mode(
      (const specbleach_denoiser*)self->instances[channel], mode);
}

bool specbleach_stereo_has_any_profile_for_channel(
    const specbleach_stereo* instance, const uint32_t channel) {
  const SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return false;
  }

  return specbleach_denoiser_has_any_profile(
      (const specbleach_denoiser*)self->instances[channel]);
}

uint32_t specbleach_stereo_get_profile_block_count_for_channel(
    const specbleach_stereo* instance, const uint32_t channel,
    const SpecbleachProfileMode mode) {
  const SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return 0;
  }

  return specbleach_denoiser_get_noise_profile_block_count_for_mode(
      (const specbleach_denoiser*)self->instances[channel], mode);
}

uint32_t specbleach_stereo_get_noise_profile_size(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;

  if (!self || self->channels == 0) {
    return 0;
  }

  return specbleach_denoiser_get_noise_profile_size(
      (const specbleach_denoiser*)self->instances[0]);
}

bool specbleach_stereo_is_transient_detected(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (specbleach_denoiser_is_transient_detected(
            (const specbleach_denoiser*)self->instances[ch])) {
      return true;
    }
  }
  return false;
}

bool specbleach_stereo_is_transient_detected_for_channel(
    const specbleach_stereo* instance, const uint32_t channel) {
  const SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return false;
  }

  return specbleach_denoiser_is_transient_detected(
      (const specbleach_denoiser*)self->instances[channel]);
}

float specbleach_stereo_get_transient_intensity(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;

  if (!self) {
    return 0.0f;
  }

  float maximum = 0.0f;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    float intensity = specbleach_denoiser_get_transient_intensity(
        (const specbleach_denoiser*)self->instances[ch]);
    if (intensity > maximum) {
      maximum = intensity;
    }
  }
  return maximum;
}

const float* specbleach_stereo_get_tonal_mask_for_channel(
    const specbleach_stereo* instance, const uint32_t channel) {
  const SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return specbleach_denoiser_get_tonal_mask(
      (const specbleach_denoiser*)self->instances[channel]);
}

const float* specbleach_stereo_get_active_noise_profile_for_channel(
    const specbleach_stereo* instance, const uint32_t channel) {
  const SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return specbleach_denoiser_get_active_noise_profile(
      (const specbleach_denoiser*)self->instances[channel]);
}

uint32_t specbleach_stereo_get_tonal_peaks_for_channel(
    const specbleach_stereo* instance, const uint32_t channel,
    float* peak_freqs_hz, const uint32_t max_peaks) {
  const SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return 0;
  }

  return specbleach_denoiser_get_tonal_peaks(
      (const specbleach_denoiser*)self->instances[channel], peak_freqs_hz,
      max_peaks);
}

SpecbleachError specbleach_stereo_get_last_error(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  return self ? self->last_error : SPECBLEACH_ERR_NULL_ARG;
}
