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
} SpecbleachStereoState;

specbleach_stereo* specbleach_stereo_initialize(const uint32_t sample_rate,
                                                const float frame_size,
                                                const uint32_t channels) {
  if (sample_rate == 0 || frame_size <= 0.0f || channels == 0) {
    return NULL;
  }

  SpecbleachStereoState* self =
      (SpecbleachStereoState*)calloc(1U, sizeof(SpecbleachStereoState));
  if (!self) {
    return NULL;
  }

  self->channels = channels;
  self->instances = calloc((size_t)channels, sizeof(void*));
  if (!self->instances) {
    specbleach_stereo_free(self);
    return NULL;
  }

  for (uint32_t ch = 0; ch < channels; ++ch) {
    self->instances[ch] =
        specbleach_denoiser_initialize(sample_rate, frame_size);
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
    const uint32_t parameters_size) {
  SpecbleachStereoState* self = instance;

  if (!self || !parameters ||
      parameters_size != sizeof(SpecbleachDenoiserParameters)) {
    return false;
  }

  bool result = true;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (!specbleach_denoiser_load_parameters(
            (specbleach_denoiser*)self->instances[ch], parameters,
            parameters_size)) {
      result = false;
    }
  }
  return result;
}

bool specbleach_stereo_process(specbleach_stereo* instance,
                               const uint32_t number_of_samples,
                               const float** input, float** output) {
  SpecbleachStereoState* self = instance;

  if (!self || number_of_samples == 0 || !input || !output) {
    return false;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (!specbleach_denoiser_process((specbleach_denoiser*)self->instances[ch],
                                     number_of_samples, input[ch],
                                     output[ch])) {
      return false;
    }
  }

  return true;
}

uint32_t specbleach_stereo_get_channel_count(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  return self ? self->channels : 0;
}

uint32_t specbleach_stereo_get_latency(const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;

  if (!self || self->channels == 0) {
    return 0;
  }

  return specbleach_denoiser_get_latency(
      (specbleach_denoiser*)self->instances[0]);
}

bool specbleach_stereo_sync_profiles(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  bool complete = true;

  for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
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
    float* profile = specbleach_denoiser_get_noise_profile_for_mode(
        (specbleach_denoiser*)self->instances[reference], mode);

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

bool specbleach_stereo_migrate_profiles_from(specbleach_stereo* instance,
                                             const specbleach_stereo* source) {
  SpecbleachStereoState* self = instance;
  const SpecbleachStereoState* origin = source;

  if (!self || !origin || self->channels != origin->channels) {
    return false;
  }

  if (specbleach_stereo_get_noise_profile_size(self) !=
      specbleach_stereo_get_noise_profile_size(origin)) {
    return false;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    bool channel_migrated = false;
    for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
         mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
      if (!specbleach_stereo_profile_available_for_channel(
              (specbleach_stereo*)origin, ch, mode)) {
        continue;
      }
      float* profile = specbleach_stereo_get_noise_profile_for_channel(
          (specbleach_stereo*)origin, ch, mode);
      if (!profile) {
        continue;
      }
      if (specbleach_stereo_load_noise_profile_for_channel(
              instance, ch, profile,
              specbleach_stereo_get_noise_profile_size(instance),
              specbleach_stereo_get_profile_block_count_for_channel(
                  (specbleach_stereo*)origin, ch, mode),
              mode)) {
        channel_migrated = true;
      }
    }
    if (!channel_migrated) {
      return false;
    }
  }

  return true;
}

bool specbleach_stereo_reset_profiles(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  bool result = true;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (!specbleach_denoiser_reset_noise_profile(
            (specbleach_denoiser*)self->instances[ch])) {
      result = false;
    }
  }
  return result;
}

float* specbleach_stereo_get_noise_profile_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const int mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return specbleach_denoiser_get_noise_profile_for_mode(
      (specbleach_denoiser*)self->instances[channel], mode);
}

bool specbleach_stereo_load_noise_profile_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const float* profile,
    const uint32_t profile_size, const uint32_t block_count, const int mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return false;
  }

  return specbleach_denoiser_load_noise_profile_for_mode(
      (specbleach_denoiser*)self->instances[channel], profile, profile_size,
      block_count, mode);
}

bool specbleach_stereo_profile_available_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const int mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return false;
  }

  return specbleach_denoiser_noise_profile_available_for_mode(
      (specbleach_denoiser*)self->instances[channel], mode);
}

uint32_t specbleach_stereo_get_profile_block_count_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const int mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return 0;
  }

  return specbleach_denoiser_get_noise_profile_block_count_for_mode(
      (specbleach_denoiser*)self->instances[channel], mode);
}

uint32_t specbleach_stereo_get_noise_profile_size(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;

  if (!self || self->channels == 0) {
    return 0;
  }

  return specbleach_denoiser_get_noise_profile_size(
      (specbleach_denoiser*)self->instances[0]);
}

bool specbleach_stereo_is_transient_detected(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (specbleach_denoiser_is_transient_detected(
            (specbleach_denoiser*)self->instances[ch])) {
      return true;
    }
  }
  return false;
}

float specbleach_stereo_get_transient_intensity(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return 0.0f;
  }

  float maximum = 0.0f;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    float intensity = specbleach_denoiser_get_transient_intensity(
        (specbleach_denoiser*)self->instances[ch]);
    if (intensity > maximum) {
      maximum = intensity;
    }
  }
  return maximum;
}

const float* specbleach_stereo_get_tonal_mask_for_channel(
    specbleach_stereo* instance, const uint32_t channel) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return specbleach_denoiser_get_tonal_mask(
      (specbleach_denoiser*)self->instances[channel]);
}

const float* specbleach_stereo_get_active_noise_profile_for_channel(
    specbleach_stereo* instance, const uint32_t channel) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return specbleach_denoiser_get_active_noise_profile(
      (specbleach_denoiser*)self->instances[channel]);
}

uint32_t specbleach_stereo_get_tonal_peaks_for_channel(
    specbleach_stereo* instance, const uint32_t channel, float* peak_freqs_hz,
    const uint32_t max_peaks) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return 0;
  }

  return specbleach_denoiser_get_tonal_peaks(
      (specbleach_denoiser*)self->instances[channel], peak_freqs_hz, max_peaks);
}
