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
  SpecbleachStereoEngine engine;
  void** instances; /* one opaque engine handle per channel */
} SpecbleachStereoState;

/**
 * Applies spectral denoiser parameters to every channel.
 *
 * @param instance Stereo denoiser instance to update.
 * @param parameters Parameters to apply to each channel.
 * @return `true` if all channels accept the parameters, `false` otherwise.
 */
static bool stereo_load_1d(specbleach_stereo* instance,
                           const SpecbleachDenoiserParameters* parameters) {
  SpecbleachStereoState* self = instance;
  bool result = true;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (!specbleach_denoiser_load_parameters(
            (specbleach_denoiser*)self->instances[ch], parameters,
            sizeof(*parameters))) {
      result = false;
    }
  }
  return result;
}

/**
 * Loads 2D denoiser parameters for every channel.
 *
 * @param instance Stereo denoiser instance.
 * @param parameters 2D denoiser parameters to apply.
 * @return `true` if the parameters are loaded successfully for all channels, `false` otherwise.
 */
static bool stereo_load_2d(specbleach_stereo* instance,
                           const Specbleach2DDenoiserParameters* parameters) {
  SpecbleachStereoState* self = instance;
  bool result = true;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    if (!specbleach_2d_load_parameters(
            (specbleach_2d_denoiser*)self->instances[ch], parameters,
            sizeof(*parameters))) {
      result = false;
    }
  }
  return result;
}

/**
 * Creates a stereo denoiser with one processing instance per channel.
 *
 * @param sample_rate Audio sample rate in hertz.
 * @param frame_size Processing frame size.
 * @param channels Number of audio channels.
 * @param engine Denoising engine to use for each channel.
 * @return A newly initialized stereo denoiser, or NULL if the arguments are
 *         invalid or initialization fails.
 */
specbleach_stereo* specbleach_stereo_initialize(
    const uint32_t sample_rate, const float frame_size, const uint32_t channels,
    const SpecbleachStereoEngine engine) {
  if (sample_rate == 0 || frame_size <= 0.0f || channels == 0 ||
      (engine != SPECBLEACH_STEREO_ENGINE_SPECTRAL &&
       engine != SPECBLEACH_STEREO_ENGINE_NLM_2D)) {
    return NULL;
  }

  SpecbleachStereoState* self =
      (SpecbleachStereoState*)calloc(1U, sizeof(SpecbleachStereoState));
  if (!self) {
    return NULL;
  }

  self->channels = channels;
  self->engine = engine;
  self->instances = calloc((size_t)channels, sizeof(void*));
  if (!self->instances) {
    specbleach_stereo_free(self);
    return NULL;
  }

  for (uint32_t ch = 0; ch < channels; ++ch) {
    if (engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL) {
      self->instances[ch] =
          specbleach_denoiser_initialize(sample_rate, frame_size);
    } else {
      self->instances[ch] = specbleach_2d_initialize(sample_rate, frame_size);
    }
    if (!self->instances[ch]) {
      specbleach_stereo_free(self);
      return NULL;
    }
  }

  return self;
}

/**
 * Releases a stereo denoiser instance and all associated channel engines.
 *
 * @param instance Stereo denoiser instance to release. May be NULL.
 */
void specbleach_stereo_free(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return;
  }

  if (self->instances) {
    for (uint32_t ch = 0; ch < self->channels; ++ch) {
      if (self->instances[ch]) {
        if (self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL) {
          specbleach_denoiser_free((specbleach_denoiser*)self->instances[ch]);
        } else {
          specbleach_2d_free((specbleach_2d_denoiser*)self->instances[ch]);
        }
      }
    }
    free(self->instances);
  }

  free(self);
}

/**
 * Loads spectral denoiser parameters for all channels.
 *
 * @param instance Stereo denoiser instance using the spectral engine.
 * @param parameters Spectral denoiser parameters to apply.
 * @param parameters_size Size of the parameters structure.
 * @returns true if the parameters are loaded successfully, false otherwise.
 */
bool specbleach_stereo_load_parameters_1d(
    specbleach_stereo* instance, const SpecbleachDenoiserParameters* parameters,
    const uint32_t parameters_size) {
  SpecbleachStereoState* self = instance;

  if (!self || !parameters ||
      parameters_size != sizeof(SpecbleachDenoiserParameters) ||
      self->engine != SPECBLEACH_STEREO_ENGINE_SPECTRAL) {
    return false;
  }

  return stereo_load_1d(self, parameters);
}

/**
 * Applies 2D NLM denoiser parameters to all channels.
 *
 * @param parameters 2D denoiser configuration.
 * @param parameters_size Size of the parameter structure in bytes.
 * @return `true` if the parameters are applied successfully, `false` otherwise.
 */
bool specbleach_stereo_load_parameters_2d(
    specbleach_stereo* instance,
    const Specbleach2DDenoiserParameters* parameters,
    const uint32_t parameters_size) {
  SpecbleachStereoState* self = instance;

  if (!self || !parameters ||
      parameters_size != sizeof(Specbleach2DDenoiserParameters) ||
      self->engine != SPECBLEACH_STEREO_ENGINE_NLM_2D) {
    return false;
  }

  return stereo_load_2d(self, parameters);
}

/**
 * Processes audio samples for all channels.
 * @param instance Stereo denoiser instance.
 * @param number_of_samples Number of samples to process per channel.
 * @param input Input audio buffers indexed by channel.
 * @param output Output audio buffers indexed by channel.
 * @return true if all channels are processed successfully, false otherwise.
 */
bool specbleach_stereo_process(specbleach_stereo* instance,
                               const uint32_t number_of_samples,
                               const float** input, float** output) {
  SpecbleachStereoState* self = instance;

  if (!self || number_of_samples == 0 || !input || !output) {
    return false;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    bool ok = false;
    if (self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL) {
      ok =
          specbleach_denoiser_process((specbleach_denoiser*)self->instances[ch],
                                      number_of_samples, input[ch], output[ch]);
    } else {
      ok = specbleach_2d_process((specbleach_2d_denoiser*)self->instances[ch],
                                 number_of_samples, input[ch], output[ch]);
    }
    if (!ok) {
      return false;
    }
  }

  return true;
}

/**
 * Gets the number of audio channels configured for a stereo denoiser instance.
 *
 * @param instance Denoiser instance to query.
 * @return The configured channel count, or 0 if instance is NULL.
 */
uint32_t specbleach_stereo_get_channel_count(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  return self ? self->channels : 0;
}

/**
 * Gets the denoising engine configured for a stereo instance.
 *
 * @param instance Stereo denoising instance.
 * @return The configured engine, or the spectral engine when instance is NULL.
 */
SpecbleachStereoEngine specbleach_stereo_get_engine(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;
  return self ? self->engine : SPECBLEACH_STEREO_ENGINE_SPECTRAL;
}

/**
 * Gets the processing latency of the stereo denoiser.
 *
 * @param instance Stereo denoiser instance.
 * @return Processing latency in samples, or 0 if the instance is invalid.
 */
uint32_t specbleach_stereo_get_latency(const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;

  if (!self || self->channels == 0) {
    return 0;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_get_latency(
                   (specbleach_denoiser*)self->instances[0])
             : specbleach_2d_get_latency(
                   (specbleach_2d_denoiser*)self->instances[0]);
}

/**
 * Synchronizes available noise profiles across all channels for each profile mode.
 *
 * @param instance Stereo denoiser instance whose channel profiles are synchronized.
 * @return `true` if every profile mode is available on at least one channel and is
 *         successfully copied to channels without that profile, `false` otherwise.
 */
static bool stereo_sync_profiles_1d(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;
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

/**
 * Synchronizes available 2D noise profiles across all channels for each profile mode.
 *
 * @param instance Stereo denoiser instance whose channel profiles are synchronized.
 * @return `true` if every profile mode is available or successfully copied to all channels, `false` otherwise.
 */
static bool stereo_sync_profiles_2d(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;
  bool complete = true;

  for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
    int reference = -1;
    for (uint32_t ch = 0; ch < self->channels && reference < 0; ++ch) {
      if (specbleach_2d_noise_profile_available_for_mode(
              (specbleach_2d_denoiser*)self->instances[ch], mode)) {
        reference = (int)ch;
      }
    }
    if (reference < 0) {
      complete = false;
      continue;
    }

    uint32_t size = specbleach_2d_get_noise_profile_size(
        (specbleach_2d_denoiser*)self->instances[reference]);
    uint32_t blocks = specbleach_2d_get_noise_profile_block_count_for_mode(
        (specbleach_2d_denoiser*)self->instances[reference], mode);
    float* profile = specbleach_2d_get_noise_profile_for_mode(
        (specbleach_2d_denoiser*)self->instances[reference], mode);

    for (uint32_t ch = 0; ch < self->channels; ++ch) {
      if ((int)ch == reference ||
          specbleach_2d_noise_profile_available_for_mode(
              (specbleach_2d_denoiser*)self->instances[ch], mode)) {
        continue;
      }
      if (!profile || !specbleach_2d_load_noise_profile_for_mode(
                          (specbleach_2d_denoiser*)self->instances[ch], profile,
                          size, blocks, mode)) {
        complete = false;
      }
    }
  }

  return complete;
}

/**
 * Synchronizes available noise profiles across all channels.
 * @param instance Stereo denoiser instance.
 * @return `true` if synchronization succeeds, `false` otherwise.
 */
bool specbleach_stereo_sync_profiles(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? stereo_sync_profiles_1d(self)
             : stereo_sync_profiles_2d(self);
}

/**
 * Migrates available noise profiles from a source stereo instance.
 * @param instance Destination stereo instance.
 * @param source Stereo instance containing the profiles to migrate.
 * @returns `true` if every channel receives at least one profile, `false` if migration cannot be completed.
 */
bool specbleach_stereo_migrate_profiles_from(specbleach_stereo* instance,
                                             const specbleach_stereo* source) {
  SpecbleachStereoState* self = instance;
  const SpecbleachStereoState* origin = source;

  if (!self || !origin || self->channels != origin->channels ||
      self->engine == origin->engine) {
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

/**
 * Resets the noise profile for every channel.
 *
 * @param instance Stereo denoiser instance whose channel profiles are reset.
 * @return `true` if every channel profile is reset successfully, `false` otherwise.
 */
bool specbleach_stereo_reset_profiles(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  bool result = true;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    bool ok = false;
    if (self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL) {
      ok = specbleach_denoiser_reset_noise_profile(
          (specbleach_denoiser*)self->instances[ch]);
    } else {
      ok = specbleach_2d_reset_noise_profile(
          (specbleach_2d_denoiser*)self->instances[ch]);
    }
    if (!ok) {
      result = false;
    }
  }
  return result;
}

/**
 * Retrieves a noise profile for a channel and profile mode.
 * @param instance Stereo denoiser instance.
 * @param channel Zero-based channel index.
 * @param mode Noise profile mode.
 * @return Pointer to the noise profile, or NULL if the instance or channel is invalid or the profile is unavailable.
 */
float* specbleach_stereo_get_noise_profile_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const int mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_get_noise_profile_for_mode(
                   (specbleach_denoiser*)self->instances[channel], mode)
             : specbleach_2d_get_noise_profile_for_mode(
                   (specbleach_2d_denoiser*)self->instances[channel], mode);
}

/**
 * Loads a noise profile for a specific channel and profile mode.
 *
 * @param instance Stereo denoiser instance.
 * @param channel Channel index.
 * @param profile Noise profile data.
 * @param profile_size Number of values in the profile.
 * @param block_count Number of blocks represented by the profile.
 * @param mode Profile mode.
 * @return `true` if the profile is loaded successfully, `false` otherwise.
 */
bool specbleach_stereo_load_noise_profile_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const float* profile,
    const uint32_t profile_size, const uint32_t block_count, const int mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return false;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_load_noise_profile_for_mode(
                   (specbleach_denoiser*)self->instances[channel], profile,
                   profile_size, block_count, mode)
             : specbleach_2d_load_noise_profile_for_mode(
                   (specbleach_2d_denoiser*)self->instances[channel], profile,
                   profile_size, block_count, mode);
}

/**
 * Determines whether a noise profile is available for a channel and mode.
 * @param instance Stereo denoiser instance.
 * @param channel Channel index.
 * @param mode Noise profile mode.
 * @return true if a noise profile is available, false otherwise.
 */
bool specbleach_stereo_profile_available_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const int mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return false;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_noise_profile_available_for_mode(
                   (specbleach_denoiser*)self->instances[channel], mode)
             : specbleach_2d_noise_profile_available_for_mode(
                   (specbleach_2d_denoiser*)self->instances[channel], mode);
}

/**
 * Gets the number of blocks available for a noise profile on a channel.
 *
 * @param channel Channel whose profile block count is requested.
 * @param mode Noise profile mode.
 * @returns The number of profile blocks, or 0 for an invalid instance or channel.
 */
uint32_t specbleach_stereo_get_profile_block_count_for_channel(
    specbleach_stereo* instance, const uint32_t channel, const int mode) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return 0;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_get_noise_profile_block_count_for_mode(
                   (specbleach_denoiser*)self->instances[channel], mode)
             : specbleach_2d_get_noise_profile_block_count_for_mode(
                   (specbleach_2d_denoiser*)self->instances[channel], mode);
}

/**
 * Gets the noise profile size for the selected denoising engine.
 * @param instance Stereo denoising instance.
 * @return The noise profile size in bytes, or 0 for an invalid instance.
 */
uint32_t specbleach_stereo_get_noise_profile_size(
    const specbleach_stereo* instance) {
  const SpecbleachStereoState* self = instance;

  if (!self || self->channels == 0) {
    return 0;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_get_noise_profile_size(
                   (specbleach_denoiser*)self->instances[0])
             : specbleach_2d_get_noise_profile_size(
                   (specbleach_2d_denoiser*)self->instances[0]);
}

/**
 * Determines whether any stereo channel has detected a transient.
 * @param instance Stereo denoiser instance to query.
 * @return `true` if a channel has detected a transient, `false` otherwise.
 */
bool specbleach_stereo_is_transient_detected(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return false;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    bool detected = false;
    if (self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL) {
      detected = specbleach_denoiser_is_transient_detected(
          (specbleach_denoiser*)self->instances[ch]);
    } else {
      detected = specbleach_2d_is_transient_detected(
          (specbleach_2d_denoiser*)self->instances[ch]);
    }
    if (detected) {
      return true;
    }
  }
  return false;
}

/**
 * Gets the highest transient intensity detected across all channels.
 *
 * @param instance Stereo denoiser instance.
 * @return Maximum transient intensity, or 0.0 if the instance is invalid.
 */
float specbleach_stereo_get_transient_intensity(specbleach_stereo* instance) {
  SpecbleachStereoState* self = instance;

  if (!self) {
    return 0.0f;
  }

  float maximum = 0.0f;
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    float intensity = 0.0f;
    if (self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL) {
      intensity = specbleach_denoiser_get_transient_intensity(
          (specbleach_denoiser*)self->instances[ch]);
    } else {
      intensity = specbleach_2d_get_transient_intensity(
          (specbleach_2d_denoiser*)self->instances[ch]);
    }
    if (intensity > maximum) {
      maximum = intensity;
    }
  }
  return maximum;
}

/**
 * Retrieves the tonal mask for a channel.
 *
 * @param channel Channel whose tonal mask is requested.
 * @return Pointer to the channel's tonal mask, or NULL if the instance or channel is invalid.
 */
const float* specbleach_stereo_get_tonal_mask_for_channel(
    specbleach_stereo* instance, const uint32_t channel) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_get_tonal_mask(
                   (specbleach_denoiser*)self->instances[channel])
             : specbleach_2d_get_tonal_mask(
                   (specbleach_2d_denoiser*)self->instances[channel]);
}

/**
 * Retrieves the active noise profile for a channel.
 * @param instance Stereo denoiser instance.
 * @param channel Zero-based channel index.
 * @return Pointer to the active noise profile, or NULL if the instance or channel is invalid.
 */
const float* specbleach_stereo_get_active_noise_profile_for_channel(
    specbleach_stereo* instance, const uint32_t channel) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return NULL;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_get_active_noise_profile(
                   (specbleach_denoiser*)self->instances[channel])
             : specbleach_2d_get_active_noise_profile(
                   (specbleach_2d_denoiser*)self->instances[channel]);
}

/**
 * Retrieves tonal peak frequencies for a channel.
 *
 * @param instance Stereo denoiser instance.
 * @param channel Channel index.
 * @param peak_freqs_hz Buffer to receive peak frequencies in hertz.
 * @param max_peaks Maximum number of peak frequencies to retrieve.
 * @return Number of tonal peak frequencies written, or 0 for an invalid instance or channel.
 */
uint32_t specbleach_stereo_get_tonal_peaks_for_channel(
    specbleach_stereo* instance, const uint32_t channel, float* peak_freqs_hz,
    const uint32_t max_peaks) {
  SpecbleachStereoState* self = instance;

  if (!self || channel >= self->channels) {
    return 0;
  }

  return self->engine == SPECBLEACH_STEREO_ENGINE_SPECTRAL
             ? specbleach_denoiser_get_tonal_peaks(
                   (specbleach_denoiser*)self->instances[channel],
                   peak_freqs_hz, max_peaks)
             : specbleach_2d_get_tonal_peaks(
                   (specbleach_2d_denoiser*)self->instances[channel],
                   peak_freqs_hz, max_peaks);
}
