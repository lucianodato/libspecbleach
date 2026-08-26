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

/**
 * @file specbleach_delay_line.c
 * @brief Multi-channel single-tap delay line implementation.
 */

#include "specbleach_delay_line.h"

#include <stdint.h>
#include <stdlib.h>

struct specbleach_delay_line {
  float* buffer; /* channels * ring_capacity samples */
  uint32_t ring_capacity;
  uint32_t max_delay_samples;
  uint32_t delay_samples;
  uint32_t channels;
  uint32_t write_pos; /* absolute sample counter per stream */
};

/**
 * Initializes a multi-channel delay line.
 * @param max_delay_samples Maximum supported delay in samples.
 * @param channels Number of independent audio channels.
 * @return A configured delay line, or NULL if the input is invalid or allocation fails.
 */
specbleach_delay_line* specbleach_delay_line_initialize(
    const uint32_t max_delay_samples, const uint32_t channels) {
  if (channels == 0) {
    return NULL;
  }

  specbleach_delay_line* instance = calloc(1, sizeof(specbleach_delay_line));
  if (!instance) {
    return NULL;
  }

  if (max_delay_samples == UINT32_MAX) {
    free(instance);
    return NULL;
  }

  /* +1 slot so read/write pointers never collide at max delay */
  instance->ring_capacity = max_delay_samples + 1U;
  if ((size_t)channels > SIZE_MAX / sizeof(float) / instance->ring_capacity) {
    free(instance);
    return NULL;
  }
  instance->buffer =
      calloc((size_t)channels * instance->ring_capacity, sizeof(float));
  if (!instance->buffer) {
    free(instance);
    return NULL;
  }

  instance->max_delay_samples = max_delay_samples;
  instance->delay_samples = 0U;
  instance->channels = channels;

  return instance;
}

/**
 * Releases a delay-line instance and its allocated sample buffer.
 * @param instance Delay-line instance to release; may be NULL.
 */
void specbleach_delay_line_free(specbleach_delay_line* instance) {
  if (!instance) {
    return;
  }
  free(instance->buffer);
  free(instance);
}

/**
 * Sets the delay applied by the delay line.
 *
 * @param instance Delay-line instance to configure.
 * @param delay_samples Number of samples to delay.
 * @return true if the delay was set, false if the instance is NULL or the delay exceeds the maximum.
 */
bool specbleach_delay_line_set_delay(specbleach_delay_line* instance,
                                     const uint32_t delay_samples) {
  if (!instance || delay_samples > instance->max_delay_samples) {
    return false;
  }
  instance->delay_samples = delay_samples;
  return true;
}

/**
 * Gets the currently configured delay in samples.
 *
 * @param instance Delay-line instance, or NULL.
 * @return Configured delay in samples, or 0 if instance is NULL.
 */
uint32_t specbleach_delay_line_get_delay(specbleach_delay_line* instance) {
  return instance ? instance->delay_samples : 0U;
}

/**
 * Gets the maximum configurable delay.
 * @param instance Delay-line instance to query.
 * @return Maximum delay in samples, or 0 if instance is NULL.
 */
uint32_t specbleach_delay_line_get_max_delay(specbleach_delay_line* instance) {
  return instance ? instance->max_delay_samples : 0U;
}

/**
 * Gets the number of audio channels configured for a delay line.
 * @param instance Delay-line instance.
 * @returns The configured channel count, or 0 if instance is NULL.
 */
uint32_t specbleach_delay_line_get_channels(specbleach_delay_line* instance) {
  return instance ? instance->channels : 0U;
}

/**
 * Processes audio samples through the configured delay line.
 *
 * @param num_samples Number of samples to process per channel.
 * @param in Input sample buffers, one per channel.
 * @param out Output sample buffers, one per channel.
 * @return `true` if processing succeeds, `false` if the instance or any buffer pointer is NULL.
 */
bool specbleach_delay_line_process(specbleach_delay_line* instance,
                                   const uint32_t num_samples,
                                   const float* const* in, float* const* out) {
  if (!instance || !in || !out) {
    return false;
  }

  const uint32_t delay = instance->delay_samples;
  const uint32_t capacity = instance->ring_capacity;

  for (uint32_t ch = 0U; ch < instance->channels; ++ch) {
    if (!in[ch] || !out[ch]) {
      return false;
    }
    float* ring = &instance->buffer[(size_t)ch * capacity];
    uint32_t pos = instance->write_pos;
    for (uint32_t s = 0U; s < num_samples; ++s) {
      /* Write first so delay == 0 degrades to a well-defined pass-through */
      ring[pos] = in[ch][s];
      const uint32_t read_pos =
          pos >= delay ? pos - delay : pos + capacity - delay;
      out[ch][s] = ring[read_pos];
      pos = pos + 1U == capacity ? 0U : pos + 1U;
    }
  }

  instance->write_pos += num_samples;
  if (instance->write_pos >= capacity) {
    instance->write_pos %= capacity;
  }

  return true;
}
