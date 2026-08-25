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

#include "specbleach_transition.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef enum TransitionPhase {
  TRANSITION_IDLE = 0,
  TRANSITION_FADING = 1,
  TRANSITION_SLEWING = 2,
} TransitionPhase;

typedef struct specbleach_transition { // NOLINT(readability-identifier-naming)
  uint32_t sample_rate;
  uint32_t channels;
  uint32_t max_alignment_delay;

  TransitionPhase phase;
  uint32_t latency_from;
  uint32_t latency_to;
  uint32_t alignment_delay;
  float fade_progress;
  float fade_step;
  float slew_progress;
  float slew_step;

  /*
   * History of the lower-latency family's live output. Both fade directions
   * blend through an aligned (delayed) copy of that family so the emitted
   * timeline never jumps; after landing on the shorter engine, the slew
   * phase slides the remaining delay back out between two versions of the
   * same signal.
   */
  float** history;
  size_t history_capacity;
  size_t history_write_pos;
} SpecbleachTransitionState;

/* Half circle in radians; equal-power fades traverse a quarter circle. */
#define TRANSITION_HALF_PI 1.5707963267948966f

specbleach_transition* specbleach_transition_initialize(
    const uint32_t sample_rate, const uint32_t max_block_size,
    const uint32_t channels, const uint32_t max_alignment_delay) {
  if (sample_rate == 0 || max_block_size == 0 || channels == 0) {
    return NULL;
  }

  SpecbleachTransitionState* self =
      (SpecbleachTransitionState*)calloc(1U, sizeof(SpecbleachTransitionState));
  if (!self) {
    return NULL;
  }

  self->sample_rate = sample_rate;
  self->channels = channels;
  self->max_alignment_delay = max_alignment_delay;
  self->history_capacity = (size_t)max_block_size + (size_t)max_alignment_delay;
  self->phase = TRANSITION_IDLE;

  self->history = (float**)calloc((size_t)channels * sizeof(float*), 1U);
  if (!self->history) {
    free(self);
    return NULL;
  }
  for (uint32_t ch = 0; ch < channels; ++ch) {
    self->history[ch] =
        (float*)calloc(self->history_capacity * sizeof(float), 1U);
    if (!self->history[ch]) {
      for (uint32_t free_ch = 0; free_ch < ch; ++free_ch) {
        free(self->history[free_ch]);
      }
      free((void*)self->history);
      free(self);
      return NULL;
    }
  }

  return self;
}

void specbleach_transition_free(specbleach_transition* instance) {
  SpecbleachTransitionState* self = instance;

  if (!self) {
    return;
  }

  if (self->history) {
    for (uint32_t ch = 0; ch < self->channels; ++ch) {
      free(self->history[ch]);
    }
    free((void*)self->history);
  }

  free(self);
}

bool specbleach_transition_begin(specbleach_transition* instance,
                                 const uint32_t latency_from,
                                 const uint32_t latency_to) {
  SpecbleachTransitionState* self = instance;

  if (!self) {
    return false;
  }

  const uint32_t diff = latency_to > latency_from ? latency_to - latency_from
                                                  : latency_from - latency_to;
  if (diff > self->max_alignment_delay) {
    return false; /* Initialize with sufficient headroom. */
  }

  const float fade_samples =
      (float)self->sample_rate * (SPECBLEACH_TRANSITION_FADE_TIME_MS / 1000.0f);
  /* Sliding the alignment delay back out needs at least its own length. */
  const float slew_floor_ms = (float)diff / (float)self->sample_rate * 1000.0f;
  const float slew_ms = slew_floor_ms > SPECBLEACH_TRANSITION_FADE_TIME_MS
                            ? slew_floor_ms
                            : SPECBLEACH_TRANSITION_FADE_TIME_MS;
  const float slew_samples = (float)self->sample_rate * (slew_ms / 1000.0f);

  /* History persists across transitions on purpose: later fades read real
   * audio from it instead of silence. */
  self->latency_from = latency_from;
  self->latency_to = latency_to;
  self->alignment_delay = diff;
  self->fade_step = fade_samples >= 1.0f ? 1.0f / fade_samples : 1.0f;
  self->slew_step = slew_samples >= 1.0f ? 1.0f / slew_samples : 1.0f;
  self->fade_progress = 0.0f;
  self->slew_progress = 0.0f;
  self->phase =
      latency_to != latency_from ? TRANSITION_FADING : TRANSITION_IDLE;
  return true;
}

bool specbleach_transition_process(specbleach_transition* instance,
                                   const uint32_t number_of_samples,
                                   const float** from_out, const float** to_out,
                                   float** blended) {
  SpecbleachTransitionState* self = instance;

  if (!self || number_of_samples == 0 || !to_out || !blended ||
      (!from_out && self->phase != TRANSITION_IDLE)) {
    return false;
  }
  assert(self->history != NULL);

  const bool to_higher = self->latency_to > self->latency_from;

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    float progress = self->fade_progress;
    float slew = self->slew_progress;

    for (uint32_t s = 0; s < number_of_samples; ++s) {
      float out_sample;

      if (self->phase == TRANSITION_FADING) {
        /* Whichever side has the shorter latency is read through the
         * alignment tap so both blend sides share the longer time origin. */
        const size_t tap = self->alignment_delay;
        const float angle = progress > 1.0f ? TRANSITION_HALF_PI
                                            : progress * TRANSITION_HALF_PI;

        if (to_higher) {
          const float source_sample = from_out[ch][s];
          self->history[ch][self->history_write_pos] = source_sample;
          const size_t read_pos =
              (self->history_write_pos + self->history_capacity - tap) %
              self->history_capacity;
          out_sample = (cosf(angle) * self->history[ch][read_pos]) +
                       (sinf(angle) * to_out[ch][s]);
        } else {
          const float target_sample = to_out[ch][s];
          self->history[ch][self->history_write_pos] = target_sample;
          const size_t read_pos =
              (self->history_write_pos + self->history_capacity - tap) %
              self->history_capacity;
          out_sample = (cosf(angle) * from_out[ch][s]) +
                       (sinf(angle) * self->history[ch][read_pos]);
        }
        progress += self->fade_step;
      } else if (self->phase == TRANSITION_SLEWING) {
        /* Landed on the shorter-latency engine: crossfade from its aligned
         * copy to its direct output, sliding the timeline smoothly. */
        const float target_sample = to_out[ch][s];
        self->history[ch][self->history_write_pos] = target_sample;
        const size_t read_pos =
            (self->history_write_pos + self->history_capacity -
             (size_t)self->alignment_delay) %
            self->history_capacity;
        const float angle =
            slew > 1.0f ? TRANSITION_HALF_PI : slew * TRANSITION_HALF_PI;
        out_sample = (cosf(angle) * self->history[ch][read_pos]) +
                     (sinf(angle) * target_sample);
        slew += self->slew_step;
      } else {
        out_sample = to_out[ch][s];
      }

      blended[ch][s] = out_sample;
      self->history_write_pos =
          (self->history_write_pos + 1U) % self->history_capacity;
    }
  }

  if (self->phase == TRANSITION_FADING) {
    self->fade_progress += (float)number_of_samples * self->fade_step;
    if (self->fade_progress >= 1.0f) {
      self->fade_progress = 1.0f;
      self->phase = to_higher ? TRANSITION_IDLE : TRANSITION_SLEWING;
    }
  } else if (self->phase == TRANSITION_SLEWING) {
    self->slew_progress += (float)number_of_samples * self->slew_step;
    if (self->slew_progress >= 1.0f) {
      self->slew_progress = 1.0f;
      self->phase = TRANSITION_IDLE;
    }
  }

  return true;
}

bool specbleach_transition_active(const specbleach_transition* instance) {
  const SpecbleachTransitionState* self = instance;

  return self ? self->phase != TRANSITION_IDLE : false;
}

uint32_t specbleach_transition_get_latency(
    const specbleach_transition* instance) {
  const SpecbleachTransitionState* self = instance;

  if (!self) {
    return 0;
  }

  return self->latency_to;
}
