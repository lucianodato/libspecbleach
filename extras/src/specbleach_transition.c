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
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef enum TransitionPhase {
  TRANSITION_IDLE = 0,
  TRANSITION_FADING = 1,
} TransitionPhase;

typedef struct specbleach_transition { // NOLINT(readability-identifier-naming)
  uint32_t sample_rate;
  uint32_t channels;

  TransitionPhase phase;
  uint32_t latency_to;
  float fade_progress; /* 0..1 within FADING */
  float fade_step;     /* progress increment per sample */
} SpecbleachTransitionState;

/* Half circle in radians; equal-power fades traverse a quarter circle. */
#define TRANSITION_HALF_PI 1.5707963267948966f

specbleach_transition* specbleach_transition_initialize(
    const uint32_t sample_rate, const uint32_t max_block_size,
    const uint32_t channels) {
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
  self->phase = TRANSITION_IDLE;

  return self;
}

void specbleach_transition_free(specbleach_transition* instance) {
  SpecbleachTransitionState* self = instance;

  if (!self) {
    return;
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

  const float fade_samples =
      (float)self->sample_rate * (SPECBLEACH_TRANSITION_FADE_TIME_MS / 1000.0f);

  self->latency_to = latency_to;
  self->fade_step = fade_samples >= 1.0f ? 1.0f / fade_samples : 1.0f;
  self->fade_progress = 0.0f;
  self->phase =
      latency_to != latency_from ? TRANSITION_FADING : TRANSITION_IDLE;
  return true;
}

bool specbleach_transition_process(specbleach_transition* instance,
                                   const uint32_t number_of_samples,
                                   const float** from_out, const float** to_out,
                                   float** blended) {
  SpecbleachTransitionState* self = instance;

  if (!self || number_of_samples == 0 || !to_out || !blended) {
    return false;
  }
  if (!from_out && self->phase == TRANSITION_FADING) {
    return false;
  }

  if (self->phase == TRANSITION_IDLE) {
    for (uint32_t ch = 0; ch < self->channels; ++ch) {
      if (blended[ch] != to_out[ch]) {
        memcpy(blended[ch], to_out[ch],
               (size_t)number_of_samples * sizeof(float));
      }
    }
    return true;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    float progress = self->fade_progress;

    for (uint32_t s = 0; s < number_of_samples; ++s) {
      const float clamped = progress > 1.0f ? 1.0f : progress;
      const float angle = clamped * TRANSITION_HALF_PI;
      const float w_from = cosf(angle);
      const float w_to = sinf(angle);

      blended[ch][s] = (w_from * from_out[ch][s]) + (w_to * to_out[ch][s]);
      progress += self->fade_step;
    }
  }

  self->fade_progress += (float)number_of_samples * self->fade_step;
  if (self->fade_progress >= 1.0f) {
    self->fade_progress = 1.0f;
    self->phase = TRANSITION_IDLE;
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
