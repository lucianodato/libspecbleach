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
  TRANSITION_SLEWING = 2,
} TransitionPhase;

typedef struct specbleach_transition {
  uint32_t sample_rate;
  uint32_t channels;
  uint32_t max_block_size;

  TransitionPhase phase;
  uint32_t latency_from;
  uint32_t latency_to;
  uint32_t latency_diff; /* |to - from|, alignment delay amount */
  float fade_progress;   /* 0..1 within FADING */
  float fade_step;       /* progress increment per sample */
  float slew_progress;   /* 0..1 within SLEWING */

  /* Per-channel ring buffers delaying the lower-latency stream. */
  float** delay_lines;
  size_t delay_capacity; /* samples per channel */
  size_t delay_write_pos;
} SpecbleachTransitionState;

/* Half circle in radians; equal-power fades traverse a quarter circle. */
#define TRANSITION_HALF_PI 1.5707963267948966f

static void transition_reset_delay(specbleach_transition* instance) {
  SpecbleachTransitionState* self = instance;

  self->delay_write_pos = 0;
  if (self->delay_lines) {
    for (uint32_t ch = 0; ch < self->channels; ++ch) {
      if (self->delay_lines[ch]) {
        memset(self->delay_lines[ch], 0, self->delay_capacity * sizeof(float));
      }
    }
  }
}

static bool transition_ensure_delay_capacity(specbleach_transition* instance,
                                             const size_t required_samples) {
  SpecbleachTransitionState* self = instance;

  if (self->delay_lines && self->delay_capacity >= required_samples) {
    return true;
  }

  const size_t capacity = required_samples > self->max_block_size
                              ? required_samples + self->max_block_size
                              : self->max_block_size * 2U;

  float** lines = (float**)calloc(self->channels, sizeof(float*));
  if (!lines) {
    return false;
  }

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    lines[ch] = (float*)calloc(capacity, sizeof(float));
    if (!lines[ch]) {
      for (uint32_t free_ch = 0; free_ch < ch; ++free_ch) {
        free(lines[free_ch]);
      }
      free(lines);
      return false;
    }
  }
  /* Freshly allocated lines are already zeroed; only the write position
   * needs resetting. Swap AFTER releasing the old buffers. */
  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    free(self->delay_lines ? self->delay_lines[ch] : NULL);
  }
  free(self->delay_lines);

  self->delay_lines = lines;
  self->delay_capacity = capacity;
  self->delay_write_pos = 0;
  return true;
}

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
  self->max_block_size = max_block_size;
  self->phase = TRANSITION_IDLE;

  return self;
}

void specbleach_transition_free(specbleach_transition* instance) {
  SpecbleachTransitionState* self = instance;

  if (!self) {
    return;
  }

  if (self->delay_lines) {
    for (uint32_t ch = 0; ch < self->channels; ++ch) {
      free(self->delay_lines[ch]);
    }
    free(self->delay_lines);
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

  if (diff > 0 && !transition_ensure_delay_capacity(self, (size_t)diff + 1U)) {
    return false;
  }

  const float fade_samples =
      self->sample_rate * (SPECBLEACH_TRANSITION_FADE_TIME_MS / 1000.0f);

  self->latency_from = latency_from;
  self->latency_to = latency_to;
  self->latency_diff = diff;
  self->fade_step = fade_samples >= 1.0f ? 1.0f / fade_samples : 1.0f;
  self->fade_progress = 0.0f;
  self->slew_progress = 0.0f;
  self->phase = diff > 0 || latency_to != latency_from ? TRANSITION_FADING
                                                       : TRANSITION_IDLE;

  transition_reset_delay(self);
  return true;
}

static void transition_write_and_read_delay(
    const SpecbleachTransitionState* self, const uint32_t channel,
    const float input_sample, const size_t delay_amount, float* output_sample) {
  float* line = self->delay_lines[channel];
  const size_t capacity = self->delay_capacity;

  line[self->delay_write_pos] = input_sample;
  const size_t read_pos =
      (self->delay_write_pos + capacity - delay_amount) % capacity;
  *output_sample = line[read_pos];
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
        memcpy(blended[ch], to_out[ch], number_of_samples * sizeof(float));
      }
    }
    return true;
  }

  const size_t delay_amount = (size_t)self->latency_diff;

  for (uint32_t ch = 0; ch < self->channels; ++ch) {
    float progress = self->fade_progress;
    float slew = self->slew_progress;
    size_t write_pos = self->delay_write_pos;

    for (uint32_t s = 0; s < number_of_samples; ++s) {
      if (self->phase == TRANSITION_FADING) {
        const float source_sample = from_out[ch][s];

        float aligned_source = source_sample;
        if (delay_amount > 0) {
          transition_write_and_read_delay(self, ch, source_sample, delay_amount,
                                          &aligned_source);
          write_pos = (write_pos + 1U) % self->delay_capacity;
        }

        const float clamped = progress > 1.0f ? 1.0f : progress;
        const float angle = clamped * TRANSITION_HALF_PI;
        const float w_from = cosf(angle);
        const float w_to = sinf(angle);

        blended[ch][s] = w_from * aligned_source + w_to * to_out[ch][s];
        progress += self->fade_step;
      } else {
        /* SLEWING: ramp the compensation delay on the target stream back
         * to zero so reported and actual latencies converge smoothly. */
        const float target_sample = to_out[ch][s];

        float aligned_target = target_sample;
        if (delay_amount > 0) {
          transition_write_and_read_delay(self, ch, target_sample, delay_amount,
                                          &aligned_target);
          write_pos = (write_pos + 1U) % self->delay_capacity;
        }

        const float clamped = slew > 1.0f ? 1.0f : slew;
        const float angle = clamped * TRANSITION_HALF_PI;
        const float w_aligned = cosf(angle);
        const float w_direct = sinf(angle);

        blended[ch][s] = w_aligned * aligned_target + w_direct * target_sample;
        slew += self->fade_step;
      }
    }
  }

  self->delay_write_pos =
      (self->delay_write_pos + number_of_samples) % self->delay_capacity;

  if (self->phase == TRANSITION_FADING) {
    self->fade_progress += (float)number_of_samples * self->fade_step;
    if (self->fade_progress >= 1.0f) {
      self->fade_progress = 1.0f;
      if (self->latency_to < self->latency_from) {
        /* Landed on the shorter-latency engine: ramp its compensation
         * delay away over the same duration. */
        self->phase = TRANSITION_SLEWING;
        self->slew_progress = 0.0f;
      } else {
        self->phase = TRANSITION_IDLE;
      }
    }
  } else if (self->phase == TRANSITION_SLEWING) {
    self->slew_progress += (float)number_of_samples * self->fade_step;
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
