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

#include "louizou_noise_estimator.h"
#include "../configurations.h"
#include "../utils/general_utils.h"
#include "../utils/spectral_utils.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct FrameSpectrum {
  float* smoothed_spectrum;
  float* local_minimum_spectrum;
  float* speech_present_probability_spectrum;
} FrameSpectrum;

struct LouizouNoiseEstimator {
  uint32_t noise_spectrum_size;
  float noisy_speech_ratio;

  FrameSpectrum* current;
  FrameSpectrum* previous;

  float* minimum_detection_thresholds;
  float* previous_noise_spectrum;
  float* time_frequency_smoothing_constant;
  uint32_t* speech_presence_detection;
  bool is_first_frame;
};

static FrameSpectrum* frame_spectrum_initialize(uint32_t frame_size) {
  FrameSpectrum* self = (FrameSpectrum*)calloc(1U, sizeof(FrameSpectrum));
  if (!self) {
    return NULL;
  }

  self->smoothed_spectrum = (float*)calloc(frame_size, sizeof(float));
  self->local_minimum_spectrum = (float*)calloc(frame_size, sizeof(float));
  self->speech_present_probability_spectrum =
      (float*)calloc(frame_size, sizeof(float));

  if (!self->smoothed_spectrum || !self->local_minimum_spectrum ||
      !self->speech_present_probability_spectrum) {
    if (self->smoothed_spectrum) {
      free(self->smoothed_spectrum);
    }
    if (self->local_minimum_spectrum) {
      free(self->local_minimum_spectrum);
    }
    if (self->speech_present_probability_spectrum) {
      free(self->speech_present_probability_spectrum);
    }
    free(self);
    return NULL;
  }

  (void)initialize_spectrum_with_value(self->local_minimum_spectrum, frame_size,
                                       FLT_MIN);

  return self;
}

static void frame_spectrum_free(FrameSpectrum* self) {
  if (!self) {
    return;
  }
  free(self->smoothed_spectrum);
  free(self->local_minimum_spectrum);
  free(self->speech_present_probability_spectrum);

  free(self);
}

static void compute_auto_thresholds(LouizouNoiseEstimator* self,
                                    uint32_t sample_rate,
                                    uint32_t noise_spectrum_size,
                                    uint32_t fft_size) {
  uint32_t crossover_bin1 =
      freq_to_fft_bin(CROSSOVER_POINT1, sample_rate, fft_size);
  uint32_t crossover_bin2 =
      freq_to_fft_bin(CROSSOVER_POINT2, sample_rate, fft_size);
  for (uint32_t k = 0U; k < noise_spectrum_size; k++) {
    if (k <= crossover_bin1) {
      self->minimum_detection_thresholds[k] = BAND_1_LEVEL;
    }
    if (k > crossover_bin1 && k < crossover_bin2) {
      self->minimum_detection_thresholds[k] = BAND_2_LEVEL;
    }
    if (k >= crossover_bin2) {
      self->minimum_detection_thresholds[k] = BAND_3_LEVEL;
    }
  }
}

static void update_frame_spectums(LouizouNoiseEstimator* self,
                                  const float* noise_spectrum) {
  memcpy(self->previous_noise_spectrum, noise_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->local_minimum_spectrum,
         self->current->local_minimum_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->smoothed_spectrum, self->current->smoothed_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->speech_present_probability_spectrum,
         self->current->speech_present_probability_spectrum,
         sizeof(float) * self->noise_spectrum_size);
}

LouizouNoiseEstimator* louizou_noise_estimator_initialize(
    uint32_t noise_spectrum_size, uint32_t sample_rate, uint32_t fft_size) {
  LouizouNoiseEstimator* self =
      (LouizouNoiseEstimator*)calloc(1U, sizeof(LouizouNoiseEstimator));
  if (!self) {
    return NULL;
  }

  self->noise_spectrum_size = noise_spectrum_size;

  self->minimum_detection_thresholds =
      (float*)calloc(self->noise_spectrum_size, sizeof(float));
  self->time_frequency_smoothing_constant =
      (float*)calloc(self->noise_spectrum_size, sizeof(float));
  self->speech_presence_detection =
      (uint32_t*)calloc(self->noise_spectrum_size, sizeof(uint32_t));
  self->previous_noise_spectrum =
      (float*)calloc(self->noise_spectrum_size, sizeof(float));

  if (!self->minimum_detection_thresholds ||
      !self->time_frequency_smoothing_constant ||
      !self->speech_presence_detection || !self->previous_noise_spectrum) {
    louizou_noise_estimator_free(self);
    return NULL;
  }

  compute_auto_thresholds(self, sample_rate, noise_spectrum_size, fft_size);
  self->current = frame_spectrum_initialize(noise_spectrum_size);
  self->previous = frame_spectrum_initialize(noise_spectrum_size);

  if (!self->current || !self->previous) {
    louizou_noise_estimator_free(self);
    return NULL;
  }

  self->noisy_speech_ratio = 0.F;
  self->is_first_frame = true;

  return self;
}

void louizou_noise_estimator_free(LouizouNoiseEstimator* self) {
  if (!self) {
    return;
  }
  free(self->minimum_detection_thresholds);
  free(self->time_frequency_smoothing_constant);
  free(self->speech_presence_detection);
  free(self->previous_noise_spectrum);

  frame_spectrum_free(self->current);
  frame_spectrum_free(self->previous);

  free(self);
}

bool louizou_noise_estimator_run(LouizouNoiseEstimator* self,
                                 const float* spectrum, float* noise_spectrum) {
  if (!self || !spectrum || !noise_spectrum) {
    return false;
  }

  if (self->is_first_frame) {
    for (uint32_t k = 0U; k < self->noise_spectrum_size; k++) {
      self->current->smoothed_spectrum[k] = spectrum[k];
      self->current->local_minimum_spectrum[k] = spectrum[k];
      noise_spectrum[k] = spectrum[k];
    }
    self->is_first_frame = false;
  } else {
    for (uint32_t k = 0U; k < self->noise_spectrum_size; k++) {
      self->current->smoothed_spectrum[k] =
          (N_SMOOTH * self->previous->smoothed_spectrum[k]) +
          ((1.F - N_SMOOTH) * spectrum[k]);

      if (self->previous->local_minimum_spectrum[k] <
          self->current->smoothed_spectrum[k]) {
        self->current->local_minimum_spectrum[k] =
            (GAMMA * self->previous->local_minimum_spectrum[k]) +
            (((1.F - GAMMA) / (1.F - BETA_AT)) *
             (self->current->smoothed_spectrum[k] -
              (BETA_AT * self->previous->smoothed_spectrum[k])));
      } else {
        self->current->local_minimum_spectrum[k] =
            self->current->smoothed_spectrum[k];
      }

      self->noisy_speech_ratio = sanitize_denormal(
          self->current->smoothed_spectrum[k] /
          (self->current->local_minimum_spectrum[k] + SPECTRAL_EPSILON));

      if (self->noisy_speech_ratio > self->minimum_detection_thresholds[k]) {
        self->speech_presence_detection[k] = 1U;
      } else {
        self->speech_presence_detection[k] = 0U;
      }

      self->current->speech_present_probability_spectrum[k] =
          (ALPHA_P * self->previous->speech_present_probability_spectrum[k]) +
          ((1.F - ALPHA_P) * (float)self->speech_presence_detection[k]);

      self->time_frequency_smoothing_constant[k] =
          ALPHA_D + ((1.F - ALPHA_D) *
                     self->current->speech_present_probability_spectrum[k]);

      noise_spectrum[k] =
          (self->time_frequency_smoothing_constant[k] *
           self->previous_noise_spectrum[k]) +
          ((1.F - self->time_frequency_smoothing_constant[k]) * spectrum[k]);
    }
  }

  update_frame_spectums(self, noise_spectrum);

  return true;
}

void louizou_noise_estimator_set_state(LouizouNoiseEstimator* self,
                                       const float* initial_profile) {
  if (!self || !initial_profile) {
    return;
  }

  for (uint32_t k = 0U; k < self->noise_spectrum_size; k++) {
    float val = fmaxf(initial_profile[k], FLT_MIN);

    self->previous_noise_spectrum[k] = val;
    self->current->smoothed_spectrum[k] = val;
    self->current->local_minimum_spectrum[k] = val;
    self->previous->smoothed_spectrum[k] = val;
    self->previous->local_minimum_spectrum[k] = val;
  }

  self->is_first_frame = false;
}

void louizou_noise_estimator_update_seed(LouizouNoiseEstimator* self,
                                         const float* seed_profile) {
  if (!self || !seed_profile) {
    return;
  }

  for (uint32_t k = 0U; k < self->noise_spectrum_size; k++) {
    float val = fmaxf(seed_profile[k], FLT_MIN);
    self->previous_noise_spectrum[k] = val;
    self->current->smoothed_spectrum[k] = val;
    self->current->local_minimum_spectrum[k] = val;
    self->previous->smoothed_spectrum[k] = val;
    self->previous->local_minimum_spectrum[k] = val;
  }

  self->is_first_frame = false;
}

void louizou_noise_estimator_apply_floor(LouizouNoiseEstimator* self,
                                         const float* floor_profile) {
  if (!self || !floor_profile) {
    return;
  }

  for (uint32_t k = 0U; k < self->noise_spectrum_size; k++) {
    float floor_val = floor_profile[k];
    if (self->previous_noise_spectrum[k] < floor_val) {
      self->previous_noise_spectrum[k] = floor_val;
    }
    if (self->current->local_minimum_spectrum[k] < floor_val) {
      self->current->local_minimum_spectrum[k] = floor_val;
    }
    if (self->previous->local_minimum_spectrum[k] < floor_val) {
      self->previous->local_minimum_spectrum[k] = floor_val;
    }
  }
}
