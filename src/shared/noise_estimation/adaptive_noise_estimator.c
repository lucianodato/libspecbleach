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

#include "adaptive_noise_estimator.h"
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

static FrameSpectrum* frame_spectrum_initialize(uint32_t frame_size);
static void frame_spectrum_free(FrameSpectrum* self);
static void compute_auto_thresholds(AdaptiveNoiseEstimator* self,
                                    uint32_t sample_rate,
                                    uint32_t noise_spectrum_size,
                                    uint32_t fft_size);
static void update_frame_spectums(AdaptiveNoiseEstimator* self,
                                  const float* noise_spectrum);

// SPP-MMSE helper functions
static float compute_spp_probability(float observation_power,
                                     float previous_noise_psd);
static float compute_mmse_noise_estimate(float spp_h1, float spp_h0,
                                         float observation_power,
                                         float previous_noise_psd);

struct AdaptiveNoiseEstimator {
  uint32_t noise_spectrum_size;
  float noisy_speech_ratio;

  FrameSpectrum* current;
  FrameSpectrum* previous;

  float* minimum_detection_thresholds;
  float* previous_noise_spectrum;
  float* time_frequency_smoothing_constant;
  uint32_t* speech_presence_detection;
  bool is_first_frame;

  // SPP-MMSE specific fields (optional, used when SPP method is selected)
  float* spp_previous_noise_psd; // Previous noise PSD estimate
  float* spp_smoothed_spp;       // Smoothed SPP for stagnation control
};

// SPP-MMSE helper function implementations
static float compute_spp_probability(float observation_power,
                                     float previous_noise_psd) {
  // Avoid division by zero and ensure numerical stability
  if (previous_noise_psd < 1e-12F) {
    previous_noise_psd = 1e-12F;
  }

  // Compute the exponent: -(|y|^2 / σ_N²(l-1)) * (ξ_H1 / (1 + ξ_H1))
  float ratio = observation_power / previous_noise_psd;
  float exponent = -ratio * (SPP_FIXED_XI_H1 / (1.F + SPP_FIXED_XI_H1));

  // Compute exp(exponent) with numerical stability check
  float exp_term = expf(exponent);
  if (!isfinite(exp_term)) {
    exp_term = (exponent > 0.F) ? FLT_MAX : 0.F;
  }

  // Compute the ratio: P(H0)/P(H1) * (1 + ξ_H1) * exp(...)
  // Since P(H0) = P(H1) = 0.5, P(H0)/P(H1) = 1
  float denominator_ratio = (1.F + SPP_FIXED_XI_H1) * exp_term;

  // Compute SPP: 1 / (1 + denominator_ratio)
  float spp = 1.F / (1.F + denominator_ratio);

  // Ensure SPP is in valid range [0, 1]
  spp = fmaxf(0.F, fminf(1.F, spp));

  return spp;
}

static float compute_mmse_noise_estimate(float spp_h1, float spp_h0,
                                         float observation_power,
                                         float previous_noise_psd) {
  // MMSE estimate: E{|N|²|y} = P(H0|y) * |y|² + P(H1|y) * σ_N²(l-1)
  return spp_h0 * observation_power + spp_h1 * previous_noise_psd;
}

AdaptiveNoiseEstimator* louizou_estimator_initialize(
    const uint32_t noise_spectrum_size, const uint32_t sample_rate,
    const uint32_t fft_size) {
  AdaptiveNoiseEstimator* self =
      (AdaptiveNoiseEstimator*)calloc(1U, sizeof(AdaptiveNoiseEstimator));
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
  self->spp_previous_noise_psd =
      (float*)calloc(self->noise_spectrum_size, sizeof(float));
  self->spp_smoothed_spp =
      (float*)calloc(self->noise_spectrum_size, sizeof(float));

  if (!self->minimum_detection_thresholds ||
      !self->time_frequency_smoothing_constant ||
      !self->speech_presence_detection || !self->previous_noise_spectrum ||
      !self->spp_previous_noise_psd || !self->spp_smoothed_spp) {
    louizou_estimator_free(self);
    return NULL;
  }

  compute_auto_thresholds(self, sample_rate, noise_spectrum_size, fft_size);
  self->current = frame_spectrum_initialize(noise_spectrum_size);
  self->previous = frame_spectrum_initialize(noise_spectrum_size);

  if (!self->current || !self->previous) {
    louizou_estimator_free(self);
    return NULL;
  }

  self->noisy_speech_ratio = 0.F;
  self->is_first_frame = true;

  return self;
}

AdaptiveNoiseEstimator* spp_mmse_estimator_initialize(
    const uint32_t noise_spectrum_size, const uint32_t sample_rate,
    const uint32_t fft_size) {
  AdaptiveNoiseEstimator* self =
      (AdaptiveNoiseEstimator*)calloc(1U, sizeof(AdaptiveNoiseEstimator));
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
  self->spp_previous_noise_psd =
      (float*)calloc(self->noise_spectrum_size, sizeof(float));
  self->spp_smoothed_spp =
      (float*)calloc(self->noise_spectrum_size, sizeof(float));

  if (!self->minimum_detection_thresholds ||
      !self->time_frequency_smoothing_constant ||
      !self->speech_presence_detection || !self->previous_noise_spectrum ||
      !self->spp_previous_noise_psd || !self->spp_smoothed_spp) {
    spp_mmse_estimator_free(self);
    return NULL;
  }

  compute_auto_thresholds(self, sample_rate, noise_spectrum_size, fft_size);
  self->current = frame_spectrum_initialize(noise_spectrum_size);
  self->previous = frame_spectrum_initialize(noise_spectrum_size);

  if (!self->current || !self->previous) {
    spp_mmse_estimator_free(self);
    return NULL;
  }

  self->noisy_speech_ratio = 0.F;
  self->is_first_frame = true;

  return self;
}

void spp_mmse_estimator_free(AdaptiveNoiseEstimator* self) {
  louizou_estimator_free(self); // Reuse the same cleanup logic
}

void louizou_estimator_free(AdaptiveNoiseEstimator* self) {
  free(self->minimum_detection_thresholds);
  free(self->time_frequency_smoothing_constant);
  free(self->speech_presence_detection);
  free(self->previous_noise_spectrum);
  free(self->spp_previous_noise_psd);
  free(self->spp_smoothed_spp);

  frame_spectrum_free(self->current);
  frame_spectrum_free(self->previous);

  free(self);
}

bool louizou_estimator_run(AdaptiveNoiseEstimator* self, const float* spectrum,
                           float* noise_spectrum) {
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
          N_SMOOTH * self->previous->smoothed_spectrum[k] +
          (1.F - N_SMOOTH) * spectrum[k];

      if (self->previous->local_minimum_spectrum[k] <
          self->current->smoothed_spectrum[k]) {
        self->current->local_minimum_spectrum[k] =
            GAMMA * self->previous->local_minimum_spectrum[k] +
            ((1.F - GAMMA) / (1.F - BETA_AT)) *
                (self->current->smoothed_spectrum[k] -
                 BETA_AT * self->previous->smoothed_spectrum[k]);
      } else {
        self->current->local_minimum_spectrum[k] =
            self->current->smoothed_spectrum[k];
      }

      self->noisy_speech_ratio = sanitize_denormal(
          self->current->smoothed_spectrum[k] /
          (self->current->local_minimum_spectrum[k] + 1e-12F));

      if (self->noisy_speech_ratio > self->minimum_detection_thresholds[k]) {
        self->speech_presence_detection[k] = 1U;
      } else {
        self->speech_presence_detection[k] = 0U;
      }

      self->current->speech_present_probability_spectrum[k] =
          ALPHA_P * self->previous->speech_present_probability_spectrum[k] +
          (1.F - ALPHA_P) * (float)self->speech_presence_detection[k];

      self->time_frequency_smoothing_constant[k] =
          ALPHA_D + (1.F - ALPHA_D) *
                        self->current->speech_present_probability_spectrum[k];

      noise_spectrum[k] =
          self->time_frequency_smoothing_constant[k] *
              self->previous_noise_spectrum[k] +
          (1.F - self->time_frequency_smoothing_constant[k]) * spectrum[k];
    }
  }

  update_frame_spectums(self, noise_spectrum);

  return true;
}

bool spp_mmse_estimator_run(AdaptiveNoiseEstimator* self, const float* spectrum,
                            float* noise_spectrum) {
  if (!self || !spectrum || !noise_spectrum) {
    return false;
  }

  if (self->is_first_frame) {
    // Initialize with first frame (assume noise-only)
    for (uint32_t k = 0U; k < self->noise_spectrum_size; k++) {
      self->spp_previous_noise_psd[k] = spectrum[k];
      self->spp_smoothed_spp[k] = 0.F; // Initialize smoothed SPP to 0
      noise_spectrum[k] = spectrum[k];
    }
    self->is_first_frame = false;
  } else {
    for (uint32_t k = 0U; k < self->noise_spectrum_size; k++) {
      // Step 1: Compute A Posteriori Speech Presence Probability
      float spp_h1 =
          compute_spp_probability(spectrum[k], self->spp_previous_noise_psd[k]);

      // Step 2: Apply stagnation control
      // If smoothed SPP > 0.99, cap current SPP at 0.99 to allow noise update
      if (self->spp_smoothed_spp[k] > SPP_STAGNATION_CAP) {
        spp_h1 = fminf(spp_h1, SPP_STAGNATION_CAP);
      }
      float spp_h0 = 1.F - spp_h1;

      // Step 3: Compute MMSE noise periodogram estimate
      float mmse_noise_estimate = compute_mmse_noise_estimate(
          spp_h1, spp_h0, spectrum[k], self->spp_previous_noise_psd[k]);

      // Step 4: Temporal smoothing
      // σ_N²(l) = α_pow * σ_N²(l-1) + (1 - α_pow) * E{|N|²|y}
      noise_spectrum[k] = SPP_ALPHA_POW * self->spp_previous_noise_psd[k] +
                          (1.F - SPP_ALPHA_POW) * mmse_noise_estimate;

      // Step 5: Update smoothed SPP for next frame's stagnation control
      // P̄(l) = 0.9 * P̄(l-1) + 0.1 * P(H1|y)
      self->spp_smoothed_spp[k] =
          SPP_SMOOTH_SPP * self->spp_smoothed_spp[k] + SPP_CURRENT_SPP * spp_h1;

      // Step 6: Store current noise estimate for next frame
      self->spp_previous_noise_psd[k] = noise_spectrum[k];
    }
  }

  // Update frame spectrums (reuse existing infrastructure for compatibility)
  update_frame_spectums(self, noise_spectrum);

  return true;
}

static void update_frame_spectums(AdaptiveNoiseEstimator* self,
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

static FrameSpectrum* frame_spectrum_initialize(const uint32_t frame_size) {
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
    frame_spectrum_free(self);
    return NULL;
  }

  (void)initialize_spectrum_with_value(self->local_minimum_spectrum, frame_size,
                                       FLT_MIN);

  return self;
}

static void frame_spectrum_free(FrameSpectrum* self) {
  free(self->smoothed_spectrum);
  free(self->local_minimum_spectrum);
  free(self->speech_present_probability_spectrum);

  free(self);
}

static void compute_auto_thresholds(AdaptiveNoiseEstimator* self,
                                    const uint32_t sample_rate,
                                    const uint32_t noise_spectrum_size,
                                    const uint32_t fft_size) {
  uint32_t LF = freq_to_fft_bin(CROSSOVER_POINT1, sample_rate, fft_size);
  uint32_t MF = freq_to_fft_bin(CROSSOVER_POINT2, sample_rate, fft_size);
  for (uint32_t k = 0U; k < noise_spectrum_size; k++) {
    if (k <= LF) {
      self->minimum_detection_thresholds[k] = BAND_1_LEVEL;
    }
    if (k > LF && k < MF) {
      self->minimum_detection_thresholds[k] = BAND_2_LEVEL;
    }
    if (k >= MF) {
      self->minimum_detection_thresholds[k] = BAND_3_LEVEL;
    }
  }
}
