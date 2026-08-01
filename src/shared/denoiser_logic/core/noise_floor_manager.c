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

#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/processing/spectral_whitening.h"
#include "shared/utils/spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseFloorManager {
  SpectralWhitening* whitening;
  float* whitening_weights;
  float* prev_gain_spectrum;
  float* gain_mean;
  float* gain_variance;
  uint32_t real_spectrum_size;
  uint32_t fft_size;
};

NoiseFloorManager* noise_floor_manager_initialize(const uint32_t fft_size) {
  if (fft_size == 0U) {
    return NULL;
  }
  NoiseFloorManager* self = 
      (NoiseFloorManager*)calloc(1U, sizeof(NoiseFloorManager));
  if (!self) {
    return NULL;
  }

  self->real_spectrum_size = (fft_size / 2U) + 1U;
  self->fft_size = fft_size;

  self->whitening = spectral_whitening_initialize(fft_size);
  if (!self->whitening) {
    free(self);
    return NULL;
  }

  self->whitening_weights = 
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->whitening_weights) {
    spectral_whitening_free(self->whitening);
    free(self);
    return NULL;
  }

  self->prev_gain_spectrum = 
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->prev_gain_spectrum) {
    free(self->whitening_weights);
    spectral_whitening_free(self->whitening);
    free(self);
    return NULL;
  }

  self->gain_mean = 
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->gain_mean) {
    free(self->prev_gain_spectrum);
    free(self->whitening_weights);
    spectral_whitening_free(self->whitening);
    free(self);
    return NULL;
  }

  self->gain_variance = 
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->gain_variance) {
    free(self->gain_mean);
    free(self->prev_gain_spectrum);
    free(self->whitening_weights);
    spectral_whitening_free(self->whitening);
    free(self);
    return NULL;
  }

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    self->prev_gain_spectrum[k] = 1.0f;
    self->gain_mean[k] = 1.0f;
    self->gain_variance[k] = 0.0f;
  }

  return self;
}

void noise_floor_manager_free(NoiseFloorManager* self) {
  if (!self) {
    return;
  }
  if (self->whitening) {
    spectral_whitening_free(self->whitening);
  }
  if (self->whitening_weights) {
    free(self->whitening_weights);
  }
  if (self->prev_gain_spectrum) {
    free(self->prev_gain_spectrum);
  }
  if (self->gain_mean) {
    free(self->gain_mean);
  }
  if (self->gain_variance) {
    free(self->gain_variance);
  }
  free(self);
}

void noise_floor_manager_apply(NoiseFloorManager* self,
                               uint32_t real_spectrum_size, uint32_t fft_size,
                               float* gain_spectrum, const float* noise_profile,
                               float reduction_amount,
                               float tonal_reduction_amount,
                               const float* tonal_mask,
                               float whitening_factor) {
  if (!self || !gain_spectrum || !noise_profile) {
    return;
  }

  if (real_spectrum_size != self->real_spectrum_size ||
      fft_size != self->fft_size) {
    real_spectrum_size = self->real_spectrum_size;
    fft_size = self->fft_size;
  }

  if (reduction_amount >= 0.999f && tonal_reduction_amount >= 0.999f) {
    // Transparency Guard: If both reduction paths are at 0dB (1.0f),
    // we force unity gain to ensure BIT TRANSPARENCY and skip processing.
    for (uint32_t k = 0U; k < real_spectrum_size; k++) {
      gain_spectrum[k] = 1.0f;
      self->prev_gain_spectrum[k] = 1.0f;
      self->gain_mean[k] = 1.0f;
      self->gain_variance[k] = 0.0f;
    }
    sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
    return;
  }

  // 1. Calculate whitening weights (including tapering)
  spectral_whitening_get_weights(self->whitening, whitening_factor,
                                 noise_profile, self->whitening_weights);

  // 2. Apply biasing + frequency-dependent floor
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float mask = (tonal_mask) ? tonal_mask[k] : 0.0f;
    if (mask > 0.0f) {
      mask = sqrtf(sqrtf(mask));
    }
    // Proportional interpolation between regular reduction and tonal reduction
    // Note: These are already linear gains (0.0 to 1.0) passed from the loading
    // layer.
    float dual_path_reduction = 
        (reduction_amount * (1.0f - mask)) + (tonal_reduction_amount * mask);

    /*
     * TONAL HAND-OFF:
     * As whitening increases, we want to smooth out the "notches" created by
     * tonal reduction to achieve a truly flat white residual.
     * We interpolate between the dual-path reduction (with notches) and the
     * broadband-only reduction level based on the whitening strength.
     */
    float target_reduction = 
        dual_path_reduction +
        (whitening_factor * (reduction_amount - dual_path_reduction));

    /*
     * DYNAMIC WHITENING DEPTH & BIT TRANSPARENCY:
     * To achieve bit transparency at 0dB reduction, we damp the whitening
     * weights toward 1.0 as the target reduction level approaches 1.0 (no
     * reduction). This also makes the effect feel more natural and dynamic.
     */
    float reduction_depth = 1.0f - target_reduction;
    float effective_whitening_weight = 
        1.0f + (reduction_depth * (self->whitening_weights[k] - 1.0f));

    float whitened_floor = target_reduction * effective_whitening_weight;

    // Safety clamp (Floor cannot exceed Unity Gain)
    if (whitened_floor > 1.0f) {
      whitened_floor = 1.0f;
    }

    /*
     * CLAMPING FLOOR:
     * We use fmaxf to ensure the gain never drops below the whitened floor.
     * This "replaces" the suppressed noise with the floor level,
     * effectively shaving off the spectral shape of the noise residual.
     */
    gain_spectrum[k] = fmaxf(whitened_floor, gain_spectrum[k]);
  }

  // Update temporal mean and variance of the raw gain spectrum
  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    float current_val = gain_spectrum[k];
    float diff = current_val - self->gain_mean[k];
    self->gain_mean[k] = 0.75f * self->gain_mean[k] + 0.25f * current_val;
    self->gain_variance[k] = 0.75f * self->gain_variance[k] + 0.25f * (diff * diff);
  }

  // Apply Neighborhood-Coherence Filter to gain_spectrum
  float* temp_coherent = (float*)malloc(real_spectrum_size * sizeof(float));
  if (temp_coherent) {
    for (uint32_t k = 0U; k < real_spectrum_size; k++) {
      float fluctuation_strength = sqrtf(self->gain_variance[k]);
      // Map fluctuation strength to coherence blend factor [0.0, 0.80]
      float g_coherence = fluctuation_strength * 2.5f;
      if (g_coherence > 0.80f) {
        g_coherence = 0.80f;
      }
      if (g_coherence < 0.0f) {
        g_coherence = 0.0f;
      }

      // Compute neighborhood mean using a 5-bin window
      float neighborhood_sum = 0.0f;
      float weight_sum = 0.0f;
      int32_t start_offset = -2;
      int32_t end_offset = 2;
      for (int32_t offset = start_offset; offset <= end_offset; offset++) {
        int32_t idx = (int32_t)k + offset;
        if (idx >= 0 && idx < (int32_t)real_spectrum_size) {
          float w = 1.0f;
          if (offset == -2 || offset == 2) {
            w = 0.5f;
          }
          neighborhood_sum += gain_spectrum[idx] * w;
          weight_sum += w;
        }
      }
      float neighborhood_mean = neighborhood_sum / weight_sum;

      // Blend raw gain with neighborhood mean based on local fluctuation
      temp_coherent[k] = (1.0f - g_coherence) * gain_spectrum[k] + g_coherence * neighborhood_mean;
    }
    memcpy(gain_spectrum, temp_coherent, real_spectrum_size * sizeof(float));
    free(temp_coherent);
  }

  // 3. Stateful Adaptive Spectral-Temporal Gain Smoothing
  float* temp_smoothed = (float*)malloc(real_spectrum_size * sizeof(float));
  if (temp_smoothed) {
    float total_gain_increase = 0.0f;
    for (uint32_t k = 0U; k < real_spectrum_size; k++) {
      float diff = gain_spectrum[k] - self->prev_gain_spectrum[k];
      if (diff > 0.0f) {
        total_gain_increase += diff;
      }
    }
    float mean_gain_increase = total_gain_increase / (float)real_spectrum_size;

    // If there is a sudden global gain increase, we are likely at a transient onset.
    // We scale up the attack rate to be almost instantaneous.
    float transient_boost = mean_gain_increase * 5.0f;
    if (transient_boost > 0.15f) {
      transient_boost = 0.15f;
    }

    float alpha_attack_base = 0.85f + transient_boost;
    float alpha_release_base = 0.15f;

    for (uint32_t k = 0U; k < real_spectrum_size; k++) {
      float mask = (tonal_mask) ? tonal_mask[k] : 0.0f;
      if (mask > 1.0f) {
        mask = 1.0f;
      }
      if (mask < 0.0f) { 
        mask = 0.0f;
      }

      float current_gain = gain_spectrum[k];
      float prev_gain = self->prev_gain_spectrum[k];
      float alpha_k;

      if (current_gain > prev_gain) {
        // Attack phase
        alpha_k = alpha_attack_base + (1.0f - alpha_attack_base) * mask;
      } else {
        // Release phase
        alpha_k = alpha_release_base + (0.50f - alpha_release_base) * mask;
      }

      // Dynamic stabilization based on temporal variance
      // If the bin is highly volatile, we slow down the smoothing (reduce alpha_k)
      // but only if we are not in a global transient.
      if (mean_gain_increase < 0.05f) {
        float fluctuation_strength = sqrtf(self->gain_variance[k]);
        float stabilization = 1.0f - fluctuation_strength * 2.0f;
        if (stabilization < 0.3f) {
          stabilization = 0.3f;
        }
        alpha_k *= stabilization;
      }

      if (alpha_k > 1.0f) {
        alpha_k = 1.0f;
      }
      if (alpha_k < 0.0f) {
        alpha_k = 0.0f;
      }

      temp_smoothed[k] = alpha_k * current_gain + (1.0f - alpha_k) * prev_gain;
      self->prev_gain_spectrum[k] = temp_smoothed[k];
    }

    // Tonal-guided spectral smoothing
    for (uint32_t k = 0U; k < real_spectrum_size; k++) {
      float mask = (tonal_mask) ? tonal_mask[k] : 0.0f;
      if (mask > 1.0f) {
        mask = 1.0f;
      }
      if (mask < 0.0f) {
        mask = 0.0f;
      }

      float spec_smoothed_gain;
      if (k == 0U) {
        spec_smoothed_gain = temp_smoothed[0U];
      } else if (k == real_spectrum_size - 1U) {
        spec_smoothed_gain = temp_smoothed[real_spectrum_size - 1U];
      } else {
        spec_smoothed_gain = 0.25f * temp_smoothed[k - 1U] +
                             0.50f * temp_smoothed[k] +
                             0.25f * temp_smoothed[k + 1U];
      }

      gain_spectrum[k] = mask * temp_smoothed[k] + (1.0f - mask) * spec_smoothed_gain;
    }

    free(temp_smoothed);
  }

  sb_apply_spectral_symmetry(gain_spectrum, real_spectrum_size, fft_size);
}