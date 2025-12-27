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

#include "postfilter.h"
#include "../configurations.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct PostFilter {
  float* intermediate_gains;

  uint32_t fft_size;
  uint32_t real_spectrum_size;
  bool preserve_minimum;
  float default_postfilter_scale;
  float min_gain_coefficient;
};

PostFilter* postfilter_initialize(const uint32_t fft_size) {
  PostFilter* self = (PostFilter*)calloc(1U, sizeof(PostFilter));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->preserve_minimum = (bool)PRESERVE_MINIMUN_GAIN;
  self->default_postfilter_scale = POSTFILTER_SCALE;
  self->min_gain_coefficient = powf(10.F, (float)POSTFILTER_MIN_GAIN_DB / 20.F);

  self->intermediate_gains =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  if (!self->intermediate_gains) {
    free(self);
    return NULL;
  }

  return self;
}

void postfilter_free(PostFilter* self) {
  if (!self) {
    return;
  }
  free(self->intermediate_gains);
  free(self);
}

static uint32_t get_adaptive_window_size(const PostFilter* self,
                                         const float* spectrum,
                                         const float snr_threshold,
                                         const float* gain_spectrum) {
  float clean_energy = 0.F;
  float noisy_energy = 0.F;

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float noisy = spectrum[k];
    const float clean = noisy * gain_spectrum[k];
    clean_energy += clean * clean;
    noisy_energy += noisy * noisy;
  }

  if (noisy_energy <= 1e-12F) {
    return 1U;
  }

  const float zeta = clean_energy / noisy_energy;
  const float zeta_t = (zeta >= snr_threshold) ? 1.F : zeta;

  if (zeta_t >= 1.F) {
    return 1U;
  }

  const float n = 2.F * roundf(self->default_postfilter_scale *
                               (1.F - zeta_t / snr_threshold)) +
                  1.F;

  return (uint32_t)n;
}

static void moving_average(const float* in, float* out, uint32_t size,
                           uint32_t n) {
  if (n <= 1U || n > size) {
    memcpy(out, in, size * sizeof(float));
    return;
  }

  const uint32_t half = n / 2U;
  double current_sum = 0.0;

  // Initial window sum (boundary handling: use clamping for start)
  for (int i = -(int)half; i <= (int)half; i++) {
    int idx = i;
    if (idx < 0)
      idx = 0;
    if (idx >= (int)size)
      idx = (int)size - 1;
    current_sum += (double)in[idx];
  }

  for (uint32_t i = 0U; i < size; i++) {
    out[i] = (float)(current_sum / (double)n);

    if (i + 1U < size) {
      // Move window: subtract oldest, add newest
      int old_idx = (int)i - (int)half;
      int new_idx = (int)i + (int)half + 1;

      if (old_idx < 0)
        old_idx = 0;
      if (new_idx >= (int)size)
        new_idx = (int)size - 1;

      current_sum -= (double)in[old_idx];
      current_sum += (double)in[new_idx];
    }
  }
}

bool postfilter_apply(PostFilter* self, const float* spectrum,
                      float* gain_spectrum,
                      const PostFiltersParameters parameters) {
  if (!self || !spectrum || !gain_spectrum) {
    return false;
  }

  const uint32_t n = get_adaptive_window_size(
      self, spectrum, parameters.snr_threshold, gain_spectrum);

  if (n > 1U) {
    moving_average(gain_spectrum, self->intermediate_gains,
                   self->real_spectrum_size, n);

    if (self->preserve_minimum) {
      for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
        gain_spectrum[k] = fminf(gain_spectrum[k], self->intermediate_gains[k]);
      }
    } else {
      memcpy(gain_spectrum, self->intermediate_gains,
             self->real_spectrum_size * sizeof(float));
    }
  }

  // Apply gain floor
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    if (gain_spectrum[k] < parameters.gain_floor) {
      gain_spectrum[k] = parameters.gain_floor;
    }
  }

  return true;
}
