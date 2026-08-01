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

#include "shared/denoiser_logic/processing/spectral_whitening.h"
#include "shared/configurations.h"
#include "shared/utils/spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int float_comparator(const void* a, const void* b) {
  float x = *(const float*)a;
  float y = *(const float*)b;
  return (x > y) - (x < y);
}

struct SpectralWhitening {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
};

SpectralWhitening* spectral_whitening_initialize(const uint32_t fft_size) {
  SpectralWhitening* self =
      (SpectralWhitening*)calloc(1U, sizeof(SpectralWhitening));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = (self->fft_size / 2U) + 1U;

  return self;
}

void spectral_whitening_free(SpectralWhitening* self) {
  if (!self) {
    return;
  }
  free(self);
}

#define WHITENING_NUM_BANDS 24U

void spectral_whitening_get_weights(SpectralWhitening* self,
                                    float whitening_factor,
                                    const float* noise_profile,
                                    float* weights_out) {
  if (!self || !weights_out || !noise_profile) {
    return;
  }

  // 1. Create a smoothed copy of the noise profile for stable weight calculation
  float* smoothed_profile =
      (float*)malloc(self->real_spectrum_size * sizeof(float));
  if (!smoothed_profile) {
    return;
  }
  memcpy(smoothed_profile, noise_profile,
         self->real_spectrum_size * sizeof(float));
  smooth_spectrum(smoothed_profile, self->real_spectrum_size, 0.5f);

  // 2. Calculate Localized Band Anchor Magnitudes using Bark-spaced sub-bands
  // Multi-band localized anchor tracking prevents low-frequency noise rumble
  // from dominating high-frequency whitening weights or vice-versa.
  float* local_anchors =
      (float*)malloc(self->real_spectrum_size * sizeof(float));
  if (!local_anchors) {
    free(smoothed_profile);
    return;
  }

  uint32_t num_bands = WHITENING_NUM_BANDS;
  if (self->real_spectrum_size < WHITENING_NUM_BANDS) {
    num_bands = self->real_spectrum_size;
  }

  float band_centers[WHITENING_NUM_BANDS];
  float band_medians[WHITENING_NUM_BANDS];

  float* band_sort_buf =
      (float*)malloc(self->real_spectrum_size * sizeof(float));
  if (!band_sort_buf) {
    free(local_anchors);
    free(smoothed_profile);
    return;
  }

  // Define band boundaries logarithmically/quadratically across frequency bins
  for (uint32_t b = 0U; b < num_bands; b++) {
    float frac_start = (float)b / (float)num_bands;
    float frac_end = (float)(b + 1U) / (float)num_bands;

    // Use squared warping for quasi-Bark band spacing
    uint32_t start_bin =
        (uint32_t)(frac_start * frac_start * (float)self->real_spectrum_size);
    uint32_t end_bin =
        (uint32_t)(frac_end * frac_end * (float)self->real_spectrum_size);

    if (end_bin <= start_bin) {
      end_bin = start_bin + 1U;
    }
    if (end_bin > self->real_spectrum_size) {
      end_bin = self->real_spectrum_size;
    }

    band_centers[b] = 0.5f * ((float)start_bin + (float)end_bin);

    uint32_t count = end_bin - start_bin;
    if (count > 0U) {
      memcpy(band_sort_buf, &smoothed_profile[start_bin], count * sizeof(float));
      qsort(band_sort_buf, count, sizeof(float), float_comparator);
      float med = band_sort_buf[count / 2U];
      band_medians[b] = (med < SPECTRAL_EPSILON) ? SPECTRAL_EPSILON : med;
    } else {
      band_medians[b] = SPECTRAL_EPSILON;
    }
  }

  free(band_sort_buf);

  // Interpolate localized anchor envelope across all frequency bins
  uint32_t current_band = 0U;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    while (current_band + 1U < num_bands &&
           (float)k > band_centers[current_band + 1U]) {
      current_band++;
    }

    if (current_band + 1U < num_bands) {
      float x0 = band_centers[current_band];
      float x1 = band_centers[current_band + 1U];
      float y0 = band_medians[current_band];
      float y1 = band_medians[current_band + 1U];

      if (x1 > x0) {
        float t = ((float)k - x0) / (x1 - x0);
        if (t < 0.0f) {
          t = 0.0f;
        }
        if (t > 1.0f) {
          t = 1.0f;
        }
        local_anchors[k] = y0 + t * (y1 - y0);
      } else {
        local_anchors[k] = y0;
      }
    } else {
      local_anchors[k] = band_medians[num_bands - 1U];
    }

    if (local_anchors[k] < SPECTRAL_EPSILON) {
      local_anchors[k] = SPECTRAL_EPSILON;
    }
  }

  // 3. Compute power-law whitening weights relative to local band anchors
  float normalized_factor = whitening_factor;

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    float weight = 1.0f;
    if (normalized_factor > 0.0f && smoothed_profile[k] > SPECTRAL_EPSILON) {
      weight = powf(local_anchors[k] / smoothed_profile[k], normalized_factor);
    }
    weights_out[k] = weight;
  }

  free(local_anchors);
  free(smoothed_profile);
}
