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

#include "transient_detector.h"
#include "../configurations.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct TransientDetector {
  uint32_t num_items;
  float* smoothed_items;
  float alpha;
  bool initialized;
};

TransientDetector* transient_detector_initialize(const uint32_t num_items) {
  TransientDetector* self =
      (TransientDetector*)calloc(1U, sizeof(TransientDetector));

  if (!self) {
    return NULL;
  }

  self->num_items = num_items;
  self->smoothed_items = (float*)calloc(num_items, sizeof(float));

  if (!self->smoothed_items) {
    transient_detector_free(self);
    return NULL;
  }

  self->alpha = TRANSIENT_SMOOTH_ALPHA;
  self->initialized = false;

  return self;
}

void transient_detector_free(TransientDetector* self) {
  if (!self) {
    return;
  }

  free(self->smoothed_items);

  free(self);
}

bool transient_detector_process(TransientDetector* self,
                                const float* band_energies,
                                float* onset_weights, float* out_intensity) {
  if (!self || !band_energies) {
    if (out_intensity) {
      *out_intensity = 0.0f;
    }
    return false;
  }

  if (!self->initialized) {
    memcpy(self->smoothed_items, band_energies,
           sizeof(float) * self->num_items);
    self->initialized = true;
    if (onset_weights) {
      memset(onset_weights, 0, sizeof(float) * self->num_items);
    }
    if (out_intensity) {
      *out_intensity = 0.0f;
    }
    return false;
  }

  float max_weight = 0.0f;
  float sum_energy = 0.0f;
  float transient_energy = 0.0f;
  uint32_t active_bands_count = 0;
  uint32_t high_bands_count = 0;
  uint32_t mid_high_start_band =
      self->num_items / 3U; // upper 2/3 of critical spectrum

  for (uint32_t j = 0; j < self->num_items; j++) {
    const float current = band_energies[j];
    const float smoothed = self->smoothed_items[j];

    // Background tracking with asymmetric attack/decay:
    // Minimum energy floor per band (0.05 linear amplitude sum over critical
    // band bins)
    const float min_band_energy = 0.05F;
    const float base_floor = fmaxf(smoothed, min_band_energy);
    const float ratio = (current - base_floor) / base_floor;

    // Onset weight:
    // Requires energy to be significantly above floor and at least +100%
    // innovation jump
    float weight = 0.0F;
    if (current >= (min_band_energy * 2.0F) && ratio >= 1.0F) {
      weight = (ratio - 1.0F) / 2.5F;
      weight = fminf(fmaxf(weight, 0.0F), 1.0F);
    }

    if (weight >= 0.25F) {
      active_bands_count++;
      if (j >= mid_high_start_band) {
        high_bands_count++;
      }
    }

    if (onset_weights) {
      onset_weights[j] = weight;
    }

    if (weight > max_weight) {
      max_weight = weight;
    }
    sum_energy += current;
    transient_energy += weight * current;

    // Asymmetric background floor adaptation:
    // - On attack (current > smoothed): adapt floor slower (alpha = 0.85) to
    // preserve onset contrast
    // - On decay/sustain (current <= smoothed): track smoothly downwards (alpha
    // = 0.50) without lagging
    float adapt_alpha = (current > smoothed) ? self->alpha : 0.50F;
    self->smoothed_items[j] = (self->smoothed_items[j] * adapt_alpha) +
                              (current * (1.0F - adapt_alpha));
  }

  // Reject noise-only or silent frames:
  if (sum_energy < (0.05F * (float)self->num_items * 2.0F)) {
    if (onset_weights) {
      memset(onset_weights, 0, sizeof(float) * self->num_items);
    }
    if (out_intensity) {
      *out_intensity = 0.0f;
    }
    return false;
  }

  // Discriminate against speech vowels / voiced phonemes vs genuine
  // plosives/clicks/plucks: Voiced speech vowels have strong low-frequency
  // formant surges with little to no high-frequency onset innovation. Genuine
  // transients (mouth clicks, plucks, plosives, snaps) exhibit high-frequency
  // onset energy across mid/high bands.
  float hf_ratio = (self->num_items > mid_high_start_band)
                       ? ((float)high_bands_count /
                          (float)(self->num_items - mid_high_start_band))
                       : 0.0f;
  float band_fraction =
      (self->num_items > 0U)
          ? ((float)active_bands_count / (float)self->num_items)
          : 0.0f;
  float energy_fraction = (sum_energy > MIN_INNOVATION_ENERGY)
                              ? (transient_energy / sum_energy)
                              : 0.0f;

  // True transients: broadband onset OR high-frequency localized snap/click
  float broadband_score = 0.55f * band_fraction + 0.45f * energy_fraction;
  // Localized HF click score: NO additive offset — requires genuine hf_ratio >
  // 0 to score Fricative onsets hit only 1 upper band slowly; mouth clicks hit
  // 2+ upper bands simultaneously
  float localized_hf_score =
      (high_bands_count >= 2U) ? (max_weight * hf_ratio) : 0.0f;
  float intensity = fmaxf(broadband_score, localized_hf_score);

  // If high frequency onset is completely absent and bandwidth is narrow
  // (voiced vowel onset), suppress
  if (high_bands_count == 0 && band_fraction < 0.25f) {
    intensity = 0.0f;
  }

  // Strict zeroing of weak/ambiguous background fluctuations
  if (intensity < 0.25f && max_weight < 0.70f) {
    intensity = 0.0f;
  }

  intensity = fminf(fmaxf(intensity, 0.0f), 1.0f);

  if (out_intensity) {
    *out_intensity = intensity;
  }

  // Localized HF path requires 2+ upper bands with strong weight to avoid
  // fricative false triggers
  return (intensity >= 0.25f ||
          (max_weight >= 0.70f && high_bands_count >= 2U));
}
