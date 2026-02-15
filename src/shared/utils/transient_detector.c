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
#include "spectral_utils.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

// Note: Transient detection constants are now imported from configurations.h

struct TransientDetector {
  uint32_t num_items;
  float* smoothed_items;
  float* previous_items;
  float alpha;
  float rolling_mean;
  uint32_t window_count;
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
  self->previous_items = (float*)calloc(num_items, sizeof(float));

  if (!self->smoothed_items || !self->previous_items) {
    transient_detector_free(self);
    return NULL;
  }

  self->alpha = TRANSIENT_SMOOTH_ALPHA;
  self->rolling_mean = 0.F;
  self->window_count = 0U;
  self->initialized = false;

  return self;
}

void transient_detector_free(TransientDetector* self) {
  if (!self) {
    return;
  }
  free(self->smoothed_items);
  free(self->previous_items);

  free(self);
}

void transient_detector_process(TransientDetector* self,
                                const float* band_energies,
                                float* onset_weights) {
  if (!self || !band_energies || !onset_weights) {
    return;
  }

  if (!self->initialized) {
    memcpy(self->smoothed_items, band_energies,
           sizeof(float) * self->num_items);
    self->initialized = true;
  }

  for (uint32_t j = 0; j < self->num_items; j++) {
    const float current = band_energies[j];
    const float smoothed = self->smoothed_items[j];

    // Calculate Ratio (Innovation)
    const float ratio = current / (smoothed + 1e-9F);

    // Onset weight: 0.0 at ratio 1.0, 1.0 at ratio 1.25 (25% increase)
    // We only trigger weight if the signal energy is significant
    float weight = (ratio - 1.0F) / ONSET_RATIO_SENSITIVITY;
    weight = (current < MIN_INNOVATION_ENERGY) ? 0.0F : weight;
    onset_weights[j] = fminf(fmaxf(weight, 0.0F), 1.0F);

    // Update smoothing reference
    float adapt_alpha = (weight > 0.0F) ? (self->alpha * 0.5F) : self->alpha;
    self->smoothed_items[j] =
        (current * (1.0F - adapt_alpha)) + (smoothed * adapt_alpha);
  }
}

bool transient_detector_run(TransientDetector* self, const float* spectrum) {
  if (!self || !spectrum) {
    return false;
  }

  const float flux =
      spectral_flux(spectrum, self->previous_items, self->num_items);

  self->window_count += 1U;

  if (self->window_count > 1U) {
    self->rolling_mean +=
        ((flux - self->rolling_mean) / (float)self->window_count);
  } else {
    self->rolling_mean = flux;
  }

  // Threshold logic matching original implementation
  const float adapted_threshold =
      ((UPPER_LIMIT - DEFAULT_TRANSIENT_THRESHOLD) * self->rolling_mean) +
      1e-6F;

  memcpy(self->previous_items, spectrum, sizeof(float) * self->num_items);

  return (flux > adapted_threshold);
}
