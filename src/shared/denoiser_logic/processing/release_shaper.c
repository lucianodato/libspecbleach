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

#include "release_shaper.h"
#include "shared/configurations.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/simd_utils.h"
#include "shared/utils/spectral_utils.h"
#include <math.h>
#include <stdlib.h>

struct ReleaseShaper {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  float hop_sec; // 0 = derive from fft/overlap
  bool seeded;

  CriticalBands* critical_bands;
  uint32_t number_bands;
  float* band_envelope; // Recent-energy memory per band (fast attack / slow
                        // release)
  float* band_current;  // Current frame band energy per band
};

ReleaseShaper* release_shaper_initialize(const uint32_t sample_rate,
                                         const uint32_t fft_size) {
  if (fft_size == 0U || sample_rate == 0U) {
    return NULL;
  }

  ReleaseShaper* self = (ReleaseShaper*)calloc(1U, sizeof(ReleaseShaper));
  if (!self) {
    return NULL;
  }

  self->fft_size = fft_size;
  self->real_spectrum_size = (fft_size / 2U) + 1U;
  self->sample_rate = sample_rate;
  self->hop_sec = 0.0F;
  self->seeded = false;

  self->critical_bands =
      critical_bands_initialize(sample_rate, fft_size, BARK_SCALE);
  if (!self->critical_bands) {
    release_shaper_free(self);
    return NULL;
  }
  self->number_bands = get_number_of_critical_bands(self->critical_bands);

  self->band_envelope = (float*)calloc(self->number_bands, sizeof(float));
  self->band_current = (float*)calloc(self->number_bands, sizeof(float));

  if (!self->band_envelope || !self->band_current) {
    release_shaper_free(self);
    return NULL;
  }

  return self;
}

void release_shaper_free(ReleaseShaper* self) {
  if (!self) {
    return;
  }
  if (self->critical_bands) {
    critical_bands_free(self->critical_bands);
  }
  free(self->band_envelope);
  free(self->band_current);
  free(self);
}

void release_shaper_reset(ReleaseShaper* self) {
  if (!self) {
    return;
  }
  self->seeded = false;
}

void release_shaper_set_hop_sec(ReleaseShaper* self, const float hop_sec) {
  if (!self) {
    return;
  }
  self->hop_sec = hop_sec;
}

static float release_shaper_get_hop_sec(const ReleaseShaper* self) {
  if (self->hop_sec > 0.0F) {
    return self->hop_sec;
  }
  const float hop_samples = (float)self->fft_size / (float)OVERLAP_FACTOR;
  return hop_samples / (float)self->sample_rate;
}

bool release_shaper_compute(ReleaseShaper* self, const float* magnitude,
                            float* release_scale) {
  if (!self || !magnitude || !release_scale) {
    return false;
  }

  // Per-frame envelope decays are a denormal hazard after a signal ends;
  // guard the whole pass like the other per-frame DSP paths.
  sb_simd_state_t old_simd_state = sb_simd_enable_ftz_daz();

  // Default: full slider release everywhere (bins outside any band, first
  // frame, and bands without collapse evidence)
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    release_scale[k] = 1.0F;
  }

  // Aggregate current band energies and scatter per-bin scales
  for (uint32_t j = 0U; j < self->number_bands; j++) {
    const CriticalBandIndexes band_indexes =
        get_band_indexes(self->critical_bands, j);

    float band_energy = 0.0F;
    for (uint32_t k = band_indexes.start_position;
         k < band_indexes.end_position; k++) {
      band_energy += magnitude[k];
    }
    self->band_current[j] = band_energy;

    if (!self->seeded) {
      self->band_envelope[j] = band_energy;
      continue;
    }

    // Collapse evidence: recent-energy envelope vs current energy. A band
    // that hard-cuts to (near) zero after a loud frame is the strongest
    // collapse case, so a nontrivial envelope with floored current energy
    // still yields the fast close. Only a trivial envelope (no past energy)
    // carries no evidence and keeps full release.
    float weight = 1.0F;
    if (self->band_envelope[j] > SPECTRAL_EPSILON) {
      const float current_energy = fmaxf(band_energy, SPECTRAL_EPSILON);
      const float gap_db =
          10.0F * log10f(self->band_envelope[j] / current_energy);
      const float normalized_gap = (gap_db - RELEASE_SHAPING_GAP_ONSET_DB) /
                                   RELEASE_SHAPING_GAP_RANGE_DB;
      weight = 1.0F - fmaxf(0.0F, fminf(1.0F, normalized_gap));
    }
    for (uint32_t k = band_indexes.start_position;
         k < band_indexes.end_position; k++) {
      release_scale[k] = weight;
    }
  }

  if (!self->seeded) {
    self->seeded = true;
    sb_simd_restore_state(old_simd_state);
    return true;
  }

  // Envelope update: fast attack (RELEASE_SHAPING_ATTACK_HOPS hops,
  // frame-rate invariant by construction), slow release (named constant)
  const float hop_sec = release_shaper_get_hop_sec(self);
  const float alpha_attack = expf(-RELEASE_SHAPING_ATTACK_HOPS);
  const float alpha_release =
      expf(-hop_sec / RELEASE_SHAPING_ENVELOPE_RELEASE_SEC);

  for (uint32_t j = 0U; j < self->number_bands; j++) {
    const float current = self->band_current[j];
    float envelope = self->band_envelope[j];
    if (current >= envelope) {
      envelope = (alpha_attack * envelope) + ((1.0F - alpha_attack) * current);
    } else {
      envelope =
          (alpha_release * envelope) + ((1.0F - alpha_release) * current);
    }
    self->band_envelope[j] = envelope;
  }

  sb_simd_restore_state(old_simd_state);

  return true;
}
