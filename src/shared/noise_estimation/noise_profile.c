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

#include "noise_profile.h"
#include "../configurations.h"
#include "../utils/spectral_utils.h"
#include <stdlib.h>
#include <string.h>

struct NoiseProfile {
  uint32_t noise_profile_size;
  uint32_t noise_profile_blocks_averaged[NOISE_PROFILE_MODES];
  float* noise_profiles[NOISE_PROFILE_MODES];
  bool noise_spectrum_available[NOISE_PROFILE_MODES];
};

NoiseProfile* noise_profile_initialize(const uint32_t size) {
  NoiseProfile* self = (NoiseProfile*)calloc(1U, sizeof(NoiseProfile));
  if (!self) {
    return NULL;
  }
  self->noise_profile_size = size;

  for (int i = 0; i < NOISE_PROFILE_MODES; i++) {
    self->noise_profile_blocks_averaged[i] = 0U;
    self->noise_spectrum_available[i] = false;
    self->noise_profiles[i] = (float*)calloc(size, sizeof(float));
    if (!self->noise_profiles[i]) {
      noise_profile_free(self);
      return NULL;
    }
  }

  return self;
}

void noise_profile_free(NoiseProfile* self) {
  if (self) {
    for (int i = 0; i < NOISE_PROFILE_MODES; i++) {
      if (self->noise_profiles[i]) {
        free(self->noise_profiles[i]);
      }
    }
    free(self);
  }
}

bool is_noise_estimation_available(NoiseProfile* self, int mode) {
  if (mode < 1 || mode > 3)
    return false;
  return self->noise_spectrum_available[mode - 1];
}

float* get_noise_profile(NoiseProfile* self, int mode) {
  if (mode < 1 || mode > 3)
    return NULL;
  return self->noise_profiles[mode - 1];
}

uint32_t get_noise_profile_size(NoiseProfile* self) {
  return self->noise_profile_size;
}

uint32_t get_noise_profile_blocks_averaged(NoiseProfile* self, int mode) {
  if (mode < 1 || mode > 3)
    return 0;
  return self->noise_profile_blocks_averaged[mode - 1];
}
void set_noise_profile_available(NoiseProfile* self, int mode) {
  if (mode >= 1 && mode <= 3) {
    self->noise_spectrum_available[mode - 1] = true;
  }
}

bool set_noise_profile(NoiseProfile* self, int mode, const float* noise_profile,
                       const uint32_t noise_profile_size,
                       const uint32_t noise_profile_blocks_averaged) {
  if (!self || mode < 1 || mode > 3 || !noise_profile ||
      noise_profile_size != self->noise_profile_size) {
    return false;
  }
  int index = mode - 1;
  memcpy(self->noise_profiles[index], noise_profile,
         noise_profile_size * sizeof(float));

  self->noise_profile_blocks_averaged[index] = noise_profile_blocks_averaged;
  self->noise_spectrum_available[index] = true;

  return true;
}

bool increment_blocks_averaged(NoiseProfile* self, int mode) {
  if (!self || mode < 1 || mode > 3) {
    return false;
  }

  int index = mode - 1;
  self->noise_profile_blocks_averaged[index]++;

  if (self->noise_profile_blocks_averaged[index] >
          MIN_NUMBER_OF_WINDOWS_NOISE_AVERAGED &&
      !self->noise_spectrum_available[index]) {
    self->noise_spectrum_available[index] = true;
  }

  return true;
}

bool reset_noise_profile(NoiseProfile* self) {
  if (!self) {
    return false;
  }

  for (int i = 0; i < NOISE_PROFILE_MODES; i++) {
    (void)initialize_spectrum_with_value(self->noise_profiles[i],
                                         self->noise_profile_size, 0.F);
    self->noise_profile_blocks_averaged[i] = 0U;
    self->noise_spectrum_available[i] = false;
  }

  return true;
}
