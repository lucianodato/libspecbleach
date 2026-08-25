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

#include "specbleach_profile_migration.h"
#include <stddef.h>

static bool migrate_profiles_generic(
    void* source, void* target, const uint32_t profile_size,
    bool (*available)(void*, int), float* (*get_profile)(void*, int),
    uint32_t (*get_blocks)(void*, int),
    bool (*load_profile)(void*, const float*, uint32_t, uint32_t, int)) {
  if (!source || !target || profile_size == 0) {
    return false;
  }

  float* fallback_profile = NULL;
  uint32_t fallback_blocks = 0;
  bool any_available = false;

  for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
    if (!available(source, mode)) {
      continue;
    }
    float* profile = get_profile(source, mode);
    if (!profile) {
      continue;
    }
    any_available = true;
    fallback_profile = profile;
    fallback_blocks = get_blocks(source, mode);
    break;
  }

  if (!any_available) {
    return false;
  }

  for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
    float* profile = NULL;
    uint32_t blocks = 0;

    if (available(source, mode)) {
      profile = get_profile(source, mode);
      blocks = get_blocks(source, mode);
    } else {
      profile = fallback_profile;
      blocks = fallback_blocks;
    }

    if (!profile) {
      continue;
    }
    if (!load_profile(target, profile, profile_size, blocks, mode)) {
      return false;
    }
  }

  return true;
}

bool specbleach_migrate_profiles_1d_to_2d(specbleach_denoiser* source,
                                          specbleach_2d_denoiser* target) {
  if (!source || !target ||
      specbleach_denoiser_get_noise_profile_size(source) !=
          specbleach_2d_get_noise_profile_size(target)) {
    return false;
  }

  return migrate_profiles_generic(
      source, target, specbleach_denoiser_get_noise_profile_size(source),
      (bool (*)(void*,
                int))specbleach_denoiser_noise_profile_available_for_mode,
      (float* (*)(void*, int))specbleach_denoiser_get_noise_profile_for_mode,
      (uint32_t (*)(void*, int))
          specbleach_denoiser_get_noise_profile_block_count_for_mode,
      (bool (*)(void*, const float*, uint32_t, uint32_t,
                int))specbleach_2d_load_noise_profile_for_mode);
}

bool specbleach_migrate_profiles_2d_to_1d(specbleach_2d_denoiser* source,
                                          specbleach_denoiser* target) {
  if (!source || !target ||
      specbleach_2d_get_noise_profile_size(source) !=
          specbleach_denoiser_get_noise_profile_size(target)) {
    return false;
  }

  return migrate_profiles_generic(
      source, target, specbleach_2d_get_noise_profile_size(source),
      (bool (*)(void*, int))specbleach_2d_noise_profile_available_for_mode,
      (float* (*)(void*, int))specbleach_2d_get_noise_profile_for_mode,
      (uint32_t (*)(void*,
                    int))specbleach_2d_get_noise_profile_block_count_for_mode,
      (bool (*)(void*, const float*, uint32_t, uint32_t,
                int))specbleach_denoiser_load_noise_profile_for_mode);
}
