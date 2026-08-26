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

/**
 * Migrates noise profiles from a source object to a target object.
 *
 * Uses each mode's available profile and applies an available profile as a
 * fallback for modes without one. Fails when no valid source profile exists
 * or when a profile cannot be loaded into the target.
 *
 * @param source Source object containing the profiles.
 * @param target Target object receiving the profiles.
 * @param profile_size Number of samples in each profile.
 * @param available Callback that checks whether a profile is available for a mode.
 * @param get_profile Callback that retrieves a profile for a mode.
 * @param get_blocks Callback that retrieves the block count for a mode.
 * @param load_profile Callback that loads a profile into the target.
 * @return true if all applicable profiles are migrated successfully, false otherwise.
 */
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

/**
 * Determines whether a noise profile is available for a mode in a 1D denoiser.
 *
 * @param ctx 1D denoiser to query.
 * @param mode Mode whose profile availability is checked.
 * @return true if a noise profile is available for the mode, false otherwise.
 */
static bool wrap_1d_available(void* ctx, int mode) {
  return specbleach_denoiser_noise_profile_available_for_mode(
      (specbleach_denoiser*)ctx, mode);
}

static float* wrap_1d_get_profile(void* ctx, int mode) {
  return specbleach_denoiser_get_noise_profile_for_mode(
      (specbleach_denoiser*)ctx, mode);
}

/**
 * Gets the number of blocks in a 1D denoiser's noise profile for a mode.
 * @param ctx 1D denoiser instance.
 * @param mode Processing mode.
 * @returns The noise profile block count.
 */
static uint32_t wrap_1d_get_blocks(void* ctx, int mode) {
  return specbleach_denoiser_get_noise_profile_block_count_for_mode(
      (specbleach_denoiser*)ctx, mode);
}

/**
 * Loads a noise profile for a mode in a 1D denoiser.
 * @param ctx 1D denoiser receiving the profile.
 * @param profile Noise profile data.
 * @param size Number of values in the profile.
 * @param blocks Number of blocks represented by the profile.
 * @param mode Processing mode associated with the profile.
 * @return `true` if the profile is loaded successfully, `false` otherwise.
 */
static bool wrap_1d_load(void* ctx, const float* profile, uint32_t size,
                         uint32_t blocks, int mode) {
  return specbleach_denoiser_load_noise_profile_for_mode(
      (specbleach_denoiser*)ctx, profile, size, blocks, mode);
}

/**
 * Determines whether a 2D denoiser has a noise profile available for a mode.
 * @param ctx 2D denoiser to query.
 * @param mode Mode whose profile availability is checked.
 * @return `true` if a profile is available for the mode, `false` otherwise.
 */
static bool wrap_2d_available(void* ctx, int mode) {
  return specbleach_2d_noise_profile_available_for_mode(
      (specbleach_2d_denoiser*)ctx, mode);
}

/**
 * Retrieves the noise profile for a mode from a 2D denoiser.
 * @param ctx 2D denoiser instance.
 * @param mode Mode whose noise profile is requested.
 * @returns The noise profile for the specified mode.
 */
static float* wrap_2d_get_profile(void* ctx, int mode) {
  return specbleach_2d_get_noise_profile_for_mode(
      (specbleach_2d_denoiser*)ctx, mode);
}

/**
 * Gets the number of noise profile blocks for a 2D denoiser mode.
 *
 * @param ctx 2D denoiser context.
 * @param mode Mode whose noise profile block count is requested.
 * @return Number of noise profile blocks for the mode.
 */
static uint32_t wrap_2d_get_blocks(void* ctx, int mode) {
  return specbleach_2d_get_noise_profile_block_count_for_mode(
      (specbleach_2d_denoiser*)ctx, mode);
}

/**
 * Loads a noise profile for a mode on a 2D denoiser.
 * @param ctx 2D denoiser that receives the profile.
 * @param profile Noise profile data.
 * @param size Number of profile values.
 * @param blocks Number of blocks associated with the profile.
 * @param mode Mode for which to load the profile.
 * @returns `true` if the profile is loaded successfully, `false` otherwise.
 */
static bool wrap_2d_load(void* ctx, const float* profile, uint32_t size,
                         uint32_t blocks, int mode) {
  return specbleach_2d_load_noise_profile_for_mode(
      (specbleach_2d_denoiser*)ctx, profile, size, blocks, mode);
}

/**
 * Migrates noise profiles from a 1D denoiser to a 2D denoiser.
 *
 * @param source Source 1D denoiser.
 * @param target Target 2D denoiser.
 * @return `true` if all profiles are migrated successfully, `false` if either
 *         denoiser is null, their profile sizes differ, or migration fails.
 */
bool specbleach_migrate_profiles_1d_to_2d(specbleach_denoiser* source,
                                          specbleach_2d_denoiser* target) {
  if (!source || !target ||
      specbleach_denoiser_get_noise_profile_size(source) !=
          specbleach_2d_get_noise_profile_size(target)) {
    return false;
  }

  return migrate_profiles_generic(
      source, target, specbleach_denoiser_get_noise_profile_size(source),
      wrap_1d_available, wrap_1d_get_profile, wrap_1d_get_blocks, wrap_2d_load);
}

/**
 * Migrates noise profiles from a 2D denoiser to a 1D denoiser.
 *
 * @param source Source 2D denoiser.
 * @param target Target 1D denoiser.
 * @return true if all profiles are migrated successfully, false if either
 *         denoiser is null, their profile sizes differ, or migration fails.
 */
bool specbleach_migrate_profiles_2d_to_1d(specbleach_2d_denoiser* source,
                                          specbleach_denoiser* target) {
  if (!source || !target ||
      specbleach_2d_get_noise_profile_size(source) !=
          specbleach_denoiser_get_noise_profile_size(target)) {
    return false;
  }

  return migrate_profiles_generic(
      source, target, specbleach_2d_get_noise_profile_size(source),
      wrap_2d_available, wrap_2d_get_profile, wrap_2d_get_blocks,
      wrap_1d_load);
}
