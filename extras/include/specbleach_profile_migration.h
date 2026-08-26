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

#ifndef SPECBLEACH_PROFILE_MIGRATION_H_INCLUDED
#define SPECBLEACH_PROFILE_MIGRATION_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "specbleach_2d_denoiser.h"
#include "specbleach_denoiser.h"
#include "specbleach_export.h"

/**
 * Copies noise profiles between different engine families so an engine
 * switch does not lose the learned noise profile.
 *
 * For every profile mode: modes available on the source are copied with
 * their block counts; missing modes are fallback-filled from the first
 * available mode so the target never ends up partially empty.
 *
 * Not real-time safe. Call before switching engines, before the target
 * group renders (see examples/README.md "Switching between processors").
 */

/**
 * Migrates all noise profiles from a spectral (1D) engine to a 2D NLM
 * engine.
 * @return false if either handle is NULL, the profile sizes differ, or the
 * source has no profile for any mode yet.
 */
SPECBLEACH_API bool specbleach_migrate_profiles_1d_to_2d(
    specbleach_denoiser* source, specbleach_2d_denoiser* target);

/**
 * Migrates all noise profiles from a 2D NLM engine to a spectral (1D)
 * engine.
 * @return false if either handle is NULL, the profile sizes differ, or the
 * source has no profile for any mode yet.
 */
SPECBLEACH_API bool specbleach_migrate_profiles_2d_to_1d(
    specbleach_2d_denoiser* source, specbleach_denoiser* target);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_PROFILE_MIGRATION_H_INCLUDED */
