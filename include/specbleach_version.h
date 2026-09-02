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

/**
 * @file specbleach_version.h
 * @brief Runtime version information for libspecbleach.
 *
 * The macros below are the single source of truth mirrored from the CMake
 * project version. CMake verifies the match at configure time, so a build
 * cannot silently drift out of sync. The accessor exists so applications
 * linking against a system-installed library (or deciding between a bundled
 * and a system copy at runtime) can query what they actually got.
 *
 * Thread-safety: all functions in this header are thread-safe and may be
 * called before any instance is created.
 */

#ifndef SPECBLEACH_VERSION_H_INCLUDED
#define SPECBLEACH_VERSION_H_INCLUDED

#include "specbleach_export.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Major version component of the library release.
 */
#define SPECBLEACH_VERSION_MAJOR 0

/**
 * Minor version component of the library release.
 */
#define SPECBLEACH_VERSION_MINOR 4

/**
 * Patch version component of the library release.
 */
#define SPECBLEACH_VERSION_PATCH 0

/**
 * Release version as a string literal (e.g., "0.4.0").
 */
#define SPECBLEACH_VERSION_STRING "0.4.0"

/**
 * Composite numeric version for #if comparisons, e.g.
 * #if SPECBLEACH_VERSION_NUMBER >= 400.
 */
#define SPECBLEACH_VERSION_NUMBER                                              \
  (SPECBLEACH_VERSION_MAJOR * 10000 + SPECBLEACH_VERSION_MINOR * 100 +         \
   SPECBLEACH_VERSION_PATCH)

/**
 * Returns the runtime version banner ("libspecbleach " +
 * SPECBLEACH_VERSION_STRING). The returned pointer is static and owned by
 * the library; callers must not free or modify it.
 */
SPECBLEACH_API const char* specbleach_get_version_string(void);

#ifdef __cplusplus
}
#endif

#endif /* SPECBLEACH_VERSION_H_INCLUDED */
