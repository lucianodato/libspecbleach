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

#ifndef SPECBLEACH_COMMON_H_INCLUDED
#define SPECBLEACH_COMMON_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Noise profile learning mode.
 *
 * Numeric values are part of the public contract. Do not change them:
 * persisted state in downstream applications depends on them.
 */
typedef enum SpecbleachLearnMode {
  SPECBLEACH_LEARN_OFF =
      0, /**< Learning disabled, reduce using stored profile */
  SPECBLEACH_LEARN_ALL = 1, /**< Learn all profile types simultaneously */
} SpecbleachLearnMode;

/**
 * Adaptive noise estimation method.
 *
 * Numeric values are part of the public contract. Do not change them:
 * persisted state in downstream applications depends on them.
 */
typedef enum SpecbleachNoiseEstimationMethod {
  SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE = 0, /**< SPP-MMSE (unbiased) */
  SPECBLEACH_NOISE_ESTIMATION_BRANDT = 1,   /**< Brandt trimmed mean */
  SPECBLEACH_NOISE_ESTIMATION_MARTIN = 2,   /**< Martin minimum statistics */
} SpecbleachNoiseEstimationMethod;

/**
 * Smoothing strategy applied by the unified spectral denoiser.
 *
 * Numeric values are part of the public contract.
 */
typedef enum SbSmoothingMode {
  SB_SMOOTHING_TEMPORAL = 0, /**< 1D temporal/spatial gain smoothing */
  SB_SMOOTHING_NLM_2D = 1,   /**< 2D Non-Local Means patch smoothing */
} SbSmoothingMode;

/**
 * Inclusive range of valid noise profile mode indexes used by the
 * profile accessors across all denoiser APIs.
 */
#define SPECBLEACH_PROFILE_MODE_FIRST 1
#define SPECBLEACH_PROFILE_MODE_LAST 4

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_COMMON_H_INCLUDED */
