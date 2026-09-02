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
 * @file specbleach_error.h
 * @brief Shared error codes for all libspecbleach public APIs.
 *
 * Both the single-channel denoiser (specbleach_denoiser.h) and the
 * multi-channel group (specbleach_stereo.h) report failures through this
 * one enum, so it lives here rather than in either API header. Each API
 * still owns its own get_last_error() reader; the human-readable form of
 * any code comes from specbleach_denoiser_get_last_error_string().
 */

#ifndef SPECBLEACH_ERROR_H_INCLUDED
#define SPECBLEACH_ERROR_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Machine-readable failure reason for the last fallible call on an instance.
 * Query per API with its get_last_error() reader; human-readable form via
 * specbleach_denoiser_get_last_error_string(). RT-safe to read.
 */
typedef enum SpecbleachError {
  SPECBLEACH_OK = 0,               /**< Last call succeeded */
  SPECBLEACH_ERR_NULL_ARG = 1,     /**< NULL instance/buffer argument */
  SPECBLEACH_ERR_ABI_MISMATCH = 2, /**< parameters_size != sizeof(parameters) */
  SPECBLEACH_ERR_SIZE_MISMATCH = 3,   /**< profile/curve size mismatch */
  SPECBLEACH_ERR_INVALID_MODE = 4,    /**< mode outside FIRST..LAST */
  SPECBLEACH_ERR_INVALID_CHANNEL = 5, /**< channel >= channel count (stereo) */
  SPECBLEACH_ERR_NO_MEMORY = 6,       /**< internal (re)allocation failed */
  SPECBLEACH_ERR_EMPTY = 7, /**< Empty block (number_of_samples == 0) */
} SpecbleachError;

#ifdef __cplusplus
}
#endif

#endif /* SPECBLEACH_ERROR_H_INCLUDED */
