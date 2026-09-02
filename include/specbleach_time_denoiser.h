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

#ifndef SPECBLEACH_TIME_DENOISER_H_INCLUDED
#define SPECBLEACH_TIME_DENOISER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_export.h"

/**
 * Adaptive noise estimation method (shared contract with the spectral
 * denoiser API).
 */
typedef enum SpecbleachNoiseEstimationMethod {
  SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE = 0, /**< SPP-MMSE (unbiased) */
  SPECBLEACH_NOISE_ESTIMATION_BRANDT = 1,   /**< Brandt trimmed mean */
  SPECBLEACH_NOISE_ESTIMATION_MARTIN = 2,   /**< Martin minimum statistics */
} SpecbleachNoiseEstimationMethod;

/**
 * Opaque handle to a single-channel zero-latency time-domain denoiser
 * (minimum-phase delayless subband adaptive filter).
 *
 * Handles are created by specbleach_time_denoiser_initialize() and must
 * only be passed to the specbleach_time_denoiser_* functions declared in
 * this header. They are NOT interchangeable with handles from other
 * libspecbleach APIs.
 */
typedef struct specbleach_time_denoiser specbleach_time_denoiser;

/**
 * Parameters for the time-domain denoiser.
 *
 * Real-time contract for the whole header: every function is classified as
 * either "RT-safe" (no allocations, no locks, no I/O — callable from a
 * real-time audio callback thread) or "setup-only" (may allocate, block,
 * or copy — call it from a control/setup thread only). Unless stated
 * otherwise, calls on the SAME instance must never run concurrently with
 * each other or with specbleach_time_denoiser_process().
 */
typedef struct SpecbleachTimeDenoiserParameters {
  /**
   * Linear gain floor for noise reduction (0.0 to 1.0, where 1.0 = 0 dB /
   * no reduction). Values outside the range are clamped.
   */
  float reduction_gain;

  /**
   * Normalized temporal smoothing of the ERB band gains (0.0 to 1.0).
   * Higher values adapt slower and sound steadier. Clamped.
   */
  float smoothing_factor;

  /**
   * Enables adaptive noise estimation, which continuously updates the
   * noise floor based on the current input signal.
   */
  bool adaptive_noise;

  /**
   * Sets the method used for adaptive noise estimation.
   */
  SpecbleachNoiseEstimationMethod noise_estimation_method;

  /**
   * Suppression aggressiveness / oversubtraction (0.0 - 1.0). Clamped.
   */
  float suppression_strength;
} SpecbleachTimeDenoiserParameters;

/**
 * Creates a single-channel zero-latency time-domain denoiser instance.
 *
 * Spawns an internal analysis worker thread. Algorithmic latency is always
 * zero regardless of configuration.
 *
 * Thread safety: setup-only (allocates and spawns a thread). Never call
 * from an audio thread.
 *
 * @return A new instance or NULL on failure. Free it with
 * specbleach_time_denoiser_free().
 */
SPECBLEACH_API specbleach_time_denoiser* specbleach_time_denoiser_initialize(
    uint32_t sample_rate);

/**
 * Frees an instance created by specbleach_time_denoiser_initialize()
 * and joins the internal worker thread. Passing NULL is a no-op. The
 * handle is invalid after this call.
 *
 * Thread safety: setup-only. Never call from an audio thread.
 */
SPECBLEACH_API void specbleach_time_denoiser_free(
    specbleach_time_denoiser* instance);

/**
 * Loads parameters for the reduction.
 *
 * @param parameters Pointer to the parameter block to load. Must not be
 * NULL. The library copies all data before returning.
 * @param parameters_size Must be exactly
 * sizeof(SpecbleachTimeDenoiserParameters). Any other value fails cleanly.
 *
 * Thread safety: setup-only.
 *
 * @return true if the parameters were loaded, false on NULL arguments or a
 * mismatched parameters_size.
 */
SPECBLEACH_API bool specbleach_time_denoiser_load_parameters(
    specbleach_time_denoiser* instance,
    const SpecbleachTimeDenoiserParameters* parameters,
    uint32_t parameters_size);

/**
 * Processes a buffer of samples with zero algorithmic latency.
 *
 * Thread safety: RT-safe (no allocations, locks, or I/O). Safe for the
 * real-time audio callback thread.
 *
 * Buffer contract: input/output are plain float arrays of
 * number_of_samples length (mono, non-interleaved). Any block size is
 * accepted. Output may alias input.
 *
 * @return true on success, false on NULL arguments or an empty block.
 */
SPECBLEACH_API bool specbleach_time_denoiser_process(
    specbleach_time_denoiser* instance, uint32_t number_of_samples,
    const float* input, float* output);

/**
 * Returns the algorithmic latency in samples. Always 0 for this
 * processor (delayless by construction); provided so hosts can query it
 * uniformly with the other denoiser APIs.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API uint32_t
specbleach_time_denoiser_get_latency(specbleach_time_denoiser* instance);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_TIME_DENOISER_H_INCLUDED */
