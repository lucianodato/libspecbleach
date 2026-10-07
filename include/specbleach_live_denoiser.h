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

#ifndef SPECBLEACH_LIVE_DENOISER_H_INCLUDED
#define SPECBLEACH_LIVE_DENOISER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_export.h"

/**
 * Adaptive speed preset for the per-band noise floor tracker.
 */
typedef enum SpecbleachLiveEstimationMethod {
  SPECBLEACH_LIVE_SPP_MMSE = 0, /**< fast adaptation */
  SPECBLEACH_LIVE_BRANDT = 1,   /**< medium adaptation */
  SPECBLEACH_LIVE_MARTIN = 2,   /**< slow adaptation */
} SpecbleachLiveEstimationMethod;

/**
 * Opaque handle to a single-channel zero-latency time-domain denoiser
 * (256-band Bark-spaced multiband gate).
 */
typedef struct specbleach_live_denoiser specbleach_live_denoiser;

/**
 * Band count of the live filterbank. Sizes the level arrays used by
 * specbleach_live_denoiser_get_band_levels().
 */
#define SPECBLEACH_LIVE_NUM_BANDS 256U

/**
 * Parameters for the live denoiser.
 *
 * Real-time contract: initialize/free/load_parameters are setup-only
 * (may allocate, block, or copy); process/get_latency are RT-safe.
 * Calls on the SAME instance must never run concurrently with each other
 * or with specbleach_live_denoiser_process().
 */
typedef struct SpecbleachLiveDenoiserParameters {
  /**
   * Linear gain floor for attenuated bands (0.0 to 1.0, 1.0 = no
   * reduction). Values outside the range are clamped.
   */
  float reduction_gain;

  /**
   * Gate opening time in seconds (attack). Clamped to
   * [LIVE_GATE_ATTACK_MIN_SEC, LIVE_GATE_ATTACK_MAX_SEC].
   */
  float attack_time;

  /**
   * Gate closing time in seconds (release). Higher values close
   * slower and sound steadier. Clamped to
   * [LIVE_GATE_RELEASE_MIN_SEC, LIVE_GATE_RELEASE_MAX_SEC].
   */
  float release_time;

  /**
   * Enables continuous noise-floor tracking. When false the floor
   * learned so far is frozen.
   */
  bool adaptive_noise;

  /**
   * Adaptation speed preset for the noise floor tracker.
   */
  SpecbleachLiveEstimationMethod noise_estimation_method;

  /**
   * Gate threshold offset in dB (like the full denoiser profile
   * offset). Higher values gate more. Clamped to
   * [LIVE_THRESHOLD_DB_MIN, LIVE_THRESHOLD_DB_MAX].
   */
  float threshold_db;

  /**
   * Soft-knee width in dB below the threshold. 0 is a hard gate.
   * Clamped to [LIVE_GATE_KNEE_MIN_DB, LIVE_GATE_KNEE_MAX_DB].
   */
  float knee_db;
} SpecbleachLiveDenoiserParameters;

/**
 * Creates a single-channel zero-latency live denoiser instance.
 * Fully synchronous: no worker thread. Algorithmic latency is always
 * zero regardless of configuration.
 *
 * Thread safety: setup-only. Never call from an audio thread.
 */
SPECBLEACH_API specbleach_live_denoiser* specbleach_live_denoiser_initialize(
    uint32_t sample_rate);

/**
 * Frees an instance. Passing NULL is a no-op.
 *
 * Thread safety: setup-only. Never call from an audio thread.
 */
SPECBLEACH_API void specbleach_live_denoiser_free(
    specbleach_live_denoiser* instance);

/**
 * Loads parameters. The library copies all data before returning.
 *
 * @param parameters_size Must be exactly
 * sizeof(SpecbleachLiveDenoiserParameters).
 *
 * Thread safety: setup-only.
 */
SPECBLEACH_API bool specbleach_live_denoiser_load_parameters(
    specbleach_live_denoiser* instance,
    const SpecbleachLiveDenoiserParameters* parameters,
    uint32_t parameters_size);

/**
 * Processes a buffer of samples with zero algorithmic latency.
 *
 * Buffer contract: mono non-interleaved float arrays of
 * number_of_samples length. Any block size is accepted. Output may
 * alias input.
 *
 * Thread safety: RT-safe (no allocations, locks, or I/O).
 */
SPECBLEACH_API bool specbleach_live_denoiser_process(
    specbleach_live_denoiser* instance, uint32_t number_of_samples,
    const float* input, float* output);

/**
 * Returns the algorithmic latency in samples. Always 0 (causal by
 * construction).
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API uint32_t
specbleach_live_denoiser_get_latency(specbleach_live_denoiser* instance);

/**
 * Copies the current per-band levels for spectrum display (linear
 * amplitudes, SPECBLEACH_LIVE_NUM_BANDS entries each; any array may be
 * NULL to skip it):
 * - input: band input energy (envelope follower state),
 * - output: post-gate band energy (gain times envelope),
 * - threshold: the actual gate threshold (learned noise floor times
 *   the aggressiveness multiplier).
 *
 * Thread safety: RT-safe (plain state copy). Call from the same thread
 * as process, or while the instance is idle.
 */
SPECBLEACH_API bool specbleach_live_denoiser_get_band_levels(
    specbleach_live_denoiser* instance, float* input, float* output,
    float* threshold);

/**
 * Re-arms the noise floor learner: the per-band floor decays back to
 * zero and re-converges on the current input (the "Learn" action).
 *
 * Thread safety: RT-safe (sets a flag honored by the next process call).
 */
SPECBLEACH_API bool specbleach_live_denoiser_reset_noise_floor(
    specbleach_live_denoiser* instance);

/**
 * Copies the filterbank band edges in Hz (SPECBLEACH_LIVE_NUM_BANDS
 * lower/upper pairs) so displays can use the exact psychoacoustic scale
 * the engine processes on. Inactive bands (above the design range) report
 * 0/0. Either array may be NULL to skip it.
 *
 * Thread safety: RT-safe (read-only query of setup state).
 */
SPECBLEACH_API bool specbleach_live_denoiser_get_band_edges(
    specbleach_live_denoiser* instance, float* lower_hz, float* upper_hz);

/**
 * Enables delta monitoring: when true, process() outputs the removed
 * noise (the attenuated part of the signal) instead of the denoised
 * signal, i.e. sum (1 - gain) * band / norm. When false, it outputs
 * the normal denoised signal. Useful for auditioning what is being
 * removed.
 *
 * Thread safety: RT-safe (atomic flag).
 */
SPECBLEACH_API bool specbleach_live_denoiser_set_delta_monitoring(
    specbleach_live_denoiser* instance, bool enabled);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_LIVE_DENOISER_H_INCLUDED */
