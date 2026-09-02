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

#ifndef SPECBLEACH_STEREO_H_INCLUDED
#define SPECBLEACH_STEREO_H_INCLUDED

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_denoiser.h"
#include "specbleach_export.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Multi-channel orchestration over single-channel libspecbleach engines.
 *
 * This module belongs to the optional extras layer: it manages one engine
 * instance per channel, fans parameter loads out to every channel, and
 * provides the profile bookkeeping multi-channel applications otherwise
 * reimplement. It contains no DSP of its own.
 *
 * Real-time contract: same classification as specbleach_denoiser.h.
 * "RT-safe" functions (process, get_latency, read-only queries) carry no
 * allocations, locks, or I/O and are safe for the audio callback thread;
 * all other functions (create, free, parameter/profile loads, sync/reset)
 * are setup-only and must be called from a control thread. Calls on the
 * same instance must never run concurrently with process().
 */

/**
 * Opaque handle to a multi-channel engine group.
 */
typedef struct specbleach_stereo
    specbleach_stereo; // NOLINT(readability-identifier-naming)

/**
 * Creates one unified spectral denoiser instance per channel.
 *
 * @param channels Number of channels (1 supported too; useful for
 * symmetric code paths).
 * @return New instance or NULL on failure. Free with
 * specbleach_stereo_free().
 *
 * Thread safety: setup-only (allocates one engine per channel).
 */
SPECBLEACH_API specbleach_stereo* specbleach_stereo_initialize(
    uint32_t sample_rate, float frame_size, uint32_t channels);

/**
 * Frees the instance and every underlying engine. NULL is a no-op.
 *
 * Thread safety: setup-only (deallocates).
 */
SPECBLEACH_API void specbleach_stereo_free(specbleach_stereo* instance);

/**
 * Loads parameters on every channel.
 *
 * parameters_size must be exactly sizeof(SpecbleachDenoiserParameters) and
 * guards against ABI drift between separately compiled binaries.
 *
 * See the engine header for allocation caveats related to
 * reduction_curve_bias (each channel owns its copy; first enabled load may
 * allocate). The smoothing_mode field can be changed at runtime; the
 * library performs the internal seamless crossfade.
 *
 * Thread safety: setup-only (copies parameters to every channel; first
 * curve-enabled load may allocate).
 */
SPECBLEACH_API bool specbleach_stereo_load_parameters(
    specbleach_stereo* instance, const SpecbleachDenoiserParameters* parameters,
    uint32_t parameters_size);

/**
 * Processes all channels. Buffers are deinterleaved (planar): input[i] /
 * output[i] address channel i, each holding number_of_samples frames. Any
 * block size is accepted.
 *
 * Thread safety: RT-safe (no allocations, locks, or I/O). Safe for the
 * real-time audio callback thread.
 */
SPECBLEACH_API bool specbleach_stereo_process(specbleach_stereo* instance,
                                              uint32_t number_of_samples,
                                              const float** input,
                                              float** output);

/**
 * Returns the channel count passed at initialization.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_channel_count(const specbleach_stereo* instance);

/**
 * Returns the algorithmic latency in samples (identical for every channel).
 * Constant across smoothing mode changes; query once at prepare time for
 * host delay compensation.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_latency(const specbleach_stereo* instance);

/**
 * Fallback-fills noise profiles across channels: any channel missing a
 * profile for a mode copies it from the first channel that has one, so
 * every mode on every channel produces consistent reduction. Channels that
 * already have a profile keep their own.
 *
 * Not real-time safe: synchronizes profile state. Call from setup, state
 * restore, or a deferred control thread.
 *
 * @return true if every channel ended up with every mode populated.
 */
SPECBLEACH_API bool specbleach_stereo_sync_profiles(
    specbleach_stereo* instance);

/**
 * Resets stored noise profiles on every channel.
 *
 * Thread safety: setup-only.
 */
SPECBLEACH_API bool specbleach_stereo_reset_profiles(
    specbleach_stereo* instance);

/**
 * Per-channel profile accessors forwarding to the underlying engine.
 * Pointers stay valid until the instance is freed. mode must be within
 * [SPECBLEACH_PROFILE_MODE_FIRST, SPECBLEACH_PROFILE_MODE_LAST].
 *
 * Thread safety: getters are RT-safe (read-only); the load function is
 * setup-only (copies profile data).
 */
SPECBLEACH_API float* specbleach_stereo_get_noise_profile_for_channel(
    specbleach_stereo* instance, uint32_t channel, int mode);
SPECBLEACH_API bool specbleach_stereo_load_noise_profile_for_channel(
    specbleach_stereo* instance, uint32_t channel, const float* profile,
    uint32_t profile_size, uint32_t block_count, int mode);
SPECBLEACH_API bool specbleach_stereo_profile_available_for_channel(
    specbleach_stereo* instance, uint32_t channel, int mode);
SPECBLEACH_API uint32_t specbleach_stereo_get_profile_block_count_for_channel(
    specbleach_stereo* instance, uint32_t channel, int mode);

/**
 * Returns the size of the noise profile spectrum in bins.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_noise_profile_size(const specbleach_stereo* instance);

/**
 * Returns the maximum transient intensity [0.0, 1.0] across channels from
 * the last processed block.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API float specbleach_stereo_get_transient_intensity(
    specbleach_stereo* instance);

/**
 * Per-channel visualization/query forwarders. See the engine headers for
 * semantics and threading notes (tonal peak queries are offline-only,
 * setup-only; the mask/active-profile pointers are RT-safe read-only
 * queries into live internal state).
 */
SPECBLEACH_API const float* specbleach_stereo_get_tonal_mask_for_channel(
    specbleach_stereo* instance, uint32_t channel);
SPECBLEACH_API const float*
specbleach_stereo_get_active_noise_profile_for_channel(
    specbleach_stereo* instance, uint32_t channel);
SPECBLEACH_API uint32_t specbleach_stereo_get_tonal_peaks_for_channel(
    specbleach_stereo* instance, uint32_t channel, float* peak_freqs_hz,
    uint32_t max_peaks);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_STEREO_H_INCLUDED */
