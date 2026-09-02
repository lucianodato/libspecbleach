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
 */
SPECBLEACH_API specbleach_stereo* specbleach_stereo_initialize(
    uint32_t sample_rate, float frame_size, uint32_t channels);

/**
 * Frees the instance and every underlying engine. NULL is a no-op.
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
 */
SPECBLEACH_API bool specbleach_stereo_load_parameters(
    specbleach_stereo* instance, const SpecbleachDenoiserParameters* parameters,
    uint32_t parameters_size);

/**
 * Processes all channels. Buffers are deinterleaved: input[i] / output[i]
 * address channel i, each holding number_of_samples frames. Real-time safe;
 * may allocate-free on first curve-enabled load as documented above.
 */
SPECBLEACH_API bool specbleach_stereo_process(specbleach_stereo* instance,
                                              uint32_t number_of_samples,
                                              const float** input,
                                              float** output);

/**
 * Returns the channel count passed at initialization.
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_channel_count(const specbleach_stereo* instance);

/**
 * Returns the algorithmic latency in samples (identical for every channel).
 * Constant across smoothing mode changes.
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
 */
SPECBLEACH_API bool specbleach_stereo_reset_profiles(
    specbleach_stereo* instance);

/**
 * Per-channel profile accessors forwarding to the underlying engine.
 * Pointers stay valid until the instance is freed. mode must be within
 * [SPECBLEACH_PROFILE_MODE_FIRST, SPECBLEACH_PROFILE_MODE_LAST].
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
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_noise_profile_size(const specbleach_stereo* instance);

/**
 * Returns the maximum transient intensity [0.0, 1.0] across channels from
 * the last processed block.
 */
SPECBLEACH_API float specbleach_stereo_get_transient_intensity(
    specbleach_stereo* instance);

/**
 * Per-channel visualization/query forwarders. See the engine headers for
 * semantics and threading notes (tonal peak queries are offline-only).
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
