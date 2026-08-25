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

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_2d_denoiser.h"
#include "specbleach_denoiser.h"
#include "specbleach_export.h"

/**
 * Multi-channel orchestration over single-channel libspecbleach engines.
 *
 * This module belongs to the optional extras layer: it manages one engine
 * instance per channel, fans parameter loads out to every channel, and
 * provides the profile bookkeeping multi-channel applications otherwise
 * reimplement. It contains no DSP of its own.
 */

/**
 * Selects which engine family the instances wrap.
 */
typedef enum SpecbleachStereoEngine {
  SPECBLEACH_STEREO_ENGINE_SPECTRAL = 0, /**< STFT spectral denoiser (1D) */
  SPECBLEACH_STEREO_ENGINE_NLM_2D = 1,   /**< 2D Non-Local Means denoiser */
} SpecbleachStereoEngine;

/**
 * Opaque handle to a multi-channel engine group.
 */
typedef struct specbleach_stereo specbleach_stereo;

/**
 * Creates one engine instance per channel.
 *
 * @param channels Number of channels (1 supported too; useful for
 * symmetric code paths).
 * @return New instance or NULL on failure. Free with
 * specbleach_stereo_free().
 */
SPECBLEACH_API specbleach_stereo* specbleach_stereo_initialize(
    uint32_t sample_rate, float frame_size, uint32_t channels,
    SpecbleachStereoEngine engine);

/**
 * Frees the instance and every underlying engine. NULL is a no-op.
 */
SPECBLEACH_API void specbleach_stereo_free(specbleach_stereo* instance);

/**
 * Loads parameters on every channel.
 *
 * These are strongly typed per engine family: calling the 1D loader with a
 * 2D group (or vice versa) is a compile-time error. parameters_size must be
 * exactly the matching sizeof and guards against ABI drift between
 * separately compiled binaries.
 *
 * See the engine headers for allocation caveats related to
 * reduction_curve_bias (each channel owns its copy; first enabled load may
 * allocate).
 */
SPECBLEACH_API bool specbleach_stereo_load_parameters_1d(
    specbleach_stereo* instance, const SpecbleachDenoiserParameters* parameters,
    uint32_t parameters_size);
SPECBLEACH_API bool specbleach_stereo_load_parameters_2d(
    specbleach_stereo* instance,
    const Specbleach2DDenoiserParameters* parameters, uint32_t parameters_size);

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
 * Returns the engine family passed at initialization.
 */
SPECBLEACH_API SpecbleachStereoEngine
specbleach_stereo_get_engine(const specbleach_stereo* instance);

/**
 * Returns the algorithmic latency in samples (identical for every channel).
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
 * Copies noise profiles from another engine family's group into this one,
 * channel by channel (e.g., before switching a spectral group for a 2D NLM
 * group). Both groups must share the same channel count and profile size.
 *
 * Not real-time safe.
 *
 * @return true if every channel was migrated successfully. Fails cleanly
 * when the source group has not learned anything yet.
 */
SPECBLEACH_API bool specbleach_stereo_migrate_profiles_from(
    specbleach_stereo* instance, const specbleach_stereo* source);

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
 * Aggregated transient detection across channels (true if ANY channel
 * detected a transient in the last processed block).
 */
SPECBLEACH_API bool specbleach_stereo_is_transient_detected(
    specbleach_stereo* instance);

/**
 * Returns the maximum transient intensity [0.0, 1.0] across channels from
 * the last processed block.
 */
SPECBLEACH_API float specbleach_stereo_get_transient_intensity(
    specbleach_stereo* instance);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_STEREO_H_INCLUDED */
