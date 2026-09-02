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
 * @file specbleach_stereo.h
 * @brief Multi-channel orchestration over single-channel engines.
 */

#ifndef SPECBLEACH_STEREO_H_INCLUDED
#define SPECBLEACH_STEREO_H_INCLUDED

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_denoiser.h"
#include "specbleach_error.h"
#include "specbleach_export.h"

/* NOTE: this header depends on the core include dir via the PUBLIC link
 * libspecbleach_extras -> libspecbleach (see CMakeLists.txt); the flat
 * #include "specbleach_denoiser.h" above resolves through that transitive
 * include path and the flattened install layout. */

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
 *
 * Channel convention: out-of-range channel indexes return the same "empty"
 * sentinel as valid-but-unpopulated state (NULL / false / 0). Callers
 * should bounds-check with channel < specbleach_stereo_get_channel_count().
 */

/**
 * Opaque handle to a multi-channel engine group.
 */
typedef struct specbleach_stereo
    specbleach_stereo; // NOLINT(readability-identifier-naming)

/**
 * Creates one unified spectral denoiser instance per channel.
 *
 * @param sample_rate Sample rate in Hz (e.g. 48000).
 * @param frame_size_ms STFT frame size in milliseconds (e.g. 46.0f).
 * @param channels Number of channels (1 supported too; useful for
 * symmetric code paths).
 * @return New instance or NULL on failure. Free with
 * specbleach_stereo_free().
 *
 * Thread safety: setup-only (allocates one engine per channel).
 */
SPECBLEACH_API specbleach_stereo* specbleach_stereo_initialize(
    uint32_t sample_rate, float frame_size_ms, uint32_t channels);

/**
 * Frees the instance and every underlying engine. NULL is a no-op.
 *
 * Thread safety: setup-only (deallocates).
 *
 * @param instance Instance to free, or NULL.
 */
SPECBLEACH_API void specbleach_stereo_free(specbleach_stereo* instance);

/**
 * Loads parameters on every channel.
 *
 * parameters_size must be exactly sizeof(SpecbleachDenoiserParameters)
 * (see SPECBLEACH_PARAMETERS_SIZE) and guards against ABI drift between
 * separately compiled binaries.
 *
 * See the engine header for allocation caveats related to
 * reduction_curve_bias (each channel owns its copy; first enabled load may
 * allocate). The smoothing_mode field can be changed at runtime; the
 * library performs the internal seamless crossfade.
 *
 * Thread safety: setup-only (copies parameters to every channel; first
 * curve-enabled load may allocate).
 *
 * @param instance Instance handle.
 * @param parameters Parameter block (copied to every channel).
 * @param parameters_size Must equal sizeof(SpecbleachDenoiserParameters).
 * @return true if loaded on every channel.
 */
SPECBLEACH_API bool specbleach_stereo_load_parameters(
    specbleach_stereo* instance, const SpecbleachDenoiserParameters* parameters,
    uint32_t parameters_size);

/**
 * Processes all channels. Buffers are deinterleaved (planar): input[i] /
 * output[i] address channel i, each holding number_of_samples frames. Any
 * block size is accepted. Every row pointer must be non-NULL (NULL rows
 * fail fast with false instead of crashing).
 *
 * Thread safety: RT-safe (no allocations, locks, or I/O). Safe for the
 * real-time audio callback thread.
 *
 * @param instance Instance handle.
 * @param number_of_samples Frames per channel (must be > 0).
 * @param input Array of channel pointers (size = channel count).
 * @param output Array of channel pointers (size = channel count).
 * @return true if every channel processed successfully.
 */
SPECBLEACH_API bool specbleach_stereo_process(specbleach_stereo* instance,
                                              uint32_t number_of_samples,
                                              const float** input,
                                              float** output);

/**
 * Returns the channel count passed at initialization.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Channel count, or 0 for NULL.
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_channel_count(const specbleach_stereo* instance);

/**
 * Returns the sample rate passed at initialization (identical for every
 * channel).
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Sample rate in Hz, or 0 for NULL.
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_sample_rate(const specbleach_stereo* instance);

/**
 * Returns the effective STFT frame size in samples (identical for every
 * channel). See specbleach_denoiser_initialize() for the ms-to-samples
 * mapping.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Frame size in samples, or 0 for NULL.
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_frame_size(const specbleach_stereo* instance);

/**
 * Returns the effective FFT size in samples (identical for every
 * channel).
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return FFT size in samples, or 0 for NULL.
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_fft_size(const specbleach_stereo* instance);

/**
 * Returns the STFT hop size in samples (identical for every channel).
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Hop size in samples, or 0 for NULL.
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_hop_size(const specbleach_stereo* instance);

/**
 * Returns the algorithmic latency in samples (identical for every channel).
 * Constant across smoothing mode changes; query once at prepare time for
 * host delay compensation.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Latency in samples, or 0 for NULL.
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_latency(const specbleach_stereo* instance);

/**
 * Fallback-fills noise profiles across channels: any channel missing a
 * profile for a mode copies it from the first channel that has one, so
 * every mode on every channel produces consistent reduction. Channels that
 * already have a profile keep their own.
 *
 * Thread safety: setup-only (copies profile state). Call from setup, state
 * restore, or a deferred control thread.
 *
 * @param instance Instance handle.
 * @return true if every channel ended up with every mode populated; false
 * means at least one mode had no source on any channel (nothing was
 * available to copy) — the already-populated modes are still synced.
 */
SPECBLEACH_API bool specbleach_stereo_sync_profiles(
    specbleach_stereo* instance);

/**
 * Resets stored noise profiles on every channel.
 *
 * Thread safety: setup-only.
 *
 * @param instance Instance handle (NULL is a no-op).
 */
SPECBLEACH_API void specbleach_stereo_reset_profiles(
    specbleach_stereo* instance);

/**
 * Flushes DSP history on every channel (see
 * specbleach_denoiser_reset_dsp_state()). Preserves loaded parameters.
 *
 * Thread safety: setup-only. Never call from an audio thread.
 *
 * @param instance Instance handle.
 * @return true if every channel reset, false on NULL/rebuild failure.
 */
SPECBLEACH_API bool specbleach_stereo_reset_dsp_state(
    specbleach_stereo* instance);

/**
 * Returns a read-only pointer to the stored noise profile for a channel
 * and mode. NULL means out-of-range channel/mode or NULL instance.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API const float* specbleach_stereo_get_noise_profile_for_channel(
    const specbleach_stereo* instance, uint32_t channel,
    SpecbleachProfileMode mode);

/**
 * Copies a profile into one channel. See the mono header for size rules.
 *
 * Thread safety: setup-only (copies profile data).
 *
 * @return false on out-of-range channel/mode, size mismatch, or NULL args.
 */
SPECBLEACH_API bool specbleach_stereo_load_noise_profile_for_channel(
    specbleach_stereo* instance, uint32_t channel, const float* profile,
    uint32_t profile_size, uint32_t block_count, SpecbleachProfileMode mode);

/**
 * Copies a profile into one channel, linearly resampling to the native
 * profile size when sizes differ.
 *
 * Thread safety: setup-only (copies profile data).
 *
 * @return false on out-of-range channel/mode, NULL args, or OOM.
 */
SPECBLEACH_API bool specbleach_stereo_load_noise_profile_resampled_for_channel(
    specbleach_stereo* instance, uint32_t channel, const float* profile,
    uint32_t source_size, uint32_t block_count, SpecbleachProfileMode mode);

/**
 * Returns whether a channel has a profile for a mode. False also covers
 * out-of-range channel/mode — bounds-check first if you need to tell them
 * apart.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API bool specbleach_stereo_profile_available_for_channel(
    const specbleach_stereo* instance, uint32_t channel,
    SpecbleachProfileMode mode);

/**
 * Returns true if ANY mode on the channel is populated (convenience over
 * looping FIRST..LAST).
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API bool specbleach_stereo_has_any_profile_for_channel(
    const specbleach_stereo* instance, uint32_t channel);

/**
 * Returns the block count for a channel/mode, or 0 on error/OOB.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API uint32_t specbleach_stereo_get_profile_block_count_for_channel(
    const specbleach_stereo* instance, uint32_t channel,
    SpecbleachProfileMode mode);

/**
 * Returns the size of the noise profile spectrum in bins.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Profile size in bins, or 0 for NULL.
 */
SPECBLEACH_API uint32_t
specbleach_stereo_get_noise_profile_size(const specbleach_stereo* instance);

/**
 * Returns true if any channel detected a transient in the last block.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API bool specbleach_stereo_is_transient_detected(
    const specbleach_stereo* instance);

/**
 * Returns true if the given channel detected a transient in the last
 * block. False also covers out-of-range channels.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API bool specbleach_stereo_is_transient_detected_for_channel(
    const specbleach_stereo* instance, uint32_t channel);

/**
 * Returns the maximum transient intensity [0.0, 1.0] across channels from
 * the last processed block.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Maximum intensity, or 0.0f for NULL.
 */
SPECBLEACH_API float specbleach_stereo_get_transient_intensity(
    const specbleach_stereo* instance);

/**
 * Returns the read-only tonal mask for a channel (NULL on OOB). Array size
 * is specbleach_stereo_get_noise_profile_size(); values 0.0..1.0.
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API const float* specbleach_stereo_get_tonal_mask_for_channel(
    const specbleach_stereo* instance, uint32_t channel);

/**
 * Returns the read-only active (morphed) noise profile for a channel
 * (NULL on OOB).
 *
 * Thread safety: RT-safe (read-only query).
 */
SPECBLEACH_API const float*
specbleach_stereo_get_active_noise_profile_for_channel(
    const specbleach_stereo* instance, uint32_t channel);

/**
 * Writes up to max_peaks tonal peak frequencies (Hz) for a channel.
 *
 * Thread safety: RT-safe (read-only query into the published
 * double-buffer). Safe for the audio thread at decimated rates.
 *
 * @return Number of peaks written (0 = none found or OOB/NULL).
 */
SPECBLEACH_API uint32_t specbleach_stereo_get_tonal_peaks_for_channel(
    const specbleach_stereo* instance, uint32_t channel, float* peak_freqs_hz,
    uint32_t max_peaks);

/**
 * Failure reason for the last fallible stereo call (channel 0's error when
 * the failure is per-channel). RT-safe to read.
 */
SPECBLEACH_API SpecbleachError
specbleach_stereo_get_last_error(const specbleach_stereo* instance);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_STEREO_H_INCLUDED */
