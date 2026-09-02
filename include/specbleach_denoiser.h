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
 * @file specbleach_denoiser.h
 * @brief Single-channel spectral denoiser public API.
 *
 * Naming convention: opaque handles are lowercase (specbleach_denoiser);
 * parameter structs and enums are CamelCase (Specbleach*); enum values and
 * macros are SPECBLEACH_*.
 */

#ifndef SPECBLEACH_DENOISER_H_INCLUDED
#define SPECBLEACH_DENOISER_H_INCLUDED

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_export.h"

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
typedef enum SpecbleachSmoothingMode {
  SPECBLEACH_SMOOTHING_TEMPORAL = 0, /**< 1D temporal/spatial gain smoothing */
  SPECBLEACH_SMOOTHING_NLM_2D = 1,   /**< 2D Non-Local Means patch smoothing */
} SpecbleachSmoothingMode;

/**
 * Noise profile mode index.
 *
 * Numeric values are part of the public contract (1..4). They match the
 * internal ROLLING_MEAN/MEDIAN/STD_DEV/CV_MASK slots.
 */
typedef enum SpecbleachProfileMode {
  SPECBLEACH_PROFILE_ROLLING_MEAN = 1, /**< Rolling mean profile */
  SPECBLEACH_PROFILE_MEDIAN = 2,       /**< Median profile */
  SPECBLEACH_PROFILE_STD_DEV = 3,      /**< Standard deviation profile */
  SPECBLEACH_PROFILE_CV_MASK = 4,      /**< Coefficient-of-variation mask */
} SpecbleachProfileMode;

/**
 * Inclusive range of valid noise profile mode indexes used by the
 * profile accessors across all denoiser APIs.
 */
#define SPECBLEACH_PROFILE_MODE_FIRST 1
#define SPECBLEACH_PROFILE_MODE_LAST 4

/**
 * Machine-readable failure reason for the last fallible call on an instance.
 * Query with specbleach_denoiser_get_last_error(); human-readable form via
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

/**
 * Opaque handle to a single-channel spectral denoiser instance.
 *
 * The underlying struct is private. Handles are created by
 * specbleach_denoiser_initialize() and must only be passed to the
 * specbleach_denoiser_* functions declared in this header. They are NOT
 * interchangeable with handles from other libspecbleach APIs; mixing them is
 * a compile-time error.
 */
typedef struct specbleach_denoiser specbleach_denoiser;

/**
 * Parameters for the spectral denoiser.
 *
 * Real-time contract for the whole header: every function is classified as
 * either "RT-safe" (no allocations, no locks, no I/O — callable from a
 * real-time audio callback thread) or "setup-only" (may allocate, block, or
 * copy — call it from a control/setup thread only). Unless stated otherwise,
 * calls on the SAME instance must never run concurrently with each other or
 * with specbleach_denoiser_process().
 *
 * All sizes are bounded by UINT32_MAX samples; size_t is intentionally not
 * used in this audio-domain API.
 */
typedef struct SpecbleachDenoiserParameters {
  /**
   * Sets the processor in listening mode to capture the noise profile.
   */
  SpecbleachLearnMode learn_noise;

  /**
   * Outputs the residue of the reduction processing instead of the denoised
   * signal.
   */
  bool residual_listen;

  /**
   * Linear gain floor for broadband noise reduction (0.0 to 1.0, where 1.0 =
   * 0 dB / no reduction).
   */
  float reduction_gain;

  /**
   * Normalized smoothing factor across frames (0.0 to 1.0). Values outside
   * the range are clamped.
   *
   * In temporal mode this controls temporal/spatial gain smoothing. In NLM
   * 2D mode it controls the NLM h parameter and 2D time-frequency artifact
   * smoothing.
   */
  float smoothing_factor;

  /**
   * Smoothing strategy used by the processor. Can be changed at runtime;
   * the library performs an internal allocation-free crossfade and both
   * modes share the same latency.
   */
  SpecbleachSmoothingMode smoothing_mode;

  /**
   * Normalized whitening factor for residue noise floor (0.0 to 1.0). Values
   * outside the range are clamped.
   */
  float whitening_factor;

  /**
   * Enables adaptive noise estimation, which continuously updates the noise
   * profile based on the current input signal.
   */
  bool adaptive_noise;

  /**
   * Sets the method used for adaptive noise estimation.
   */
  SpecbleachNoiseEstimationMethod noise_estimation_method;

  /**
   * Masking veto/protection depth (0.0 - 1.0): depth of signal energy
   * preservation. Applies to both smoothing modes. Values outside the range
   * are clamped.
   */
  float masking_depth;

  /**
   * Suppression aggressiveness (0.0 - 1.0): Berouti oversubtraction factor.
   * Values outside the range are clamped.
   */
  float suppression_strength;

  /**
   * Intelligent steering: -1.0 (Median/Min) to 1.0 (Max), 0.0 (Mean).
   */
  float aggressiveness;

  /**
   * Linear gain floor for tonal noise components (0.0 to 1.0, where 1.0 = 0
   * dB / no reduction). Values outside the range are clamped.
   */
  float tonal_reduction_gain;

  /**
   * Enables HPSS transient protection.
   */
  bool hpss_enable;

  /**
   * Noise Profile Linear Scale — multiplier for the noise power spectrum.
   * Positive inputs are clamped to [0.01f, 100.0f]. Non-positive inputs use
   * the 1.0f default. Values > 1.0 shift threshold higher (more reduction).
   */
  float noise_profile_scale;

  /**
   * Frequency-dependent reduction bias curve.
   * Array of dB offsets per frequency bin, or NULL if disabled.
   * Positive values = more reduction at that frequency.
   *
   * OWNERSHIP: the library copies this array during
   * specbleach_denoiser_load_parameters(). The caller may free or reuse the
   * buffer as soon as that call returns. Array length must match
   * specbleach_denoiser_get_noise_profile_size().
   */
  const float* reduction_curve_bias;

  /**
   * Enables the reduction curve bias. When false, reduction_curve_bias is
   * ignored even if non-NULL.
   */
  bool reduction_curve_enabled;

  /**
   * Length of the reduction_curve_bias array in bins. Must match
   * specbleach_denoiser_get_noise_profile_size() when
   * reduction_curve_enabled is true; ignored otherwise.
   */
  uint32_t reduction_curve_size;

  /**
   * Tonal Noise Profile Linear Scale — multiplier applied to the noise power
   * spectrum only at detected tonal components. Positive inputs are clamped
   * to [0.01f, 100.0f]. Non-positive inputs use the 1.0f default. Values >
   * 1.0 shift threshold higher (more reduction) at tonal bins.
   */
  float tonal_noise_profile_scale;
} SpecbleachDenoiserParameters;

/**
 * Compile-time ABI guard value: always sizeof(SpecbleachDenoiserParameters).
 * Pass SPECBLEACH_PARAMETERS_SIZE (or sizeof your local copy) as
 * parameters_size to specbleach_denoiser_load_parameters().
 */
#define SPECBLEACH_PARAMETERS_SIZE (sizeof(SpecbleachDenoiserParameters))

/**
 * Returns a safe baseline parameter block. Prefer this over "= {0}":
 * zero-initialized reduction_gain means maximum reduction, not "no change".
 *
 * Intended use: call once at init, override the fields you care about, then
 * cache the struct and mutate it (e.g. per UI change). It is a plain struct
 * copy with no allocation, but there is no reason to call it per audio
 * block — keep your own copy and reload it with load_parameters().
 *
 * Thread safety: thread-safe, callable before any instance exists.
 *
 * @return Baseline parameters (learn off, moderate reduction, temporal
 * smoothing, scales at 1.0, curve disabled).
 */
SPECBLEACH_API SpecbleachDenoiserParameters
specbleach_denoiser_get_default_parameters(void);

/**
 * Creates a single-channel spectral denoiser instance.
 *
 * Sample rate can be anything from 4000 Hz to 192 kHz. Recommended frame
 * size range is between 20 ms and 100 ms.
 *
 * Frame-to-FFT mapping: frame_size_ms is converted to samples
 * (frame_ms / 1000 * sample_rate, truncated) and the FFT size is the
 * smallest PFFFT-valid size >= that frame. PFFFT does NOT require a power
 * of two — a size is valid when it is a multiple of 32 whose remaining
 * prime factors are only 2, 3, or 5 (e.g. 50 ms at 48 kHz gives exactly
 * 2400 samples, used as-is). The frame is centered and zero-padded when
 * the FFT size exceeds it. Padding guarantee: FFT buffer regions outside
 * the centered frame are always exactly zero, independent of the window
 * type, and padding never affects latency. Reported latency is frame_size
 * plus the NLM look-ahead (NLM_SEARCH_RANGE_TIME_FUTURE frames x hop),
 * constant across smoothing modes; query it at runtime with
 * specbleach_denoiser_get_latency(). Query the
 * effective values with specbleach_denoiser_get_frame_size(),
 * specbleach_denoiser_get_fft_size(), and
 * specbleach_denoiser_get_hop_size().
 *
 * Thread safety: setup-only (allocates). Never call from an audio thread.
 *
 * @param sample_rate Sample rate in Hz (e.g. 48000).
 * @param frame_size_ms STFT frame size in milliseconds (e.g. 46.0f). Must
 * yield at least OVERLAP (4) samples, otherwise initialization fails.
 * @return A new instance or NULL on allocation failure. Free it with
 * specbleach_denoiser_free().
 */
SPECBLEACH_API specbleach_denoiser* specbleach_denoiser_initialize(
    uint32_t sample_rate, float frame_size_ms);

/**
 * Frees an instance created by specbleach_denoiser_initialize().
 * Passing NULL is a no-op. The handle is invalid after this call.
 *
 * Thread safety: setup-only (deallocates). Never call from an audio thread.
 *
 * @param instance Instance to free, or NULL.
 */
SPECBLEACH_API void specbleach_denoiser_free(specbleach_denoiser* instance);

/**
 * Loads parameters for the reduction.
 *
 * @param instance Instance handle. Must not be NULL.
 * @param parameters Pointer to the parameter block to load. Must not be
 * NULL. The library copies all data it needs, including the reduction curve,
 * before returning.
 * @param parameters_size Must be exactly
 * sizeof(SpecbleachDenoiserParameters) (see SPECBLEACH_PARAMETERS_SIZE).
 * This guards against ABI mismatches between separately compiled binaries;
 * any other value fails cleanly instead of reading out-of-bounds memory.
 *
 * Thread safety: setup-only. May allocate on the FIRST call after enabling
 * reduction_curve_enabled (the internal copy buffer is allocated then);
 * subsequent calls reuse that buffer. Load parameters once during setup if
 * strict audio-thread allocation freedom is required.
 *
 * @return true if the parameters were loaded, false on NULL arguments, a
 * mismatched parameters_size, a reduction_curve_size mismatch, or an
 * allocation failure. See specbleach_denoiser_get_last_error().
 */
SPECBLEACH_API bool specbleach_denoiser_load_parameters(
    specbleach_denoiser* instance,
    const SpecbleachDenoiserParameters* parameters, uint32_t parameters_size);

/**
 * Processes a buffer of samples in place-capable buffers.
 *
 * Thread safety: RT-safe (no allocations, locks, or I/O). Safe for the
 * real-time audio callback thread.
 *
 * Buffer contract: input/output are plain float arrays of
 * number_of_samples length (mono, non-interleaved). Any block size is
 * accepted - internal buffering handles alignment. Output may alias input.
 * Returns false only for invalid arguments; internal processing cannot fail.
 *
 * @param instance Instance handle.
 * @param number_of_samples Block length in samples (must be > 0).
 * @param input Mono input buffer (must not be NULL).
 * @param output Mono output buffer (must not be NULL, may equal input).
 * @return true on success, false on NULL arguments or an empty block.
 */
SPECBLEACH_API bool specbleach_denoiser_process(specbleach_denoiser* instance,
                                                uint32_t number_of_samples,
                                                const float* input,
                                                float* output);

/**
 * Returns the algorithmic latency in samples introduced by the instance.
 * Hosts should report this value for delay compensation.
 *
 * Thread safety: RT-safe (read-only query). The value is stable for the
 * lifetime of the instance; query it once at prepare time. Smoothing mode
 * switches never change it.
 *
 * @param instance Instance handle.
 * @return Latency in samples, or 0 for a NULL instance.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_latency(const specbleach_denoiser* instance);

/**
 * Returns the sample rate passed at initialization.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Sample rate in Hz, or 0 for a NULL instance.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_sample_rate(const specbleach_denoiser* instance);

/**
 * Returns the effective STFT frame size in samples derived from
 * frame_size_ms at initialization (frame_ms / 1000 * sample_rate,
 * truncated). The hop is frame / 4 (4x overlap).
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Frame size in samples, or 0 for a NULL instance.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_frame_size(const specbleach_denoiser* instance);

/**
 * Returns the effective FFT size in samples: the smallest PFFFT-valid
 * size >= the frame size (see specbleach_denoiser_initialize() for the
 * validity rule). The noise profile size is fft_size / 2 + 1.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return FFT size in samples, or 0 for a NULL instance.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_fft_size(const specbleach_denoiser* instance);

/**
 * Returns the STFT hop size in samples (frame / 4).
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Hop size in samples, or 0 for a NULL instance.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_hop_size(const specbleach_denoiser* instance);

/**
 * Returns the size of the noise profile spectrum in bins.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Profile size in bins, or 0 for a NULL instance.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_noise_profile_size(const specbleach_denoiser* instance);

/**
 * Loads a custom noise profile for a specific mode.
 *
 * Thread safety: setup-only (copies profile data). Never call from an
 * audio thread.
 *
 * @param instance Instance handle.
 * @param restored_profile Source profile array (copied before returning).
 * @param profile_size Must equal get_noise_profile_size().
 * @param block_count Number of blocks the profile was estimated from.
 * @param mode One of [SPECBLEACH_PROFILE_MODE_FIRST,
 * SPECBLEACH_PROFILE_MODE_LAST].
 * @return false on invalid mode, size mismatch, or NULL arguments.
 */
SPECBLEACH_API bool specbleach_denoiser_load_noise_profile_for_mode(
    specbleach_denoiser* instance, const float* restored_profile,
    uint32_t profile_size, uint32_t block_count, SpecbleachProfileMode mode);

/**
 * Loads a noise profile, linearly resampling it to the instance's native
 * profile size when sizes differ. This absorbs the sample-rate/version
 * migration helper downstream apps otherwise reimplement.
 *
 * Thread safety: setup-only (copies profile data). Never call from an
 * audio thread.
 *
 * @param instance Instance handle.
 * @param restored_profile Source profile array of source_size bins.
 * @param source_size Source array length in bins (must be > 0).
 * @param block_count Number of blocks the profile was estimated from.
 * @param mode One of [SPECBLEACH_PROFILE_MODE_FIRST,
 * SPECBLEACH_PROFILE_MODE_LAST].
 * @return true on success, false on NULL arguments, invalid mode, or
 * allocation failure.
 */
SPECBLEACH_API bool specbleach_denoiser_load_noise_profile_resampled(
    specbleach_denoiser* instance, const float* restored_profile,
    uint32_t source_size, uint32_t block_count, SpecbleachProfileMode mode);

/**
 * Resets the internal noise profile of the instance (profile data only;
 * DSP history is preserved — see specbleach_denoiser_reset_dsp_state()
 * for a full flush, e.g. on transport seek).
 *
 * Thread safety: setup-only. Never call from an audio thread.
 *
 * @param instance Instance handle (NULL is a no-op).
 */
SPECBLEACH_API void specbleach_denoiser_reset_noise_profile(
    specbleach_denoiser* instance);

/**
 * Flushes all DSP history (STFT overlap buffers, NLM history, smoothers,
 * transient detector) and clears stored noise profiles. Output after this
 * call is as if the instance had just been created. Preserves loaded
 * parameters and the reduction-curve copy buffer (no allocation when the
 * curve is already sized).
 *
 * Thread safety: setup-only (tears down and rebuilds internal engines).
 * Never call from an audio thread.
 *
 * @param instance Instance handle.
 * @return true on success, false on NULL instance or rebuild failure.
 */
SPECBLEACH_API bool specbleach_denoiser_reset_dsp_state(
    specbleach_denoiser* instance);

/**
 * Returns the number of blocks used for the noise profile calculation for a
 * specific mode.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @param mode Profile mode index.
 * @return Block count, or 0 for NULL instance / out-of-range mode.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_noise_profile_block_count_for_mode(
    const specbleach_denoiser* instance, SpecbleachProfileMode mode);

/**
 * Returns a pointer to the internally owned noise profile for a specific
 * mode. The pointer stays valid until the instance is freed or the DSP
 * state is reset. Returns NULL for a NULL instance or an out-of-range
 * mode. Use specbleach_denoiser_noise_profile_available_for_mode() to check
 * whether learning has completed and the profile has been populated.
 *
 * Thread safety: RT-safe (read-only query). The returned pointer points to
 * live internal state and must not be written to or retained past free().
 *
 * @param instance Instance handle.
 * @param mode Profile mode index.
 * @return Read-only profile pointer, or NULL on error.
 */
SPECBLEACH_API const float* specbleach_denoiser_get_noise_profile_for_mode(
    const specbleach_denoiser* instance, SpecbleachProfileMode mode);

/**
 * Returns whether a noise profile has been calculated for a specific mode.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @param mode Profile mode index.
 * @return true if the profile for mode is populated.
 */
SPECBLEACH_API bool specbleach_denoiser_noise_profile_available_for_mode(
    const specbleach_denoiser* instance, SpecbleachProfileMode mode);

/**
 * Returns true if ANY profile mode is populated (convenience over looping
 * FIRST..LAST with noise_profile_available_for_mode()).
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return true if at least one mode has a profile.
 */
SPECBLEACH_API bool specbleach_denoiser_has_any_profile(
    const specbleach_denoiser* instance);

/**
 * Returns a pointer to the tonal mask array detected during spectral
 * processing. Array size matches
 * specbleach_denoiser_get_noise_profile_size(instance). Values range from
 * 0.0 (broadband) to 1.0 (pure tone).
 *
 * Thread safety: RT-safe (read-only query). Pointer points to live internal
 * state; valid until the instance is freed.
 *
 * @param instance Instance handle.
 * @return Read-only mask pointer, or NULL for a NULL instance.
 */
SPECBLEACH_API const float* specbleach_denoiser_get_tonal_mask(
    const specbleach_denoiser* instance);

/**
 * Returns the detected tonal peak frequencies in Hz.
 * Writes up to max_peaks peak frequencies into peak_freqs_hz.
 *
 * Thread safety: RT-safe (read-only query into the published double-buffer;
 * no allocation, no locks). Safe for the audio thread; prefer calling at a
 * decimated visualization rate.
 *
 * @param instance Instance handle.
 * @param peak_freqs_hz Output array (may be NULL only if max_peaks is 0).
 * @param max_peaks Capacity of peak_freqs_hz.
 * @return Number of peak frequencies written (0 = none found or NULL
 * instance).
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_tonal_peaks(const specbleach_denoiser* instance,
                                    float* peak_freqs_hz, uint32_t max_peaks);

/**
 * Returns a pointer to the active morphed noise profile array.
 *
 * Thread safety: RT-safe (read-only query). Pointer points to live internal
 * state; valid until the instance is freed.
 *
 * @param instance Instance handle.
 * @return Read-only profile pointer, or NULL for a NULL instance.
 */
SPECBLEACH_API const float* specbleach_denoiser_get_active_noise_profile(
    const specbleach_denoiser* instance);

/**
 * Returns true if a transient was detected in the last processed frame.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Transient flag for the last block.
 */
SPECBLEACH_API bool specbleach_denoiser_is_transient_detected(
    const specbleach_denoiser* instance);

/**
 * Returns the detected transient intensity [0.0, 1.0] from the last
 * processed frame.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle.
 * @return Intensity, or 0.0f for a NULL instance.
 */
SPECBLEACH_API float specbleach_denoiser_get_transient_intensity(
    const specbleach_denoiser* instance);

/**
 * Returns the failure reason for the last fallible call on this instance
 * (load/process/profile APIs). Cleared to SPECBLEACH_OK on every success.
 *
 * Thread safety: RT-safe (read-only query).
 *
 * @param instance Instance handle (NULL yields SPECBLEACH_ERR_NULL_ARG).
 * @return Last error code.
 */
SPECBLEACH_API SpecbleachError
specbleach_denoiser_get_last_error(const specbleach_denoiser* instance);

/**
 * Human-readable form of a SpecbleachError code. The pointer is static;
 * do not free it.
 *
 * Thread safety: thread-safe, callable without an instance.
 *
 * @param error Error code.
 * @return Short English description ("ok", "null argument", ...).
 */
SPECBLEACH_API const char* specbleach_denoiser_get_last_error_string(
    SpecbleachError error);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_DENOISER_H_INCLUDED */
