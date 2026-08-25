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

#ifndef SPECBLEACH_DENOISER_H_INCLUDED
#define SPECBLEACH_DENOISER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_common.h"
#include "specbleach_export.h"

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
 * Thread safety contract for the whole header: unless explicitly stated
 * otherwise, functions taking a handle may be called from any thread, but
 * calls on the SAME instance must never run concurrently with each other or
 * with specbleach_denoiser_process().
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
   */
  float smoothing_factor;

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
   * Masking veto depth (0.0 - 1.0): depth of signal energy preservation.
   * Values outside the range are clamped.
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
   * Tonal Noise Profile Linear Scale — multiplier applied to the noise power
   * spectrum only at detected tonal components. Positive inputs are clamped
   * to [0.01f, 100.0f]. Non-positive inputs use the 1.0f default. Values >
   * 1.0 shift threshold higher (more reduction) at tonal bins.
   */
  float tonal_noise_profile_scale;
} SpecbleachDenoiserParameters;

/**
 * Creates a single-channel spectral denoiser instance.
 *
 * Sample rate can be anything from 4000 Hz to 192 kHz. Recommended frame
 * size range is between 20 ms and 100 ms.
 *
 * @return A new instance or NULL on allocation failure. Free it with
 * specbleach_denoiser_free().
 */
SPECBLEACH_API specbleach_denoiser* specbleach_denoiser_initialize(
    uint32_t sample_rate, float frame_size);

/**
 * Frees an instance created by specbleach_denoiser_initialize().
 * Passing NULL is a no-op. The handle is invalid after this call.
 */
SPECBLEACH_API void specbleach_denoiser_free(specbleach_denoiser* instance);

/**
 * Loads parameters for the reduction.
 *
 * @param parameters Pointer to the parameter block to load. Must not be
 * NULL. The library copies all data it needs, including the reduction curve,
 * before returning.
 * @param parameters_size Must be exactly
 * sizeof(SpecbleachDenoiserParameters). This guards against ABI mismatches
 * between separately compiled binaries; any other value fails cleanly
 * instead of reading out-of-bounds memory.
 *
 * Not guaranteed to be allocation-free on the FIRST call after enabling
 * reduction_curve_enabled (the internal copy buffer is allocated then);
 * subsequent calls reuse that buffer. Load parameters once during setup if
 * strict audio-thread allocation freedom is required.
 *
 * @return true if the parameters were loaded, false on NULL arguments or a
 * mismatched parameters_size.
 */
SPECBLEACH_API bool specbleach_denoiser_load_parameters(
    specbleach_denoiser* instance,
    const SpecbleachDenoiserParameters* parameters, uint32_t parameters_size);

/**
 * Processes a buffer of samples in place-capable buffers.
 *
 * Real-time safe: contains no allocations, locks, or I/O.
 *
 * @return true on success, false on NULL arguments or an empty block.
 */
SPECBLEACH_API bool specbleach_denoiser_process(specbleach_denoiser* instance,
                                                uint32_t number_of_samples,
                                                const float* input,
                                                float* output);

/**
 * Returns the algorithmic latency in samples introduced by the instance.
 * Hosts should report this value for delay compensation.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_latency(specbleach_denoiser* instance);

/**
 * Returns the size of the noise profile spectrum in bins.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_noise_profile_size(specbleach_denoiser* instance);

/**
 * Loads a custom noise profile for a specific mode.
 *
 * @param mode One of [SPECBLEACH_PROFILE_MODE_FIRST,
 * SPECBLEACH_PROFILE_MODE_LAST].
 * @return false on invalid mode, size mismatch, or NULL arguments.
 */
SPECBLEACH_API bool specbleach_denoiser_load_noise_profile_for_mode(
    specbleach_denoiser* instance, const float* restored_profile,
    uint32_t profile_size, uint32_t block_count, int mode);

/**
 * Resets the internal noise profile of the instance.
 */
SPECBLEACH_API bool specbleach_denoiser_reset_noise_profile(
    specbleach_denoiser* instance);

/**
 * Returns the number of blocks used for the noise profile calculation for a
 * specific mode.
 */
SPECBLEACH_API uint32_t
specbleach_denoiser_get_noise_profile_block_count_for_mode(
    specbleach_denoiser* instance, int mode);

/**
 * Returns a pointer to the internally owned noise profile for a specific
 * mode. The pointer stays valid until the instance is freed. May return NULL
 * if the profile was never populated.
 */
SPECBLEACH_API float* specbleach_denoiser_get_noise_profile_for_mode(
    specbleach_denoiser* instance, int mode);

/**
 * Returns whether a noise profile has been calculated for a specific mode.
 */
SPECBLEACH_API bool specbleach_denoiser_noise_profile_available_for_mode(
    specbleach_denoiser* instance, int mode);

/**
 * Returns a pointer to the tonal mask array detected during spectral
 * processing. Array size matches
 * specbleach_denoiser_get_noise_profile_size(instance). Values range from
 * 0.0 (broadband) to 1.0 (pure tone).
 */
SPECBLEACH_API const float* specbleach_denoiser_get_tonal_mask(
    specbleach_denoiser* instance);

/**
 * Returns the detected tonal peak frequencies in Hz.
 * Writes up to max_peaks peak frequencies into peak_freqs_hz.
 *
 * Offline/query-only API; must NOT be called from the real-time audio
 * thread.
 *
 * @return Number of peak frequencies written.
 */
SPECBLEACH_API uint32_t specbleach_denoiser_get_tonal_peaks(
    specbleach_denoiser* instance, float* peak_freqs_hz, uint32_t max_peaks);

/**
 * Returns peak frequencies in Hz computed directly from a caller-provided
 * noise profile array.
 *
 * Offline/query-only API; must NOT be called from the real-time audio
 * thread.
 */
SPECBLEACH_API uint32_t specbleach_denoiser_get_tonal_peaks_for_profile(
    specbleach_denoiser* instance, const float* profile, uint32_t profile_size,
    float* peak_freqs_hz, uint32_t max_peaks);

/**
 * Returns a pointer to the active morphed noise profile array.
 */
SPECBLEACH_API const float* specbleach_denoiser_get_active_noise_profile(
    specbleach_denoiser* instance);

/**
 * Returns true if a transient was detected in the last processed frame.
 */
SPECBLEACH_API bool specbleach_denoiser_is_transient_detected(
    specbleach_denoiser* instance);

/**
 * Returns the detected transient intensity [0.0, 1.0] from the last
 * processed frame.
 */
SPECBLEACH_API float specbleach_denoiser_get_transient_intensity(
    specbleach_denoiser* instance);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_DENOISER_H_INCLUDED */
