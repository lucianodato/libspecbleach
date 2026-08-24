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

typedef void* SpectralBleachHandle;

typedef struct SpectralBleachDenoiserParameters {
  /* Sets the processor in listening mode to capture the noise profile.
   * 0 is disabled, 1 will learn all profile types simultaneously. */
  int learn_noise;

  /* Enables outputting the residue of the reduction processing. It's either
   * true or false */
  bool residual_listen;

  /* Linear gain floor for broadband noise reduction (0.0 to 1.0, where 1.0 = 0
   * dB / no reduction). */
  float reduction_gain;

  /* Normalized smoothing factor across frames (0.0 to 1.0). */
  float smoothing_factor;

  /* Normalized whitening factor for residue noise floor (0.0 to 1.0). */
  float whitening_factor;

  /* Enables the adaptive noise estimation, which will continuously update the
   * noise profile based on the current input signal. */
  int adaptive_noise;

  /* Sets the method used for adaptive noise estimation.
   * 0: SPP-MMSE method (unbiased estimation)
   * 1: Brandt (Trimmed Mean)
   * 2: Martin Minimum Statistics */
  int noise_estimation_method;

  /** Masking Veto depth (0.0 - 1.0: Depth of signal energy preservation) */
  float masking_depth;

  /** Suppression aggressiveness (0.0 - 1.0: Berouti oversubtraction factor) */
  float suppression_strength;

  /* Intelligent Steering */
  float aggressiveness; /**< -1.0 (Median/Min) to 1.0 (Max), 0.0 (Mean) */

  /* Linear gain floor for tonal noise components (0.0 to 1.0, where 1.0 = 0 dB
   * / no reduction). */
  float tonal_reduction_gain;

  /* Enable HPSS transient protection (0 = disabled, 1 = enabled). Default: 1 */
  int hpss_enable;

  /* Noise Profile Linear Scale — multiplier for the noise power spectrum.
   * Positive inputs are clamped to [0.01f, 100.0f]. Non-positive inputs use
   * the 1.0f default. Values > 1.0 shift threshold higher (more reduction). */
  float noise_profile_scale;

  /* Tonal Noise Profile Linear Scale — multiplier applied to the noise power
   * spectrum only at detected tonal components. Positive inputs are clamped to
   * [0.01f, 100.0f]. Non-positive inputs use the 1.0f default. Values > 1.0
   * shift threshold higher (more reduction) at tonal bins. */
  float tonal_noise_profile_scale;

  /* Frequency-dependent reduction bias curve.
   * Array of dB offsets per frequency bin, or NULL if disabled.
   * Positive values = more reduction at that frequency.
   * Array length must match specbleach_get_noise_profile_size(). */
  const float* reduction_curve_bias;

  /* Enables the reduction curve bias. When false, reduction_curve_bias
   * is ignored even if non-NULL. */
  bool reduction_curve_enabled;
} SpectralBleachDenoiserParameters;

/**
 * Returns a handle to an instance of the library for the spectral based
 * noise reduction. Sample rate could be anything from 4000hz to 192khz.
 * Recommended frame size range is between 20ms and 100ms
 */
SpectralBleachHandle specbleach_initialize(uint32_t sample_rate,
                                           float frame_size);
/**
 * Free instance associated to the handle passed
 */
void specbleach_free(SpectralBleachHandle instance);
/**
 * Loads the parameters for the reduction.
 * This has to be called before processing
 */
bool specbleach_load_parameters(SpectralBleachHandle instance,
                                SpectralBleachDenoiserParameters parameters);
/**
 * Process buffer of a number of samples
 */
bool specbleach_process(SpectralBleachHandle instance,
                        uint32_t number_of_samples, const float* input,
                        float* output);
/**
 * Returns the latency in samples associated with the library instance
 */
uint32_t specbleach_get_latency(SpectralBleachHandle instance);
/**
 * Returns the size of the noise profile spectrum
 */
uint32_t specbleach_get_noise_profile_size(SpectralBleachHandle instance);
/**
 * Allows to load a custom noise profile for a specific mode
 */
bool specbleach_load_noise_profile_for_mode(SpectralBleachHandle instance,
                                            const float* restored_profile,
                                            uint32_t profile_size,
                                            uint32_t block_count, int mode);
/**
 * Resets the internal noise profile of the library instance
 */
bool specbleach_reset_noise_profile(SpectralBleachHandle instance);

/**
 * Returns the number of blocks used for the noise profile calculation for a
 * specific mode
 */
uint32_t specbleach_get_noise_profile_block_count_for_mode(
    SpectralBleachHandle instance, int mode);

/**
 * Returns a pointer to the noise profile for a specific mode
 */
float* specbleach_get_noise_profile_for_mode(SpectralBleachHandle instance,
                                             int mode);

bool specbleach_noise_profile_available_for_mode(SpectralBleachHandle instance,
                                                 int mode);

/**
 * Returns a pointer to the tonal mask array detected during spectral
 * processing. Array size matches specbleach_get_noise_profile_size(instance).
 * Values range from 0.0 (broadband) to 1.0 (pure tone).
 */
const float* specbleach_get_tonal_mask(SpectralBleachHandle instance);

/**
 * Returns the detected tonal peak frequencies in Hz.
 * Writes up to max_peaks peak frequencies into peak_freqs_hz.
 * Note: This is an offline/query-only API and must NOT be called from the
 * real-time audio thread.
 * @return Number of peak frequencies written.
 */
uint32_t specbleach_get_tonal_peaks(SpectralBleachHandle instance,
                                    float* peak_freqs_hz, uint32_t max_peaks);

/**
 * Returns peak frequencies in Hz directly for a given noise profile array.
 * Note: This is an offline/query-only API and must NOT be called from the
 * real-time audio thread.
 */
uint32_t specbleach_get_tonal_peaks_for_profile(SpectralBleachHandle instance,
                                                const float* profile,
                                                uint32_t profile_size,
                                                float* peak_freqs_hz,
                                                uint32_t max_peaks);

/**
 * Returns a pointer to the active morphed noise profile array.
 */
const float* specbleach_get_active_noise_profile(SpectralBleachHandle instance);

/**
 * Returns true if a transient was detected in the last processed frame.
 */
bool specbleach_is_transient_detected(SpectralBleachHandle instance);

/**
 * Returns the detected transient intensity [0.0, 1.0] from the last processed
 * frame.
 */
float specbleach_get_transient_intensity(SpectralBleachHandle instance);

#ifdef __cplusplus
}
#endif
#endif
