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

#ifndef SPECBLEACH_2D_DENOISER_H_INCLUDED
#define SPECBLEACH_2D_DENOISER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

typedef void* SpectralBleachHandle;

/**
 * Parameters for the 2D Non-Local Means denoiser.
 * Uses 2D pattern matching on spectrograms to suppress musical noise.
 */
typedef struct SpectralBleach2DDenoiserParameters {
  /**
   * Sets the processor in listening mode to capture the noise profile.
   * 0 is disabled, 1 will learn all profile types simultaneously.
   */
  int learn_noise;
  bool residual_listen;

  /**
   * Linear gain floor for broadband noise reduction (0.0 to 1.0, where 1.0 = 0
   * dB / no reduction).
   */
  float reduction_gain;

  /**
   * Normalized 2D smoothing factor (0.0 to 1.0).
   * Controls 2D time-frequency gain-mask artifact smoothing across time and
   * frequency. Also controls the NLM h parameter: non-zero values activate and
   * scale NLM patch similarity smoothing up to NLM_MAX_H_PARAMETER.
   */
  float smoothing_factor;

  /**
   * Normalized whitening factor for residue noise floor (0.0 to 1.0).
   */
  float whitening_factor;

  /**
   * Enabled the adaptive noise estimation.
   */
  int adaptive_noise;

  /**
   * Sets the method used for adaptive noise estimation.
   * 0: SPP-MMSE method
   * 1: Brandt (Trimmed Mean)
   * 2: Martin Minimum Statistics
   */
  int noise_estimation_method;

  /**
   * Sets the masking protection depth for the NLM filter (0.0 to 1.0).
   */
  float nlm_masking_protection;

  /**
   * Sets the suppression aggressiveness (0.0 to 1.0).
   */
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

  /* Frequency-dependent reduction bias curve.
   * Array of dB offsets per frequency bin, or NULL if disabled.
   * Positive values = more reduction at that frequency.
   * Array length must match specbleach_2d_get_noise_profile_size(). */
  const float* reduction_curve_bias;

  /* Enables the reduction curve bias. When false, reduction_curve_bias
   * is ignored even if non-NULL. */
  bool reduction_curve_enabled;
} SpectralBleach2DDenoiserParameters;

/**
 * Returns a handle to an instance of the library for the 2D NLM based
 * noise reduction. Sample rate could be anything from 4000hz to 192khz.
 * Recommended frame size range is between 20ms and 100ms.
 *
 * Note: This processor has additional latency due to NLM look-ahead.
 * Use specbleach_2d_get_latency() to query the total latency.
 */
SpectralBleachHandle specbleach_2d_initialize(uint32_t sample_rate,
                                              float frame_size);

/**
 * Free instance associated to the handle passed.
 */
void specbleach_2d_free(SpectralBleachHandle instance);

/**
 * Loads the parameters for the reduction.
 * This has to be called before processing.
 */
bool specbleach_2d_load_parameters(
    SpectralBleachHandle instance,
    SpectralBleach2DDenoiserParameters parameters);

/**
 * Process buffer of a number of samples.
 */
bool specbleach_2d_process(SpectralBleachHandle instance,
                           uint32_t number_of_samples, const float* input,
                           float* output);

/**
 * Returns the latency in samples associated with the library instance.
 * Includes both STFT latency and NLM look-ahead latency.
 */
uint32_t specbleach_2d_get_latency(SpectralBleachHandle instance);

/**
 * Returns the size of the noise profile spectrum.
 */
uint32_t specbleach_2d_get_noise_profile_size(SpectralBleachHandle instance);

/**
 * Allows to load a custom noise profile for a specific mode.
 */
bool specbleach_2d_load_noise_profile_for_mode(SpectralBleachHandle instance,
                                               const float* restored_profile,
                                               uint32_t profile_size,
                                               uint32_t block_count, int mode);

/**
 * Resets the internal noise profile of the library instance.
 */
bool specbleach_2d_reset_noise_profile(SpectralBleachHandle instance);

/**
 * Returns the number of blocks used for the noise profile calculation for a
 * specific mode.
 */
uint32_t specbleach_2d_get_noise_profile_block_count_for_mode(
    SpectralBleachHandle instance, int mode);

/**
 * Returns a pointer to the noise profile for a specific mode.
 */
float* specbleach_2d_get_noise_profile_for_mode(SpectralBleachHandle instance,
                                                int mode);

/**
 * Returns if the instance has a noise profile calculated for a specific mode.
 */
bool specbleach_2d_noise_profile_available_for_mode(
    SpectralBleachHandle instance, int mode);

/**
 * Returns a pointer to the tonal mask array detected during spectral
 * processing. Array size matches
 * specbleach_2d_get_noise_profile_size(instance). Values range from 0.0
 * (broadband) to 1.0 (pure tone).
 */
const float* specbleach_2d_get_tonal_mask(SpectralBleachHandle instance);

/**
 * Returns the detected tonal peak frequencies in Hz.
 * Writes up to max_peaks peak frequencies into peak_freqs_hz.
 * Note: This is an offline/query-only API and must NOT be called from the
 * real-time audio thread.
 * @return Number of peak frequencies written.
 */
uint32_t specbleach_2d_get_tonal_peaks(SpectralBleachHandle instance,
                                       float* peak_freqs_hz,
                                       uint32_t max_peaks);

/**
 * Returns peak frequencies in Hz directly for a given noise profile array.
 * Note: This is an offline/query-only API and must NOT be called from the
 * real-time audio thread.
 */
uint32_t specbleach_2d_get_tonal_peaks_for_profile(
    SpectralBleachHandle instance, const float* profile, uint32_t profile_size,
    float* peak_freqs_hz, uint32_t max_peaks);

/**
 * Returns a pointer to the active morphed noise profile array.
 */
const float* specbleach_2d_get_active_noise_profile(
    SpectralBleachHandle instance);

/**
 * Returns true if a transient was detected and protected by HPSS in the last
 * processed frame.
 */
bool specbleach_2d_is_transient_detected(SpectralBleachHandle instance);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_2D_DENOISER_H_INCLUDED */
