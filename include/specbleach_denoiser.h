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

  /* Sets the amount of dBs that the noise will be attenuated. It goes from 0 dB
   * to 40 dB. This controls both the gain attenuation strength and the residual
   * noise mixing level for optimal noise reduction. */
  float reduction_amount;

  /* Percentage of smoothing to apply. Averages the reduction calculation frame
   * per frame so the rate of change is less resulting in less musical noise but
   * if too strong it can blur transient and reduce high frequencies. It goes
   * from 0 to 100 percent */
  float smoothing_factor;

  /* Percentage of whitening that is going to be applied to the residue of the
   * reduction. It modifies the noise floor to be more like white noise. This
   * can help hide musical noise when the noise is colored. It goes from 0 to
   * 100 percent */
  float whitening_factor;

  /* Enables the adaptive noise estimation, which will continuously update the
   * noise profile based on the current input signal. */
  int adaptive_noise;

  /* Sets the method used for adaptive noise estimation.
   * 0: SPP-MMSE method (unbiased estimation)
   * 1: Brandt (Trimmed Mean)
   * 2: Martin Minimum Statistics */
  int noise_estimation_method;

  /** Masking Veto depth (0.0 - 1.0) */
  float masking_depth; // 0.0 - 1.0: Depth of signal energy preservation

  /** Suppression aggressiveness (0.0 - 1.0) */
  float suppression_strength; // 0.0 - 1.0: Berouti oversubtraction factor

  /* Intelligent Steering */
  float aggressiveness; /**< -1.0 (Median/Min) to 1.0 (Max), 0.0 (Mean) */
  /* Tonal Separation */
  float tonal_reduction; // 0.0 to 1.0: Independent reduction for tones

  /* HPSS quality mode (0 = Off, 1 = Low, 2 = Mid, 3 = High). Default: 0 (Off)
   */
  int hpss_quality_mode;

  /* HPSS transient detection sensitivity (0.0 to 1.0). Default: 0.5 */
  float hpss_sensitivity;

  /* Noise Profile Offset — shifts the noise threshold up/down in dB.
   * Positive values make detection more aggressive (more noise removed).
   * Default: 2.0, Range: [-6.0, +6.0] */
  float noise_profile_offset_db;

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
 * Returns true if a transient was detected and protected by HPSS in the last
 * processed frame.
 */
bool specbleach_is_transient_detected(SpectralBleachHandle instance);

#ifdef __cplusplus
}
#endif
#endif
