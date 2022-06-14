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

typedef void *SpectralBleachHandle;

typedef struct SpectralBleachParameters {

  /* Sets the processor in listening mode to capture the noise profile. 0 is
   * disabled, 1 will learn the average profile, 2 will learn the maximun median
   * profile and 3 will learn the max profile. For the average and median
   * profile you need at least 5 frames of audio */
  int learn_noise;

  /* Enables outputting the residue of the reduction processing. It's either
   * true or false */
  bool residual_listen;

  /* Sets the amount of dBs that the noise will be attenuated. It goes from 0 dB
   * to 40 dB */
  float reduction_amount;

  /* Percentage of smoothing to apply. Averages the reduction calculation frame
   * per frame so the rate of change is less resulting in less musical noise but
   * if too strong it can blur transient and reduce high frequencies. It goes
   * from 0 to 100 percent */
  float smoothing_factor;

  /* Enables or disables the transient protection when smoothing_factor is being
   * used. This can help to preserve transient content when smoothing is strong.
   */
  bool transient_protection;

  /* Percentage of whitening that is going to be applied to the residue of the
   * reduction. It modifies the noise floor to be more like white noise. This
   * can help hide musical noise when the noise is colored. It goes from 0 to
   * 100 percent */
  float whitening_factor;

  /* Type of algorithm used to scale noise in order to apply over or under
   * subtraction in different parts of the spectrum while calculating the
   * reduction. 0 is a-posteriori snr scaling using the complete spectrum, 1 is
   * a-posteriori using critical bands and 3 is using masking thresholds
   */
  int noise_scaling_type;

  /* Strength in which the reduction will be applied. It uses the masking
   * thresholds of the signal to determine where in the spectrum the reduction
   * needs to be stronger. This parameter scales how much in each of the
   * frequencies the reduction is going to be applied. It can be a positive dB
   * value in between 0 dB and 12 dB */
  float noise_rescale;

  /* Sets the SNR threshold in dB in which the post-filter will start to blur
   * musical noise. It can be a positive or negative dB value in between -10 dB
   * and 10 dB */
  float post_filter_threshold;
} SpectralBleachParameters;

/**
 * Returns a handle to an instance of the library for the adaptive based
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
                                SpectralBleachParameters parameters);
/**
 * Process buffer of a number of samples
 */
bool specbleach_process(SpectralBleachHandle instance,
                        uint32_t number_of_samples, const float *input,
                        float *output);
/**
 * Returns the latency in samples associated with the library instance
 */
uint32_t specbleach_get_latency(SpectralBleachHandle instance);
/**
 * Returns the size of the noise profile spectrum
 */
uint32_t specbleach_get_noise_profile_size(SpectralBleachHandle instance);
/**
 * Returns a pointer to the noise profile calculated inside the instance
 */
float *specbleach_get_noise_profile(SpectralBleachHandle instance);
/**
 * Allows to load a custom noise profile
 */
bool specbleach_load_noise_profile(SpectralBleachHandle instance,
                                   const float *restored_profile,
                                   uint32_t profile_size,
                                   uint32_t profile_blocks);
/**
 * Resets the internal noise profile of the library instance
 */
bool specbleach_reset_noise_profile(SpectralBleachHandle instance);
/**
 * Returns if the instance has a noise profile calculated internally
 */
bool specbleach_noise_profile_available(SpectralBleachHandle instance);
/**
 * Returns the number of blocks used for the noise profile calculation
 */
uint32_t
specbleach_get_noise_profile_blocks_averaged(SpectralBleachHandle instance);

#ifdef __cplusplus
}
#endif
#endif