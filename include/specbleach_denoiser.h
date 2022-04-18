/*
libspecbleach - A spectral processing library

Copyright 2021 Luciano Dato <lucianodato@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/
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
  bool learn_noise;
  bool residual_listen;
  bool transient_protection;
  float reduction_amount;
  float smoothing_factor;
  float whitening_factor;
  float noise_rescale;
} SpectralBleachParameters;

/**
 * Returns a handle to an instance of the library for the manual based
 * noise reduction
 */
SpectralBleachHandle specbleach_initialize(uint32_t sample_rate);
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