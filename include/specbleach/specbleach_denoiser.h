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

#include <specbleach/specbleach_parameters.h>
#include <stdbool.h>
#include <stdint.h>

typedef void* SpectralBleachHandle;

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
 * Returns a pointer to the noise profile calculated inside the instance
 */
float* specbleach_get_noise_profile(SpectralBleachHandle instance);
/**
 * Allows to load a custom noise profile
 */
bool specbleach_load_noise_profile(SpectralBleachHandle instance,
                                   const float* restored_profile,
                                   uint32_t profile_size,
                                   uint32_t profile_blocks);
/**
 * Allows to load a custom noise profile for a specific mode
 */
bool specbleach_load_noise_profile_for_mode(SpectralBleachHandle instance,
                                            const float* restored_profile,
                                            uint32_t profile_size,
                                            uint32_t profile_blocks, int mode);
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
uint32_t specbleach_get_noise_profile_blocks_averaged(
    SpectralBleachHandle instance);

/**
 * Returns the number of blocks used for the noise profile calculation for a
 * specific mode
 */
uint32_t specbleach_get_noise_profile_blocks_averaged_for_mode(
    SpectralBleachHandle instance, int mode);

/**
 * Returns a pointer to the noise profile for a specific mode
 */
float* specbleach_get_noise_profile_for_mode(SpectralBleachHandle instance,
                                             int mode);

/**
 * Returns if the instance has a noise profile calculated for a specific mode
 */
bool specbleach_noise_profile_available_for_mode(SpectralBleachHandle instance,
                                                 int mode);

#ifdef __cplusplus
}
#endif
#endif
