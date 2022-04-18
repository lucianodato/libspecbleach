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

#ifndef SPECBLEACH_ADENOISER_H_INCLUDED
#define SPECBLEACH_ADENOISER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

typedef void *SpectralBleachHandle;

typedef struct SpectralBleachParameters {
  bool residual_listen;
  float reduction_amount;
  float smoothing_factor;
  float noise_rescale;
} SpectralBleachParameters;

/**
 * Returns a handle to an instance of the library for the adaptive based
 * noise reduction
 */
SpectralBleachHandle specbleach_adaptive_initialize(uint32_t sample_rate);

/**
 * Free instance associated to the handle passed
 */
void specbleach_adaptive_free(SpectralBleachHandle instance);
/**
 * Loads the parameters for the reduction.
 * This has to be called before processing
 */
bool specbleach_adaptive_load_parameters(SpectralBleachHandle instance,
                                         SpectralBleachParameters parameters);
/**
 * Returns the latency in samples associated with the library instance
 */
uint32_t specbleach_adaptive_get_latency(SpectralBleachHandle instance);
/**
 * Process buffer of a number of samples
 */
bool specbleach_adaptive_process(SpectralBleachHandle instance,
                                 uint32_t number_of_samples, const float *input,
                                 float *output);

#ifdef __cplusplus
}
#endif
#endif