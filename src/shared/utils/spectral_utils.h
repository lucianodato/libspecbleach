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

#ifndef SPECTRAL_UTILS_H
#define SPECTRAL_UTILS_H

#include <stdbool.h>
#include <stdint.h>

typedef enum WindowTypes {
  HANN_WINDOW = 0,
  HAMMING_WINDOW = 1,
  BLACKMAN_WINDOW = 2,
  VORBIS_WINDOW = 3
} WindowTypes;

bool get_fft_window(float *window, uint32_t fft_size, WindowTypes window_type);
bool initialize_spectrum_with_value(float *spectrum, uint32_t spectrum_size,
                                    float value);
bool direct_matrix_to_vector_spectral_convolution(const float *matrix_spectum,
                                                  const float *spectrum,
                                                  float *out_spectrum,
                                                  uint32_t spectrum_size);
float max_spectral_value(const float *spectrum, uint32_t real_spectrum_size);
float min_spectral_value(const float *spectrum, uint32_t real_spectrum_size);
bool min_spectrum(float *spectrum_one, const float *spectrum_two,
                  uint32_t spectrum_size);
bool max_spectrum(float *spectrum_one, const float *spectrum_two,
                  uint32_t spectrum_size);
float fft_bin_to_freq(uint32_t bin_index, uint32_t sample_rate,
                      uint32_t fft_size);
uint32_t freq_to_fft_bin(float freq, uint32_t sample_rate, uint32_t fft_size);
float spectral_flux(const float *spectrum, const float *previous_spectrum,
                    uint32_t spectrum_size);
bool get_rolling_mean_spectrum(float *averaged_spectrum,
                               const float *current_spectrum,
                               uint32_t number_of_blocks,
                               uint32_t spectrum_size);
bool get_rolling_median_spectrum(float *median_spectrum,
                               const float *current_spectrum_buffer,
                               uint32_t number_of_blocks,
                               uint32_t spectrum_size);

#endif