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

#ifndef FFT_TRANSFORM_H
#define FFT_TRANSFORM_H

#include <stdbool.h>
#include <stdint.h>

// C17 enum validation
enum ZeroPaddingType {
  NEXT_POWER_OF_TWO = 0,
  FIXED_AMOUNT = 1,
  // Smallest PFFFT-valid size >= frame, frame centered, the rest zero-padded.
  // This IS padding (up to the next valid size) — the name states the policy,
  // not the absence of it. Padding never affects latency (latency =
  // frame_size + NLM look-ahead, constant across smoothing modes).
  PAD_TO_VALID_SIZE = 2,
};

// Compile-time validation of enum values
_Static_assert(NEXT_POWER_OF_TWO == 0, "NEXT_POWER_OF_TWO must be 0");
_Static_assert(FIXED_AMOUNT == 1, "FIXED_AMOUNT must be 1");
_Static_assert(PAD_TO_VALID_SIZE == 2, "PAD_TO_VALID_SIZE must be 2");

typedef enum ZeroPaddingType ZeroPaddingType;

typedef struct FftTransform FftTransform;

FftTransform* fft_transform_initialize(uint32_t frame_size,
                                       ZeroPaddingType padding_type,
                                       uint32_t zeropadding_amount);
FftTransform* fft_transform_initialize_bins(uint32_t fft_size);
void fft_transform_free(FftTransform* self);
bool fft_load_input_samples(FftTransform* self, const float* input);
bool fft_get_output_samples(FftTransform* self, float* output);
/**
 * Accumulates FFT output samples into an accumulator buffer.
 * @param accumulator Target buffer containing at least frame_size float
 * elements.
 */
bool fft_accumulate_output_samples(FftTransform* self, float* accumulator);
uint32_t get_fft_size(FftTransform* self);
uint32_t get_fft_real_spectrum_size(FftTransform* self);
bool compute_forward_fft(FftTransform* self);
bool compute_backward_fft(FftTransform* self);
float* get_fft_input_buffer(FftTransform* self);
float* get_fft_output_buffer(FftTransform* self);

#endif
