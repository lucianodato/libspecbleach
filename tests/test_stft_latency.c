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

#include "../src/shared/stft/stft_processor.h"
#include <assert.h>
#include <stdio.h>

int main(void) {
  const uint32_t sample_rate = 44100;
  const float stft_frame_size = 20.0f; // 20ms
  const uint32_t overlap_factor = 4;

  // Initialize STFT processor
  StftProcessor* processor =
      stft_processor_initialize(sample_rate, stft_frame_size, overlap_factor,
                                NO_PADDING, 0, HANN_WINDOW, HANN_WINDOW);

  assert(processor != NULL);

  // Get the latency
  uint32_t latency = get_stft_latency(processor);

  // Calculate expected values
  uint32_t expected_frame_size =
      (uint32_t)((stft_frame_size / 1000.0f) * (float)sample_rate);
  uint32_t expected_hop = expected_frame_size / overlap_factor;

  // The latency should be equal to the frame size (one full frame)
  // NOT frame_size - hop as it was before
  uint32_t expected_latency = expected_frame_size;

  printf("Frame size: %u samples\n", expected_frame_size);
  printf("Hop size: %u samples\n", expected_hop);
  printf("Expected latency: %u samples (%.2f ms)\n", expected_latency,
         (float)expected_latency * 1000.0f / (float)sample_rate);
  printf("Actual latency: %u samples (%.2f ms)\n", latency,
         (float)latency * 1000.0f / (float)sample_rate);

  // Test: latency should equal frame_size
  assert(latency == expected_latency);

  // Test: latency should NOT be frame_size - hop (the old incorrect value)
  uint32_t old_incorrect_latency = expected_frame_size - expected_hop;
  assert(latency != old_incorrect_latency);

  printf("\n✓ STFT latency calculation is correct!\n");
  printf("  Latency = frame_size (%u samples)\n", latency);
  printf("  (Previously was incorrectly: frame_size - hop = %u samples)\n",
         old_incorrect_latency);

  stft_processor_free(processor);

  return 0;
}
