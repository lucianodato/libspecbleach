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

/**
 * @file simple_embed.c
 * @brief Zero-dependency embedding example: denoises a synthetic signal.
 *
 * Demonstrates the complete processing lifecycle with no external libraries
 * (no libsndfile), so it builds and runs anywhere the library builds:
 * generate noisy audio -> learn a noise profile -> denoise -> verify the
 * noise actually dropped -> clean up. Every call here maps 1:1 to a
 * real-time audio callback integration.
 *
 * Build: cmake -B build -DENABLE_EXAMPLES=ON && cmake --build build
 * Run:   ./build/simple_embed
 */

#include <math.h>
#include <specbleach_denoiser.h>
#include <stdint.h>
#include <stdio.h>

#define SAMPLE_RATE 48000U
#define FRAME_SIZE_MS 46.0f
#define BLOCK_SIZE 512U
#define LEARN_BLOCKS 8U
#define TOTAL_BLOCKS 64U

static float make_noisy_sample(uint32_t sample_rate, uint32_t index) {
  // 220 Hz tone plus white noise floor
  float tone = 0.25f * sinf(2.0f * (float)M_PI * 220.0f * (float)index /
                            (float)sample_rate);
  float noise =
      ((float)((index * 1103515245U + 12345U) >> 16) / (float)(1U << 15)) -
      1.0f;
  return tone + 0.05f * noise;
}

int main(void) {
  float input[BLOCK_SIZE];
  float output[BLOCK_SIZE];

  // 1. CREATE
  specbleach_denoiser* denoiser =
      specbleach_denoiser_initialize(SAMPLE_RATE, FRAME_SIZE_MS);
  if (!denoiser) {
    fprintf(stderr, "Error: failed to initialize denoiser\n");
    return 1;
  }

  // Report algorithmic latency to your host for delay compensation
  printf("Latency: %u samples\n", specbleach_denoiser_get_latency(denoiser));

  // 2. CONFIGURE (learning on): start from the safe baseline and
  // override only what this example needs
  SpecbleachDenoiserParameters parameters =
      specbleach_denoiser_get_default_parameters();
  parameters.learn_noise = SPECBLEACH_LEARN_ALL;
  parameters.reduction_gain = 0.1f; // -20 dB
  if (!specbleach_denoiser_load_parameters(denoiser, &parameters,
                                           sizeof(parameters))) {
    fprintf(stderr, "Error: failed to load parameters\n");
    specbleach_denoiser_free(denoiser);
    return 1;
  }

  // 3. LEARN — feed blocks containing ONLY noise
  for (uint32_t i = 0; i < LEARN_BLOCKS; i++) {
    for (uint32_t n = 0; n < BLOCK_SIZE; n++) {
      input[n] = make_noisy_sample(SAMPLE_RATE, i * BLOCK_SIZE + n);
    }
    if (!specbleach_denoiser_process(denoiser, BLOCK_SIZE, input, output)) {
      fprintf(stderr, "Error: learn stage failed\n");
      specbleach_denoiser_free(denoiser);
      return 1;
    }
  }

  // 4. FINALIZE — profiles are only usable after learning turns OFF
  parameters.learn_noise = SPECBLEACH_LEARN_OFF;
  if (!specbleach_denoiser_load_parameters(denoiser, &parameters,
                                           sizeof(parameters))) {
    fprintf(stderr, "Error: failed to finalize parameters\n");
    specbleach_denoiser_free(denoiser);
    return 1;
  }

  // 5. REDUCE — same call shape for the rest of the stream
  float in_rms = 0.0f;
  float out_rms = 0.0f;
  for (uint32_t i = LEARN_BLOCKS; i < TOTAL_BLOCKS; i++) {
    for (uint32_t n = 0; n < BLOCK_SIZE; n++) {
      input[n] = make_noisy_sample(SAMPLE_RATE, i * BLOCK_SIZE + n);
    }
    if (!specbleach_denoiser_process(denoiser, BLOCK_SIZE, input, output)) {
      fprintf(stderr, "Error: reduction stage failed\n");
      specbleach_denoiser_free(denoiser);
      return 1;
    }
    for (uint32_t n = 0; n < BLOCK_SIZE; n++) {
      in_rms += input[n] * input[n];
      out_rms += output[n] * output[n];
    }
  }
  in_rms = sqrtf(in_rms / (float)((TOTAL_BLOCKS - LEARN_BLOCKS) * BLOCK_SIZE));
  out_rms =
      sqrtf(out_rms / (float)((TOTAL_BLOCKS - LEARN_BLOCKS) * BLOCK_SIZE));
  printf("Input RMS: %f  Output RMS: %f\n", (double)in_rms, (double)out_rms);

  // 6. CLEANUP
  specbleach_denoiser_free(denoiser);

  if (out_rms >= in_rms) {
    fprintf(stderr, "FAIL: output level did not decrease (no reduction)\n");
    return 1;
  }
  printf("PASS: denoiser reduced signal level\n");
  return 0;
}
