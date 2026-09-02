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

#include "specbleach_time_denoiser.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SAMPLE_RATE 48000U

static void test_public_contract(void) {
  assert(specbleach_time_denoiser_initialize(0U) == NULL);
  assert(specbleach_time_denoiser_get_latency(NULL) == 0U);
  assert(!specbleach_time_denoiser_load_parameters(NULL, NULL, 0U));
  assert(!specbleach_time_denoiser_process(NULL, 1U, NULL, NULL));

  specbleach_time_denoiser_free(NULL);

  specbleach_time_denoiser* instance =
      specbleach_time_denoiser_initialize(SAMPLE_RATE);
  assert(instance != NULL);
  assert(specbleach_time_denoiser_get_latency(instance) == 0U);

  SpecbleachTimeDenoiserParameters params;
  params.reduction_gain = 0.063F;
  params.smoothing_factor = 0.5F;
  params.adaptive_noise = true;
  params.noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE;
  params.suppression_strength = 0.5F;

  // ABI guard: wrong size must fail cleanly
  assert(!specbleach_time_denoiser_load_parameters(instance, &params, 0U));
  assert(specbleach_time_denoiser_load_parameters(
      instance, &params, sizeof(SpecbleachTimeDenoiserParameters)));

  specbleach_time_denoiser_free(instance);
}

static void test_impulse_zero_latency_all_block_sizes(void) {
  // An impulse in a block must emerge at the same position with unity
  // gains (delta filter before convergence) regardless of block size
  const uint32_t block_sizes[] = {1U, 13U, 64U, 256U, 512U};

  for (uint32_t s = 0U; s < sizeof(block_sizes) / sizeof(uint32_t); s++) {
    const uint32_t block = block_sizes[s];
    specbleach_time_denoiser* instance =
        specbleach_time_denoiser_initialize(SAMPLE_RATE);
    assert(instance != NULL);

    SpecbleachTimeDenoiserParameters params;
    params.reduction_gain = 1.0F;   // 0 dB floor: unity gain everywhere
    params.smoothing_factor = 0.0F; // no temporal lag: filter converges now
    params.adaptive_noise = false;
    params.noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE;
    params.suppression_strength = 0.0F;
    assert(specbleach_time_denoiser_load_parameters(
        instance, &params, sizeof(SpecbleachTimeDenoiserParameters)));

    float* input = (float*)calloc(block, sizeof(float));
    float* output = (float*)calloc(block, sizeof(float));

    // Warm up with enough silence for multiple 512-sample analysis frames so
    // the worker publishes an initial filter
    for (uint32_t b = 0U; b < 8U * (512U / block + 1U); b++) {
      assert(specbleach_time_denoiser_process(instance, block, input, output));
    }

    // Single impulse at a mid-block position must appear at that position
    const uint32_t impulse_pos = block / 2U;
    for (uint32_t b = 0U; b < 4U; b++) {
      memset(input, 0, block * sizeof(float));
      memset(output, 0, block * sizeof(float));
      if (b == 0U) {
        input[impulse_pos] = 1.0F;
      }
      assert(specbleach_time_denoiser_process(instance, block, input, output));
      for (uint32_t i = 0U; i < block; i++) {
        assert(isfinite(output[i]));
        if (b == 0U && i == impulse_pos) {
          assert(output[i] > 0.9F); // still ~delta near n=0 (min-phase)
        } else if (b == 0U) {
          assert(fabsf(output[i]) < 0.1F); // minimal smearing
        }
      }
    }

    free(input);
    free(output);
    specbleach_time_denoiser_free(instance);
  }
}

static void test_sweep_no_artifacts(void) {
  specbleach_time_denoiser* instance =
      specbleach_time_denoiser_initialize(SAMPLE_RATE);
  assert(instance != NULL);

  SpecbleachTimeDenoiserParameters params;
  params.reduction_gain = 0.5F;
  params.smoothing_factor = 0.5F;
  params.adaptive_noise = true;
  params.noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE;
  params.suppression_strength = 0.5F;
  assert(specbleach_time_denoiser_load_parameters(
      instance, &params, sizeof(SpecbleachTimeDenoiserParameters)));

  const uint32_t block = 64U;
  const uint32_t total = 48000U; // 1 second
  float input[block];
  float output[block];

  // Logarithmic sine sweep 20 Hz -> 20 kHz; assert bounded finite output
  for (uint32_t n = 0U; n < total; n += block) {
    for (uint32_t i = 0U; i < block; i++) {
      const float t = (float)(n + i) / (float)SAMPLE_RATE;
      const float freq = 20.0F * powf(1000.0F, t); // 20 Hz -> 20 kHz over 1 s
      input[i] = 0.5F * sinf(2.0F * M_PI * freq * t);
    }
    assert(specbleach_time_denoiser_process(instance, block, input, output));
    for (uint32_t i = 0U; i < block; i++) {
      assert(isfinite(output[i]));
      assert(fabsf(output[i]) <= 2.0F);
    }
  }

  specbleach_time_denoiser_free(instance);
}

static float lcg_rand(void) {
  static uint32_t state = 12345U;
  state = state * 1664525U + 1013904223U;
  return (float)((double)state / 4294967295.0) * 2.0F - 1.0F;
}

static void test_adaptive_noise_reduction_all_block_sizes(void) {
  // Stationary noise + adaptive estimation must yield real reduction for
  // ANY host block size, including sizes larger than the internal ring
  // capacity (regression: whole-block ring push used to drop everything)
  const uint32_t block_sizes[] = {64U, 512U, 1024U, 4096U};

  for (uint32_t s = 0U; s < sizeof(block_sizes) / sizeof(uint32_t); s++) {
    const uint32_t block = block_sizes[s];
    specbleach_time_denoiser* instance =
        specbleach_time_denoiser_initialize(SAMPLE_RATE);
    assert(instance != NULL);

    SpecbleachTimeDenoiserParameters params;
    params.reduction_gain = 0.01F; // -40 dB gain floor
    params.smoothing_factor = 0.0F;
    params.adaptive_noise = true;
    params.noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE;
    params.suppression_strength = 1.0F;
    assert(specbleach_time_denoiser_load_parameters(
        instance, &params, sizeof(SpecbleachTimeDenoiserParameters)));

    float* input = (float*)malloc(block * sizeof(float));
    float* output = (float*)malloc(block * sizeof(float));

    double in_energy = 0.0;
    double out_energy = 0.0;
    const uint32_t warmup_blocks = (SAMPLE_RATE * 2U) / block;
    const uint32_t measure_blocks = (SAMPLE_RATE * 2U) / block;

    for (uint32_t b = 0U; b < warmup_blocks + measure_blocks; b++) {
      for (uint32_t i = 0U; i < block; i++) {
        input[i] = 0.1F * lcg_rand();
      }
      assert(specbleach_time_denoiser_process(instance, block, input, output));
      if (b >= warmup_blocks) {
        for (uint32_t i = 0U; i < block; i++) {
          in_energy += (double)input[i] * (double)input[i];
          out_energy += (double)output[i] * (double)output[i];
        }
      }
    }

    const double in_rms = sqrt(in_energy / (double)(measure_blocks * block));
    const double out_rms = sqrt(out_energy / (double)(measure_blocks * block));
    assert(isfinite(in_rms) && in_rms > 0.0);
    assert(isfinite(out_rms));
    assert(out_rms < in_rms * 0.25); // expect close to the -40 dB floor

    free(input);
    free(output);
    specbleach_time_denoiser_free(instance);
  }
}

int main(void) {
  test_public_contract();
  test_impulse_zero_latency_all_block_sizes();
  test_sweep_no_artifacts();
  test_adaptive_noise_reduction_all_block_sizes();

  printf("✓ Time denoiser public API tests passed\n");
  return 0;
}
