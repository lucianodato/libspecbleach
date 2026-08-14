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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "specbleach_2d_denoiser.h"

#define ROLLING_MEAN 1

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define SAMPLE_RATE 44100
#define FRAME_SIZE 20.0f // ms

void test_initialization_failure(void) {
  printf("Testing initialization failure cases...\n");

  // Invalid parameters
  TEST_ASSERT(specbleach_2d_initialize(0, FRAME_SIZE) == NULL,
              "Should fail with 0 sample rate");
  TEST_ASSERT(specbleach_2d_initialize(SAMPLE_RATE, 0.0f) == NULL,
              "Should fail with 0 frame size");
}

void test_null_handling(void) {
  printf("Testing NULL handling...\n");

  specbleach_2d_free(NULL); // Should not crash

  TEST_ASSERT(specbleach_2d_get_latency(NULL) == 0,
              "NULL handle should return 0 latency");

  TEST_ASSERT(specbleach_2d_process(NULL, 100, NULL, NULL) == false,
              "NULL handle/buffers process should fail");

  float buffer[100];
  SpectralBleachHandle h = specbleach_2d_initialize(SAMPLE_RATE, FRAME_SIZE);
  TEST_ASSERT(specbleach_2d_process(h, 0, buffer, buffer) == false,
              "Zero samples process should fail");
  TEST_ASSERT(specbleach_2d_process(h, 100, NULL, buffer) == false,
              "NULL input buffer should fail");
  TEST_ASSERT(specbleach_2d_process(h, 100, buffer, NULL) == false,
              "NULL output buffer should fail");

  TEST_ASSERT(specbleach_2d_get_noise_profile_size(NULL) == 0,
              "NULL handle should return 0 profile size");
  TEST_ASSERT(specbleach_2d_get_noise_profile_block_count_for_mode(
                  NULL, ROLLING_MEAN) == 0,
              "NULL handle should return 0 blocks averaged");
  TEST_ASSERT(
      specbleach_2d_get_noise_profile_for_mode(NULL, ROLLING_MEAN) == NULL,
      "NULL handle should return NULL profile");

  TEST_ASSERT(specbleach_2d_load_noise_profile_for_mode(NULL, NULL, 0, 0,
                                                        ROLLING_MEAN) == false,
              "NULL handle loading profile should fail");
  TEST_ASSERT(specbleach_2d_load_noise_profile_for_mode(h, NULL, 0, 0,
                                                        ROLLING_MEAN) == false,
              "NULL profile data loading should fail");

  TEST_ASSERT(specbleach_2d_reset_noise_profile(NULL) == false,
              "NULL handle reset should fail");
  TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(
                  NULL, ROLLING_MEAN) == false,
              "NULL handle available check should fail");

  specbleach_2d_free(h);
}

void test_noise_profile_api(void) {
  printf("Testing Noise Profile API...\n");

  SpectralBleachHandle h = specbleach_2d_initialize(SAMPLE_RATE, FRAME_SIZE);

  SpectralBleach2DDenoiserParameters params = {
      .learn_noise = 0,
      .reduction_amount = 20.0f,
      .smoothing_factor = 1.0f,
      .nlm_masking_protection = 0.5f,
      .tonal_reduction = 0.0f,
      .aggressiveness = 0.0f,
  };
  specbleach_2d_load_parameters(h, params);

  uint32_t profile_size = specbleach_2d_get_noise_profile_size(h);
  TEST_ASSERT(profile_size > 0, "Profile size should be valid");

  // Test loading a fake profile
  float* fake_profile = calloc(profile_size, sizeof(float));
  for (uint32_t i = 0; i < profile_size; i++) {
    fake_profile[i] = 1.0f;
  }

  TEST_ASSERT(specbleach_2d_load_noise_profile_for_mode(
                  h, fake_profile, profile_size, 10, ROLLING_MEAN),
              "Load profile success");
  TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(h, ROLLING_MEAN),
              "Profile available true");
  TEST_ASSERT(specbleach_2d_get_noise_profile_block_count_for_mode(
                  h, ROLLING_MEAN) == 10,
              "Blocks averaged 10");

  // Test fetching profile
  float* retrieved = specbleach_2d_get_noise_profile_for_mode(h, ROLLING_MEAN);
  TEST_ASSERT(retrieved != NULL, "Should be able to retrieve profile");
  TEST_ASSERT(retrieved[0] == 1.0f, "Retrieved data should match");

  // Test reset
  specbleach_2d_reset_noise_profile(h);
  TEST_ASSERT(!specbleach_2d_noise_profile_available_for_mode(h, ROLLING_MEAN),
              "Profile should be unavailable after reset");

  // Test per-mode API
  TEST_ASSERT(specbleach_2d_load_noise_profile_for_mode(h, fake_profile,
                                                        profile_size, 5, 1),
              "Loading for mode 1 should succeed");

  TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(h, 1),
              "Mode 1 should be available");
  TEST_ASSERT(!specbleach_2d_noise_profile_available_for_mode(h, 2),
              "Mode 2 should be unavailable");

  TEST_ASSERT(specbleach_2d_get_noise_profile_block_count_for_mode(h, 1) == 5,
              "Mode 1 blocks should match");

  TEST_ASSERT(specbleach_2d_get_noise_profile_for_mode(h, 1) != NULL,
              "Mode 1 profile should be retrievable");
  TEST_ASSERT(specbleach_2d_get_noise_profile_for_mode(NULL, 1) == NULL,
              "NULL handle get mode profile should return NULL");

  // Test invalid load params
  TEST_ASSERT(specbleach_2d_load_noise_profile_for_mode(
                  h, fake_profile, profile_size + 1, 10, ROLLING_MEAN) == false,
              "Size mismatch load should fail");
  TEST_ASSERT(specbleach_2d_load_noise_profile_for_mode(
                  h, fake_profile, profile_size, 10, 99) == false,
              "Invalid mode load should fail");

  // Test invalid mode validation (modes outside 1-4)
  TEST_ASSERT(specbleach_2d_get_noise_profile_block_count_for_mode(h, 0) == 0,
              "Invalid mode 0 should return 0 blocks");
  TEST_ASSERT(specbleach_2d_get_noise_profile_block_count_for_mode(h, 5) == 0,
              "Invalid mode 5 should return 0 blocks");
  TEST_ASSERT(specbleach_2d_get_noise_profile_for_mode(h, 0) == NULL,
              "Invalid mode 0 should return NULL profile");
  TEST_ASSERT(specbleach_2d_get_noise_profile_for_mode(h, 5) == NULL,
              "Invalid mode 5 should return NULL profile");
  TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(h, 0) == false,
              "Invalid mode 0 should return false availability");
  TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(h, 5) == false,
              "Invalid mode 5 should return false availability");

  free(fake_profile);
  specbleach_2d_free(h);
}

void test_2d_parameter_switching(void) {
  printf("Testing 2D parameter switching and adaptive methods...\n");
  SpectralBleachHandle h = specbleach_2d_initialize(SAMPLE_RATE, FRAME_SIZE);

  SpectralBleach2DDenoiserParameters params = {
      .learn_noise = 0,
      .reduction_amount = 20.0f,
      .smoothing_factor = 1.0f,
      .adaptive_noise = 1,
      .noise_estimation_method = 0, // Louizou
      .nlm_masking_protection = 0.5f,
      .tonal_reduction = 0.0f,
      .aggressiveness = 0.0f,
  };

  // 1. Load Louizou adaptive
  TEST_ASSERT(specbleach_2d_load_parameters(h, params),
              "Load Louizou adaptive should succeed");

  // 2. Switch to SPP-MMSE adaptive
  params.noise_estimation_method = 1;
  TEST_ASSERT(specbleach_2d_load_parameters(h, params),
              "Switch to SPP-MMSE should succeed");

  // 3. Switch adaptive off
  params.adaptive_noise = 0;
  TEST_ASSERT(specbleach_2d_load_parameters(h, params),
              "Switch adaptive off should succeed");

  // 4. Test different reduction modes
  for (int mode = 2; mode <= 3; mode++) {
    // No longer noise_reduction_mode
    TEST_ASSERT(specbleach_2d_load_parameters(h, params),
                "Switch reduction mode should succeed");
  }

  uint32_t profile_size = specbleach_2d_get_noise_profile_size(h);
  float* input = calloc(1024, sizeof(float));
  float* output = calloc(1024, sizeof(float));
  float* profile = calloc(profile_size, sizeof(float));
  for (uint32_t i = 0; i < profile_size; i++) {
    profile[i] = 0.001f;
  }

  specbleach_2d_load_noise_profile_for_mode(h, profile, profile_size, 1,
                                            ROLLING_MEAN);

  // Run with SPP-MMSE adaptive (method 1)
  params.adaptive_noise = 1;
  params.noise_estimation_method = 1;
  specbleach_2d_load_parameters(h, params);
  specbleach_2d_process(h, 1024, input, output);

  // Run with Louizou adaptive (method 0)
  params.noise_estimation_method = 0;
  specbleach_2d_load_parameters(h, params);
  specbleach_2d_process(h, 1024, input, output);

  free(input);
  free(output);
  free(profile);
  specbleach_2d_free(h);
}

void test_process_loop(void) {
  printf("Testing process loop (happy/unhappy paths)...\n");
  SpectralBleachHandle h = specbleach_2d_initialize(SAMPLE_RATE, FRAME_SIZE);

  float input[1024] = {0};
  float output[1024] = {0};

  // 1. Process without noise profile (unhappy path for reduction)
  SpectralBleach2DDenoiserParameters params = {
      .learn_noise = 0,
      .reduction_amount = 20.0f,
      .smoothing_factor = 1.0f,
      .tonal_reduction = 0.0f,
      .aggressiveness = 0.0f,
  };
  specbleach_2d_load_parameters(h, params);

  // Should succeed (likely silence or passthrough depending on impl)
  TEST_ASSERT(specbleach_2d_process(h, 1024, input, output),
              "Process without profile should succeed safely");

  // 2. Process with learn mode
  params.learn_noise = 1;
  specbleach_2d_load_parameters(h, params);

  // Process enough samples to fill latency buffers and accumulate profile stats
  // We need at least ~5 frames + latency.
  // 1 second of audio (44100 samples) should be plenty.
  float* big_buffer = calloc(44100, sizeof(float));
  TEST_ASSERT(specbleach_2d_process(h, 44100, big_buffer, big_buffer),
              "Process in learn mode should succeed");
  free(big_buffer);

  TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(h, ROLLING_MEAN),
              "Profile should be available after learning");

  // New PR #99 2D getters
  specbleach_2d_get_tonal_mask(h);
  specbleach_2d_get_active_noise_profile(h);
  float peak_freqs[10];
  specbleach_2d_get_tonal_peaks(h, peak_freqs, 10);
  uint32_t prof_size = specbleach_2d_get_noise_profile_size(h);
  float* dummy_prof = calloc(prof_size, sizeof(float));
  specbleach_2d_get_tonal_peaks_for_profile(h, dummy_prof, prof_size,
                                            peak_freqs, 10);
  free(dummy_prof);

  // Process with tonal reduction and HPSS
  float in_buf[1024] = {0};
  float out_buf[1024] = {0};
  params.learn_noise = 0;
  params.tonal_reduction = 0.5f;
  params.reduction_amount = 20.0f;
  params.hpss_quality_mode = 2;
  specbleach_2d_load_parameters(h, params);
  specbleach_2d_process(h, 1024, in_buf, out_buf);
  specbleach_2d_get_tonal_mask(h);
  specbleach_2d_get_tonal_peaks(h, peak_freqs, 10);

  // Switch HPSS modes
  params.hpss_quality_mode = 1;
  specbleach_2d_load_parameters(h, params);
  specbleach_2d_process(h, 1024, in_buf, out_buf);

  params.hpss_quality_mode = 3;
  specbleach_2d_load_parameters(h, params);
  specbleach_2d_process(h, 1024, in_buf, out_buf);

  params.hpss_quality_mode = 0;
  specbleach_2d_load_parameters(h, params);
  specbleach_2d_process(h, 1024, in_buf, out_buf);

  // NULL checks
  TEST_ASSERT(specbleach_2d_get_tonal_mask(NULL) == NULL, "NULL 2D tonal mask");
  TEST_ASSERT(specbleach_2d_get_active_noise_profile(NULL) == NULL,
              "NULL 2D active profile");
  TEST_ASSERT(specbleach_2d_get_tonal_peaks(NULL, peak_freqs, 10) == 0,
              "NULL 2D tonal peaks");
  TEST_ASSERT(specbleach_2d_get_tonal_peaks_for_profile(NULL, NULL, 0,
                                                        peak_freqs, 10) == 0,
              "NULL 2D tonal peaks for profile");

  specbleach_2d_free(h);
}

void test_2d_smoothing_factor_responsiveness(void) {
  printf("Testing 2D smoothing factor responsiveness (0%% vs 100%%)...\n");

  SpectralBleachHandle h0 = specbleach_2d_initialize(SAMPLE_RATE, FRAME_SIZE);
  TEST_ASSERT(h0 != NULL, "Initialization of h0 should succeed");

  SpectralBleachHandle h100 = specbleach_2d_initialize(SAMPLE_RATE, FRAME_SIZE);
  TEST_ASSERT(h100 != NULL, "Initialization of h100 should succeed");

  // Create noisy input signal
  const int n_samples = 44100;
  float* noise_input = calloc(n_samples, sizeof(float));
  float* temp_out = calloc(n_samples, sizeof(float));
  float* out_smooth_0 = calloc(n_samples, sizeof(float));
  float* out_smooth_100 = calloc(n_samples, sizeof(float));

  for (int i = 0; i < n_samples; i++) {
    noise_input[i] = 0.1f * ((float)(i % 100) / 100.0f - 0.5f);
  }

  // 1. Learn noise profile for both handles
  SpectralBleach2DDenoiserParameters params_learn = {
      .learn_noise = 1,
      .reduction_amount = 20.0f,
      .smoothing_factor = 0.0f,
  };

  TEST_ASSERT(specbleach_2d_load_parameters(h0, params_learn),
              "Loading learn params for h0 should succeed");
  TEST_ASSERT(specbleach_2d_process(h0, n_samples, noise_input, temp_out),
              "Processing learn mode for h0 should succeed");
  TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(h0, ROLLING_MEAN),
              "Profile for h0 should be available after learning");

  TEST_ASSERT(specbleach_2d_load_parameters(h100, params_learn),
              "Loading learn params for h100 should succeed");
  TEST_ASSERT(specbleach_2d_process(h100, n_samples, noise_input, temp_out),
              "Processing learn mode for h100 should succeed");
  TEST_ASSERT(
      specbleach_2d_noise_profile_available_for_mode(h100, ROLLING_MEAN),
      "Profile for h100 should be available after learning");

  // 2. Process with 0% smoothing on h0
  SpectralBleach2DDenoiserParameters params_0 = {
      .learn_noise = 0,
      .reduction_amount = 20.0f,
      .smoothing_factor = 0.0f,
  };
  TEST_ASSERT(specbleach_2d_load_parameters(h0, params_0),
              "Loading 0% smoothing params for h0 should succeed");
  TEST_ASSERT(specbleach_2d_process(h0, n_samples, noise_input, out_smooth_0),
              "Processing 0% smoothing for h0 should succeed");

  // 3. Process with 100% smoothing on h100
  SpectralBleach2DDenoiserParameters params_100 = {
      .learn_noise = 0,
      .reduction_amount = 20.0f,
      .smoothing_factor = 100.0f,
  };
  TEST_ASSERT(specbleach_2d_load_parameters(h100, params_100),
              "Loading 100% smoothing params for h100 should succeed");
  TEST_ASSERT(
      specbleach_2d_process(h100, n_samples, noise_input, out_smooth_100),
      "Processing 100% smoothing for h100 should succeed");

  // Compare output buffers: 0% vs 100% smoothing should produce different
  // results
  float diff_sum = 0.0f;
  for (int i = 2000; i < n_samples; i++) {
    diff_sum += fabsf(out_smooth_0[i] - out_smooth_100[i]);
  }
  TEST_ASSERT(diff_sum > 1e-4f,
              "0% and 100% smoothing factor should produce different output");

  free(noise_input);
  free(temp_out);
  free(out_smooth_0);
  free(out_smooth_100);
  specbleach_2d_free(h0);
  specbleach_2d_free(h100);
  printf("✓ 2D smoothing factor responsiveness tests passed\n");
}

int main(void) {
  printf("Running specbleach_2d_denoiser wrapper tests...\n");

  test_initialization_failure();
  test_null_handling();
  test_noise_profile_api();
  test_2d_parameter_switching();
  test_process_loop();
  test_2d_smoothing_factor_responsiveness();

  printf("✅ All specbleach_2d_denoiser tests passed!\n");
  return 0;
}
