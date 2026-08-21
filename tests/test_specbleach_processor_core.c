/*
libspecbleach - A spectral processing library

Copyright 2022-2026 Luciano Dato <lucianodato@gmail.com>

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

#include "processors/specbleach_processor_core.h"
#include "shared/configurations.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_processor_core_init_and_free(void) {
  printf("Testing processor core init and free...\n");

  // Invalid sample rate / frame size
  TEST_ASSERT(sb_processor_core_initialize(
                  0, 20.0f, OVERLAP_FACTOR_1D, PADDING_CONFIGURATION_1D,
                  ZEROPADDING_AMOUNT_1D, INPUT_WINDOW_TYPE_1D,
                  OUTPUT_WINDOW_TYPE_1D, 0) == NULL,
              "sample_rate == 0 should return NULL");
  TEST_ASSERT(sb_processor_core_initialize(
                  44100, 0.0f, OVERLAP_FACTOR_1D, PADDING_CONFIGURATION_1D,
                  ZEROPADDING_AMOUNT_1D, INPUT_WINDOW_TYPE_1D,
                  OUTPUT_WINDOW_TYPE_1D, 0) == NULL,
              "frame_size <= 0 should return NULL");

  // NULL safety for free
  sb_processor_core_free(NULL);

  // Real spectrum default init
  SbProcessorCore* core_default = sb_processor_core_initialize(
      44100, 20.0f, OVERLAP_FACTOR_1D, PADDING_CONFIGURATION_1D,
      ZEROPADDING_AMOUNT_1D, INPUT_WINDOW_TYPE_1D, OUTPUT_WINDOW_TYPE_1D,
      SB_PROCESSOR_CORE_DEFAULT_REAL_SPECTRUM);
  TEST_ASSERT(core_default != NULL, "Init with default spectrum size");
  TEST_ASSERT(sb_processor_core_get_noise_profile_size(core_default) > 0,
              "Valid noise profile size");

  // Full FFT spectrum init
  SbProcessorCore* core_full = sb_processor_core_initialize(
      44100, 20.0f, OVERLAP_FACTOR_2D, PADDING_CONFIGURATION_2D,
      ZEROPADDING_AMOUNT_2D, INPUT_WINDOW_TYPE_2D, OUTPUT_WINDOW_TYPE_2D,
      SB_PROCESSOR_CORE_FULL_FFT_SPECTRUM);
  TEST_ASSERT(core_full != NULL, "Init with full FFT spectrum size");
  TEST_ASSERT(sb_processor_core_get_noise_profile_size(core_full) >
                  sb_processor_core_get_noise_profile_size(core_default),
              "Full FFT spectrum size should be larger than real spectrum");

  sb_processor_core_free(core_default);
  sb_processor_core_free(core_full);
  printf("✓ Processor core init and free tests passed\n");
}

void test_processor_core_noise_profile_accessors(void) {
  printf("Testing processor core noise profile accessors...\n");

  // NULL safety tests
  TEST_ASSERT(sb_processor_core_get_noise_profile_size(NULL) == 0,
              "NULL get profile size");
  TEST_ASSERT(sb_processor_core_load_noise_profile_for_mode(NULL, NULL, 0, 0,
                                                            1) == false,
              "NULL load noise profile");
  TEST_ASSERT(sb_processor_core_reset_noise_profile(NULL) == false,
              "NULL reset profile");
  TEST_ASSERT(
      sb_processor_core_get_noise_profile_block_count_for_mode(NULL, 1) == 0,
      "NULL block count");
  TEST_ASSERT(sb_processor_core_get_noise_profile_for_mode(NULL, 1) == NULL,
              "NULL get profile for mode");
  TEST_ASSERT(
      sb_processor_core_noise_profile_available_for_mode(NULL, 1) == false,
      "NULL noise profile available");

  SbProcessorCore* core = sb_processor_core_initialize(
      44100, 20.0f, OVERLAP_FACTOR_1D, PADDING_CONFIGURATION_1D,
      ZEROPADDING_AMOUNT_1D, INPUT_WINDOW_TYPE_1D, OUTPUT_WINDOW_TYPE_1D, 0);
  TEST_ASSERT(core != NULL, "Valid core");

  uint32_t profile_size = sb_processor_core_get_noise_profile_size(core);
  float* test_profile = (float*)malloc(profile_size * sizeof(float));
  for (uint32_t i = 0; i < profile_size; i++) {
    test_profile[i] = 1.0f;
  }

  // Test mode out of bounds
  TEST_ASSERT(sb_processor_core_load_noise_profile_for_mode(
                  core, test_profile, profile_size, 5, 0) == false,
              "Mode 0 load should fail");
  TEST_ASSERT(sb_processor_core_load_noise_profile_for_mode(
                  core, test_profile, profile_size, 5, 5) == false,
              "Mode 5 load should fail");
  TEST_ASSERT(sb_processor_core_load_noise_profile_for_mode(
                  core, test_profile, profile_size + 1, 5, 1) == false,
              "Wrong profile size should fail");

  // Valid load
  TEST_ASSERT(sb_processor_core_load_noise_profile_for_mode(
                  core, test_profile, profile_size, 10, 1) == true,
              "Valid profile load mode 1");
  TEST_ASSERT(
      sb_processor_core_noise_profile_available_for_mode(core, 1) == true,
      "Profile 1 available");
  TEST_ASSERT(
      sb_processor_core_get_noise_profile_block_count_for_mode(core, 1) == 10,
      "Block count 10");
  TEST_ASSERT(sb_processor_core_get_noise_profile_for_mode(core, 1) != NULL,
              "Profile ptr valid");

  // Invalid mode getters
  TEST_ASSERT(
      sb_processor_core_get_noise_profile_block_count_for_mode(core, 0) == 0,
      "Mode 0 count");
  TEST_ASSERT(
      sb_processor_core_get_noise_profile_block_count_for_mode(core, 5) == 0,
      "Mode 5 count");
  TEST_ASSERT(sb_processor_core_get_noise_profile_for_mode(core, 0) == NULL,
              "Mode 0 ptr");
  TEST_ASSERT(sb_processor_core_get_noise_profile_for_mode(core, 5) == NULL,
              "Mode 5 ptr");
  TEST_ASSERT(
      sb_processor_core_noise_profile_available_for_mode(core, 0) == false,
      "Mode 0 avail");
  TEST_ASSERT(
      sb_processor_core_noise_profile_available_for_mode(core, 5) == false,
      "Mode 5 avail");

  // Reset
  TEST_ASSERT(sb_processor_core_reset_noise_profile(core) == true, "Reset");
  TEST_ASSERT(
      sb_processor_core_noise_profile_available_for_mode(core, 1) == false,
      "Profile 1 unavailable after reset");

  free(test_profile);
  sb_processor_core_free(core);
  printf("✓ Processor core noise profile accessors tests passed\n");
}

void test_params_sanitize_functions(void) {
  printf("Testing params sanitize functions...\n");

  SpectralBleachDenoiserParameters p1 = {
      .learn_noise = true,
      .residual_listen = false,
      .reduction_amount = -12.0f,
      .smoothing_factor = 50.0f,
      .whitening_factor = 100.0f,
      .adaptive_noise = false,
      .noise_estimation_method = 1,
      .masking_depth = 0.5f,
      .suppression_strength = 80.0f,
      .aggressiveness = 1.0f,
      .tonal_reduction = -6.0f,
      .hpss_quality_mode = 1,
      .hpss_sensitivity = 0.7f,
      .noise_profile_offset_db = 2.0f,
      .reduction_curve_bias = NULL,
      .reduction_curve_enabled = false,
  };

  DenoiserParameters sp1 = sb_denoiser_params_sanitize(p1);
  TEST_ASSERT(sp1.learn_noise == true, "Sanitize learn_noise");
  TEST_FLOAT_CLOSE(sp1.whitening_factor, 1.0f, 1e-4f);
  TEST_FLOAT_CLOSE(sp1.suppression_strength, 0.8f, 1e-4f);
  TEST_FLOAT_CLOSE(sp1.hpss_sensitivity, 0.7f, 1e-4f);
  TEST_FLOAT_CLOSE(sp1.noise_profile_offset_linear, 1.2589f, 1e-3f);
  TEST_ASSERT(sp1.reduction_curve_bias == NULL, "Sanitize disabled curve bias");

  SpectralBleach2DDenoiserParameters p2 = {
      .learn_noise = false,
      .residual_listen = true,
      .reduction_amount = -18.0f,
      .smoothing_factor = 0.3f,
      .whitening_factor = 50.0f,
      .adaptive_noise = true,
      .noise_estimation_method = 2,
      .nlm_masking_protection = 0.8f,
      .suppression_strength = 100.0f,
      .aggressiveness = 2.0f,
      .tonal_reduction = -3.0f,
      .hpss_quality_mode = 2,
      .hpss_sensitivity = 0.85f,
      .noise_profile_offset_db = -6.0f,
      .reduction_curve_bias = (const float*)0x1234,
      .reduction_curve_enabled = true,
  };

  Denoiser2DParameters sp2 = sb_denoiser_2d_params_sanitize(p2);
  TEST_ASSERT(sp2.residual_listen == true, "Sanitize 2d residual_listen");
  TEST_FLOAT_CLOSE(sp2.whitening_factor, 0.5f, 1e-4f);
  TEST_FLOAT_CLOSE(sp2.suppression_strength, 1.0f, 1e-4f);
  TEST_FLOAT_CLOSE(sp2.hpss_sensitivity, 0.85f, 1e-4f);
  TEST_FLOAT_CLOSE(sp2.noise_profile_offset_linear, 0.5011f, 1e-3f);
  TEST_ASSERT(sp2.reduction_curve_bias == (const float*)0x1234,
              "Sanitize enabled curve bias ptr");

  // Bounds clamping checks for noise_profile_offset_db (-12dB to +12dB)
  SpectralBleachDenoiserParameters p_min = p1;
  p_min.noise_profile_offset_db = -20.0f;
  DenoiserParameters sp_min = sb_denoiser_params_sanitize(p_min);
  TEST_FLOAT_CLOSE(sp_min.noise_profile_offset_linear, 0.2512f, 1e-3f);

  SpectralBleachDenoiserParameters p_max = p1;
  p_max.noise_profile_offset_db = 20.0f;
  DenoiserParameters sp_max = sb_denoiser_params_sanitize(p_max);
  TEST_FLOAT_CLOSE(sp_max.noise_profile_offset_linear, 3.9811f, 1e-3f);

  SpectralBleach2DDenoiserParameters p2_min = p2;
  p2_min.noise_profile_offset_db = -20.0f;
  Denoiser2DParameters sp2_min = sb_denoiser_2d_params_sanitize(p2_min);
  TEST_FLOAT_CLOSE(sp2_min.noise_profile_offset_linear, 0.2512f, 1e-3f);

  SpectralBleach2DDenoiserParameters p2_max = p2;
  p2_max.noise_profile_offset_db = 20.0f;
  Denoiser2DParameters sp2_max = sb_denoiser_2d_params_sanitize(p2_max);
  TEST_FLOAT_CLOSE(sp2_max.noise_profile_offset_linear, 3.9811f, 1e-3f);

  printf("✓ Params sanitize functions tests passed\n");
}

int main(void) {
  printf("Running specbleach processor core tests...\n");

  test_processor_core_init_and_free();
  test_processor_core_noise_profile_accessors();
  test_params_sanitize_functions();

  printf("✅ All specbleach processor core tests passed!\n");
  return 0;
}
