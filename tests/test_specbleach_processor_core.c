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

  SpecbleachDenoiserParameters p1 = {
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .residual_listen = false,
      .reduction_gain = 0.25f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 1.0f,
      .adaptive_noise = false,
      .noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_BRANDT,
      .masking_depth = 0.5f,
      .suppression_strength = 0.8f,
      .aggressiveness = 1.0f,
      .tonal_reduction_gain = 0.5f,
      .hpss_enable = true,
      .noise_profile_scale = 1.5849f,
      .tonal_noise_profile_scale = 3.9811f,
      .reduction_curve_bias = NULL,
      .reduction_curve_enabled = false,
  };

  DenoiserParameters sp1 = sb_denoiser_params_sanitize(&p1);
  TEST_ASSERT(sp1.learn_noise == true, "Sanitize learn_noise");
  TEST_FLOAT_CLOSE(sp1.reduction_amount, 0.25f, 1e-4f);
  TEST_FLOAT_CLOSE(sp1.whitening_factor, 1.0f, 1e-4f);
  TEST_FLOAT_CLOSE(sp1.suppression_strength, 0.8f, 1e-4f);
  TEST_ASSERT(sp1.hpss_enable == 1, "Sanitize hpss_enable");
  TEST_FLOAT_CLOSE(sp1.noise_profile_offset_linear, 1.5849f, 1e-3f);
  TEST_FLOAT_CLOSE(sp1.tonal_noise_profile_offset_linear, 3.9811f, 1e-3f);
  TEST_ASSERT(sp1.reduction_curve_bias == NULL, "Sanitize disabled curve bias");

  Specbleach2DDenoiserParameters p2 = {
      .learn_noise = SPECBLEACH_LEARN_OFF,
      .residual_listen = true,
      .reduction_gain = 0.125f,
      .smoothing_factor = 0.3f,
      .whitening_factor = 0.5f,
      .adaptive_noise = true,
      .noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_MARTIN,
      .nlm_masking_protection = 0.8f,
      .suppression_strength = 1.0f,
      .aggressiveness = 2.0f,
      .tonal_reduction_gain = 0.7f,
      .hpss_enable = true,
      .noise_profile_scale = 0.2512f,
      .tonal_noise_profile_scale = 0.631f,
      .reduction_curve_bias = (const float*)0x1234,
      .reduction_curve_enabled = true,
  };

  Denoiser2DParameters sp2 = sb_denoiser_2d_params_sanitize(&p2);
  TEST_ASSERT(sp2.residual_listen == true, "Sanitize 2d residual_listen");
  TEST_FLOAT_CLOSE(sp2.reduction_amount, 0.125f, 1e-4f);
  TEST_FLOAT_CLOSE(sp2.whitening_factor, 0.5f, 1e-4f);
  TEST_FLOAT_CLOSE(sp2.suppression_strength, 1.0f, 1e-4f);
  TEST_ASSERT(sp2.hpss_enable == 1, "Sanitize 2d hpss_enable");
  TEST_FLOAT_CLOSE(sp2.noise_profile_offset_linear, 0.2512f, 1e-3f);
  TEST_FLOAT_CLOSE(sp2.tonal_noise_profile_offset_linear, 0.631f, 1e-3f);
  TEST_ASSERT(sp2.reduction_curve_bias == NULL,
              "Sanitize never passes caller curve ptr");

  // Bounds clamping checks for linear parameters
  SpecbleachDenoiserParameters p_min = p1;
  p_min.reduction_gain = -5.0f;
  p_min.noise_profile_scale = -2.0f;
  p_min.tonal_noise_profile_scale = -2.0f;
  DenoiserParameters sp_min = sb_denoiser_params_sanitize(&p_min);
  TEST_FLOAT_CLOSE(sp_min.reduction_amount, 0.0f, 1e-3f);
  TEST_FLOAT_CLOSE(sp_min.noise_profile_offset_linear, 1.0f, 1e-3f);
  TEST_FLOAT_CLOSE(sp_min.tonal_noise_profile_offset_linear, 1.0f, 1e-3f);

  SpecbleachDenoiserParameters p_scale_low = p1;
  p_scale_low.noise_profile_scale = 0.001f;
  p_scale_low.tonal_noise_profile_scale = 0.001f;
  DenoiserParameters sp_scale_low = sb_denoiser_params_sanitize(&p_scale_low);
  TEST_FLOAT_CLOSE(sp_scale_low.noise_profile_offset_linear, 0.01f, 1e-3f);
  TEST_FLOAT_CLOSE(sp_scale_low.tonal_noise_profile_offset_linear, 0.01f,
                   1e-3f);

  SpecbleachDenoiserParameters p_max = p1;
  p_max.reduction_gain = 5.0f;
  p_max.noise_profile_scale = 200.0f;
  p_max.tonal_noise_profile_scale = 200.0f;
  DenoiserParameters sp_max = sb_denoiser_params_sanitize(&p_max);
  TEST_FLOAT_CLOSE(sp_max.reduction_amount, 1.0f, 1e-3f);
  TEST_FLOAT_CLOSE(sp_max.noise_profile_offset_linear, 100.0f, 1e-3f);
  TEST_FLOAT_CLOSE(sp_max.tonal_noise_profile_offset_linear, 100.0f, 1e-3f);

  Specbleach2DDenoiserParameters p2_min = p2;
  p2_min.reduction_gain = -5.0f;
  p2_min.noise_profile_scale = -2.0f;
  p2_min.tonal_noise_profile_scale = -2.0f;
  Denoiser2DParameters sp2_min = sb_denoiser_2d_params_sanitize(&p2_min);
  TEST_FLOAT_CLOSE(sp2_min.reduction_amount, 0.0f, 1e-3f);
  TEST_FLOAT_CLOSE(sp2_min.noise_profile_offset_linear, 1.0f, 1e-3f);
  TEST_FLOAT_CLOSE(sp2_min.tonal_noise_profile_offset_linear, 1.0f, 1e-3f);

  Specbleach2DDenoiserParameters p2_scale_low = p2;
  p2_scale_low.noise_profile_scale = 0.001f;
  p2_scale_low.tonal_noise_profile_scale = 0.001f;
  Denoiser2DParameters sp2_scale_low =
      sb_denoiser_2d_params_sanitize(&p2_scale_low);
  TEST_FLOAT_CLOSE(sp2_scale_low.noise_profile_offset_linear, 0.01f, 1e-3f);
  TEST_FLOAT_CLOSE(sp2_scale_low.tonal_noise_profile_offset_linear, 0.01f,
                   1e-3f);

  Specbleach2DDenoiserParameters p2_max = p2;
  p2_max.reduction_gain = 5.0f;
  p2_max.noise_profile_scale = 200.0f;
  p2_max.tonal_noise_profile_scale = 200.0f;
  Denoiser2DParameters sp2_max = sb_denoiser_2d_params_sanitize(&p2_max);
  TEST_FLOAT_CLOSE(sp2_max.reduction_amount, 1.0f, 1e-3f);
  TEST_FLOAT_CLOSE(sp2_max.noise_profile_offset_linear, 100.0f, 1e-3f);
  TEST_FLOAT_CLOSE(sp2_max.tonal_noise_profile_offset_linear, 100.0f, 1e-3f);

  printf("✓ Params sanitize functions tests passed\n");
}

void test_curve_bias_copy(void) {
  printf("Testing curve bias copy helper...\n");

  float* buffer = NULL;
  uint32_t capacity = 0;

  // Disabled curve yields NULL and allocates nothing
  const float source_disabled[3] = {1.0f, 2.0f, 3.0f};
  TEST_ASSERT(
      sb_curve_bias_copy(&buffer, &capacity, 3, false, source_disabled) == NULL,
      "Copy disabled returns NULL");
  TEST_ASSERT(buffer == NULL, "Copy disabled does not allocate");
  TEST_ASSERT(capacity == 0, "Copy disabled leaves capacity zero");

  // NULL source yields NULL
  TEST_ASSERT(sb_curve_bias_copy(&buffer, &capacity, 3, true, NULL) == NULL,
              "Copy NULL source returns NULL");

  // Zero size yields NULL
  TEST_ASSERT(
      sb_curve_bias_copy(&buffer, &capacity, 0, true, source_disabled) == NULL,
      "Copy zero size returns NULL");

  // Enabled curve copies contents into owned memory
  float source[4] = {0.5f, -1.5f, 2.25f, 0.0f};
  const float* copied = sb_curve_bias_copy(&buffer, &capacity, 4, true, source);
  TEST_ASSERT(copied != NULL, "Copy enabled returns pointer");
  TEST_ASSERT(copied != source, "Copy returns owned buffer, not source");
  TEST_ASSERT(capacity == 4, "Capacity tracks required size");
  for (uint32_t i = 0; i < 4; ++i) {
    TEST_FLOAT_CLOSE(copied[i], source[i], 1e-6f);
  }

  // Caller may mutate/free its source afterwards without affecting the copy
  source[0] = 99.0f;
  TEST_FLOAT_CLOSE(copied[0], 0.5f, 1e-6f);

  // Same-size reload reuses the buffer (stable pointer, updated contents)
  const float* previous = buffer;
  float source_b[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  const float* recopied =
      sb_curve_bias_copy(&buffer, &capacity, 4, true, source_b);
  TEST_ASSERT(recopied == previous, "Reload reuses existing buffer");
  TEST_FLOAT_CLOSE(recopied[0], 1.0f, 1e-6f);

  // Growth path reallocates and preserves capacity semantics
  float big_source[8] = {0.0f};
  const float* grown =
      sb_curve_bias_copy(&buffer, &capacity, 8, true, big_source);
  TEST_ASSERT(grown != NULL, "Growth copy succeeds");
  TEST_ASSERT(capacity == 8, "Capacity grows to required size");
  TEST_FLOAT_CLOSE(grown[7], 0.0f, 1e-6f);

  // Shrink case reuses allocation and preserves original capacity
  float small_source[4] = {7.0f, 8.0f, 9.0f, 10.0f};
  const float* shrunk =
      sb_curve_bias_copy(&buffer, &capacity, 4, true, small_source);
  TEST_ASSERT(shrunk != NULL, "Shrink copy returns valid buffer");
  TEST_ASSERT(capacity == 8, "Capacity preserved on shrink");
  TEST_ASSERT(shrunk == grown, "Shrink reuses existing buffer");
  for (uint32_t i = 0; i < 4; ++i) {
    TEST_FLOAT_CLOSE(shrunk[i], small_source[i], 1e-6f);
  }

  free(buffer);
  printf("✓ Curve bias copy tests passed\n");
}

int main(void) {
  printf("Running specbleach processor core tests...\n");

  test_processor_core_init_and_free();
  test_processor_core_noise_profile_accessors();
  test_params_sanitize_functions();
  test_curve_bias_copy();

  printf("✅ All specbleach processor core tests passed!\n");
  return 0;
}
