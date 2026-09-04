/*
 * Unit tests for SpectralBleach Denoiser API functions
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/denoiser_logic/processing/hpss_filter.h"
#include "specbleach_denoiser.h"

#define ROLLING_MEAN 1

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_specbleach_noise_profile_mode_functions(void) {
  printf("Testing specbleach noise profile mode functions...\n");

  // Initialize denoiser
  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  SpecbleachDenoiserParameters params = {
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading parameters should succeed");

  // Test NULL handle cases
  TEST_ASSERT(
      specbleach_denoiser_get_noise_profile_block_count_for_mode(NULL, 1) == 0,
      "NULL handle should return 0");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_for_mode(NULL, 1) == NULL,
              "NULL handle should return NULL");
  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(NULL, 1) == false,
      "NULL handle should return false");

  // Test invalid modes
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_block_count_for_mode(
                  handle, 0) == 0,
              "Invalid mode should return 0");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_block_count_for_mode(
                  handle, 5) == 0,
              "Invalid mode should return 0");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_for_mode(handle, 0) == NULL,
              "Invalid mode should return NULL");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_for_mode(handle, 5) == NULL,
              "Invalid mode should return NULL");
  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(handle, 0) == false,
      "Invalid mode should return false");
  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(handle, 5) == false,
      "Invalid mode should return false");

  // Initially no profiles should be available
  for (int mode = 1; mode <= 4; mode++) {
    TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(
                    handle, mode) == false,
                "Profile should not be available initially");
    TEST_ASSERT(specbleach_denoiser_get_noise_profile_block_count_for_mode(
                    handle, mode) == 0,
                "Should have 0 blocks initially");
    // Note: get_noise_profile_for_mode returns the profile array even if not
    // available (this matches the behavior of the original single-mode API)
    TEST_ASSERT(
        specbleach_denoiser_get_noise_profile_for_mode(handle, mode) != NULL,
        "Should return profile array even if not available");
  }

  specbleach_denoiser_free(handle);
  printf("✓ Specbleach noise profile mode functions tests passed\n");
}

void test_specbleach_load_noise_profile_with_mode(void) {
  printf("Testing specbleach load noise profile with mode...\n");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  // Load
  SpecbleachDenoiserParameters params = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading parameters should succeed");

  // Get profile size
  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);
  TEST_ASSERT(profile_size > 0, "Profile size should be valid");

  // Create test profile
  float* test_profile = (float*)malloc(profile_size * sizeof(float));
  TEST_ASSERT(test_profile != NULL, "Profile allocation should succeed");

  for (uint32_t i = 0; i < profile_size; i++) {
    test_profile[i] = 0.1f + ((float)i * 0.001f);
  }

  // Test loading profile with mode
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, ROLLING_MEAN) == true,
              "Loading noise profile should succeed");

  // Verify profile is available and has correct block count
  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(handle, 1) == true,
      "Profile should be available for mode 1");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_block_count_for_mode(
                  handle, 1) == 10,
              "Should have 10 blocks averaged");

  // Get the profile back and verify it matches
  const float* retrieved_profile =
      specbleach_denoiser_get_noise_profile_for_mode(handle, 1);
  TEST_ASSERT(retrieved_profile != NULL, "Should get valid profile");

  for (uint32_t i = 0; i < profile_size; i++) {
    TEST_FLOAT_CLOSE(retrieved_profile[i], test_profile[i], 0.001f);
  }

  // Test that other modes are unaffected
  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(handle, 2) == false,
      "Mode 2 should still be unavailable");
  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(handle, 3) == false,
      "Mode 3 should still be unavailable");

  free(test_profile);
  specbleach_denoiser_free(handle);
  printf("✓ Specbleach load noise profile with mode tests passed\n");
}

void test_specbleach_mode_switching(void) {
  printf("Testing specbleach mode switching...\n");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);

  // Create different profiles for each mode
  float* profile1 = (float*)malloc(profile_size * sizeof(float));
  float* profile2 = (float*)malloc(profile_size * sizeof(float));
  float* profile3 = (float*)malloc(profile_size * sizeof(float));

  TEST_ASSERT(profile1 != NULL && profile2 != NULL && profile3 != NULL,
              "Profile allocations should succeed");

  for (uint32_t i = 0; i < profile_size; i++) {
    profile1[i] = 1.0f;
    profile2[i] = 2.0f;
    profile3[i] = 3.0f;
  }

  // Load profiles for each mode
  SpecbleachDenoiserParameters params;

  // Mode 1
  params = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading mode 1 parameters should succeed");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, profile1, profile_size, 5, ROLLING_MEAN) == true,
              "Loading mode 1 profile should succeed");

  // Mode 2
  // Explicitly loading for mode 2 since noise_reduction_mode is removed
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading mode 2 parameters should succeed");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, profile2, profile_size, 10, 2) == true,
              "Loading mode 2 profile should succeed");

  // Mode 3
  // Explicitly loading for mode 3
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading mode 3 parameters should succeed");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, profile3, profile_size, 15, 3) == true,
              "Loading mode 3 profile should succeed");

  // Verify all modes have their profiles
  for (int mode = 1; mode <= 3; mode++) {
    TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(
                    handle, mode) == true,
                "Profile should be available");
  }

  // Check block counts
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_block_count_for_mode(
                  handle, 1) == 5,
              "Mode 1 should have 5 blocks");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_block_count_for_mode(
                  handle, 2) == 10,
              "Mode 2 should have 10 blocks");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_block_count_for_mode(
                  handle, 3) == 15,
              "Mode 3 should have 15 blocks");

  // Check profile values
  const float* retrieved1 =
      specbleach_denoiser_get_noise_profile_for_mode(handle, 1);
  const float* retrieved2 =
      specbleach_denoiser_get_noise_profile_for_mode(handle, 2);
  const float* retrieved3 =
      specbleach_denoiser_get_noise_profile_for_mode(handle, 3);

  TEST_ASSERT(retrieved1 != NULL && retrieved2 != NULL && retrieved3 != NULL,
              "All profiles should be retrievable");

  for (uint32_t i = 0; i < profile_size; i++) {
    TEST_FLOAT_CLOSE(retrieved1[i], 1.0f, 0.001f);
    TEST_FLOAT_CLOSE(retrieved2[i], 2.0f, 0.001f);
    TEST_FLOAT_CLOSE(retrieved3[i], 3.0f, 0.001f);
  }

  free(profile1);
  free(profile2);
  free(profile3);
  specbleach_denoiser_free(handle);
  printf("✓ Specbleach mode switching tests passed\n");
}

void test_specbleach_reset_noise_profile(void) {
  printf("Testing specbleach reset noise profile...\n");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);
  float* test_profile = (float*)malloc(profile_size * sizeof(float));
  TEST_ASSERT(test_profile != NULL, "Profile allocation should succeed");

  for (uint32_t i = 0; i < profile_size; i++) {
    test_profile[i] = 1.0f;
  }

  // Load parameters and profile
  SpecbleachDenoiserParameters params = {
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading parameters should succeed");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, ROLLING_MEAN) == true,
              "Loading profile should succeed");

  // Verify profile is available
  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(handle, 1) == true,
      "Profile should be available before reset");

  // Reset profile (void: NULL-safe no-op, cannot fail on valid handle)
  specbleach_denoiser_reset_noise_profile(handle);

  // Verify all profiles are reset
  for (int mode = 1; mode <= 4; mode++) {
    TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(
                    handle, mode) == false,
                "Profile should not be available after reset");
    TEST_ASSERT(specbleach_denoiser_get_noise_profile_block_count_for_mode(
                    handle, mode) == 0,
                "Blocks should be 0 after reset");
  }

  free(test_profile);
  specbleach_denoiser_free(handle);
  printf("✓ Specbleach reset noise profile tests passed\n");
}

void test_specbleach_load_noise_profile_for_mode(void) {
  printf("Testing specbleach load noise profile for mode...\n");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  // Get profile size
  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);
  TEST_ASSERT(profile_size > 0, "Profile size should be valid");

  // Create test profile
  float* test_profile = (float*)malloc(profile_size * sizeof(float));
  TEST_ASSERT(test_profile != NULL, "Profile allocation should succeed");

  for (uint32_t i = 0; i < profile_size; i++) {
    test_profile[i] = 0.1f + ((float)i * 0.001f);
  }

  // Test loading profile for each mode directly
  for (int mode = 1; mode <= 4; mode++) {
    TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                    handle, test_profile, profile_size, 10, mode) == true,
                "Loading noise profile for mode should succeed");
  }

  // Test invalid mode
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, 0) == false,
              "Loading noise profile for invalid mode should fail");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, 5) == false,
              "Loading noise profile for invalid mode should fail");

  free(test_profile);
  specbleach_denoiser_free(handle);

  printf("✓ Specbleach load noise profile for mode tests passed\n");
}

void test_specbleach_run_features(void) {
  printf("Testing specbleach run features (whitening, residual listen)...\n");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  float* input = (float*)malloc(1024 * sizeof(float));
  float* output = (float*)malloc(1024 * sizeof(float));
  for (int i = 0; i < 1024; i++) {
    input[i] = 0.5f;
    output[i] = 0.0f;
  }

  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);
  float* profile = (float*)malloc(profile_size * sizeof(float));
  for (uint32_t i = 0; i < profile_size; i++) {
    profile[i] = 0.1f;
  }

  specbleach_denoiser_load_noise_profile_for_mode(handle, profile, profile_size,
                                                  1, ROLLING_MEAN);

  SpecbleachDenoiserParameters params = {
      .learn_noise = SPECBLEACH_LEARN_OFF,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 1.0f, // Test whitening
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };

  specbleach_denoiser_load_parameters(handle, &params, sizeof(params));
  specbleach_denoiser_process(handle, 1024, input, output);

  // Test residual listen
  params.residual_listen = true;
  params.whitening_factor = 0.0f;
  specbleach_denoiser_load_parameters(handle, &params, sizeof(params));
  specbleach_denoiser_process(handle, 1024, input, output);

  free(input);
  free(output);
  free(profile);
  specbleach_denoiser_free(handle);
  printf("✓ Specbleach denoiser features tests passed\n");
}

void test_specbleach_silence_bypass(void) {
  printf("Testing specbleach silence bypass with active profile...\n");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);
  float* profile = (float*)malloc(profile_size * sizeof(float));
  for (uint32_t i = 0; i < profile_size; i++) {
    profile[i] = 0.1f;
  }
  specbleach_denoiser_load_noise_profile_for_mode(handle, profile, profile_size,
                                                  1, ROLLING_MEAN);

  SpecbleachDenoiserParameters params = {
      .learn_noise = SPECBLEACH_LEARN_OFF,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.5f,
  };
  specbleach_denoiser_load_parameters(handle, &params, sizeof(params));

  // Process completely silent input buffer
  float silent_in[1024] = {0};
  float silent_out[1024] = {0};
  specbleach_denoiser_process(handle, 1024, silent_in, silent_out);

  // Change aggressiveness and process silent again
  params.aggressiveness = -0.5f;
  specbleach_denoiser_load_parameters(handle, &params, sizeof(params));
  specbleach_denoiser_process(handle, 1024, silent_in, silent_out);

  for (int i = 0; i < 1024; i++) {
    TEST_ASSERT(silent_out[i] == 0.0f,
                "Silent input must produce silent output");
  }

  free(profile);
  specbleach_denoiser_free(handle);
  printf("✓ Specbleach silence bypass test passed\n");
}

void test_specbleach_smoothing_transition_and_validation(void) {
  printf("Testing smoothing mode transitions and API validation...\n");

  // Invalid initialization arguments
  TEST_ASSERT(specbleach_denoiser_initialize(0, 20.0f) == NULL,
              "Zero sample rate must fail initialization");
  TEST_ASSERT(specbleach_denoiser_initialize(44100, 0.0f) == NULL,
              "Zero frame size must fail initialization");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  // Process with NULL buffers must be rejected
  float in_buf[1024];
  float out_buf[1024];
  float transient_like[1024];
  for (int i = 0; i < 1024; ++i) {
    in_buf[i] = 0.1f * ((float)(i % 64) / 64.0f);
    transient_like[i] = (i % 2 == 0) ? 0.9f : -0.9f; // sharp alternating onset
  }
  TEST_ASSERT(specbleach_denoiser_process(handle, 1024, NULL, out_buf) == false,
              "NULL input must be rejected");
  TEST_ASSERT(specbleach_denoiser_process(handle, 1024, in_buf, NULL) == false,
              "NULL output must be rejected");

  SpecbleachDenoiserParameters params = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_OFF,
      .reduction_gain = 0.3f,
      .smoothing_factor = 0.5f,
      .smoothing_mode = SPECBLEACH_SMOOTHING_NLM_2D,
  };
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading parameters should succeed");

  // Wrong parameter struct size must be rejected
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params) - 8) == false,
              "Wrong parameters size must be rejected");

  // Load a profile so both chains run with real data
  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);
  TEST_ASSERT(profile_size > 0, "Profile size should be valid");
  float* profile = (float*)malloc(profile_size * sizeof(float));
  TEST_ASSERT(profile != NULL, "Profile allocation should succeed");
  for (uint32_t i = 0; i < profile_size; i++) {
    profile[i] = 0.1f;
  }
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, profile, profile_size, 10, ROLLING_MEAN) == true,
              "Loading noise profile should succeed");

  // Process in NLM mode to establish the active chain
  for (int f = 0; f < 5; ++f) {
    specbleach_denoiser_process(handle, 1024, in_buf, out_buf);
  }

  // Switch to TEMPORAL at runtime: exercises the crossfade transition
  // machinery and the temporal smoothing chain
  params.smoothing_mode = SPECBLEACH_SMOOTHING_TEMPORAL;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading temporal mode should succeed");
  for (int f = 0; f < 15; ++f) {
    specbleach_denoiser_process(handle, 1024, in_buf, out_buf);
    TEST_ASSERT(specbleach_denoiser_get_latency(handle) > 0,
                "Latency must stay positive during transition");
  }
  for (int i = 0; i < 1024; ++i) {
    TEST_ASSERT(isfinite(out_buf[i]) && out_buf[i] < 100.0f,
                "Output must be finite during/after transition");
  }

  // Transient getters with a valid handle
  (void)specbleach_denoiser_is_transient_detected(handle);
  (void)specbleach_denoiser_get_transient_intensity(handle);

  // TEMPORAL -> NLM+DFTT at runtime (crossfade into the refined map), then
  // NLM+DFTT -> NLM (intra-family instant flip, no crossfade). DFTT rings
  // stay warm on every NLM pass so both directions stay finite.
  params.smoothing_mode = SPECBLEACH_SMOOTHING_NLM_2D_DFTT;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading NLM+DFTT mode should succeed");
  for (int f = 0; f < 15; ++f) {
    specbleach_denoiser_process(handle, 1024, in_buf, out_buf);
    TEST_ASSERT(specbleach_denoiser_get_latency(handle) > 0,
                "Latency must stay positive in DFTT mode");
  }
  for (int i = 0; i < 1024; ++i) {
    TEST_ASSERT(isfinite(out_buf[i]) && out_buf[i] < 100.0f,
                "Output must be finite in DFTT mode");
  }
  params.smoothing_mode = SPECBLEACH_SMOOTHING_NLM_2D;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading NLM mode should succeed");
  for (int f = 0; f < 5; ++f) {
    specbleach_denoiser_process(handle, 1024, in_buf, out_buf);
  }
  for (int i = 0; i < 1024; ++i) {
    TEST_ASSERT(isfinite(out_buf[i]), "Output must be finite after flip back");
  }

  // Transient during NLM mode: transient-mask loops in the NLM chain
  params.hpss_enable = true;
  params.smoothing_mode = SPECBLEACH_SMOOTHING_NLM_2D;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading NLM with HPSS should succeed");
  for (int f = 0; f < 15; ++f) {
    specbleach_denoiser_process(handle, 1024, transient_like, out_buf);
  }
  params.hpss_enable = false;

  // Profile with mismatched size must be rejected
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, profile, profile_size + 1, 10, ROLLING_MEAN) == false,
              "Wrong profile size must be rejected");

  // Reduction curve with mismatched size must be rejected
  uint32_t bad_size = specbleach_denoiser_get_noise_profile_size(handle);
  float* bad_curve = calloc(bad_size, sizeof(float));
  params.reduction_curve_enabled = true;
  params.reduction_curve_bias = bad_curve;
  params.reduction_curve_size = bad_size > 0 ? bad_size - 1 : 0;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == false,
              "Wrong curve size must be rejected");

  // Out-of-range scales fall back to safe defaults through the sanitizers
  params.reduction_curve_enabled = false;
  params.reduction_curve_bias = NULL;
  params.reduction_curve_size = 0;
  params.noise_profile_scale = 0.0f;
  params.tonal_noise_profile_scale = 0.0f;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Zero scales must be sanitized");
  free(bad_curve);

  // Transient accessors with NULL handle
  TEST_ASSERT(specbleach_denoiser_is_transient_detected(NULL) == false,
              "NULL transient detected must be false");
  TEST_ASSERT(specbleach_denoiser_get_transient_intensity(NULL) == 0.0f,
              "NULL transient intensity must be 0");

  free(profile);
  specbleach_denoiser_free(handle);
  printf("✓ Smoothing transition and validation test passed\n");
}

void test_specbleach_adaptive_method_switch(void) {
  printf("Testing adaptive estimation method re-initialization...\n");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  float in_buf[1024];
  float out_buf[1024];
  for (int i = 0; i < 1024; ++i) {
    in_buf[i] = 0.1f;
  }

  SpecbleachDenoiserParameters params = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_OFF,
      .adaptive_noise = true,
      .noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE,
      .reduction_gain = 0.2f,
  };
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Loading adaptive SPP-MMSE should succeed");
  for (int f = 0; f < 5; ++f) {
    specbleach_denoiser_process(handle, 1024, in_buf, out_buf);
  }

  // Switching method must re-initialize the adaptive estimator
  params.noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_MARTIN;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Switching to Martin must succeed");
  for (int f = 0; f < 5; ++f) {
    specbleach_denoiser_process(handle, 1024, in_buf, out_buf);
  }

  params.noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_BRANDT;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &params,
                                                  sizeof(params)) == true,
              "Switching to Brandt must succeed");
  for (int f = 0; f < 5; ++f) {
    specbleach_denoiser_process(handle, 1024, in_buf, out_buf);
  }

  specbleach_denoiser_free(handle);
  printf("✓ Adaptive method switch test passed\n");
}

void test_specbleach_redesigned_api_coverage(void) {
  printf("Testing redesigned API coverage...\n");

  // Init validation: zero/negative parameters and sub-hop frame must fail
  TEST_ASSERT(specbleach_denoiser_initialize(0, 20.0f) == NULL,
              "Zero sample rate should fail");
  TEST_ASSERT(specbleach_denoiser_initialize(44100, 0.0f) == NULL,
              "Zero frame size should fail");
  TEST_ASSERT(specbleach_denoiser_initialize(44100, -5.0f) == NULL,
              "Negative frame size should fail");
  TEST_ASSERT(specbleach_denoiser_initialize(44100, 0.01f) == NULL,
              "Sub-hop frame size should fail");

  // Default parameters carry safe, documented values
  SpecbleachDenoiserParameters defaults =
      specbleach_denoiser_get_default_parameters();
  TEST_ASSERT(defaults.learn_noise == SPECBLEACH_LEARN_OFF,
              "Default should not learn");
  TEST_ASSERT(defaults.smoothing_mode == SPECBLEACH_SMOOTHING_TEMPORAL,
              "Default smoothing should be temporal");
  TEST_ASSERT(defaults.reduction_curve_bias == NULL,
              "Default curve bias should be NULL");

  specbleach_denoiser* handle = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  // Geometry getters on a live instance and NULL
  TEST_ASSERT(specbleach_denoiser_get_sample_rate(handle) == 44100,
              "Sample rate getter should match");
  TEST_ASSERT(specbleach_denoiser_get_sample_rate(NULL) == 0,
              "NULL sample rate should be 0");
  uint32_t frame = specbleach_denoiser_get_frame_size(handle);
  TEST_ASSERT(frame == 882, "20 ms at 44.1 kHz should be 882 samples");
  TEST_ASSERT(specbleach_denoiser_get_frame_size(NULL) == 0,
              "NULL frame size should be 0");
  uint32_t fft = specbleach_denoiser_get_fft_size(handle);
  TEST_ASSERT(fft >= frame && fft % 32 == 0, "FFT should pad frame to 32x");
  TEST_ASSERT(specbleach_denoiser_get_fft_size(NULL) == 0,
              "NULL FFT size should be 0");
  uint32_t hop = specbleach_denoiser_get_hop_size(handle);
  TEST_ASSERT(hop == frame / 4, "Hop should be frame / overlap 4");
  TEST_ASSERT(specbleach_denoiser_get_hop_size(NULL) == 0,
              "NULL hop should be 0");
  // NLM look-ahead is ms-anchored (#152, past-heavy per #155): 20ms frame
  // -> ~5ms hop -> round(46ms/5ms) = 9 (46ms default is 4). Latency in
  // ms is sample-rate independent.
  TEST_ASSERT(specbleach_denoiser_get_latency(handle) == frame + 9 * hop,
              "Latency should be frame plus ms-anchored NLM look-ahead");

  // No profile available before learning
  TEST_ASSERT(specbleach_denoiser_has_any_profile(NULL) == false,
              "NULL should have no profile");
  TEST_ASSERT(specbleach_denoiser_has_any_profile(handle) == false,
              "Fresh instance should have no profile");

  // Process and parameter failure paths report false (no error codes;
  // fallible calls fail fast with bool)
  float in_buf[256] = {0};
  float out_buf[256] = {0};
  TEST_ASSERT(specbleach_denoiser_process(NULL, 256, in_buf, out_buf) == false,
              "NULL instance should fail");
  TEST_ASSERT(specbleach_denoiser_process(handle, 256, NULL, out_buf) == false,
              "NULL input should fail");
  TEST_ASSERT(specbleach_denoiser_process(handle, 0, in_buf, out_buf) == false,
              "Empty block should fail");
  TEST_ASSERT(specbleach_denoiser_load_parameters(NULL, &defaults,
                                                  sizeof(defaults)) == false,
              "NULL instance params should fail");
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, NULL,
                                                  sizeof(defaults)) == false,
              "NULL params should fail");
  TEST_ASSERT(specbleach_denoiser_load_parameters(
                  handle, &defaults, sizeof(defaults) - 1) == false,
              "ABI-mismatched params should fail");

  // Non-positive profile scales fall back to 1.0x
  SpecbleachDenoiserParameters flat_scales = defaults;
  flat_scales.noise_profile_scale = 0.0f;
  flat_scales.tonal_noise_profile_scale = -1.0f;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &flat_scales,
                                                  sizeof(flat_scales)) == true,
              "Non-positive scales should fall back");

  // Positive profile scales take the direct path
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &defaults,
                                                  sizeof(defaults)) == true,
              "Default params should load");

  // Enabled curve without bias buffer reports NULL_ARG
  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);
  SpecbleachDenoiserParameters bad_curve = defaults;
  bad_curve.reduction_curve_enabled = true;
  bad_curve.reduction_curve_bias = NULL;
  bad_curve.reduction_curve_size = profile_size;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &bad_curve,
                                                  sizeof(bad_curve)) == false,
              "Curve without bias should fail");

  // Resampled load: same size delegates, other sizes interpolate
  float* native_profile = (float*)malloc(profile_size * sizeof(float));
  TEST_ASSERT(native_profile != NULL, "Profile allocation should succeed");
  for (uint32_t i = 0; i < profile_size; i++) {
    native_profile[i] = 0.1f + ((float)i * 0.001f);
  }
  TEST_ASSERT(
      specbleach_denoiser_load_noise_profile_resampled(
          handle, native_profile, profile_size, 8, ROLLING_MEAN) == true,
      "Same-size resampled load should succeed");
  TEST_ASSERT(specbleach_denoiser_has_any_profile(handle) == true,
              "Profile should exist after load");
  uint32_t half_size = profile_size / 2;
  float* half_profile = (float*)malloc(half_size * sizeof(float));
  TEST_ASSERT(half_profile != NULL, "Half profile allocation should succeed");
  for (uint32_t i = 0; i < half_size; i++) {
    half_profile[i] = 0.2f + ((float)i * 0.002f);
  }
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_resampled(
                  handle, half_profile, half_size, 4, ROLLING_MEAN) == true,
              "Half-size resampled load should succeed");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_resampled(
                  NULL, native_profile, profile_size, 1, ROLLING_MEAN) == false,
              "NULL resampled load should fail");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_resampled(
                  handle, NULL, profile_size, 1, ROLLING_MEAN) == false,
              "NULL resampled source should fail");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_resampled(
                  handle, native_profile, 0, 1, ROLLING_MEAN) == false,
              "Empty resampled source should fail");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_resampled(
                  handle, native_profile, profile_size, 1, 0) == false,
              "Mode 0 resampled load should fail");
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_resampled(
                  handle, native_profile, profile_size, 1, 5) == false,
              "Mode 5 resampled load should fail");

  // NULL source with a live instance fails
  TEST_ASSERT(specbleach_denoiser_load_noise_profile_for_mode(
                  handle, NULL, profile_size, 1, ROLLING_MEAN) == false,
              "NULL for_mode source should fail");

  // DSP-state reset rebuilds engines and keeps processing
  TEST_ASSERT(specbleach_denoiser_reset_dsp_state(NULL) == false,
              "NULL reset should fail");
  TEST_ASSERT(specbleach_denoiser_reset_dsp_state(handle) == true,
              "Reset should succeed");
  TEST_ASSERT(specbleach_denoiser_has_any_profile(handle) == false,
              "Reset should clear profiles");
  TEST_ASSERT(specbleach_denoiser_process(handle, 256, in_buf, out_buf) == true,
              "Process should work after reset");

  free(native_profile);
  free(half_profile);
  specbleach_denoiser_free(handle);
  printf("✓ Redesigned API coverage tests passed\n");
}

int main(void) {
  printf("Running specbleach denoiser tests...\n");

  test_specbleach_noise_profile_mode_functions();
  test_specbleach_load_noise_profile_with_mode();
  test_specbleach_load_noise_profile_for_mode();
  test_specbleach_mode_switching();
  test_specbleach_reset_noise_profile();
  test_specbleach_run_features();
  test_specbleach_silence_bypass();
  test_specbleach_smoothing_transition_and_validation();
  test_specbleach_adaptive_method_switch();
  test_specbleach_redesigned_api_coverage();

  // Getter Coverage (Extra) and NULL Safety
  printf("Testing API Getters and NULL safety for coverage...\n");
  specbleach_denoiser* h = specbleach_denoiser_initialize(44100, 20.0f);

  // Verify getters work with valid handle
  specbleach_denoiser_get_noise_profile_for_mode(h, ROLLING_MEAN);
  uint32_t prof_size = specbleach_denoiser_get_noise_profile_size(h);
  specbleach_denoiser_get_latency(h);
  specbleach_denoiser_noise_profile_available_for_mode(h, ROLLING_MEAN);

  // New PR #99 getters
  specbleach_denoiser_get_tonal_mask(h);
  specbleach_denoiser_get_active_noise_profile(h);
  float peak_freqs[10];
  specbleach_denoiser_get_tonal_peaks(h, peak_freqs, 10);

  // Process with tonal reduction and HPSS enabled
  float in_buf[1024] = {0};
  float out_buf[1024] = {0};
  SpecbleachDenoiserParameters t_params = {
      .learn_noise = SPECBLEACH_LEARN_OFF,
      .tonal_reduction_gain = 0.5f,
      .reduction_gain = 0.1f,
      .hpss_enable = true,
  };
  specbleach_denoiser_load_parameters(h, &t_params, sizeof(t_params));
  specbleach_denoiser_process(h, 1024, in_buf, out_buf);
  specbleach_denoiser_get_tonal_mask(h);
  // Feed non-zero audio with active processing
  for (int i = 0; i < 1024; ++i) {
    in_buf[i] = ((float)(i % 100) / 100.0f) * 0.5f;
  }
  t_params.smoothing_factor = 0.5f;
  t_params.whitening_factor = 0.5f;
  t_params.residual_listen = 1;
  t_params.masking_depth = 0.5f;
  t_params.suppression_strength = 0.5f;
  specbleach_denoiser_load_parameters(h, &t_params, sizeof(t_params));
  for (int f = 0; f < 10; ++f) {
    specbleach_denoiser_process(h, 1024, in_buf, out_buf);
  }

  // Feed a sharp transient/onset to trigger onset ducking
  float transient_buf[1024];
  for (int i = 0; i < 1024; ++i) {
    transient_buf[i] = (i % 2 == 0) ? 0.9f : -0.9f;
  }
  uint32_t cb_size = specbleach_denoiser_get_noise_profile_size(h);
  float* curve_bias = calloc(cb_size, sizeof(float));
  t_params.reduction_curve_enabled = true;
  t_params.reduction_curve_bias = curve_bias;
  t_params.reduction_curve_size = cb_size;
  TEST_ASSERT(
      specbleach_denoiser_load_parameters(h, &t_params, sizeof(t_params)),
      "Load with curve bias should succeed");
  specbleach_denoiser_process(h, 1024, transient_buf, out_buf);

  // Switch HPSS modes and verify latency
  t_params.residual_listen = 0;
  t_params.reduction_curve_enabled = false;
  t_params.reduction_curve_bias = NULL;
  t_params.reduction_curve_size = 0;
  free(curve_bias);

  t_params.hpss_enable = false;
  specbleach_denoiser_load_parameters(h, &t_params, sizeof(t_params));
  specbleach_denoiser_process(h, 1024, in_buf, out_buf);
  uint32_t lat0 = specbleach_denoiser_get_latency(h);

  t_params.hpss_enable = true;
  specbleach_denoiser_load_parameters(h, &t_params, sizeof(t_params));
  specbleach_denoiser_process(h, 1024, in_buf, out_buf);
  uint32_t lat1 = specbleach_denoiser_get_latency(h);
  TEST_ASSERT(lat1 == lat0, "Sliding HPSS introduces zero lookahead latency");

  // Verify NULL handle protections
  TEST_ASSERT(specbleach_denoiser_get_latency(NULL) == 0, "NULL latency");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_size(NULL) == 0,
              "NULL size");
  TEST_ASSERT(specbleach_denoiser_get_noise_profile_for_mode(
                  NULL, ROLLING_MEAN) == NULL,
              "NULL profile");
  TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(
                  NULL, ROLLING_MEAN) == false,
              "NULL available");
  specbleach_denoiser_reset_noise_profile(NULL); // must not crash
  TEST_ASSERT(true, "NULL reset");
  SpecbleachDenoiserParameters null_params = {0};
  TEST_ASSERT(specbleach_denoiser_load_parameters(NULL, &null_params,
                                                  sizeof(null_params)) == false,
              "NULL load");
  TEST_ASSERT(specbleach_denoiser_get_tonal_mask(NULL) == NULL,
              "NULL tonal mask");
  TEST_ASSERT(specbleach_denoiser_get_active_noise_profile(NULL) == NULL,
              "NULL active profile");
  TEST_ASSERT(specbleach_denoiser_get_tonal_peaks(NULL, peak_freqs, 10) == 0,
              "NULL tonal peaks");

  specbleach_denoiser_free(h);

  // Test adaptive noise profile persistence when deactivating adaptive mode
  printf("Testing adaptive profile persistence across mode toggles...\n");
  specbleach_denoiser* h_adapt = specbleach_denoiser_initialize(44100, 20.0f);
  TEST_ASSERT(h_adapt != NULL, "Initialization should succeed");

  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(h_adapt, 1) == false,
      "Profile should not be available initially");

  SpecbleachDenoiserParameters adapt_params = {
      .adaptive_noise = true,
      .noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE,
      .reduction_gain = 0.25f,
  };
  specbleach_denoiser_load_parameters(h_adapt, &adapt_params,
                                      sizeof(adapt_params));

  float adapt_in[1024];
  float adapt_out[1024];
  for (int i = 0; i < 1024; i++)
    adapt_in[i] = 0.1f;

  for (int b = 0; b < 10; b++) {
    specbleach_denoiser_process(h_adapt, 1024, adapt_in, adapt_out);
  }

  TEST_ASSERT(specbleach_denoiser_get_active_noise_profile(h_adapt) != NULL,
              "Active profile should exist during adaptive mode");
  TEST_ASSERT(
      specbleach_denoiser_noise_profile_available_for_mode(h_adapt, 1) == false,
      "Manual profile should remain unavailable in standalone adaptive mode");

  // Deactivate adaptive mode
  adapt_params.adaptive_noise = false;
  specbleach_denoiser_load_parameters(h_adapt, &adapt_params,
                                      sizeof(adapt_params));
  specbleach_denoiser_process(h_adapt, 1024, adapt_in, adapt_out);

  specbleach_denoiser_free(h_adapt);

  printf("✅ All specbleach denoiser tests passed!\n");
  return 0;
}
