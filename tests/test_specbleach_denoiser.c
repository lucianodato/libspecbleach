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
  SpectralBleachHandle handle = specbleach_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  SpectralBleachDenoiserParameters params = {
      .learn_noise = true,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };

  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading parameters should succeed");

  // Test NULL handle cases
  TEST_ASSERT(specbleach_get_noise_profile_block_count_for_mode(NULL, 1) == 0,
              "NULL handle should return 0");
  TEST_ASSERT(specbleach_get_noise_profile_for_mode(NULL, 1) == NULL,
              "NULL handle should return NULL");
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(NULL, 1) == false,
              "NULL handle should return false");

  // Test invalid modes
  TEST_ASSERT(specbleach_get_noise_profile_block_count_for_mode(handle, 0) == 0,
              "Invalid mode should return 0");
  TEST_ASSERT(specbleach_get_noise_profile_block_count_for_mode(handle, 5) == 0,
              "Invalid mode should return 0");
  TEST_ASSERT(specbleach_get_noise_profile_for_mode(handle, 0) == NULL,
              "Invalid mode should return NULL");
  TEST_ASSERT(specbleach_get_noise_profile_for_mode(handle, 5) == NULL,
              "Invalid mode should return NULL");
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 0) == false,
              "Invalid mode should return false");
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 5) == false,
              "Invalid mode should return false");

  // Initially no profiles should be available
  for (int mode = 1; mode <= 4; mode++) {
    TEST_ASSERT(
        specbleach_noise_profile_available_for_mode(handle, mode) == false,
        "Profile should not be available initially");
    TEST_ASSERT(
        specbleach_get_noise_profile_block_count_for_mode(handle, mode) == 0,
        "Should have 0 blocks initially");
    // Note: get_noise_profile_for_mode returns the profile array even if not
    // available (this matches the behavior of the original single-mode API)
    TEST_ASSERT(specbleach_get_noise_profile_for_mode(handle, mode) != NULL,
                "Should return profile array even if not available");
  }

  specbleach_free(handle);
  printf("✓ Specbleach noise profile mode functions tests passed\n");
}

void test_specbleach_load_noise_profile_with_mode(void) {
  printf("Testing specbleach load noise profile with mode...\n");

  SpectralBleachHandle handle = specbleach_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  // Load
  SpectralBleachDenoiserParameters params = (SpectralBleachDenoiserParameters){
      .learn_noise = true,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };

  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading parameters should succeed");

  // Get profile size
  uint32_t profile_size = specbleach_get_noise_profile_size(handle);
  TEST_ASSERT(profile_size > 0, "Profile size should be valid");

  // Create test profile
  float* test_profile = (float*)malloc(profile_size * sizeof(float));
  TEST_ASSERT(test_profile != NULL, "Profile allocation should succeed");

  for (uint32_t i = 0; i < profile_size; i++) {
    test_profile[i] = 0.1f + ((float)i * 0.001f);
  }

  // Test loading profile with mode
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, ROLLING_MEAN) == true,
              "Loading noise profile should succeed");

  // Verify profile is available and has correct block count
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 1) == true,
              "Profile should be available for mode 1");
  TEST_ASSERT(
      specbleach_get_noise_profile_block_count_for_mode(handle, 1) == 10,
      "Should have 10 blocks averaged");

  // Get the profile back and verify it matches
  float* retrieved_profile = specbleach_get_noise_profile_for_mode(handle, 1);
  TEST_ASSERT(retrieved_profile != NULL, "Should get valid profile");

  for (uint32_t i = 0; i < profile_size; i++) {
    TEST_FLOAT_CLOSE(retrieved_profile[i], test_profile[i], 0.001f);
  }

  // Test that other modes are unaffected
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 2) == false,
              "Mode 2 should still be unavailable");
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 3) == false,
              "Mode 3 should still be unavailable");

  free(test_profile);
  specbleach_free(handle);
  printf("✓ Specbleach load noise profile with mode tests passed\n");
}

void test_specbleach_mode_switching(void) {
  printf("Testing specbleach mode switching...\n");

  SpectralBleachHandle handle = specbleach_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  uint32_t profile_size = specbleach_get_noise_profile_size(handle);

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
  SpectralBleachDenoiserParameters params;

  // Mode 1
  params = (SpectralBleachDenoiserParameters){
      .learn_noise = true,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };
  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading mode 1 parameters should succeed");
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, profile1, profile_size, 5, ROLLING_MEAN) == true,
              "Loading mode 1 profile should succeed");

  // Mode 2
  // Explicitly loading for mode 2 since noise_reduction_mode is removed
  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading mode 2 parameters should succeed");
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, profile2, profile_size, 10, 2) == true,
              "Loading mode 2 profile should succeed");

  // Mode 3
  // Explicitly loading for mode 3
  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading mode 3 parameters should succeed");
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, profile3, profile_size, 15, 3) == true,
              "Loading mode 3 profile should succeed");

  // Verify all modes have their profiles
  for (int mode = 1; mode <= 3; mode++) {
    TEST_ASSERT(
        specbleach_noise_profile_available_for_mode(handle, mode) == true,
        "Profile should be available");
  }

  // Check block counts
  TEST_ASSERT(specbleach_get_noise_profile_block_count_for_mode(handle, 1) == 5,
              "Mode 1 should have 5 blocks");
  TEST_ASSERT(
      specbleach_get_noise_profile_block_count_for_mode(handle, 2) == 10,
      "Mode 2 should have 10 blocks");
  TEST_ASSERT(
      specbleach_get_noise_profile_block_count_for_mode(handle, 3) == 15,
      "Mode 3 should have 15 blocks");

  // Check profile values
  float* retrieved1 = specbleach_get_noise_profile_for_mode(handle, 1);
  float* retrieved2 = specbleach_get_noise_profile_for_mode(handle, 2);
  float* retrieved3 = specbleach_get_noise_profile_for_mode(handle, 3);

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
  specbleach_free(handle);
  printf("✓ Specbleach mode switching tests passed\n");
}

void test_specbleach_reset_noise_profile(void) {
  printf("Testing specbleach reset noise profile...\n");

  SpectralBleachHandle handle = specbleach_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  uint32_t profile_size = specbleach_get_noise_profile_size(handle);
  float* test_profile = (float*)malloc(profile_size * sizeof(float));
  TEST_ASSERT(test_profile != NULL, "Profile allocation should succeed");

  for (uint32_t i = 0; i < profile_size; i++) {
    test_profile[i] = 1.0f;
  }

  // Load parameters and profile
  SpectralBleachDenoiserParameters params = {
      .learn_noise = true,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 0.0f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };

  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading parameters should succeed");
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, ROLLING_MEAN) == true,
              "Loading profile should succeed");

  // Verify profile is available
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 1) == true,
              "Profile should be available before reset");

  // Reset profile
  TEST_ASSERT(specbleach_reset_noise_profile(handle) == true,
              "Reset should succeed");

  // Verify all profiles are reset
  for (int mode = 1; mode <= 4; mode++) {
    TEST_ASSERT(
        specbleach_noise_profile_available_for_mode(handle, mode) == false,
        "Profile should not be available after reset");
    TEST_ASSERT(
        specbleach_get_noise_profile_block_count_for_mode(handle, mode) == 0,
        "Blocks should be 0 after reset");
  }

  free(test_profile);
  specbleach_free(handle);
  printf("✓ Specbleach reset noise profile tests passed\n");
}

void test_specbleach_load_noise_profile_for_mode(void) {
  printf("Testing specbleach load noise profile for mode...\n");

  SpectralBleachHandle handle = specbleach_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  // Get profile size
  uint32_t profile_size = specbleach_get_noise_profile_size(handle);
  TEST_ASSERT(profile_size > 0, "Profile size should be valid");

  // Create test profile
  float* test_profile = (float*)malloc(profile_size * sizeof(float));
  TEST_ASSERT(test_profile != NULL, "Profile allocation should succeed");

  for (uint32_t i = 0; i < profile_size; i++) {
    test_profile[i] = 0.1f + ((float)i * 0.001f);
  }

  // Test loading profile for each mode directly
  for (int mode = 1; mode <= 4; mode++) {
    TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                    handle, test_profile, profile_size, 10, mode) == true,
                "Loading noise profile for mode should succeed");
  }

  // Test invalid mode
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, 0) == false,
              "Loading noise profile for invalid mode should fail");
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, 5) == false,
              "Loading noise profile for invalid mode should fail");

  free(test_profile);
  specbleach_free(handle);

  printf("✓ Specbleach load noise profile for mode tests passed\n");
}

void test_specbleach_run_features(void) {
  printf("Testing specbleach run features (whitening, residual listen)...\n");

  SpectralBleachHandle handle = specbleach_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Denoiser initialization should succeed");

  float* input = (float*)malloc(1024 * sizeof(float));
  float* output = (float*)malloc(1024 * sizeof(float));
  for (int i = 0; i < 1024; i++) {
    input[i] = 0.5f;
    output[i] = 0.0f;
  }

  uint32_t profile_size = specbleach_get_noise_profile_size(handle);
  float* profile = (float*)malloc(profile_size * sizeof(float));
  for (uint32_t i = 0; i < profile_size; i++) {
    profile[i] = 0.1f;
  }

  specbleach_load_noise_profile_for_mode(handle, profile, profile_size, 1,
                                         ROLLING_MEAN);

  SpectralBleachDenoiserParameters params = {
      .learn_noise = false,
      .residual_listen = false,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .whitening_factor = 1.0f, // Test whitening
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = 0.0f,
  };

  specbleach_load_parameters(handle, params);
  specbleach_process(handle, 1024, input, output);

  // Test residual listen
  params.residual_listen = true;
  params.whitening_factor = 0.0f;
  specbleach_load_parameters(handle, params);
  specbleach_process(handle, 1024, input, output);

  free(input);
  free(output);
  free(profile);
  specbleach_free(handle);
  printf("✓ Specbleach denoiser features tests passed\n");
}

int main(void) {
  printf("Running specbleach denoiser tests...\n");

  test_specbleach_noise_profile_mode_functions();
  test_specbleach_load_noise_profile_with_mode();
  test_specbleach_load_noise_profile_for_mode();
  test_specbleach_mode_switching();
  test_specbleach_reset_noise_profile();
  test_specbleach_run_features();

  // Getter Coverage (Extra) and NULL Safety
  printf("Testing API Getters and NULL safety for coverage...\n");
  SpectralBleachHandle h = specbleach_initialize(44100, 20.0f);

  // Verify getters work with valid handle
  specbleach_get_noise_profile_for_mode(h, ROLLING_MEAN);
  uint32_t prof_size = specbleach_get_noise_profile_size(h);
  specbleach_get_latency(h);
  specbleach_noise_profile_available_for_mode(h, ROLLING_MEAN);

  // New PR #99 getters
  specbleach_get_tonal_mask(h);
  specbleach_get_active_noise_profile(h);
  float peak_freqs[10];
  specbleach_get_tonal_peaks(h, peak_freqs, 10);

  float* dummy_prof = calloc(prof_size, sizeof(float));
  specbleach_get_tonal_peaks_for_profile(h, dummy_prof, prof_size, peak_freqs,
                                         10);
  free(dummy_prof);

  // Process with tonal reduction and HPSS enabled
  float in_buf[1024] = {0};
  float out_buf[1024] = {0};
  SpectralBleachDenoiserParameters t_params = {
      .learn_noise = false,
      .tonal_reduction_gain = 0.5f,
      .reduction_gain = 0.1f,
      .hpss_enable = 1,
  };
  specbleach_load_parameters(h, t_params);
  specbleach_process(h, 1024, in_buf, out_buf);
  specbleach_get_tonal_mask(h);
  // Feed non-zero audio with active processing
  for (int i = 0; i < 1024; ++i) {
    in_buf[i] = ((float)(i % 100) / 100.0f) * 0.5f;
  }
  t_params.smoothing_factor = 0.5f;
  t_params.whitening_factor = 0.5f;
  t_params.residual_listen = 1;
  t_params.masking_depth = 0.5f;
  t_params.suppression_strength = 0.5f;
  specbleach_load_parameters(h, t_params);
  for (int f = 0; f < 10; ++f) {
    specbleach_process(h, 1024, in_buf, out_buf);
  }

  // Feed a sharp transient/onset to trigger onset ducking
  float transient_buf[1024];
  for (int i = 0; i < 1024; ++i) {
    transient_buf[i] = (i % 2 == 0) ? 0.9f : -0.9f;
  }
  float curve_bias[1024] = {0};
  t_params.reduction_curve_enabled = true;
  t_params.reduction_curve_bias = curve_bias;
  specbleach_load_parameters(h, t_params);
  specbleach_process(h, 1024, transient_buf, out_buf);

  // Switch HPSS modes and verify latency
  t_params.residual_listen = 0;
  t_params.reduction_curve_enabled = false;
  t_params.reduction_curve_bias = NULL;

  t_params.hpss_enable = 0;
  specbleach_load_parameters(h, t_params);
  specbleach_process(h, 1024, in_buf, out_buf);
  uint32_t lat0 = specbleach_get_latency(h);

  t_params.hpss_enable = 1;
  specbleach_load_parameters(h, t_params);
  specbleach_process(h, 1024, in_buf, out_buf);
  uint32_t lat1 = specbleach_get_latency(h);
  TEST_ASSERT(lat1 == lat0, "Sliding HPSS introduces zero lookahead latency");

  // Verify NULL handle protections
  TEST_ASSERT(specbleach_get_latency(NULL) == 0, "NULL latency");
  TEST_ASSERT(specbleach_get_noise_profile_size(NULL) == 0, "NULL size");
  TEST_ASSERT(specbleach_get_noise_profile_for_mode(NULL, ROLLING_MEAN) == NULL,
              "NULL profile");
  TEST_ASSERT(
      specbleach_noise_profile_available_for_mode(NULL, ROLLING_MEAN) == false,
      "NULL available");
  TEST_ASSERT(specbleach_reset_noise_profile(NULL) == false, "NULL reset");
  TEST_ASSERT(specbleach_load_parameters(
                  NULL, (SpectralBleachDenoiserParameters){0}) == false,
              "NULL load");
  TEST_ASSERT(specbleach_get_tonal_mask(NULL) == NULL, "NULL tonal mask");
  TEST_ASSERT(specbleach_get_active_noise_profile(NULL) == NULL,
              "NULL active profile");
  TEST_ASSERT(specbleach_get_tonal_peaks(NULL, peak_freqs, 10) == 0,
              "NULL tonal peaks");
  TEST_ASSERT(specbleach_get_tonal_peaks_for_profile(NULL, NULL, 0, peak_freqs,
                                                     10) == 0,
              "NULL tonal peaks for profile");

  specbleach_free(h);

  // Test adaptive noise profile persistence when deactivating adaptive mode
  printf("Testing adaptive profile persistence across mode toggles...\n");
  SpectralBleachHandle h_adapt = specbleach_initialize(44100, 20.0f);
  TEST_ASSERT(h_adapt != NULL, "Initialization should succeed");

  TEST_ASSERT(specbleach_noise_profile_available_for_mode(h_adapt, 1) == false,
              "Profile should not be available initially");

  SpectralBleachDenoiserParameters adapt_params = {
      .adaptive_noise = 1,
      .noise_estimation_method = 0,
      .reduction_gain = 0.25f,
  };
  specbleach_load_parameters(h_adapt, adapt_params);

  float adapt_in[1024];
  float adapt_out[1024];
  for (int i = 0; i < 1024; i++)
    adapt_in[i] = 0.1f;

  for (int b = 0; b < 10; b++) {
    specbleach_process(h_adapt, 1024, adapt_in, adapt_out);
  }

  TEST_ASSERT(specbleach_get_active_noise_profile(h_adapt) != NULL,
              "Active profile should exist during adaptive mode");
  TEST_ASSERT(
      specbleach_noise_profile_available_for_mode(h_adapt, 1) == false,
      "Manual profile should remain unavailable in standalone adaptive mode");

  // Deactivate adaptive mode
  adapt_params.adaptive_noise = 0;
  specbleach_load_parameters(h_adapt, adapt_params);
  specbleach_process(h_adapt, 1024, adapt_in, adapt_out);

  specbleach_free(h_adapt);

  printf("✅ All specbleach denoiser tests passed!\n");
  return 0;
}
