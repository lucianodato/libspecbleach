/*
 * Unit tests for SpectralBleach Denoiser API functions
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "specbleach_denoiser.h"

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

  // Load parameters with mode 2
  SpectralBleachDenoiserParameters params = {
      .learn_noise = true,
      .noise_reduction_mode = 2, // Use mode 2
      .residual_listen = false,
      .noise_scaling_type = 0,
      .reduction_amount = -20.0f,
      .noise_rescale = 1.0f,
      .smoothing_factor = 50.0f,
      .whitening_factor = 0.0f,
      .post_filter_threshold = -30.0f,
  };

  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading parameters should succeed");

  // Test NULL handle cases
  TEST_ASSERT(
      specbleach_get_noise_profile_blocks_averaged_for_mode(NULL, 1) == 0,
      "NULL handle should return 0");
  TEST_ASSERT(specbleach_get_noise_profile_for_mode(NULL, 1) == NULL,
              "NULL handle should return NULL");
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(NULL, 1) == false,
              "NULL handle should return false");

  // Test invalid modes
  TEST_ASSERT(
      specbleach_get_noise_profile_blocks_averaged_for_mode(handle, 0) == 0,
      "Invalid mode should return 0");
  TEST_ASSERT(
      specbleach_get_noise_profile_blocks_averaged_for_mode(handle, 4) == 0,
      "Invalid mode should return 0");
  TEST_ASSERT(specbleach_get_noise_profile_for_mode(handle, 0) == NULL,
              "Invalid mode should return NULL");
  TEST_ASSERT(specbleach_get_noise_profile_for_mode(handle, 4) == NULL,
              "Invalid mode should return NULL");
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 0) == false,
              "Invalid mode should return false");
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 4) == false,
              "Invalid mode should return false");

  // Initially no profiles should be available
  for (int mode = 1; mode <= 3; mode++) {
    TEST_ASSERT(
        specbleach_noise_profile_available_for_mode(handle, mode) == false,
        "Profile should not be available initially");
    TEST_ASSERT(specbleach_get_noise_profile_blocks_averaged_for_mode(
                    handle, mode) == 0,
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

  // Load parameters
  SpectralBleachDenoiserParameters params = {
      .learn_noise = true,
      .noise_reduction_mode = 1,
      .residual_listen = false,
      .noise_scaling_type = 0,
      .reduction_amount = -20.0f,
      .noise_rescale = 1.0f,
      .smoothing_factor = 50.0f,
      .whitening_factor = 0.0f,
      .post_filter_threshold = -30.0f,
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
    test_profile[i] = 0.1f + (float)i * 0.001f;
  }

  // Test loading profile with mode
  TEST_ASSERT(specbleach_load_noise_profile(handle, test_profile, profile_size,
                                            10) == true,
              "Loading noise profile should succeed");

  // Verify profile is available and has correct block count
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 1) == true,
              "Profile should be available for mode 1");
  TEST_ASSERT(
      specbleach_get_noise_profile_blocks_averaged_for_mode(handle, 1) == 10,
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
      .noise_reduction_mode = 1,
      .residual_listen = false,
      .noise_scaling_type = 0,
      .reduction_amount = -20.0f,
      .noise_rescale = 1.0f,
      .smoothing_factor = 50.0f,
      .whitening_factor = 0.0f,
      .post_filter_threshold = -30.0f,
  };
  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading mode 1 parameters should succeed");
  TEST_ASSERT(
      specbleach_load_noise_profile(handle, profile1, profile_size, 5) == true,
      "Loading mode 1 profile should succeed");

  // Mode 2
  params.noise_reduction_mode = 2;
  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading mode 2 parameters should succeed");
  TEST_ASSERT(
      specbleach_load_noise_profile(handle, profile2, profile_size, 10) == true,
      "Loading mode 2 profile should succeed");

  // Mode 3
  params.noise_reduction_mode = 3;
  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading mode 3 parameters should succeed");
  TEST_ASSERT(
      specbleach_load_noise_profile(handle, profile3, profile_size, 15) == true,
      "Loading mode 3 profile should succeed");

  // Verify all modes have their profiles
  for (int mode = 1; mode <= 3; mode++) {
    TEST_ASSERT(
        specbleach_noise_profile_available_for_mode(handle, mode) == true,
        "Profile should be available");
  }

  // Check block counts
  TEST_ASSERT(
      specbleach_get_noise_profile_blocks_averaged_for_mode(handle, 1) == 5,
      "Mode 1 should have 5 blocks");
  TEST_ASSERT(
      specbleach_get_noise_profile_blocks_averaged_for_mode(handle, 2) == 10,
      "Mode 2 should have 10 blocks");
  TEST_ASSERT(
      specbleach_get_noise_profile_blocks_averaged_for_mode(handle, 3) == 15,
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
      .noise_reduction_mode = 1,
      .residual_listen = false,
      .noise_scaling_type = 0,
      .reduction_amount = -20.0f,
      .noise_rescale = 1.0f,
      .smoothing_factor = 50.0f,
      .whitening_factor = 0.0f,
      .post_filter_threshold = -30.0f,
  };

  TEST_ASSERT(specbleach_load_parameters(handle, params) == true,
              "Loading parameters should succeed");
  TEST_ASSERT(specbleach_load_noise_profile(handle, test_profile, profile_size,
                                            10) == true,
              "Loading profile should succeed");

  // Verify profile is available
  TEST_ASSERT(specbleach_noise_profile_available_for_mode(handle, 1) == true,
              "Profile should be available before reset");

  // Reset profile
  TEST_ASSERT(specbleach_reset_noise_profile(handle) == true,
              "Reset should succeed");

  // Verify all profiles are reset
  for (int mode = 1; mode <= 3; mode++) {
    TEST_ASSERT(
        specbleach_noise_profile_available_for_mode(handle, mode) == false,
        "Profile should not be available after reset");
    TEST_ASSERT(specbleach_get_noise_profile_blocks_averaged_for_mode(
                    handle, mode) == 0,
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
    test_profile[i] = 0.1f + (float)i * 0.001f;
  }

  // Test loading profile for each mode directly
  for (int mode = 1; mode <= 3; mode++) {
    TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                    handle, test_profile, profile_size, 10, mode) == true,
                "Loading noise profile for mode should succeed");
  }

  // Test invalid mode
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, 0) == false,
              "Loading noise profile for invalid mode should fail");
  TEST_ASSERT(specbleach_load_noise_profile_for_mode(
                  handle, test_profile, profile_size, 10, 4) == false,
              "Loading noise profile for invalid mode should fail");

  free(test_profile);
  specbleach_free(handle);

  printf("✓ Specbleach load noise profile for mode tests passed\n");
}

void test_specbleach_run_features(void) {
  printf("Testing specbleach denoiser features (run)...\n");

  SpectralBleachHandle handle = specbleach_initialize(44100, 20.0f);
  TEST_ASSERT(handle != NULL, "Initialization should succeed");

  uint32_t profile_size = specbleach_get_noise_profile_size(handle);
  float* input = (float*)calloc(1024, sizeof(float));
  float* output = (float*)calloc(1024, sizeof(float));
  float* profile = (float*)calloc(profile_size, sizeof(float));

  for (int i = 0; i < 1024; i++)
    input[i] = 0.5f;
  for (uint32_t i = 0; i < profile_size; i++)
    profile[i] = 0.1f;

  specbleach_load_noise_profile(handle, profile, profile_size, 1);

  SpectralBleachDenoiserParameters params = {
      .learn_noise = false,
      .noise_reduction_mode = 1,
      .residual_listen = false,
      .noise_scaling_type = 0,
      .reduction_amount = 20.0f,
      .noise_rescale = 1.0f,
      .smoothing_factor = 50.0f,
      .whitening_factor = 1.0f, // Test whitening
      .post_filter_threshold = -30.0f,
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

  printf("✅ All specbleach denoiser tests passed!\n");
  return 0;
}
