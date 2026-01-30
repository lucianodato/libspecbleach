/*
 * Unit tests for the 2D Denoiser Wrapper (specbleach_2d_denoiser)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "specbleach_2d_denoiser.h"

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
  TEST_ASSERT(specbleach_2d_get_noise_profile_blocks_averaged(NULL) == 0,
              "NULL handle should return 0 blocks averaged");
  TEST_ASSERT(specbleach_2d_get_noise_profile(NULL) == NULL,
              "NULL handle should return NULL profile");

  TEST_ASSERT(specbleach_2d_load_noise_profile(NULL, NULL, 0, 0) == false,
              "NULL handle loading profile should fail");
  TEST_ASSERT(specbleach_2d_load_noise_profile(h, NULL, 0, 0) == false,
              "NULL profile data loading should fail");

  TEST_ASSERT(specbleach_2d_reset_noise_profile(NULL) == false,
              "NULL handle reset should fail");
  TEST_ASSERT(specbleach_2d_noise_profile_available(NULL) == false,
              "NULL handle available check should fail");

  specbleach_2d_free(h);
}

void test_noise_profile_api(void) {
  printf("Testing Noise Profile API...\n");

  SpectralBleachHandle h = specbleach_2d_initialize(SAMPLE_RATE, FRAME_SIZE);

  // Set default parameters first to ensure valid noise reduction mode
  SpectralBleach2DDenoiserParameters params = {
      .learn_noise = 0,
      .noise_reduction_mode = 1, // ROLLING_MEAN
      .reduction_amount = 20.0f,
      .smoothing_factor = 1.0f,
  };
  specbleach_2d_load_parameters(h, params);

  uint32_t profile_size = specbleach_2d_get_noise_profile_size(h);
  TEST_ASSERT(profile_size > 0, "Profile size should be valid");

  // Test loading a fake profile
  float* fake_profile = calloc(profile_size, sizeof(float));
  for (uint32_t i = 0; i < profile_size; i++) {
    fake_profile[i] = 1.0f;
  }

  TEST_ASSERT(
      specbleach_2d_load_noise_profile(h, fake_profile, profile_size, 10),
      "Loading valid profile should succeed");

  TEST_ASSERT(specbleach_2d_noise_profile_available(h),
              "Profile should be available");

  TEST_ASSERT(specbleach_2d_get_noise_profile_blocks_averaged(h) == 10,
              "Blocks averaged should match");

  // Test fetching profile
  float* retrieved = specbleach_2d_get_noise_profile(h);
  TEST_ASSERT(retrieved != NULL, "Should be able to retrieve profile");
  TEST_ASSERT(retrieved[0] == 1.0f, "Retrieved data should match");

  // Test reset
  specbleach_2d_reset_noise_profile(h);
  TEST_ASSERT(!specbleach_2d_noise_profile_available(h),
              "Profile should be unavailable after reset");

  // Test per-mode API
  TEST_ASSERT(specbleach_2d_load_noise_profile_for_mode(h, fake_profile,
                                                        profile_size, 5, 1),
              "Loading for mode 1 should succeed");

  TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(h, 1),
              "Mode 1 should be available");
  TEST_ASSERT(!specbleach_2d_noise_profile_available_for_mode(h, 2),
              "Mode 2 should be unavailable");

  TEST_ASSERT(
      specbleach_2d_get_noise_profile_blocks_averaged_for_mode(h, 1) == 5,
      "Mode 1 blocks should match");

  TEST_ASSERT(specbleach_2d_get_noise_profile_for_mode(h, 1) != NULL,
              "Mode 1 profile should be retrievable");
  TEST_ASSERT(specbleach_2d_get_noise_profile_for_mode(NULL, 1) == NULL,
              "NULL handle get mode profile should return NULL");

  // Test invalid load params
  TEST_ASSERT(specbleach_2d_load_noise_profile(h, fake_profile,
                                               profile_size + 1, 10) == false,
              "Size mismatch load should fail");
  TEST_ASSERT(specbleach_2d_load_noise_profile_for_mode(
                  h, fake_profile, profile_size, 10, 99) == false,
              "Invalid mode load should fail");

  free(fake_profile);
  specbleach_2d_free(h);
}

void test_2d_parameter_switching(void) {
  printf("Testing 2D parameter switching and adaptive methods...\n");
  SpectralBleachHandle h = specbleach_2d_initialize(SAMPLE_RATE, FRAME_SIZE);

  SpectralBleach2DDenoiserParameters params = {
      .learn_noise = 0,
      .noise_reduction_mode = 1,
      .reduction_amount = 20.0f,
      .smoothing_factor = 1.0f,
      .adaptive_noise = 1,
      .noise_estimation_method = 0 // Louizou
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
    params.noise_reduction_mode = mode;
    TEST_ASSERT(specbleach_2d_load_parameters(h, params),
                "Switch reduction mode should succeed");
  }

  uint32_t profile_size = specbleach_2d_get_noise_profile_size(h);
  float* input = calloc(1024, sizeof(float));
  float* output = calloc(1024, sizeof(float));
  float* profile = calloc(profile_size, sizeof(float));
  for (uint32_t i = 0; i < profile_size; i++)
    profile[i] = 0.001f;

  specbleach_2d_load_noise_profile(h, profile, profile_size, 1);

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
      .noise_reduction_mode = 1,
      .reduction_amount = 20.0f,
      .smoothing_factor = 1.0f,
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

  TEST_ASSERT(specbleach_2d_noise_profile_available(h),
              "Profile should be available after learning");

  specbleach_2d_free(h);
}

int main(void) {
  printf("Running specbleach_2d_denoiser wrapper tests...\n");

  test_initialization_failure();
  test_null_handling();
  test_noise_profile_api();
  test_2d_parameter_switching();
  test_process_loop();

  printf("✅ All specbleach_2d_denoiser tests passed!\n");
  return 0;
}
