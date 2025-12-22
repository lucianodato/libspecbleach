/*
 * Unit tests for Noise Profile functionality
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/noise_estimation/noise_profile.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_noise_profile_lifecycle(void) {
  printf("Testing Noise Profile lifecycle...\n");

  uint32_t profile_size = 513;
  NoiseProfile* np = noise_profile_initialize(profile_size);
  TEST_ASSERT(np != NULL, "Noise profile initialization should succeed");

  TEST_ASSERT(get_noise_profile_size(np) == profile_size,
              "Profile size should match");

  noise_profile_free(np);
  printf("✓ Noise Profile lifecycle tests passed\n");
}

void test_noise_profile_operations(void) {
  printf("Testing Noise Profile operations...\n");

  uint32_t profile_size = 513;
  NoiseProfile* np = noise_profile_initialize(profile_size);
  TEST_ASSERT(np != NULL, "Noise profile initialization should succeed");

  // Test invalid mode handling
  TEST_ASSERT(is_noise_estimation_available(np, 0) == false,
              "Mode 0 should be invalid");
  TEST_ASSERT(is_noise_estimation_available(np, 4) == false,
              "Mode 4 should be invalid");
  TEST_ASSERT(get_noise_profile(np, 0) == NULL, "Mode 0 should return NULL");
  TEST_ASSERT(get_noise_profile(np, 4) == NULL, "Mode 4 should return NULL");
  TEST_ASSERT(get_noise_profile_blocks_averaged(np, 0) == 0,
              "Mode 0 should return 0 blocks");

  // Initially no profiles should be available
  for (int mode = 1; mode <= 3; mode++) {
    TEST_ASSERT(is_noise_estimation_available(np, mode) == false,
                "Profile should not be available initially");
    TEST_ASSERT(get_noise_profile_blocks_averaged(np, mode) == 0,
                "Should have 0 blocks averaged initially");
  }

  noise_profile_free(np);
  printf("✓ Noise Profile operations tests passed\n");
}

void test_noise_profile_set_and_get(void) {
  printf("Testing Noise Profile set and get...\n");

  uint32_t profile_size = 513;
  NoiseProfile* np = noise_profile_initialize(profile_size);
  TEST_ASSERT(np != NULL, "Noise profile initialization should succeed");

  // Create test profile data
  float test_profile[513];
  for (int i = 0; i < 513; i++) {
    test_profile[i] = (float)i * 0.01f;
  }

  // Test setting profile for mode 1
  TEST_ASSERT(set_noise_profile(np, 1, test_profile, profile_size, 10) == true,
              "Setting profile for mode 1 should succeed");

  TEST_ASSERT(is_noise_estimation_available(np, 1) == true,
              "Profile should be available after setting");
  TEST_ASSERT(get_noise_profile_blocks_averaged(np, 1) == 10,
              "Should have 10 blocks averaged");

  // Test getting the profile back
  float* retrieved_profile = get_noise_profile(np, 1);
  TEST_ASSERT(retrieved_profile != NULL, "Should get valid profile pointer");

  for (int i = 0; i < 513; i++) {
    TEST_FLOAT_CLOSE(retrieved_profile[i], test_profile[i], 0.001f);
  }

  // Test that other modes are still unavailable
  TEST_ASSERT(is_noise_estimation_available(np, 2) == false,
              "Mode 2 should still be unavailable");
  TEST_ASSERT(is_noise_estimation_available(np, 3) == false,
              "Mode 3 should still be unavailable");

  noise_profile_free(np);
  printf("✓ Noise Profile set and get tests passed\n");
}

void test_noise_profile_increment_blocks(void) {
  printf("Testing Noise Profile increment blocks...\n");

  uint32_t profile_size = 513;
  NoiseProfile* np = noise_profile_initialize(profile_size);
  TEST_ASSERT(np != NULL, "Noise profile initialization should succeed");

  // Test incrementing blocks for mode 1
  for (int i = 1; i <= 6; i++) {
    TEST_ASSERT(increment_blocks_averaged(np, 1) == true,
                "Increment should succeed");
    TEST_ASSERT(get_noise_profile_blocks_averaged(np, 1) == (uint32_t)i,
                "Blocks should match expected count");
  }

  // Profile should become available after more than
  // MIN_NUMBER_OF_WINDOWS_NOISE_AVERAGED blocks
  TEST_ASSERT(is_noise_estimation_available(np, 1) == true,
              "Profile should become available after enough blocks");

  // Test invalid mode
  TEST_ASSERT(increment_blocks_averaged(np, 0) == false,
              "Increment with invalid mode should fail");

  noise_profile_free(np);
  printf("✓ Noise Profile increment blocks tests passed\n");
}

void test_noise_profile_reset(void) {
  printf("Testing Noise Profile reset...\n");

  uint32_t profile_size = 513;
  NoiseProfile* np = noise_profile_initialize(profile_size);
  TEST_ASSERT(np != NULL, "Noise profile initialization should succeed");

  // Set up some profiles
  float test_profile[513] = {1.0f};
  TEST_ASSERT(set_noise_profile(np, 1, test_profile, profile_size, 10) == true,
              "Setting profile should succeed");
  TEST_ASSERT(set_noise_profile(np, 2, test_profile, profile_size, 20) == true,
              "Setting profile should succeed");

  // Verify they're set
  TEST_ASSERT(is_noise_estimation_available(np, 1) == true,
              "Profile 1 should be available");
  TEST_ASSERT(is_noise_estimation_available(np, 2) == true,
              "Profile 2 should be available");

  // Reset
  TEST_ASSERT(reset_noise_profile(np) == true, "Reset should succeed");

  // Verify all profiles are reset
  for (int mode = 1; mode <= 3; mode++) {
    TEST_ASSERT(is_noise_estimation_available(np, mode) == false,
                "Profile should not be available after reset");
    TEST_ASSERT(get_noise_profile_blocks_averaged(np, mode) == 0,
                "Blocks should be 0 after reset");
  }

  noise_profile_free(np);
  printf("✓ Noise Profile reset tests passed\n");
}

void test_noise_profile_multiple_modes(void) {
  printf("Testing Noise Profile multiple modes...\n");

  uint32_t profile_size = 513;
  NoiseProfile* np = noise_profile_initialize(profile_size);
  TEST_ASSERT(np != NULL, "Noise profile initialization should succeed");

  // Create different profiles for each mode
  float profile1[513], profile2[513], profile3[513];
  for (int i = 0; i < 513; i++) {
    profile1[i] = 1.0f;
    profile2[i] = 2.0f;
    profile3[i] = 3.0f;
  }

  // Set profiles for all modes
  TEST_ASSERT(set_noise_profile(np, 1, profile1, profile_size, 5) == true,
              "Setting mode 1 should succeed");
  TEST_ASSERT(set_noise_profile(np, 2, profile2, profile_size, 10) == true,
              "Setting mode 2 should succeed");
  TEST_ASSERT(set_noise_profile(np, 3, profile3, profile_size, 15) == true,
              "Setting mode 3 should succeed");

  // Verify each mode has correct data
  float* retrieved1 = get_noise_profile(np, 1);
  float* retrieved2 = get_noise_profile(np, 2);
  float* retrieved3 = get_noise_profile(np, 3);

  TEST_ASSERT(retrieved1 != NULL && retrieved2 != NULL && retrieved3 != NULL,
              "All profiles should be retrievable");

  for (int i = 0; i < 513; i++) {
    TEST_FLOAT_CLOSE(retrieved1[i], 1.0f, 0.001f);
    TEST_FLOAT_CLOSE(retrieved2[i], 2.0f, 0.001f);
    TEST_FLOAT_CLOSE(retrieved3[i], 3.0f, 0.001f);
  }

  // Verify block counts
  TEST_ASSERT(get_noise_profile_blocks_averaged(np, 1) == 5,
              "Mode 1 should have 5 blocks");
  TEST_ASSERT(get_noise_profile_blocks_averaged(np, 2) == 10,
              "Mode 2 should have 10 blocks");
  TEST_ASSERT(get_noise_profile_blocks_averaged(np, 3) == 15,
              "Mode 3 should have 15 blocks");

  noise_profile_free(np);
  printf("✓ Noise Profile multiple modes tests passed\n");
}

int main(void) {
  printf("Running noise profile tests...\n");

  test_noise_profile_lifecycle();
  test_noise_profile_operations();
  test_noise_profile_set_and_get();
  test_noise_profile_increment_blocks();
  test_noise_profile_reset();
  test_noise_profile_multiple_modes();

  printf("✅ All noise profile tests passed!\n");
  return 0;
}
