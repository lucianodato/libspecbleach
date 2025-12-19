/*
 * Unit tests for utility functions
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Include internal headers for testing
#include "../src/shared/configurations.h"
#include "../src/shared/utils/general_utils.h"
#include "../src/shared/utils/spectral_utils.h"

// Function prototypes
void test_sanitize_denormal(void);
void test_from_db_to_coefficient(void);
void test_get_next_divisible_two(void);
void test_get_next_power_two(void);
void test_remap_percentage_log_like_unity(void);

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_sanitize_denormal() {
  printf("Testing sanitize_denormal...\n");

  // Test normal values
  TEST_FLOAT_CLOSE(sanitize_denormal(1.0f), 1.0f, 1e-6f);
  TEST_FLOAT_CLOSE(sanitize_denormal(-1.0f), -1.0f, 1e-6f);
  TEST_FLOAT_CLOSE(sanitize_denormal(0.5f), 0.5f, 1e-6f);

  // Test zero
  TEST_FLOAT_CLOSE(sanitize_denormal(0.0f), 0.0f, 1e-6f);

  // Test very small values (denormals)
  float denormal = 1e-40f; // This should be treated as denormal
  TEST_FLOAT_CLOSE(sanitize_denormal(denormal), 0.0f, 1e-6f);

  printf("✓ sanitize_denormal tests passed\n");
}

void test_from_db_to_coefficient() {
  printf("Testing from_db_to_coefficient...\n");

  // Test 0 dB = 1.0 coefficient
  TEST_FLOAT_CLOSE(from_db_to_coefficient(0.0f), 1.0f, 1e-6f);

  // Test -6 dB ≈ 0.501 coefficient
  TEST_FLOAT_CLOSE(from_db_to_coefficient(-6.0f), 0.501f, 0.01f);

  // Test -12 dB ≈ 0.251 coefficient
  TEST_FLOAT_CLOSE(from_db_to_coefficient(-12.0f), 0.251f, 0.01f);

  // Test positive dB
  TEST_ASSERT(from_db_to_coefficient(6.0f) > 1.5f,
              "Positive dB should give coefficient > 1.5");

  printf("✓ from_db_to_coefficient tests passed\n");
}

void test_get_next_divisible_two() {
  printf("Testing get_next_divisible_two...\n");

  TEST_ASSERT(get_next_divisible_two(1) == 2, "1 -> 2");
  TEST_ASSERT(get_next_divisible_two(2) == 2, "2 -> 2");
  TEST_ASSERT(get_next_divisible_two(3) == 4, "3 -> 4");
  TEST_ASSERT(get_next_divisible_two(4) == 4, "4 -> 4");
  TEST_ASSERT(get_next_divisible_two(5) == 6, "5 -> 6");
  TEST_ASSERT(get_next_divisible_two(7) == 8, "7 -> 8");
  TEST_ASSERT(get_next_divisible_two(8) == 8, "8 -> 8");
  TEST_ASSERT(get_next_divisible_two(9) == 10, "9 -> 10");

  printf("✓ get_next_divisible_two tests passed\n");
}

void test_get_next_power_two() {
  printf("Testing get_next_power_two...\n");

  TEST_ASSERT(get_next_power_two(1) == 1, "1 -> 1 (2^0)");
  TEST_ASSERT(get_next_power_two(2) == 2, "2 -> 2 (2^1)");
  TEST_ASSERT(get_next_power_two(3) == 4, "3 -> 4 (2^2)");
  TEST_ASSERT(get_next_power_two(4) == 4, "4 -> 4 (2^2)");
  TEST_ASSERT(get_next_power_two(5) == 8, "5 -> 8 (2^3)");
  TEST_ASSERT(get_next_power_two(7) == 8, "7 -> 8 (2^3)");
  TEST_ASSERT(get_next_power_two(8) == 8, "8 -> 8 (2^3)");
  TEST_ASSERT(get_next_power_two(9) == 16, "9 -> 16 (2^4)");
  TEST_ASSERT(get_next_power_two(15) == 16, "15 -> 16 (2^4)");
  TEST_ASSERT(get_next_power_two(16) == 16, "16 -> 16 (2^4)");
  TEST_ASSERT(get_next_power_two(17) == 32, "17 -> 32 (2^5)");

  printf("✓ get_next_power_two tests passed\n");
}

void test_remap_percentage_log_like_unity() {
  printf("Testing remap_percentage_log_like_unity...\n");

  // Test edge cases
  TEST_FLOAT_CLOSE(remap_percentage_log_like_unity(0.0f), 0.0f, 1e-6f);
  TEST_FLOAT_CLOSE(remap_percentage_log_like_unity(1.0f), 0.9502f, 0.001f);

  // Test monotonicity - function should be increasing
  float val1 = remap_percentage_log_like_unity(0.1f);
  float val2 = remap_percentage_log_like_unity(0.5f);
  float val3 = remap_percentage_log_like_unity(0.9f);
  TEST_ASSERT(val1 < val2 && val2 < val3,
              "Function should be monotonically increasing");

  printf("✓ remap_percentage_log_like_unity tests passed\n");
}

int main() {
  printf("Running utility function tests...\n\n");

  test_sanitize_denormal();
  test_from_db_to_coefficient();
  test_get_next_divisible_two();
  test_get_next_power_two();
  test_remap_percentage_log_like_unity();

  printf("\n✅ All utility function tests passed!\n");
  return 0;
}
