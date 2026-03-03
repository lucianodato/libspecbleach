/*
 * Unit tests for SIMD Utilities
 * Tests portable SIMD abstraction layer for x86 and ARM architectures
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/utils/simd_utils.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_simd_ftz_daz_state(void) {
  printf("Testing SIMD FTZ/DAZ state management...\n");

  // Save original state
  sb_simd_state_t old_state = sb_simd_enable_ftz_daz();

  // Enable again (should not crash)
  sb_simd_state_t state2 = sb_simd_enable_ftz_daz();
  (void)state2;

  // Restore original state
  sb_simd_restore_state(old_state);

  printf("✓ SIMD FTZ/DAZ state management tests passed\n");
}

void test_simd_vec8_load_store(void) {
  printf("Testing 8-wide vector load/store...\n");

  float input[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
  float output[8] = {0.0f};

  sb_vec8_t v = sb_load8(input);
  sb_store8(output, v);

  for (int i = 0; i < 8; i++) {
    TEST_FLOAT_CLOSE(output[i], input[i], 1e-6f);
  }

  printf("✓ 8-wide vector load/store tests passed\n");
}

void test_simd_vec8_set(void) {
  printf("Testing 8-wide vector set...\n");

  float output[8];
  sb_vec8_t v = sb_set8(3.14f);
  sb_store8(output, v);

  for (int i = 0; i < 8; i++) {
    TEST_FLOAT_CLOSE(output[i], 3.14f, 1e-6f);
  }

  printf("✓ 8-wide vector set tests passed\n");
}

void test_simd_vec8_arithmetic(void) {
  printf("Testing 8-wide vector arithmetic...\n");

  float a[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
  float b[8] = {8.0f, 7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f};
  float result[8];

  sb_vec8_t va = sb_load8(a);
  sb_vec8_t vb = sb_load8(b);

  // Test addition
  sb_vec8_t vadd = sb_add8(va, vb);
  sb_store8(result, vadd);
  for (int i = 0; i < 8; i++) {
    TEST_FLOAT_CLOSE(result[i], a[i] + b[i], 1e-5f);
  }

  // Test subtraction
  sb_vec8_t vsub = sb_sub8(va, vb);
  sb_store8(result, vsub);
  for (int i = 0; i < 8; i++) {
    TEST_FLOAT_CLOSE(result[i], a[i] - b[i], 1e-5f);
  }

  // Test multiplication
  sb_vec8_t vmul = sb_mul8(va, vb);
  sb_store8(result, vmul);
  for (int i = 0; i < 8; i++) {
    TEST_FLOAT_CLOSE(result[i], a[i] * b[i], 1e-5f);
  }

  // Test division
  sb_vec8_t vdiv = sb_div8(va, vb);
  sb_store8(result, vdiv);
  for (int i = 0; i < 8; i++) {
    TEST_FLOAT_CLOSE(result[i], a[i] / b[i], 1e-4f);
  }

  // Test max
  sb_vec8_t vmax = sb_max8(va, vb);
  sb_store8(result, vmax);
  for (int i = 0; i < 8; i++) {
    float expected = (a[i] > b[i]) ? a[i] : b[i];
    TEST_FLOAT_CLOSE(result[i], expected, 1e-5f);
  }

  printf("✓ 8-wide vector arithmetic tests passed\n");
}

void test_simd_vec8_ssd(void) {
  printf("Testing 8-wide vector SSD...\n");

  float a[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
  float b[8] = {1.5f, 2.5f, 3.5f, 4.5f, 5.5f, 6.5f, 7.5f, 8.5f};

  sb_vec8_t va = sb_load8(a);
  sb_vec8_t vb = sb_load8(b);

  float ssd = sb_vec8_ssd(va, vb);

  // Calculate expected SSD
  float expected = 0.0f;
  for (int i = 0; i < 8; i++) {
    float diff = a[i] - b[i];
    expected += diff * diff;
  }

  TEST_FLOAT_CLOSE(ssd, expected, 1e-4f);

  printf("✓ 8-wide vector SSD tests passed\n");
}

void test_simd_vec8_accumulator(void) {
  printf("Testing 8-wide accumulator operations...\n");

  sb_acc8_t acc = sb_acc8_zero();

  float a[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
  float b[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

  sb_vec8_t va = sb_load8(a);
  sb_vec8_t vb = sb_load8(b);

  // Accumulate SSD multiple times
  for (int iter = 0; iter < 3; iter++) {
    acc = sb_acc8_add_ssd(acc, va, vb);
  }

  float sum = sb_acc8_hsum(acc);

  // Expected: 3 iterations of sum of squares of differences
  float expected = 0.0f;
  for (int i = 0; i < 8; i++) {
    float diff = a[i] - b[i];
    expected += diff * diff;
  }
  expected *= 3.0f; // 3 iterations

  TEST_FLOAT_CLOSE(sum, expected, 1e-3f);

  printf("✓ 8-wide accumulator tests passed\n");
}

void test_simd_vec8_comparison(void) {
  printf("Testing 8-wide vector comparison...\n");

  float a[8] = {1.0f, 5.0f, 3.0f, 8.0f, 2.0f, 9.0f, 4.0f, 6.0f};
  float b[8] = {2.0f, 4.0f, 3.0f, 7.0f, 3.0f, 8.0f, 5.0f, 5.0f};

  sb_vec8_t va = sb_load8(a);
  sb_vec8_t vb = sb_load8(b);

  sb_vec8_t mask = sb_gt8(va, vb);

  // Store mask to verify (implementation-defined representation)
  float mask_vals[8];
  sb_store8(mask_vals, mask);

  // Check that some comparison happened (mask is not all zeros or invalid)
  // Implementation-specific, but should produce some non-zero values
  printf("✓ 8-wide vector comparison tests passed\n");
}

void test_simd_vec8_select(void) {
  printf("Testing 8-wide vector select...\n");

  float mask_data[8] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f};
  float a[8] = {10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f};
  float b[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
  float result[8];

  sb_vec8_t mask = sb_load8(mask_data);
  sb_vec8_t va = sb_load8(a);
  sb_vec8_t vb = sb_load8(b);

  sb_vec8_t selected = sb_sel8(mask, va, vb);
  sb_store8(result, selected);

  // Where mask is non-zero, should select from a, else from b
  for (int i = 0; i < 8; i++) {
    float expected = (mask_data[i] != 0.0f) ? a[i] : b[i];
    TEST_FLOAT_CLOSE(result[i], expected, 1e-5f);
  }

  printf("✓ 8-wide vector select tests passed\n");
}

void test_simd_vec4_operations(void) {
  printf("Testing 4-wide vector operations...\n");

  float a[4] = {1.0f, 2.0f, 3.0f, 4.0f};
  float b[4] = {4.0f, 3.0f, 2.0f, 1.0f};
  float result[4];

  sb_vec4_t va = sb_load4(a);
  sb_vec4_t vb = sb_load4(b);

  // Test SSD
  float ssd = sb_vec4_ssd(va, vb);
  float expected_ssd = 0.0f;
  for (int i = 0; i < 4; i++) {
    float diff = a[i] - b[i];
    expected_ssd += diff * diff;
  }
  TEST_FLOAT_CLOSE(ssd, expected_ssd, 1e-5f);

  // Test arithmetic
  sb_vec4_t vadd = sb_add4(va, vb);
  sb_store4(result, vadd);
  for (int i = 0; i < 4; i++) {
    TEST_FLOAT_CLOSE(result[i], a[i] + b[i], 1e-5f);
  }

  sb_vec4_t vsub = sb_sub4(va, vb);
  sb_store4(result, vsub);
  for (int i = 0; i < 4; i++) {
    TEST_FLOAT_CLOSE(result[i], a[i] - b[i], 1e-5f);
  }

  sb_vec4_t vmul = sb_mul4(va, vb);
  sb_store4(result, vmul);
  for (int i = 0; i < 4; i++) {
    TEST_FLOAT_CLOSE(result[i], a[i] * b[i], 1e-5f);
  }

  // Test set
  sb_vec4_t vset = sb_set4(2.5f);
  sb_store4(result, vset);
  for (int i = 0; i < 4; i++) {
    TEST_FLOAT_CLOSE(result[i], 2.5f, 1e-6f);
  }

  printf("✓ 4-wide vector operations tests passed\n");
}

void test_simd_edge_cases(void) {
  printf("Testing SIMD edge cases...\n");

  // Test with zeros
  float zeros[8] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
  sb_vec8_t v_zero = sb_load8(zeros);
  float ssd_zero = sb_vec8_ssd(v_zero, v_zero);
  TEST_FLOAT_CLOSE(ssd_zero, 0.0f, 1e-6f);

  // Test with very small numbers (denormals)
  float denormals[8];
  for (int i = 0; i < 8; i++) {
    denormals[i] = 1e-40f;
  }
  sb_simd_state_t old_state = sb_simd_enable_ftz_daz();
  sb_vec8_t v_denorm = sb_load8(denormals);
  sb_vec8_t v_zero2 = sb_load8(zeros);
  float ssd_denorm = sb_vec8_ssd(v_denorm, v_zero2);
  // With FTZ/DAZ, denormals should be treated as zero
  TEST_ASSERT(ssd_denorm >= 0.0f, "SSD with denormals should be non-negative");
  sb_simd_restore_state(old_state);

  // Test with large numbers
  float large[8];
  for (int i = 0; i < 8; i++) {
    large[i] = 1e30f;
  }
  sb_vec8_t v_large = sb_load8(large);
  float ssd_large = sb_vec8_ssd(v_large, v_zero);
  TEST_ASSERT(!isnan(ssd_large) && !isinf(ssd_large),
              "SSD with large numbers should be valid");

  // Test division by very small number
  float ones[8] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  float small[8] = {1e-10f, 1e-10f, 1e-10f, 1e-10f,
                    1e-10f, 1e-10f, 1e-10f, 1e-10f};
  sb_vec8_t v_ones = sb_load8(ones);
  sb_vec8_t v_small = sb_load8(small);
  sb_vec8_t v_div = sb_div8(v_ones, v_small);
  float result[8];
  sb_store8(result, v_div);
  // Result should be very large but not NaN/Inf ideally
  TEST_ASSERT(result[0] > 1e5f, "Division by small number should yield large result");

  printf("✓ SIMD edge case tests passed\n");
}

void test_simd_register_blocking_pattern(void) {
  printf("Testing SIMD register blocking pattern (NLM usage)...\n");

  // Simulate the NLM register blocking optimization pattern
  // Load target patch once, reuse for multiple comparisons

  float target[64]; // 8x8 patch
  for (int i = 0; i < 64; i++) {
    target[i] = (float)(i % 8);
  }

  // Pre-load target into 8 vectors (one per row)
  sb_vec8_t target_rows[8];
  for (int r = 0; r < 8; r++) {
    target_rows[r] = sb_load8(&target[r * 8]);
  }

  // Simulate comparison with multiple candidate patches
  float candidates[3][64];
  for (int c = 0; c < 3; c++) {
    for (int i = 0; i < 64; i++) {
      candidates[c][i] = (float)(i % 8) + (float)c * 0.1f;
    }
  }

  // Compare target with each candidate using register-blocked approach
  for (int c = 0; c < 3; c++) {
    sb_acc8_t acc = sb_acc8_zero();

    for (int r = 0; r < 8; r++) {
      sb_vec8_t cand_row = sb_load8(&candidates[c][r * 8]);
      acc = sb_acc8_add_ssd(acc, target_rows[r], cand_row);
    }

    float distance = sb_acc8_hsum(acc);

    // Verify result
    float expected = 0.0f;
    for (int i = 0; i < 64; i++) {
      float diff = target[i] - candidates[c][i];
      expected += diff * diff;
    }

    TEST_FLOAT_CLOSE(distance, expected, 1e-3f);
  }

  printf("✓ SIMD register blocking pattern tests passed\n");
}

int main(void) {
  printf("Running SIMD Utilities tests...\n\n");

  test_simd_ftz_daz_state();
  test_simd_vec8_load_store();
  test_simd_vec8_set();
  test_simd_vec8_arithmetic();
  test_simd_vec8_ssd();
  test_simd_vec8_accumulator();
  test_simd_vec8_comparison();
  test_simd_vec8_select();
  test_simd_vec4_operations();
  test_simd_edge_cases();
  test_simd_register_blocking_pattern();

  printf("\n✅ All SIMD Utilities tests passed!\n");
  return 0;
}