/*
 * Unit tests for NLM Filter (Non-Local Means)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/post_estimation/nlm_filter.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_nlm_filter_lifecycle(void) {
  printf("Testing NLM filter lifecycle...\n");

  NlmFilterConfig config = {
      .spectrum_size = 513,
      .time_buffer_size = 21,
      .patch_size = 8,
      .paste_block_size = 4,
      .search_range_freq = 8,
      .search_range_time_past = 16,
      .search_range_time_future = 4,
      .h_parameter = 1.0f,
      .distance_threshold = 0.0f,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  nlm_filter_free(filter);

  // Test with NULL config values (should use defaults)
  NlmFilterConfig config_minimal = {
      .spectrum_size = 256,
  };

  filter = nlm_filter_initialize(config_minimal);
  TEST_ASSERT(filter != NULL,
              "NLM filter initialization with minimal config should succeed");
  nlm_filter_free(filter);

  // Test with invalid spectrum_size
  NlmFilterConfig config_invalid = {
      .spectrum_size = 0,
  };

  filter = nlm_filter_initialize(config_invalid);
  TEST_ASSERT(filter == NULL,
              "NLM filter initialization with zero spectrum_size should fail");

  printf("✓ NLM filter lifecycle tests passed\n");
}

void test_nlm_filter_push_frame(void) {
  printf("Testing NLM filter push frame...\n");

  NlmFilterConfig config = {
      .spectrum_size = 64,
      .time_buffer_size = 5,
      .patch_size = 4,
      .paste_block_size = 2,
      .search_range_freq = 4,
      .search_range_time_past = 2,
      .search_range_time_future = 2,
      .h_parameter = 1.0f,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  // Initially not ready
  TEST_ASSERT(!nlm_filter_is_ready(filter),
              "Filter should not be ready initially");

  // Push frames
  float frame[64];
  for (int i = 0; i < 64; i++) {
    frame[i] = 1.0f;
  }

  for (int f = 0; f < 4; f++) {
    nlm_filter_push_frame(filter, frame);
    TEST_ASSERT(!nlm_filter_is_ready(filter),
                "Filter should not be ready until buffer is full");
  }

  // Push one more frame to fill buffer
  nlm_filter_push_frame(filter, frame);
  TEST_ASSERT(nlm_filter_is_ready(filter),
              "Filter should be ready after buffer is filled");

  nlm_filter_free(filter);
  printf("✓ NLM filter push frame tests passed\n");
}

void test_nlm_filter_process_uniform(void) {
  printf("Testing NLM filter process with uniform input...\n");

  NlmFilterConfig config = {
      .spectrum_size = 32,
      .time_buffer_size = 5,
      .patch_size = 4,
      .paste_block_size = 2,
      .search_range_freq = 2,
      .search_range_time_past = 2,
      .search_range_time_future = 2,
      .h_parameter = 1.0f,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  // Push uniform frames (should result in identical output)
  float frame[32];
  float output[32];
  for (int i = 0; i < 32; i++) {
    frame[i] = 5.0f; // Uniform SNR value
  }

  // Fill buffer
  for (int f = 0; f < 5; f++) {
    nlm_filter_push_frame(filter, frame);
  }

  // Process
  TEST_ASSERT(nlm_filter_process(filter, output),
              "Process should succeed when buffer is full");

  // With uniform input, output should be close to input
  for (int i = 0; i < 32; i++) {
    TEST_FLOAT_CLOSE(output[i], 5.0f, 0.01f);
  }

  nlm_filter_free(filter);
  printf("✓ NLM filter uniform process tests passed\n");
}

void test_nlm_filter_process_noisy(void) {
  printf("Testing NLM filter process with noisy input...\n");

  NlmFilterConfig config = {
      .spectrum_size = 32,
      .time_buffer_size = 5,
      .patch_size = 4,
      .paste_block_size = 2,
      .search_range_freq = 2,
      .search_range_time_past = 2,
      .search_range_time_future = 2,
      .h_parameter = 5.0f, // Higher h for more smoothing
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  // Create frames with a consistent signal + some noise
  float frames[5][32];
  for (int f = 0; f < 5; f++) {
    for (int i = 0; i < 32; i++) {
      // Base pattern: harmonic-like structure
      frames[f][i] = 2.0f + sinf((float)i * 0.5f);
      // Add some "noise" (different per frame)
      frames[f][i] += 0.1f * (float)(f % 3);
    }
  }

  // Push frames
  for (int f = 0; f < 5; f++) {
    nlm_filter_push_frame(filter, frames[f]);
  }

  float output[32];
  TEST_ASSERT(nlm_filter_process(filter, output),
              "Process should succeed when buffer is full");

  // Output should be non-negative (SNR values)
  for (int i = 0; i < 32; i++) {
    TEST_ASSERT(output[i] >= 0.0f, "Output SNR should be non-negative");
  }

  nlm_filter_free(filter);
  printf("✓ NLM filter noisy process tests passed\n");
}

void test_nlm_filter_reset(void) {
  printf("Testing NLM filter reset...\n");

  NlmFilterConfig config = {
      .spectrum_size = 32,
      .time_buffer_size = 3,
      .h_parameter = 1.0f,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  float frame[32];
  for (int i = 0; i < 32; i++) {
    frame[i] = 1.0f;
  }

  // Fill buffer
  for (int f = 0; f < 3; f++) {
    nlm_filter_push_frame(filter, frame);
  }
  TEST_ASSERT(nlm_filter_is_ready(filter), "Filter should be ready");

  // Reset
  nlm_filter_reset(filter);
  TEST_ASSERT(!nlm_filter_is_ready(filter),
              "Filter should not be ready after reset");

  nlm_filter_free(filter);
  printf("✓ NLM filter reset tests passed\n");
}

void test_nlm_filter_h_parameter(void) {
  printf("Testing NLM filter h parameter update...\n");

  NlmFilterConfig config = {
      .spectrum_size = 32,
      .time_buffer_size = 3,
      .h_parameter = 1.0f,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  // Update h parameter
  nlm_filter_set_h_parameter(filter, 2.5f);

  // Fill buffer and process
  float frame[32];
  float output[32];
  for (int i = 0; i < 32; i++) {
    frame[i] = 3.0f;
  }

  for (int f = 0; f < 3; f++) {
    nlm_filter_push_frame(filter, frame);
  }

  TEST_ASSERT(nlm_filter_process(filter, output),
              "Process should succeed after h parameter update");

  nlm_filter_free(filter);
  printf("✓ NLM filter h parameter tests passed\n");
}

void test_nlm_filter_latency(void) {
  printf("Testing NLM filter latency reporting...\n");

  NlmFilterConfig config = {
      .spectrum_size = 64,
      .search_range_time_past = 16,
      .search_range_time_future = 4,
      .h_parameter = 1.0f,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  uint32_t latency = nlm_filter_get_latency_frames(filter);
  TEST_ASSERT(latency == 4, "Latency should equal search_range_time_future");

  nlm_filter_free(filter);

  // Test with NULL
  TEST_ASSERT(nlm_filter_get_latency_frames(NULL) == 0,
              "NULL filter should return 0 latency");

  printf("✓ NLM filter latency tests passed\n");
}

void test_nlm_filter_null_handling(void) {
  printf("Testing NLM filter NULL handling...\n");

  // Test NULL filter operations
  nlm_filter_free(NULL); // Should not crash

  nlm_filter_set_h_parameter(NULL, 1.0f); // Should not crash

  nlm_filter_push_frame(NULL, NULL); // Should not crash

  TEST_ASSERT(!nlm_filter_is_ready(NULL), "NULL filter should not be ready");

  float output[32];
  TEST_ASSERT(!nlm_filter_process(NULL, output),
              "NULL filter process should fail");

  NlmFilterConfig config = {.spectrum_size = 32};
  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(!nlm_filter_process(filter, NULL),
              "NULL output process should fail");

  nlm_filter_reset(NULL); // Should not crash

  nlm_filter_free(filter);
  printf("✓ NLM filter NULL handling tests passed\n");
}

int main(void) {
  printf("Running NLM filter tests...\n\n");

  test_nlm_filter_lifecycle();
  test_nlm_filter_push_frame();
  test_nlm_filter_process_uniform();
  test_nlm_filter_process_noisy();
  test_nlm_filter_reset();
  test_nlm_filter_h_parameter();
  test_nlm_filter_latency();
  test_nlm_filter_null_handling();

  printf("\n✅ All NLM filter tests passed!\n");
  return 0;
}
