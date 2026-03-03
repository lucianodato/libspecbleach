/*
 * Unit tests for NLM Filter (Non-Local Means)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/denoiser_logic/processing/nlm_filter.h"

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

  float input[32] = {0.0f};
  // Push enough frames to fill the buffer
  for (int i = 0; i < 3; i++) {
    nlm_filter_push_frame(filter, input);
  }

  float output[32];
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

void test_nlm_filter_process_patch8(void) {
  printf("Testing NLM filter process with patch_size=8...\n");

  NlmFilterConfig config = {
      .spectrum_size = 64,
      .time_buffer_size = 5,
      .patch_size = 8,
      .paste_block_size = 4,
      .search_range_freq = 4,
      .search_range_time_past = 2,
      .search_range_time_future = 2,
      .h_parameter = 2.0f,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  float frame[64];
  float output[64];
  for (int i = 0; i < 64; i++) {
    frame[i] = 5.0f + (float)(i % 2); // Alternating pattern
  }

  float* alpha = (float*)malloc(64 * sizeof(float));
  for (int i = 0; i < 64; i++) {
    alpha[i] = 1.0f;
  }

  // Fill buffer
  for (int f = 0; f < 5; f++) {
    nlm_filter_push_frame(filter, frame);
  }

  TEST_ASSERT(nlm_filter_process(filter, output),
              "Process should succeed with patch_size=8");

  // Verify output (should be smoothed but not zero)
  for (int i = 0; i < 64; i++) {
    TEST_ASSERT(output[i] > 0.0f, "Output should be positive");
  }

  free(alpha);
  nlm_filter_free(filter);
  printf("✓ NLM filter patch_size=8 process tests passed\n");
}

void test_nlm_filter_snr_calculation(void) {
  printf("Testing NLM filter SNR calculation...\n");

  NlmFilterConfig config = {.spectrum_size = 32};
  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "Filter initialization should succeed");

  float reference[32];
  float noise[32];
  float snr_output[32];

  // Test case 1: Normal SNR calculation
  for (int i = 0; i < 32; i++) {
    reference[i] = 10.0f;
    noise[i] = 2.0f;
  }

  nlm_filter_calculate_snr(filter, reference, noise, snr_output);

  for (int i = 0; i < 32; i++) {
    TEST_FLOAT_CLOSE(snr_output[i], 5.0f, 1e-5f); // 10/2 = 5
  }

  // Test case 2: Very small noise (should handle gracefully)
  for (int i = 0; i < 32; i++) {
    reference[i] = 1.0f;
    noise[i] = 1e-20f;
  }

  nlm_filter_calculate_snr(filter, reference, noise, snr_output);

  for (int i = 0; i < 32; i++) {
    TEST_ASSERT(snr_output[i] > 0.0f && !isinf(snr_output[i]),
                "SNR should be valid with very small noise");
  }

  // Test case 3: Zero reference
  for (int i = 0; i < 32; i++) {
    reference[i] = 0.0f;
    noise[i] = 1.0f;
  }

  nlm_filter_calculate_snr(filter, reference, noise, snr_output);

  for (int i = 0; i < 32; i++) {
    TEST_FLOAT_CLOSE(snr_output[i], 0.0f, 1e-5f);
  }

  // Test NULL handling
  nlm_filter_calculate_snr(NULL, reference, noise, snr_output);
  nlm_filter_calculate_snr(filter, NULL, noise, snr_output);
  nlm_filter_calculate_snr(filter, reference, NULL, snr_output);
  nlm_filter_calculate_snr(filter, reference, noise, NULL);

  nlm_filter_free(filter);
  printf("✓ NLM filter SNR calculation tests passed\n");
}

void test_nlm_filter_magnitude_reconstruction(void) {
  printf("Testing NLM filter magnitude reconstruction...\n");

  NlmFilterConfig config = {.spectrum_size = 32};
  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "Filter initialization should succeed");

  float smoothed_snr[32];
  float noise[32];
  float magnitude_output[32];

  // Test case 1: Normal reconstruction
  for (int i = 0; i < 32; i++) {
    smoothed_snr[i] = 5.0f;
    noise[i] = 2.0f;
  }

  nlm_filter_reconstruct_magnitude(filter, smoothed_snr, noise,
                                   magnitude_output);

  for (int i = 0; i < 32; i++) {
    TEST_FLOAT_CLOSE(magnitude_output[i], 10.0f, 1e-5f); // 5*2 = 10
  }

  // Test case 2: Very small noise
  for (int i = 0; i < 32; i++) {
    smoothed_snr[i] = 100.0f;
    noise[i] = 1e-20f;
  }

  nlm_filter_reconstruct_magnitude(filter, smoothed_snr, noise,
                                   magnitude_output);

  for (int i = 0; i < 32; i++) {
    TEST_ASSERT(magnitude_output[i] >= 0.0f && !isinf(magnitude_output[i]),
                "Reconstructed magnitude should be valid");
  }

  // Test case 3: Zero SNR
  for (int i = 0; i < 32; i++) {
    smoothed_snr[i] = 0.0f;
    noise[i] = 5.0f;
  }

  nlm_filter_reconstruct_magnitude(filter, smoothed_snr, noise,
                                   magnitude_output);

  for (int i = 0; i < 32; i++) {
    TEST_FLOAT_CLOSE(magnitude_output[i], 0.0f, 1e-5f);
  }

  // Test NULL handling
  nlm_filter_reconstruct_magnitude(NULL, smoothed_snr, noise, magnitude_output);
  nlm_filter_reconstruct_magnitude(filter, NULL, noise, magnitude_output);
  nlm_filter_reconstruct_magnitude(filter, smoothed_snr, NULL, magnitude_output);
  nlm_filter_reconstruct_magnitude(filter, smoothed_snr, noise, NULL);

  nlm_filter_free(filter);
  printf("✓ NLM filter magnitude reconstruction tests passed\n");
}

void test_nlm_filter_round_trip(void) {
  printf("Testing NLM filter SNR calculation and reconstruction round trip...\n");

  NlmFilterConfig config = {.spectrum_size = 32};
  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "Filter initialization should succeed");

  float original_spectrum[32];
  float noise[32];
  float snr[32];
  float reconstructed[32];

  // Generate test data
  for (int i = 0; i < 32; i++) {
    original_spectrum[i] = 10.0f + (float)i * 0.5f;
    noise[i] = 2.0f + (float)i * 0.1f;
  }

  // Calculate SNR
  nlm_filter_calculate_snr(filter, original_spectrum, noise, snr);

  // Reconstruct magnitude from SNR
  nlm_filter_reconstruct_magnitude(filter, snr, noise, reconstructed);

  // Verify round trip preserves original spectrum (approximately)
  for (int i = 0; i < 32; i++) {
    TEST_FLOAT_CLOSE(reconstructed[i], original_spectrum[i], 1e-3f);
  }

  nlm_filter_free(filter);
  printf("✓ NLM filter round trip tests passed\n");
}

void test_nlm_filter_boundary_conditions(void) {
  printf("Testing NLM filter with boundary conditions...\n");

  // Test with very small spectrum size
  NlmFilterConfig config_small = {
      .spectrum_size = 8,
      .patch_size = 4,
      .search_range_freq = 2,
      .h_parameter = 1.0f,
  };

  NlmFilter* filter_small = nlm_filter_initialize(config_small);
  TEST_ASSERT(filter_small != NULL,
              "Small spectrum initialization should succeed");

  float frame_small[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
  for (int i = 0; i < config_small.search_range_time_past +
                          config_small.search_range_time_future + 1;
       i++) {
    nlm_filter_push_frame(filter_small, frame_small);
  }

  float output_small[8];
  TEST_ASSERT(nlm_filter_process(filter_small, output_small),
              "Process should succeed with small spectrum");

  nlm_filter_free(filter_small);

  // Test with large spectrum size
  NlmFilterConfig config_large = {
      .spectrum_size = 2048,
      .patch_size = 8,
      .search_range_freq = 4,
      .h_parameter = 1.0f,
  };

  NlmFilter* filter_large = nlm_filter_initialize(config_large);
  TEST_ASSERT(filter_large != NULL,
              "Large spectrum initialization should succeed");

  float* frame_large = (float*)calloc(2048, sizeof(float));
  for (int i = 0; i < 2048; i++) {
    frame_large[i] = 1.0f + sinf((float)i * 0.01f);
  }

  for (int i = 0; i < config_large.search_range_time_past +
                          config_large.search_range_time_future + 1;
       i++) {
    nlm_filter_push_frame(filter_large, frame_large);
  }

  float* output_large = (float*)calloc(2048, sizeof(float));
  TEST_ASSERT(nlm_filter_process(filter_large, output_large),
              "Process should succeed with large spectrum");

  free(frame_large);
  free(output_large);
  nlm_filter_free(filter_large);

  printf("✓ NLM filter boundary condition tests passed\n");
}

int main(void) {
  printf("Running NLM filter tests...\n\n");

  test_nlm_filter_lifecycle();
  test_nlm_filter_push_frame();
  test_nlm_filter_process_uniform();
  test_nlm_filter_process_noisy();
  test_nlm_filter_process_patch8();
  test_nlm_filter_reset();
  test_nlm_filter_h_parameter();
  test_nlm_filter_latency();
  test_nlm_filter_null_handling();
  test_nlm_filter_snr_calculation();
  test_nlm_filter_magnitude_reconstruction();
  test_nlm_filter_round_trip();
  test_nlm_filter_boundary_conditions();

  printf("\n✅ All NLM filter tests passed!\n");
  return 0;
}