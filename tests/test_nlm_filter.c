/*
 * Unit tests for NLM Filter (Non-Local Means)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/denoiser_logic/processing/nlm_filter_internal.h"

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

  // Process generic
  filter->process_fn = nlm_filter_process_generic;
  TEST_ASSERT(nlm_filter_process(filter, output),
              "Process should succeed when buffer is full (generic)");

  // With uniform input, output should be close to input
  for (int i = 0; i < 32; i++) {
    TEST_FLOAT_CLOSE(output[i], 5.0f, 0.01f);
  }

#if defined(__x86_64__) || defined(__i386__)
  if (__builtin_cpu_supports("avx")) {
    float output_avx[32];
    filter->process_fn = nlm_filter_process_avx;
    TEST_ASSERT(nlm_filter_process(filter, output_avx),
                "Process should succeed when buffer is full (AVX)");

    for (int i = 0; i < 32; i++) {
      TEST_FLOAT_CLOSE(output_avx[i], 5.0f, 0.01f);
    }
  }
#endif

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

  // Process generic
  filter->process_fn = nlm_filter_process_generic;
  TEST_ASSERT(nlm_filter_process(filter, output),
              "Process should succeed when buffer is full (generic)");

  // Output should be non-negative (SNR values)
  for (int i = 0; i < 32; i++) {
    TEST_ASSERT(output[i] >= 0.0f, "Output SNR (generic) should be >= 0");
  }

#if defined(__x86_64__) || defined(__i386__)
  if (__builtin_cpu_supports("avx")) {
    float output_avx[32];
    filter->process_fn = nlm_filter_process_avx;
    TEST_ASSERT(nlm_filter_process(filter, output_avx),
                "Process should succeed when buffer is full (AVX)");

    for (int i = 0; i < 32; i++) {
      TEST_ASSERT(output_avx[i] >= 0.0f, "Output SNR (AVX) should be >= 0");
    }
  }
#endif

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

void test_nlm_filter_process_silent(void) {
  printf("Testing NLM filter process with silent input...\n");

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

  float frame[32] = {0.0f};
  float output[32];

  // Fill buffer
  for (int f = 0; f < 5; f++) {
    nlm_filter_push_frame(filter, frame);
  }

  // Process generic
  filter->process_fn = nlm_filter_process_generic;
  TEST_ASSERT(nlm_filter_process(filter, output),
              "Process should succeed when buffer is full (generic)");

  // With silent input, output should be zero
  for (int i = 0; i < 32; i++) {
    TEST_FLOAT_CLOSE(output[i], 0.0f, 1e-6f);
  }

#if defined(__x86_64__) || defined(__i386__)
  if (__builtin_cpu_supports("avx")) {
    float output_avx[32];
    filter->process_fn = nlm_filter_process_avx;
    TEST_ASSERT(nlm_filter_process(filter, output_avx),
                "Process should succeed when buffer is full (AVX)");
    for (int i = 0; i < 32; i++) {
      TEST_FLOAT_CLOSE(output_avx[i], 0.0f, 1e-6f);
    }
  }
#endif

  nlm_filter_free(filter);
  printf("✓ NLM filter silent process tests passed\n");
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

  // Fill buffer
  for (int f = 0; f < 5; f++) {
    nlm_filter_push_frame(filter, frame);
  }

  // Process generic
  filter->process_fn = nlm_filter_process_generic;
  TEST_ASSERT(nlm_filter_process(filter, output),
              "Process should succeed with patch_size=8 (generic)");

  // Verify output (should be smoothed but not zero)
  for (int i = 0; i < 64; i++) {
    TEST_ASSERT(output[i] > 0.0f, "Output should be positive");
  }

#if defined(__x86_64__) || defined(__i386__)
  if (__builtin_cpu_supports("avx")) {
    float output_avx[64];
    filter->process_fn = nlm_filter_process_avx;
    TEST_ASSERT(nlm_filter_process(filter, output_avx),
                "Process should succeed with patch_size=8 (AVX)");

    for (int i = 0; i < 64; i++) {
      TEST_ASSERT(output_avx[i] > 0.0f, "Output AVX should be positive");
    }
  }
#endif

  nlm_filter_free(filter);
  printf("✓ NLM filter patch_size=8 process tests passed\n");
}

void test_nlm_filter_snr_and_reconstruct(void) {
  printf("Testing NLM filter SNR calculation and reconstruction...\n");

  NlmFilterConfig config = {
      .spectrum_size = 4,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "NLM filter initialization should succeed");

  float reference[4] = {10.0f, 20.0f, 0.0f, 40.0f};
  float noise[4] = {2.0f, 10.0f, 0.0f, 0.0f}; // Test with zero noise bins
  float snr[4] = {0.0f};

  nlm_filter_calculate_snr(filter, reference, noise, snr);

  TEST_FLOAT_CLOSE(snr[0], 5.0f, 1e-4f);
  TEST_FLOAT_CLOSE(snr[1], 2.0f, 1e-4f);
  TEST_ASSERT(snr[2] == 0.0f, "Zero reference with zero noise should be 0 SNR");
  TEST_ASSERT(snr[3] > 1000.0f,
              "High reference with zero noise should produce large SNR");

  float reconstructed[4] = {0.0f};
  float smoothed_snr[4] = {5.0f, 2.0f, 1.0f, 10.0f};

  nlm_filter_reconstruct_magnitude(filter, smoothed_snr, noise, reconstructed);

  TEST_FLOAT_CLOSE(reconstructed[0], 10.0f, 1e-4f); // 5 * 2
  TEST_FLOAT_CLOSE(reconstructed[1], 20.0f, 1e-4f); // 2 * 10
  TEST_FLOAT_CLOSE(reconstructed[2], 0.0f, 1e-4f);  // 1 * epsilon ~ 0
  TEST_FLOAT_CLOSE(reconstructed[3], 0.0f, 1e-4f);  // 10 * epsilon ~ 0

  // Null testing
  nlm_filter_calculate_snr(NULL, reference, noise, snr);
  nlm_filter_calculate_snr(filter, NULL, noise, snr);
  nlm_filter_calculate_snr(filter, reference, NULL, snr);
  nlm_filter_calculate_snr(filter, reference, noise, NULL);

  nlm_filter_reconstruct_magnitude(NULL, smoothed_snr, noise, reconstructed);
  nlm_filter_reconstruct_magnitude(filter, NULL, noise, reconstructed);
  nlm_filter_reconstruct_magnitude(filter, smoothed_snr, NULL, reconstructed);
  nlm_filter_reconstruct_magnitude(filter, smoothed_snr, noise, NULL);

  nlm_filter_free(filter);
  printf("✓ NLM filter SNR calculation and reconstruction tests passed\n");
}

void test_nlm_filter_frame_cache(void) {
  printf("Testing NLM filter frame cache indexing...\n");

  NlmFilterConfig config = {
      .spectrum_size = 32,
      .time_buffer_size = 32,
      .search_range_time_past = 8,
      .search_range_time_future = 8,
  };

  NlmFilter* filter = nlm_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "Initialization failed");
  
  // Fill buffer
  float frame[32] = {0.0f};
  for (int f = 0; f < 32; f++) {
    nlm_filter_push_frame(filter, frame);
  }

  // Precompute cache
  populate_frame_ptrs(filter);

  // Check offsets in cache range (-past-4 to +future+4)
  // For past=8, future=8, range is -12 to 12.
  for (int32_t dt = -12; dt <= 12; dt++) {
    float* cached = cached_get_frame(filter, dt);
    float* direct = get_frame(filter, dt);
    TEST_ASSERT(cached == direct, "Cached pointer should match direct lookup");
  }

  nlm_filter_free(filter);
  printf("✓ NLM filter frame cache tests passed\n");
}

int main(void) {
  printf("Running NLM filter tests...\n\n");

  test_nlm_filter_lifecycle();
  test_nlm_filter_push_frame();
  test_nlm_filter_process_uniform();
  test_nlm_filter_process_noisy();
  test_nlm_filter_process_silent();
  test_nlm_filter_process_patch8();
  test_nlm_filter_reset();
  test_nlm_filter_h_parameter();
  test_nlm_filter_latency();
  test_nlm_filter_null_handling();
  test_nlm_filter_snr_and_reconstruct();
  test_nlm_filter_frame_cache();

  printf("\n✅ All NLM filter tests passed!\n");
  return 0;
}
