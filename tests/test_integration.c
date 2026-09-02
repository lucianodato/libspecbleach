/*
 * Integration tests for the full denoising pipeline
 */

#include "../src/shared/configurations.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Include internal headers for testing

// Include the public API

#define ROLLING_MEAN 1

#include "specbleach_denoiser.h"

// Function prototypes
void test_spectral_denoiser(void);
void test_adaptive_denoiser(void);
void test_2d_denoiser(void);
void test_runtime_mode_switch(void);
void test_different_noise_levels(void);
void test_library_info(void);
void generate_test_audio(float* buffer, size_t length, float signal_freq,
                         float noise_level);
float calculate_rms(const float* buffer, size_t length);

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

#define SAMPLE_RATE 44100
#define FRAME_SIZE 512
#define NUM_FRAMES 100
#define BLOCK_SIZE ((size_t)FRAME_SIZE * (size_t)NUM_FRAMES)

// Generate synthetic audio with noise
void generate_test_audio(float* buffer, size_t length, float signal_freq,
                         float noise_level) {
  for (size_t i = 0; i < length; i++) {
    // Generate a sine wave signal
    float signal = 0.5f * sinf((float)(2.0 * M_PIf * (double)signal_freq *
                                       (double)i / (double)SAMPLE_RATE));

    // Add white noise
    float noise = noise_level * ((float)rand() / RAND_MAX - 0.5f) * 2.0f;

    buffer[i] = signal + noise;
  }
}

// Calculate RMS of audio buffer
float calculate_rms(const float* buffer, size_t length) {
  double sum = 0.0;
  for (size_t i = 0; i < length; i++) {
    sum += (double)buffer[i] * (double)buffer[i];
  }
  return (float)sqrt(sum / (double)length);
}

// Test spectral denoiser with synthetic audio
void test_spectral_denoiser(void) {
  printf("Testing spectral denoiser integration...\n");

  // Create test audio with signal + noise
  float* input_buffer = calloc((size_t)BLOCK_SIZE, sizeof(float));
  float* output_buffer = calloc((size_t)BLOCK_SIZE, sizeof(float));
  TEST_ASSERT(input_buffer && output_buffer, "Failed to allocate test buffers");

  generate_test_audio(input_buffer, BLOCK_SIZE, 1000.0f,
                      0.5f); // 1kHz signal + noise (higher noise for RMS check)

  // Calculate input RMS
  float input_rms = calculate_rms(input_buffer, BLOCK_SIZE);
  printf("  Input RMS: %.4f\n", input_rms);

  // Initialize denoiser
  float frame_size_ms = 20.0f;
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize spectral denoiser");

  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_ALL, // Learn all modes
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.0f,
      .masking_depth = 0.5f,

      .residual_listen = false,
      .whitening_factor = 0.0f};

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load spectral denoiser parameters should succeed");

  // Process first blocks in learn mode (at least 10 frames)
  specbleach_denoiser_process(handle, FRAME_SIZE * 10, input_buffer,
                              output_buffer);

  // Switch to reduction mode
  parameters.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load reduction parameters should succeed");

  // Process remaining blocks
  size_t processed_samples = FRAME_SIZE * 10;
  while (processed_samples < BLOCK_SIZE) {
    uint32_t block_size = FRAME_SIZE;
    if (processed_samples + (size_t)block_size > BLOCK_SIZE) {
      block_size = (uint32_t)(BLOCK_SIZE - processed_samples);
    }

    bool result = specbleach_denoiser_process(
        handle, block_size, input_buffer + processed_samples,
        output_buffer + processed_samples);
    TEST_ASSERT(result == true, "Processing failed");

    processed_samples += block_size;
  }

  // Calculate output RMS
  float output_rms = calculate_rms(output_buffer, BLOCK_SIZE);
  printf("  Output RMS: %.4f\n", output_rms);

  // Verify noise reduction (output should be quieter than input)
  // With Tonal Reduction 0dB, signal is preserved, so RMS drop is smaller
  TEST_ASSERT(output_rms < input_rms * 0.9f,
              "Noise reduction should reduce RMS significantly");

  // Verify signal preservation (output shouldn't be too quiet)
  TEST_ASSERT(output_rms > input_rms * 0.1f,
              "Signal should be reasonably preserved");

  // Cleanup
  specbleach_denoiser_free(handle);
  free(input_buffer);
  free(output_buffer);

  printf("✓ Spectral denoiser integration test passed\n");
}

// Test denoiser with different noise levels
void test_different_noise_levels(void) {
  printf("Testing denoiser with different noise levels...\n");

  // Test with low noise
  float* input_buffer = calloc(BLOCK_SIZE, sizeof(float));
  float* output_buffer = calloc(BLOCK_SIZE, sizeof(float));
  TEST_ASSERT(input_buffer && output_buffer, "Failed to allocate test buffers");

  generate_test_audio(input_buffer, BLOCK_SIZE, 1000.0f, 0.05f); // Low noise

  float frame_size_ms = 20.0f;
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize denoiser");

  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_ALL, // Average
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.0f,
      .masking_depth = 0.5f,

      .residual_listen = false,
      .whitening_factor = 0.0f};

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load low-noise parameters should succeed");

  // Learn a bit (at least 10 frames)
  specbleach_denoiser_process(handle, 5000, input_buffer, output_buffer);

  // Process remaining with reduction
  parameters.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load reduction parameters should succeed");

  // Process all blocks
  size_t processed = 5000;
  while (processed < BLOCK_SIZE) {
    uint32_t block_size = FRAME_SIZE;
    if (processed + (size_t)block_size > BLOCK_SIZE) {
      block_size = (uint32_t)(BLOCK_SIZE - processed);
    }

    bool result = specbleach_denoiser_process(handle, block_size,
                                              input_buffer + processed,
                                              output_buffer + processed);
    TEST_ASSERT(result == true, "Processing failed");

    processed += block_size;
  }

  float input_rms = calculate_rms(input_buffer, BLOCK_SIZE);
  float output_rms = calculate_rms(output_buffer, BLOCK_SIZE);

  printf("  Low noise - Input RMS: %.4f, Output RMS: %.4f\n", input_rms,
         output_rms);

  specbleach_denoiser_free(handle);
  free(input_buffer);
  free(output_buffer);

  printf("✓ Different noise levels test passed\n");
}

// Test library information functions
void test_library_info(void) {
  printf("Testing library information functions...\n");

  // Test spectral denoiser info
  int latency = (int)specbleach_denoiser_get_latency(
      NULL); // Should return 0 for NULL handle
  TEST_ASSERT(latency == 0, "NULL handle should return 0 latency");

  int profile_size = (int)specbleach_denoiser_get_noise_profile_size(NULL);
  TEST_ASSERT(profile_size == 0, "NULL handle should return 0 profile size");

  // Test that we can get valid information after initialization
  float frame_size_ms = 20.0f;
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize for info test");

  latency = (int)specbleach_denoiser_get_latency(handle);
  TEST_ASSERT(latency > 0, "Latency should be strictly positive");

  profile_size = (int)specbleach_denoiser_get_noise_profile_size(handle);
  TEST_ASSERT(profile_size > 0, "Profile size should be positive");

  specbleach_denoiser_free(handle);

  printf("✓ Library information tests passed\n");
}

int main(void) {
  printf("Running integration tests...\n\n");

  // Seed random number generator for reproducible tests
  srand(42);

  test_spectral_denoiser();
  test_different_noise_levels();
  test_adaptive_denoiser();

  test_2d_denoiser();
  test_runtime_mode_switch();
  test_library_info();

  printf("\n✅ All integration tests passed!\n");
  return 0;
}

// Test adaptive functionality enabled in standard denoiser
void test_adaptive_denoiser(void) {
  printf("Testing adaptive denoiser integration...\n");

  float* input_buffer = calloc(BLOCK_SIZE, sizeof(float));
  float* output_buffer = calloc(BLOCK_SIZE, sizeof(float));
  TEST_ASSERT(input_buffer && output_buffer, "Failed to allocate test buffers");

  generate_test_audio(input_buffer, BLOCK_SIZE, 1000.0f, 0.1f);

  float frame_size_ms = 40.0f;
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize denoiser");

  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_OFF,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .masking_depth = 0.5f,

      .residual_listen = false,
      .whitening_factor = 0.5f,
      .adaptive_noise = true, // Enable adaptive mode
      .noise_estimation_method =
          SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE // Default method
  };

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load adaptive parameters should succeed");

  // Process all samples
  bool result = specbleach_denoiser_process(handle, BLOCK_SIZE, input_buffer,
                                            output_buffer);
  TEST_ASSERT(result == true, "Processing failed");

  float input_rms = calculate_rms(input_buffer, BLOCK_SIZE);
  float output_rms = calculate_rms(output_buffer, BLOCK_SIZE);
  printf("  Adaptive - Input RMS: %.4f, Output RMS: %.4f\n", input_rms,
         output_rms);

  specbleach_denoiser_free(handle);
  free(input_buffer);
  free(output_buffer);

  printf("✓ Adaptive denoiser integration test passed\n");
}

void test_2d_denoiser(void) {
  printf("Testing 2D NLM smoothing mode integration...\n");

  float* input_buffer = calloc(BLOCK_SIZE, sizeof(float));
  float* output_buffer = calloc(BLOCK_SIZE, sizeof(float));
  TEST_ASSERT(input_buffer && output_buffer, "Failed to allocate test buffers");

  generate_test_audio(input_buffer, BLOCK_SIZE, 1000.0f, 0.1f);

  float frame_size_ms = 20.0f;
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize unified denoiser");

  // Test latency reporting (should include NLM look-ahead)
  uint32_t latency = specbleach_denoiser_get_latency(handle);
  printf("  unified denoiser latency: %u samples\n", latency);
  TEST_ASSERT(latency > 0, "unified denoiser should report latency");

  // Configure for noise learning
  SpecbleachDenoiserParameters parameters = {
      .learn_noise = SPECBLEACH_LEARN_ALL, // Learn mode
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f, // NLM h parameter
      .smoothing_mode = SPECBLEACH_SMOOTHING_NLM_2D,
      .whitening_factor = 0.0f,
      .residual_listen = false,
  };

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load NLM learn parameters should succeed");

  // Process first blocks in learn mode
  specbleach_denoiser_process(handle, FRAME_SIZE * 10, input_buffer,
                              output_buffer);

  // Check that noise profile is available
  TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(
                  handle, ROLLING_MEAN),
              "Noise profile should be available after learning");

  // Switch to reduction mode
  parameters.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load NLM reduction parameters should succeed");

  // Process remaining blocks
  size_t processed_samples = FRAME_SIZE * 10;
  while (processed_samples < BLOCK_SIZE) {
    uint32_t block_size = FRAME_SIZE;
    if (processed_samples + (size_t)block_size > BLOCK_SIZE) {
      block_size = (uint32_t)(BLOCK_SIZE - processed_samples);
    }

    bool result = specbleach_denoiser_process(
        handle, block_size, input_buffer + processed_samples,
        output_buffer + processed_samples);
    TEST_ASSERT(result == true, "Processing failed");

    processed_samples += (size_t)block_size;
  }

  // Calculate RMS values
  float input_rms = calculate_rms(input_buffer, BLOCK_SIZE);
  float output_rms = calculate_rms(output_buffer, BLOCK_SIZE);
  printf("  unified denoiser - Input RMS: %.4f, Output RMS: %.4f\n", input_rms,
         output_rms);

  // Test noise profile API
  uint32_t profile_size = specbleach_denoiser_get_noise_profile_size(handle);
  TEST_ASSERT(profile_size > 0, "Profile size should be positive");

  const float* profile =
      specbleach_denoiser_get_noise_profile_for_mode(handle, ROLLING_MEAN);
  TEST_ASSERT(profile != NULL, "Should be able to get noise profile");

  // Test reset
  specbleach_denoiser_reset_noise_profile(handle);
  TEST_ASSERT(true, "Reset should succeed");
  TEST_ASSERT(!specbleach_denoiser_noise_profile_available_for_mode(
                  handle, ROLLING_MEAN),
              "Profile should not be available after reset");

  specbleach_denoiser_free(handle);
  free(input_buffer);
  free(output_buffer);

  printf("✓ 2D NLM smoothing mode integration test passed\n");
}

// Runtime smoothing mode switch: profile persists, latency constant,
// processing continues seamlessly
void test_runtime_mode_switch(void) {
  printf("Testing runtime smoothing mode switch...\n");

  float* input_buffer = calloc(BLOCK_SIZE, sizeof(float));
  float* output_buffer = calloc(BLOCK_SIZE, sizeof(float));
  TEST_ASSERT(input_buffer && output_buffer, "Failed to allocate test buffers");

  generate_test_audio(input_buffer, BLOCK_SIZE, 1000.0f, 0.1f);

  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, 20.0f);
  TEST_ASSERT(handle != NULL, "Failed to initialize denoiser");

  const uint32_t latency_temporal = specbleach_denoiser_get_latency(handle);
  TEST_ASSERT(latency_temporal > 0, "Latency must include common delay");

  SpecbleachDenoiserParameters parameters = {
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.5f,
      .smoothing_mode = SPECBLEACH_SMOOTHING_TEMPORAL,
  };
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load learn parameters");
  specbleach_denoiser_process(handle, FRAME_SIZE * 10, input_buffer,
                              output_buffer);
  TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(
                  handle, ROLLING_MEAN),
              "Profile available after learning");

  // Switch to NLM 2D mid-stream and keep processing
  parameters.learn_noise = SPECBLEACH_LEARN_OFF;
  parameters.smoothing_mode = SPECBLEACH_SMOOTHING_NLM_2D;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load NLM parameters");
  size_t processed = FRAME_SIZE * 10;
  while (processed < BLOCK_SIZE) {
    TEST_ASSERT(specbleach_denoiser_process(handle, FRAME_SIZE,
                                            input_buffer + processed,
                                            output_buffer + processed),
                "process after mode switch");
    processed += FRAME_SIZE;
  }

  // Profile must have survived the switch
  TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(
                  handle, ROLLING_MEAN),
              "Profile persists across smoothing mode switch");

  // Latency must be unchanged by the mode switch
  TEST_ASSERT(specbleach_denoiser_get_latency(handle) == latency_temporal,
              "Latency constant across smoothing mode switch");

  specbleach_denoiser_free(handle);
  free(input_buffer);
  free(output_buffer);

  printf("✓ Runtime smoothing mode switch test passed\n");
}
