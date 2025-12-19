/*
 * Integration tests for the full denoising pipeline
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Include internal headers for testing
#include "../src/shared/configurations.h"

// Include the public API
#include "specbleach/specbleach_denoiser.h"

// Function prototypes
void test_spectral_denoiser(void);
void test_different_noise_levels(void);
void test_library_info(void);
float calculate_rms(const float* buffer, int length);

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
#define BLOCK_SIZE (FRAME_SIZE * NUM_FRAMES)

// Generate synthetic audio with noise
void generate_test_audio(float* buffer, int length, float signal_freq,
                         float noise_level) {
  for (int i = 0; i < length; i++) {
    // Generate a sine wave signal
    float signal = 0.5f * sinf(2.0f * M_PI * signal_freq * i / SAMPLE_RATE);

    // Add white noise
    float noise = noise_level * ((float)rand() / RAND_MAX - 0.5f) * 2.0f;

    buffer[i] = signal + noise;
  }
}

// Calculate RMS of audio buffer
float calculate_rms(const float* buffer, int length) {
  double sum = 0.0;
  for (int i = 0; i < length; i++) {
    sum += buffer[i] * buffer[i];
  }
  return sqrtf(sum / length);
}

// Test spectral denoiser with synthetic audio
void test_spectral_denoiser() {
  printf("Testing spectral denoiser integration...\n");

  // Create test audio with signal + noise
  float* input_buffer = calloc(BLOCK_SIZE, sizeof(float));
  float* output_buffer = calloc(BLOCK_SIZE, sizeof(float));
  TEST_ASSERT(input_buffer && output_buffer, "Failed to allocate test buffers");

  generate_test_audio(input_buffer, BLOCK_SIZE, 1000.0f,
                      0.1f); // 1kHz signal + noise

  // Calculate input RMS
  float input_rms = calculate_rms(input_buffer, BLOCK_SIZE);
  printf("  Input RMS: %.4f\n", input_rms);

  // Initialize denoiser
  SpectralBleachHandle handle = specbleach_initialize(SAMPLE_RATE, FRAME_SIZE);
  TEST_ASSERT(handle != NULL, "Failed to initialize spectral denoiser");

  SpectralBleachParameters parameters =
      (SpectralBleachParameters){.learn_noise = 1, // Average
                                 .reduction_amount = 20.0f,
                                 .smoothing_factor = 0.0f,
                                 .noise_rescale = 0.0f,
                                 .noise_scaling_type = 0,
                                 .post_filter_threshold = 0.0f,
                                 .residual_listen = false,
                                 .transient_protection = false,
                                 .whitening_factor = 0.0f};

  specbleach_load_parameters(handle, parameters);

  // Process first block in learn mode
  specbleach_process(handle, FRAME_SIZE, input_buffer, output_buffer);

  // Switch to reduction mode
  parameters.learn_noise = 0;
  specbleach_load_parameters(handle, parameters);

  // Process remaining blocks
  int processed_samples = FRAME_SIZE;
  while (processed_samples < BLOCK_SIZE) {
    int block_size = FRAME_SIZE;
    if (processed_samples + block_size > BLOCK_SIZE) {
      block_size = BLOCK_SIZE - processed_samples;
    }

    bool result =
        specbleach_process(handle, block_size, input_buffer + processed_samples,
                           output_buffer + processed_samples);
    TEST_ASSERT(result == true, "Processing failed");

    processed_samples += block_size;
  }

  // Calculate output RMS
  float output_rms = calculate_rms(output_buffer, BLOCK_SIZE);
  printf("  Output RMS: %.4f\n", output_rms);

  // Verify noise reduction (output should be quieter than input)
  TEST_ASSERT(output_rms < input_rms * 0.8f,
              "Noise reduction should reduce RMS significantly");

  // Verify signal preservation (output shouldn't be too quiet)
  TEST_ASSERT(output_rms > input_rms * 0.1f,
              "Signal should be reasonably preserved");

  // Cleanup
  specbleach_free(handle);
  free(input_buffer);
  free(output_buffer);

  printf("✓ Spectral denoiser integration test passed\n");
}

// Test denoiser with different noise levels
void test_different_noise_levels() {
  printf("Testing denoiser with different noise levels...\n");

  // Test with low noise
  float* input_buffer = calloc(BLOCK_SIZE, sizeof(float));
  float* output_buffer = calloc(BLOCK_SIZE, sizeof(float));
  TEST_ASSERT(input_buffer && output_buffer, "Failed to allocate test buffers");

  generate_test_audio(input_buffer, BLOCK_SIZE, 1000.0f, 0.05f); // Low noise

  SpectralBleachHandle handle = specbleach_initialize(SAMPLE_RATE, FRAME_SIZE);
  TEST_ASSERT(handle != NULL, "Failed to initialize denoiser");

  SpectralBleachParameters parameters =
      (SpectralBleachParameters){.learn_noise = 1, // Average
                                 .reduction_amount = 20.0f,
                                 .smoothing_factor = 0.0f,
                                 .noise_rescale = 0.0f,
                                 .noise_scaling_type = 0,
                                 .post_filter_threshold = 0.0f,
                                 .residual_listen = false,
                                 .transient_protection = false,
                                 .whitening_factor = 0.0f};

  specbleach_load_parameters(handle, parameters);

  // Learn a bit
  specbleach_process(handle, 1000, input_buffer, output_buffer);

  // Process remaining with reduction
  parameters.learn_noise = 0;
  specbleach_load_parameters(handle, parameters);

  // Process all at once for simplicity
  bool result = specbleach_process(handle, BLOCK_SIZE - 1000,
                                   input_buffer + 1000, output_buffer + 1000);
  TEST_ASSERT(result == true, "Processing failed");

  float input_rms = calculate_rms(input_buffer, BLOCK_SIZE);
  float output_rms = calculate_rms(output_buffer, BLOCK_SIZE);

  printf("  Low noise - Input RMS: %.4f, Output RMS: %.4f\n", input_rms,
         output_rms);

  specbleach_free(handle);
  free(input_buffer);
  free(output_buffer);

  printf("✓ Different noise levels test passed\n");
}

// Test library information functions
void test_library_info() {
  printf("Testing library information functions...\n");

  // Test spectral denoiser info
  int latency = specbleach_get_latency(NULL); // Should return 0 for NULL handle
  TEST_ASSERT(latency == 0, "NULL handle should return 0 latency");

  int profile_size = specbleach_get_noise_profile_size(NULL);
  TEST_ASSERT(profile_size == 0, "NULL handle should return 0 profile size");

  // Test that we can get valid information after initialization
  float frame_size_ms = 20.0f;
  SpectralBleachHandle handle =
      specbleach_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize for info test");

  latency = specbleach_get_latency(handle);
  TEST_ASSERT(latency >= 0, "Latency should be non-negative");

  profile_size = specbleach_get_noise_profile_size(handle);
  TEST_ASSERT(profile_size > 0, "Profile size should be positive");

  specbleach_free(handle);

  printf("✓ Library information tests passed\n");
}

int main() {
  printf("Running integration tests...\n\n");

  // Seed random number generator for reproducible tests
  srand(42);

  test_spectral_denoiser();
  test_different_noise_levels();
  test_library_info();

  printf("\n✅ All integration tests passed!\n");
  return 0;
}
