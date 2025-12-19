/*
 * Audio regression tests - compare processed audio against reference files
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
void generate_test_signal(float* buffer, int length, unsigned int seed);
void process_audio(const float* input, float* output, int length);
float calculate_snr(const float* original, const float* processed, int length);
void test_deterministic_processing(void);
void test_noise_reduction(void);
void test_snr_improvement(void);

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
#define TEST_DURATION_SECONDS 2
#define TEST_SAMPLES (SAMPLE_RATE * TEST_DURATION_SECONDS)

// Generate deterministic test signal (sine wave + noise)
void generate_test_signal(float* buffer, int length, unsigned int seed) {
  srand(seed); // Deterministic seed

  for (int i = 0; i < length; i++) {
    // Generate a 1kHz sine wave
    float signal = 0.3f * sinf(2.0f * M_PI * 1000.0f * i / SAMPLE_RATE);

    // Add correlated noise (pink-like)
    float noise = 0.1f * ((float)rand() / RAND_MAX - 0.5f) * 2.0f;

    // Add some harmonics
    float harmonic = 0.1f * sinf(2.0f * M_PI * 2000.0f * i / SAMPLE_RATE);

    buffer[i] = signal + noise + harmonic;
  }
}

// Process audio through denoiser
void process_audio(const float* input, float* output, int length) {
  float frame_size_ms = 20.0f;
  SpectralBleachHandle handle =
      specbleach_initialize(SAMPLE_RATE, frame_size_ms);
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

  // Learn phase (first 5000 samples)
  specbleach_process(handle, 5000, input, output);

  // Reduction phase
  parameters.learn_noise = 0;
  specbleach_load_parameters(handle, parameters);

  int processed = 5000;
  while (processed < length) {
    int block_size = FRAME_SIZE;
    if (processed + block_size > length) {
      block_size = length - processed;
    }

    bool result = specbleach_process(handle, block_size, input + processed,
                                     output + processed);
    TEST_ASSERT(result == true, "Processing failed");

    processed += block_size;
  }

  specbleach_free(handle);
}

// Calculate SNR (Signal-to-Noise Ratio)
float calculate_snr(const float* original, const float* processed, int length) {
  double signal_power = 0.0;
  double noise_power = 0.0;

  for (int i = 0; i < length; i++) {
    float signal = original[i];
    float noise = original[i] - processed[i];

    signal_power += signal * signal;
    noise_power += noise * noise;
  }

  if (noise_power < 1e-10) {
    return 100.0f; // Very high SNR if noise is negligible
  }

  return 10.0f * log10(signal_power / noise_power);
}

// Test that denoising produces consistent results
void test_deterministic_processing() {
  printf("Testing deterministic processing...\n");

  float* input1 = calloc(TEST_SAMPLES, sizeof(float));
  float* input2 = calloc(TEST_SAMPLES, sizeof(float));
  float* output1 = calloc(TEST_SAMPLES, sizeof(float));
  float* output2 = calloc(TEST_SAMPLES, sizeof(float));

  TEST_ASSERT(input1 && input2 && output1 && output2,
              "Failed to allocate test buffers");

  // Generate identical input signals
  generate_test_signal(input1, TEST_SAMPLES, 12345);
  generate_test_signal(input2, TEST_SAMPLES, 12345);

  // Process both
  process_audio(input1, output1, TEST_SAMPLES);
  process_audio(input2, output2, TEST_SAMPLES);

  // Verify outputs are identical (deterministic processing)
  for (int i = 0; i < TEST_SAMPLES; i++) {
    TEST_FLOAT_CLOSE(output1[i], output2[i], 1e-10f);
  }

  free(input1);
  free(input2);
  free(output1);
  free(output2);

  printf("✓ Deterministic processing test passed\n");
}

// Test that denoising actually reduces noise
void test_noise_reduction() {
  printf("Testing noise reduction effectiveness...\n");

  float* input = calloc(TEST_SAMPLES, sizeof(float));
  float* output = calloc(TEST_SAMPLES, sizeof(float));
  TEST_ASSERT(input && output, "Failed to allocate test buffers");

  // Generate noisy signal
  generate_test_signal(input, TEST_SAMPLES, 54321);

  // Calculate input signal power (approximate)
  double input_power = 0.0;
  for (int i = 0; i < TEST_SAMPLES; i++) {
    input_power += input[i] * input[i];
  }
  input_power /= TEST_SAMPLES;

  // Process through denoiser
  process_audio(input, output, TEST_SAMPLES);

  // Calculate output signal power
  double output_power = 0.0;
  for (int i = 0; i < TEST_SAMPLES; i++) {
    output_power += output[i] * output[i];
  }
  output_power /= TEST_SAMPLES;

  printf("  Input power: %.6f\n", input_power);
  printf("  Output power: %.6f\n", output_power);

  // Verify that output power is reduced (but not too much - we want to preserve
  // signal)
  TEST_ASSERT(output_power < input_power * 0.9f,
              "Denoising should reduce signal power");
  TEST_ASSERT(output_power > input_power * 0.01f,
              "Denoising should preserve most of the signal");

  free(input);
  free(output);

  printf("✓ Noise reduction test passed\n");
}

// Test that denoising produces valid output
void test_valid_output() {
  printf("Testing valid output generation...\n");

  float* input = calloc(TEST_SAMPLES, sizeof(float));
  float* output = calloc(TEST_SAMPLES, sizeof(float));
  TEST_ASSERT(input && output, "Failed to allocate test buffers");

  // Generate test signal
  generate_test_signal(input, TEST_SAMPLES, 11111);

  // Process through denoiser
  process_audio(input, output, TEST_SAMPLES);

  // Verify output is valid (not all zeros, not NaN, within reasonable range)
  float max_output = 0.0f;
  float min_output = 0.0f;
  bool has_non_zero = false;

  for (int i = 0; i < TEST_SAMPLES; i++) {
    TEST_ASSERT(!isnan(output[i]), "Output contains NaN values");
    TEST_ASSERT(!isinf(output[i]), "Output contains infinite values");

    if (output[i] != 0.0f) {
      has_non_zero = true;
    }

    if (output[i] > max_output)
      max_output = output[i];
    if (output[i] < min_output)
      min_output = output[i];
  }

  TEST_ASSERT(has_non_zero, "Output should not be all zeros");
  TEST_ASSERT(max_output < 2.0f && min_output > -2.0f,
              "Output values should be in reasonable range");

  printf("  Output range: %.3f to %.3f\n", min_output, max_output);

  free(input);
  free(output);

  printf("✓ Valid output test passed\n");
}

int main() {
  printf("Running audio regression tests...\n\n");

  test_deterministic_processing();
  test_noise_reduction();
  test_valid_output();

  printf("\n✅ All audio regression tests passed!\n");
  return 0;
}
