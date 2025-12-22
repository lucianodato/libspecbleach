/*
 * Audio regression tests - compare processed audio against reference files
 */

#include "../src/shared/configurations.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Include internal headers for testing

// Include the public API
#include "specbleach/specbleach_adenoiser.h"
#include "specbleach/specbleach_denoiser.h"

// Function prototypes
void generate_test_signal(float* buffer, int length, unsigned int seed);
void process_audio(const float* input, float* output, int length);
void process_audio_adaptive(const float* input, float* output, int length);
float calculate_snr(const float* original, const float* processed, int length);
void test_deterministic_processing(void);
void test_noise_reduction(void);
void test_adaptive_denoising(void);
void test_noise_estimation_methods(void);
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
#define BLOCK_SIZE FRAME_SIZE
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

  SpectralBleachParameters parameters = (SpectralBleachParameters){
      .learn_noise = 1,          // Learn all modes
      .noise_reduction_mode = 1, // Use average when processing
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

// Process audio through adaptive denoiser
void process_audio_adaptive(const float* input, float* output, int length) {
  float frame_size_ms = 20.0f;
  SpectralBleachHandle handle =
      specbleach_adaptive_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize adaptive denoiser");

  SpectralBleachParameters parameters =
      (SpectralBleachParameters){.reduction_amount = 20.0f,
                                 .smoothing_factor = 0.0f,
                                 .noise_rescale = 0.0f,
                                 .noise_scaling_type = 0,
                                 .post_filter_threshold = 0.0f,
                                 .residual_listen = false,
                                 .transient_protection = false,
                                 .whitening_factor = 0.0f};

  specbleach_adaptive_load_parameters(handle, parameters);

  int processed = 0;
  while (processed < length) {
    int block_size = FRAME_SIZE;
    if (processed + block_size > length) {
      block_size = length - processed;
    }

    bool result = specbleach_adaptive_process(
        handle, block_size, input + processed, output + processed);
    TEST_ASSERT(result == true, "Adaptive processing failed");

    processed += block_size;
  }

  specbleach_adaptive_free(handle);
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
void test_deterministic_processing(void) {
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
void test_noise_reduction(void) {
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
void test_valid_output(void) {
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

// Test that adaptive denoiser works and is different from static
void test_adaptive_denoising(void) {
  printf("Testing adaptive denoiser effectiveness...\n");

  float* input = calloc(TEST_SAMPLES, sizeof(float));
  float* output_static = calloc(TEST_SAMPLES, sizeof(float));
  float* output_adaptive = calloc(TEST_SAMPLES, sizeof(float));
  TEST_ASSERT(input && output_static && output_adaptive,
              "Failed to allocate test buffers");

  // Generate noisy signal
  generate_test_signal(input, TEST_SAMPLES, 99999);

  // Process through both denoisers
  process_audio(input, output_static, TEST_SAMPLES);
  process_audio_adaptive(input, output_adaptive, TEST_SAMPLES);

  // Verify adaptive denoiser reduced noise
  double input_power = 0.0;
  double adaptive_output_power = 0.0;
  for (int i = 0; i < TEST_SAMPLES; i++) {
    input_power += input[i] * input[i];
    adaptive_output_power += output_adaptive[i] * output_adaptive[i];
  }
  input_power /= TEST_SAMPLES;
  adaptive_output_power /= TEST_SAMPLES;

  printf("  Input power: %.6f\n", input_power);
  printf("  Adaptive output power: %.6f\n", adaptive_output_power);

  TEST_ASSERT(adaptive_output_power < input_power * 0.95f,
              "Adaptive denoising should reduce signal power");

  // Verify they are NOT identical (especially since they use different
  // overlap/windows by default)
  bool identical = true;
  for (int i = 5000; i < TEST_SAMPLES; i++) { // Skip initial frames
    if (fabsf(output_static[i] - output_adaptive[i]) > 1e-4f) {
      identical = false;
      break;
    }
  }
  TEST_ASSERT(!identical,
              "Adaptive and Static denoisers should not be identical");

  free(input);
  free(output_static);
  free(output_adaptive);

  printf("✓ Adaptive denoiser test passed\n");
}

// Test that both noise estimation methods work correctly
void test_noise_estimation_methods(void) {
  printf("Testing noise estimation methods (Louizou vs SPP-MMSE)...\n");

  float frame_size_ms = 20.0f;
  float* input = calloc(TEST_SAMPLES, sizeof(float));
  float* output_louizou = calloc(TEST_SAMPLES, sizeof(float));
  float* output_spp_mmse = calloc(TEST_SAMPLES, sizeof(float));
  TEST_ASSERT(input && output_louizou && output_spp_mmse,
              "Failed to allocate test buffers");

  // Generate noisy signal
  generate_test_signal(input, TEST_SAMPLES, 12345);

  // Process with Louizou method (default)
  SpectralBleachParameters params_louizou =
      (SpectralBleachParameters){.reduction_amount = 20.0f,
                                 .smoothing_factor = 0.0f,
                                 .noise_rescale = 0.0f,
                                 .noise_scaling_type = 0,
                                 .post_filter_threshold = 0.0f,
                                 .residual_listen = false,
                                 .transient_protection = false,
                                 .whitening_factor = 0.0f,
                                 .noise_estimation_method = LOUIZOU_METHOD};

  SpectralBleachHandle handle_louizou =
      specbleach_adaptive_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle_louizou != NULL, "Failed to initialize Louizou denoiser");

  TEST_ASSERT(
      specbleach_adaptive_load_parameters(handle_louizou, params_louizou),
      "Failed to load Louizou parameters");

  for (int i = 0; i < TEST_SAMPLES; i += BLOCK_SIZE) {
    int block_size =
        (i + BLOCK_SIZE > TEST_SAMPLES) ? TEST_SAMPLES - i : BLOCK_SIZE;
    TEST_ASSERT(specbleach_adaptive_process(handle_louizou, block_size,
                                            input + i, output_louizou + i),
                "Failed to process with Louizou method");
  }

  specbleach_adaptive_free(handle_louizou);

  // Process with SPP-MMSE method
  SpectralBleachParameters params_spp_mmse =
      (SpectralBleachParameters){.reduction_amount = 20.0f,
                                 .smoothing_factor = 0.0f,
                                 .noise_rescale = 0.0f,
                                 .noise_scaling_type = 0,
                                 .post_filter_threshold = 0.0f,
                                 .residual_listen = false,
                                 .transient_protection = false,
                                 .whitening_factor = 0.0f,
                                 .noise_estimation_method = SPP_MMSE_METHOD};

  SpectralBleachHandle handle_spp_mmse =
      specbleach_adaptive_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle_spp_mmse != NULL,
              "Failed to initialize SPP-MMSE denoiser");

  TEST_ASSERT(
      specbleach_adaptive_load_parameters(handle_spp_mmse, params_spp_mmse),
      "Failed to load SPP-MMSE parameters");

  for (int i = 0; i < TEST_SAMPLES; i += BLOCK_SIZE) {
    int block_size =
        (i + BLOCK_SIZE > TEST_SAMPLES) ? TEST_SAMPLES - i : BLOCK_SIZE;
    TEST_ASSERT(specbleach_adaptive_process(handle_spp_mmse, block_size,
                                            input + i, output_spp_mmse + i),
                "Failed to process with SPP-MMSE method");
  }

  specbleach_adaptive_free(handle_spp_mmse);

  // Verify both methods produced valid output (finite values, reduced noise)
  double input_power = 0.0, louizou_power = 0.0, spp_mmse_power = 0.0;
  for (int i = 0; i < TEST_SAMPLES; i++) {
    TEST_ASSERT(isfinite(output_louizou[i]),
                "Louizou output contains non-finite values");
    TEST_ASSERT(isfinite(output_spp_mmse[i]),
                "SPP-MMSE output contains non-finite values");

    input_power += input[i] * input[i];
    louizou_power += output_louizou[i] * output_louizou[i];
    spp_mmse_power += output_spp_mmse[i] * output_spp_mmse[i];
  }

  input_power /= TEST_SAMPLES;
  louizou_power /= TEST_SAMPLES;
  spp_mmse_power /= TEST_SAMPLES;

  printf("  Input power: %.6f\n", input_power);
  printf("  Louizou output power: %.6f\n", louizou_power);
  printf("  SPP-MMSE output power: %.6f\n", spp_mmse_power);

  // Both methods should reduce noise
  TEST_ASSERT(louizou_power < input_power * 0.95f,
              "Louizou method should reduce signal power");
  TEST_ASSERT(spp_mmse_power < input_power * 0.95f,
              "SPP-MMSE method should reduce signal power");

  // Methods should produce different results (they use different algorithms)
  bool identical = true;
  for (int i = 5000; i < TEST_SAMPLES; i++) { // Skip initial frames
    if (fabsf(output_louizou[i] - output_spp_mmse[i]) > 1e-4f) {
      identical = false;
      break;
    }
  }
  TEST_ASSERT(!identical,
              "Louizou and SPP-MMSE methods should produce different results");

  free(input);
  free(output_louizou);
  free(output_spp_mmse);

  printf("✓ Noise estimation methods test passed\n");
}

int main(void) {
  printf("Running audio regression tests...\n\n");

  test_deterministic_processing();
  test_noise_reduction();
  test_valid_output();
  test_adaptive_denoising();
  test_noise_estimation_methods();

  printf("\n✅ All audio regression tests passed!\n");
  return 0;
}
