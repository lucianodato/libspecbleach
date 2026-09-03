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

#include "specbleach_denoiser.h"

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
void test_frame_size_invariance(void);

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
#define TEST_SAMPLES ((size_t)SAMPLE_RATE * (size_t)TEST_DURATION_SECONDS)

// Generate deterministic test signal (sine wave + noise)
void generate_test_signal(float* buffer, int length, unsigned int seed) {
  srand(seed); // Deterministic seed

  for (int i = 0; i < length; i++) {
    // Generate a 1kHz sine wave
    float signal =
        0.3f *
        sinf((float)(2.0 * M_PIf * 1000.0 * (double)i / (double)SAMPLE_RATE));

    // Add correlated noise (pink-like)
    float noise = 0.1f * ((float)rand() / RAND_MAX - 0.5f) * 2.0f;

    // Add some harmonics
    float harmonic =
        0.1f *
        sinf((float)(2.0 * M_PIf * 2000.0 * (double)i / (double)SAMPLE_RATE));

    buffer[i] = signal + noise + harmonic;
  }
}

// Process audio through denoiser
void process_audio(const float* input, float* output, int length) {
  float frame_size_ms = 20.0f;
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize denoiser");

  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_ALL, // Learn all modes
      .tonal_reduction_gain = 1.0f,
      .aggressiveness = -1.0f, // Use median when processing
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.0f,
      .masking_depth = 0.5f,

      .residual_listen = false,
      .whitening_factor = 0.0f};

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load spectral parameters should succeed");

  // Learn phase (first 5000 samples)
  specbleach_denoiser_process(handle, 5000, input, output);

  // Reduction phase
  parameters.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load reduction parameters should succeed");

  int processed = 5000;
  while (processed < length) {
    int block_size = FRAME_SIZE;
    if (processed + block_size > length) {
      block_size = length - processed;
    }

    bool result = specbleach_denoiser_process(
        handle, block_size, input + processed, output + processed);
    TEST_ASSERT(result == true, "Processing failed");

    processed += block_size;
  }

  specbleach_denoiser_free(handle);
}

// Process audio through adaptive denoiser
void process_audio_adaptive(const float* input, float* output, int length) {
  float frame_size_ms = 20.0f;
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize adaptive denoiser");

  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.0f,
      .masking_depth = 0.5f,

      .residual_listen = false,
      .whitening_factor = 0.0f,
      .adaptive_noise = true,
      .noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE};

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load adaptive parameters should succeed");

  int processed = 0;
  while (processed < length) {
    int block_size = FRAME_SIZE;
    if (processed + block_size > length) {
      block_size = length - processed;
    }

    bool result = specbleach_denoiser_process(
        handle, block_size, input + processed, output + processed);
    TEST_ASSERT(result == true, "Adaptive processing failed");

    processed += block_size;
  }

  specbleach_denoiser_free(handle);
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

  return (float)(10.0 * log10(signal_power / noise_power));
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
  for (size_t i = 0; i < TEST_SAMPLES; i++) {
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
  for (size_t i = 0; i < TEST_SAMPLES; i++) {
    input_power += (double)input[i] * (double)input[i];
  }
  input_power /= (double)TEST_SAMPLES;

  // Process through denoiser
  process_audio(input, output, TEST_SAMPLES);

  // Calculate output signal power
  double output_power = 0.0;
  for (size_t i = 0; i < TEST_SAMPLES; i++) {
    output_power += (double)output[i] * (double)output[i];
  }
  output_power /= (double)TEST_SAMPLES;

  printf("  Input power: %.6f\n", input_power);
  printf("  Output power: %.6f\n", output_power);

  // Verify that output power is reduced (but not too much - we want to preserve
  // signal)
  // Ratio was ~0.95 under fixed-frame tuning; ms-anchored estimators (#152)
  // shift the 20ms operating point slightly, so bound is 0.99.
  TEST_ASSERT(output_power < input_power * 0.99f,
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

  for (size_t i = 0; i < TEST_SAMPLES; i++) {
    TEST_ASSERT(!isnan(output[i]), "Output contains NaN values");
    TEST_ASSERT(!isinf(output[i]), "Output contains infinite values");

    if (output[i] != 0.0f) {
      has_non_zero = true;
    }

    if (output[i] > max_output) {
      max_output = output[i];
    }
    if (output[i] < min_output) {
      min_output = output[i];
    }
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
  for (size_t i = 0; i < TEST_SAMPLES; i++) {
    input_power += (double)input[i] * (double)input[i];
    adaptive_output_power +=
        (double)output_adaptive[i] * (double)output_adaptive[i];
  }
  input_power /= (double)TEST_SAMPLES;
  adaptive_output_power /= (double)TEST_SAMPLES;

  printf("  Input power: %.6f\n", input_power);
  printf("  Adaptive output power: %.6f\n", adaptive_output_power);

  TEST_ASSERT(adaptive_output_power < input_power * 0.98f,
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
  printf("Testing noise estimation methods (Martin MS vs SPP-MMSE)...\n");

  float frame_size_ms = 20.0f;
  float* input = calloc(TEST_SAMPLES, sizeof(float));
  float* output_martin = calloc(TEST_SAMPLES, sizeof(float));
  float* output_spp_mmse = calloc(TEST_SAMPLES, sizeof(float));
  TEST_ASSERT(input && output_martin && output_spp_mmse,
              "Failed to allocate test buffers");

  // Generate noisy signal
  generate_test_signal(input, TEST_SAMPLES, 12345);

  // Process with Martin MS method (Default)
  SpecbleachDenoiserParameters params_martin = (SpecbleachDenoiserParameters){
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.0f,
      .masking_depth = 0.5f,

      .residual_listen = false,
      .whitening_factor = 0.0f,
      .adaptive_noise = true,
      .noise_estimation_method =
          SPECBLEACH_NOISE_ESTIMATION_MARTIN}; // 2: Martin MS

  specbleach_denoiser* handle_martin =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle_martin != NULL, "Failed to initialize Martin denoiser");

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle_martin, &params_martin,
                                                  sizeof(params_martin)),
              "Load Martin parameters should succeed");

  for (size_t i = 0; i < TEST_SAMPLES; i += (size_t)BLOCK_SIZE) {
    int block_size = (i + (size_t)BLOCK_SIZE > TEST_SAMPLES)
                         ? (int)(TEST_SAMPLES - i)
                         : BLOCK_SIZE;
    TEST_ASSERT(specbleach_denoiser_process(handle_martin, block_size,
                                            input + i, output_martin + i),
                "Failed to process with Martin method");
  }

  specbleach_denoiser_free(handle_martin);

  // Process with SPP-MMSE method
  SpecbleachDenoiserParameters params_spp_mmse = (SpecbleachDenoiserParameters){
      .reduction_gain = 0.1f,
      .smoothing_factor = 0.0f,
      .masking_depth = 0.5f,

      .residual_listen = false,
      .whitening_factor = 0.0f,
      .adaptive_noise = true,
      .noise_estimation_method =
          SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE}; // 0: SPP-MMSE

  specbleach_denoiser* handle_spp_mmse =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle_spp_mmse != NULL,
              "Failed to initialize SPP-MMSE denoiser");

  TEST_ASSERT(specbleach_denoiser_load_parameters(
                  handle_spp_mmse, &params_spp_mmse, sizeof(params_spp_mmse)),
              "Load SPP-MMSE parameters should succeed");

  for (size_t i = 0; i < TEST_SAMPLES; i += (size_t)BLOCK_SIZE) {
    int block_size = (i + (size_t)BLOCK_SIZE > TEST_SAMPLES)
                         ? (int)(TEST_SAMPLES - i)
                         : BLOCK_SIZE;
    TEST_ASSERT(specbleach_denoiser_process(handle_spp_mmse, block_size,
                                            input + i, output_spp_mmse + i),
                "Failed to process with SPP-MMSE method");
  }

  specbleach_denoiser_free(handle_spp_mmse);

  // Verify both methods produced valid output (finite values, reduced noise)
  double input_power = 0.0;
  double martin_power = 0.0;
  double spp_mmse_power = 0.0;
  for (size_t i = 0; i < TEST_SAMPLES; i++) {
    TEST_ASSERT(isfinite(output_martin[i]),
                "Martin output contains non-finite values");
    TEST_ASSERT(isfinite(output_spp_mmse[i]),
                "SPP-MMSE output contains non-finite values");

    input_power += (double)input[i] * (double)input[i];
    martin_power += (double)output_martin[i] * (double)output_martin[i];
    spp_mmse_power += (double)output_spp_mmse[i] * (double)output_spp_mmse[i];
  }

  input_power /= (double)TEST_SAMPLES;
  martin_power /= (double)TEST_SAMPLES;
  spp_mmse_power /= (double)TEST_SAMPLES;

  printf("  Input power: %.6f\n", input_power);
  printf("  Martin output power: %.6f\n", martin_power);
  printf("  SPP-MMSE output power: %.6f\n", spp_mmse_power);

  // Both methods should reduce noise
  TEST_ASSERT(martin_power < input_power * 0.98f,
              "Martin method should reduce signal power");
  TEST_ASSERT(spp_mmse_power < input_power * 0.98f,
              "SPP-MMSE method should reduce signal power");

  // Methods should produce different results
  bool identical = true;
  for (int i = 5000; i < TEST_SAMPLES; i++) { // Skip initial frames
    if (fabsf(output_martin[i] - output_spp_mmse[i]) > 1e-4f) {
      identical = false;
      break;
    }
  }
  TEST_ASSERT(!identical,
              "Martin and SPP-MMSE methods should produce different results");

  free(input);
  free(output_martin);
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
  test_frame_size_invariance();

  printf("\n✅ All audio regression tests passed!\n");
  return 0;
}

/* Frame-size invariance (#152): the same transient material processed at
 * 23ms and 93ms must not diverge in time fidelity beyond hop granularity.
 * NLM geometry is fixed in ms (patch 92ms / symmetric 128ms search), so
 * both runs apply the same real-time smear; this test fails if geometry is
 * counted in raw frames. Each run is aligned by its own reported latency
 * before comparing. Pre-echo (future smear) and tail (past smear) around
 * onsets are the discriminating metrics. */
static void process_transient_nlm(const float* input, float* output, int length,
                                  float frame_size_ms, uint32_t* latency_out) {
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, frame_size_ms);
  TEST_ASSERT(handle != NULL, "Failed to initialize denoiser");

  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .tonal_reduction_gain = 0.0f,
      .aggressiveness = 0.0f,
      .reduction_gain = 0.3f,
      .smoothing_factor = 0.5f,
      .smoothing_mode = SPECBLEACH_SMOOTHING_NLM_2D,
      .masking_depth = 0.5f,
      .residual_listen = false,
      .whitening_factor = 0.0f};

  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load learn parameters should succeed");

  /* Learn on noise-only head (first 0.5s), then denoise transients. */
  const int learn_samples = SAMPLE_RATE / 2;
  specbleach_denoiser_process(handle, learn_samples, input, output);

  parameters.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "Load reduction parameters should succeed");

  int processed = learn_samples;
  while (processed < length) {
    int block_size = FRAME_SIZE;
    if (processed + block_size > length) {
      block_size = length - processed;
    }
    TEST_ASSERT(specbleach_denoiser_process(
                    handle, block_size, input + processed, output + processed),
                "Processing failed");
    processed += block_size;
  }

  if (latency_out) {
    *latency_out = specbleach_denoiser_get_latency(handle);
  }
  specbleach_denoiser_free(handle);
}

void test_frame_size_invariance(void) {
  printf("Testing frame-size invariance (23ms vs 93ms, NLM)...\n");

  /* Longer window than the shared 2s fixtures: several isolated onsets. */
  const int len = SAMPLE_RATE * 4;
  float* input = calloc(len, sizeof(float));
  float* out_23 = calloc(len, sizeof(float));
  float* out_93 = calloc(len, sizeof(float));
  TEST_ASSERT(input && out_23 && out_93, "Failed to allocate test buffers");

  /* Noise floor + transient impulses every 0.6s (alternating polarity).
   * Spacing exceeds the 400ms tail window so per-onset smear windows never
   * overlap. */
  srand(1234);
  for (int i = 0; i < len; i++) {
    input[i] = 0.05f * ((float)rand() / RAND_MAX - 0.5f) * 2.0f;
  }
  for (int t = SAMPLE_RATE / 2; t < len; t += SAMPLE_RATE * 6 / 10) {
    input[t] = (t % (SAMPLE_RATE / 2) == 0) ? 0.9f : -0.9f;
    if (t + 1 < len) {
      input[t + 1] = -input[t] * 0.5f;
    }
  }

  uint32_t lat_23 = 0, lat_93 = 0;
  process_transient_nlm(input, out_23, len, 23.0f, &lat_23);
  process_transient_nlm(input, out_93, len, 93.0f, &lat_93);

  /* Sanity: larger frame => larger latency, both finite. */
  TEST_ASSERT(lat_93 > lat_23, "93ms latency should exceed 23ms latency");
  for (int i = 0; i < len; i++) {
    TEST_ASSERT(isfinite(out_23[i]), "23ms output must be finite");
    TEST_ASSERT(isfinite(out_93[i]), "93ms output must be finite");
  }

  /* Latency-corrected, input-time-aligned comparison. Skip the longest
   * latency plus a settle margin so all smoothers are in steady state. */
  const int skip = (int)lat_93 + SAMPLE_RATE / 2;
  const int hop_23 = (int)(0.023f * SAMPLE_RATE / 4);
  const int hop_93 = (int)(0.093f * SAMPLE_RATE / 4);
  double in_pow = 0.0, p23 = 0.0, p93 = 0.0, diff_pow = 0.0;
  double peak_23 = 0.0, peak_93 = 0.0;
  /* Pre-echo (60ms before onset, excl. 1 hop) and tail (up to 150ms after
   * onset, excl. 2 hops) energies, normalized per onset by peak energy. */
  double pre_23 = 0.0, pre_93 = 0.0, tail_23 = 0.0, tail_93 = 0.0;
  double pk_e_23 = 0.0, pk_e_93 = 0.0;
  int n = 0, onsets = 0;
  /* Aligned walk: input time i maps to out_23[i + lat_23], out_93[i + lat_93].
   * Stay within bounds for both. */
  const int start = skip;
  const int end = len - (int)lat_93 - 1;
  for (int i = start; i < end; i++) {
    double a23 = out_23[i + lat_23];
    double a93 = out_93[i + lat_93];
    in_pow += (double)input[i] * input[i];
    p23 += a23 * a23;
    p93 += a93 * a93;
    double d = a23 - a93;
    diff_pow += d * d;
    if (fabs(a23) > peak_23) {
      peak_23 = fabs(a23);
    }
    if (fabs(a93) > peak_93) {
      peak_93 = fabs(a93);
    }
  }
  n = end - start;
  TEST_ASSERT(n > 0, "Comparison window must be non-empty");

  /* Per-onset smear metrics. Windows cover the widest legacy smear
   * (93ms future / 372ms past at 93ms frame under frame-counted geometry)
   * so a regression cannot hide outside the window. */
  for (int t = SAMPLE_RATE / 2; t < len - SAMPLE_RATE * 6 / 10;
       t += SAMPLE_RATE * 6 / 10) {
    if (t < start || t >= end) {
      continue;
    }
    onsets++;
    /* Peak energy near onset (within 2 hops). */
    double pk23 = 0.0, pk93 = 0.0;
    for (int j = t; j < t + 2 * hop_23 && j < end; j++) {
      double v = out_23[j + lat_23];
      if (v * v > pk23) {
        pk23 = v * v;
      }
    }
    for (int j = t; j < t + 2 * hop_93 && j < end; j++) {
      double v = out_93[j + lat_93];
      if (v * v > pk93) {
        pk93 = v * v;
      }
    }
    pk_e_23 += pk23;
    pk_e_93 += pk93;
    /* Pre-echo: 150ms before onset, excluding the hop right before it. */
    const int pre_len = (int)(0.150f * SAMPLE_RATE);
    for (int j = t - pre_len; j < t - hop_23; j++) {
      if (j >= start) {
        double v = out_23[j + lat_23];
        pre_23 += v * v;
      }
    }
    for (int j = t - pre_len; j < t - hop_93; j++) {
      if (j >= start) {
        double v = out_93[j + lat_93];
        pre_93 += v * v;
      }
    }
    /* Tail: 400ms after onset, excluding 2 hops right after it. */
    const int tail_len = (int)(0.400f * SAMPLE_RATE);
    for (int j = t + 2 * hop_23; j < t + tail_len && j < end; j++) {
      double v = out_23[j + lat_23];
      tail_23 += v * v;
    }
    for (int j = t + 2 * hop_93; j < t + tail_len && j < end; j++) {
      double v = out_93[j + lat_93];
      tail_93 += v * v;
    }
  }
  TEST_ASSERT(onsets > 0, "Must cover at least one onset");
  in_pow /= n;
  p23 /= n;
  p93 /= n;
  diff_pow /= n;

  printf("  lat23=%u lat93=%u in=%.6f out23=%.6f out93=%.6f diff=%.6f\n",
         lat_23, lat_93, in_pow, p23, p93, diff_pow);
  printf("  peak23=%.4f peak93=%.4f onsets=%d\n", peak_23, peak_93, onsets);
  printf(
      "  pre23=%.6f pre93=%.6f tail23=%.6f tail93=%.6f pk23=%.6f pk93=%.6f\n",
      pre_23, pre_93, tail_23, tail_93, pk_e_23, pk_e_93);

  /* Both runs must reduce noise but preserve transient peaks. */
  TEST_ASSERT(p23 < in_pow, "23ms run should reduce noise power");
  TEST_ASSERT(p93 < in_pow, "93ms run should reduce noise power");
  TEST_ASSERT(peak_23 > 0.1, "23ms run should preserve transient peaks");
  TEST_ASSERT(peak_93 > 0.1, "93ms run should preserve transient peaks");
  TEST_ASSERT(diff_pow < in_pow,
              "Cross-run difference should stay below input power");

  /* Smear invariance: pre-echo (future context) and tail (past context)
   * per unit peak energy must match across frame sizes. Absolute
   * suppression may differ (finer frequency resolution at larger frames is
   * the intended tradeoff), but the temporal spread in ms must not.
   * Measured: fixed-ms geometry agrees ~1.2x; frame-counted geometry
   * diverges ~2.2x. */
  TEST_ASSERT(pk_e_23 > 0.0 && pk_e_93 > 0.0, "Peak energy must be positive");
  const double pre_ratio = (pre_93 / pk_e_93) / (pre_23 / pk_e_23);
  const double tail_ratio = (tail_93 / pk_e_93) / (tail_23 / pk_e_23);
  printf("  pre_ratio=%.3f tail_ratio=%.3f\n", pre_ratio, tail_ratio);
  TEST_ASSERT(pre_ratio > 1.0 / 1.6 && pre_ratio < 1.6,
              "Pre-echo smear should match within 1.6x across frame sizes");
  TEST_ASSERT(tail_ratio > 1.0 / 1.6 && tail_ratio < 1.6,
              "Tail smear should match within 1.6x across frame sizes");

  free(input);
  free(out_23);
  free(out_93);

  printf("✓ Frame-size invariance test passed\n");
}
