/*
 * Unit tests for FFT Transform
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/shared/configurations.h"
#include "../src/shared/stft/fft_transform.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_fft_lifecycle(void) {
  printf("Testing FFT lifecycle...\n");

  uint32_t frame_size = 512;
  FftTransform* fft = fft_transform_initialize(frame_size, NO_PADDING, 0);
  TEST_ASSERT(fft != NULL, "FFT initialization should succeed");
  TEST_ASSERT(get_fft_size(fft) == 512, "FFT size should be 512");
  TEST_ASSERT(get_fft_real_spectrum_size(fft) == 257,
              "Spectrum size should be 257");

  fft_transform_free(fft);
  printf("✓ FFT lifecycle tests passed\n");
}

void test_fft_padding(void) {
  printf("Testing FFT padding modes...\n");

  // NEXT_POWER_OF_TWO
  FftTransform* fft_pow2 = fft_transform_initialize(500, NEXT_POWER_OF_TWO, 0);
  TEST_ASSERT(get_fft_size(fft_pow2) == 512, "500 -> 512 with pow2");
  fft_transform_free(fft_pow2);

  // FIXED_AMOUNT
  FftTransform* fft_fixed = fft_transform_initialize(512, FIXED_AMOUNT, 512);
  TEST_ASSERT(get_fft_size(fft_fixed) == 1024, "512 + 512 = 1024");
  fft_transform_free(fft_fixed);

  printf("✓ FFT padding tests passed\n");
}

void test_fft_load_store(void) {
  printf("Testing FFT load/store...\n");

  uint32_t frame_size = 128;
  FftTransform* fft = fft_transform_initialize(frame_size, NO_PADDING, 0);

  float input[128];
  float output[128];
  for (int i = 0; i < 128; i++)
    input[i] = (float)i;

  TEST_ASSERT(fft_load_input_samples(fft, input),
              "Load samples should succeed");
  TEST_ASSERT(fft_get_output_samples(fft, output),
              "Get samples should succeed");

  for (int i = 0; i < 128; i++) {
    TEST_FLOAT_CLOSE(input[i], output[i], 1e-6f);
  }

  fft_transform_free(fft);
  printf("✓ FFT load/store tests passed\n");
}

void test_fft_computation(void) {
  printf("Testing FFT computation (Forward/Backward)...\n");

  uint32_t frame_size = 1024;
  FftTransform* fft = fft_transform_initialize(frame_size, NO_PADDING, 0);

  float input[1024];
  float output[1024];
  // Simple sine wave
  for (int i = 0; i < 1024; i++) {
    input[i] = sinf(2.0f * M_PI * 10.0f * (float)i / 1024.0f);
  }

  fft_load_input_samples(fft, input);
  TEST_ASSERT(compute_forward_fft(fft), "Forward FFT should succeed");
  TEST_ASSERT(compute_backward_fft(fft), "Backward FFT should succeed");
  fft_get_output_samples(fft, output);

  // Result of Forward then Backward (R2HC/HC2R) should be scaled by N
  float scale = 1.0f / 1024.0f;
  for (int i = 0; i < 1024; i++) {
    TEST_FLOAT_CLOSE(input[i], output[i] * scale, 1e-4f);
  }

  fft_transform_free(fft);
  printf("✓ FFT computation tests passed\n");
}

int main(void) {
  printf("Running FFT transform tests...\n\n");

  test_fft_lifecycle();
  test_fft_padding();
  test_fft_load_store();
  test_fft_computation();

  printf("\n✅ All FFT transform tests passed!\n");
  return 0;
}
