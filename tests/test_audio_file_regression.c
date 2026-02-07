/*
 * Audio file regression tests - compare library output against reference files
 */

#include <math.h>
#include <sndfile.h>

#include <specbleach_denoiser.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define BLOCK_SIZE 512
#define NOISE_FRAMES 8
#define FRAME_SIZE_MS 46.0f
#define FRAME_SIZE_ADAPTIVE_MS 20.0f

// Canonical parameters used in generate_reference_files.sh
static const SpectralBleachDenoiserParameters CANONICAL_DENOISER_PARAMS = {
    .residual_listen = false,
    .learn_noise = 1,          // Learn all modes
    .noise_reduction_mode = 3, // Use maximum profile
    .reduction_amount = 20.0f,
    .smoothing_factor = 0.0f,
    .whitening_factor = 50.0f,
    .noise_scaling_type = 2,
    .noise_rescale = 6.0f};

static const SpectralBleachDenoiserParameters CANONICAL_ADENOISER_PARAMS = {
    .residual_listen = false,
    .reduction_amount = 20.0f,
    .smoothing_factor = 0.0f,
    .whitening_factor = 50.0f,
    .noise_scaling_type = 2,
    .noise_rescale = 6.0f,

    .adaptive_noise = 1,
    .noise_estimation_method = 0};

void test_denoiser_file_regression(void) {
  printf("Testing spectral denoiser file regression...\n");

  const char* input_path = TEST_DATA_DIR "Speech.wav";
  const char* reference_path = TEST_DATA_DIR "Speech_denoised.wav";

  SF_INFO in_info, ref_info;
  SNDFILE* in_sf = sf_open(input_path, SFM_READ, &in_info);
  SNDFILE* ref_sf = sf_open(reference_path, SFM_READ, &ref_info);

  if (!ref_sf) {
    printf("  Note: Reference file %s not found. Skipping test.\n",
           reference_path);
    if (in_sf)
      sf_close(in_sf);
    return;
  }

  TEST_ASSERT(in_sf != NULL, "Failed to open input file");
  TEST_ASSERT(in_info.channels == 1, "Input file must be mono");

  SpectralBleachHandle handle =
      specbleach_initialize((uint32_t)in_info.samplerate, FRAME_SIZE_MS);
  TEST_ASSERT(handle != NULL, "Failed to initialize denoiser");

  SpectralBleachDenoiserParameters params = CANONICAL_DENOISER_PARAMS;
  specbleach_load_parameters(handle, params);

  float* in_buf = malloc(BLOCK_SIZE * sizeof(float));
  float* out_buf = malloc(BLOCK_SIZE * sizeof(float));
  float* ref_buf = malloc(BLOCK_SIZE * sizeof(float));

  // Learn stage
  for (int i = 0; i < NOISE_FRAMES; i++) {
    sf_count_t read = sf_readf_float(in_sf, in_buf, BLOCK_SIZE);
    TEST_ASSERT(read == BLOCK_SIZE, "Failed to read learn frames");
    specbleach_process(handle, (uint32_t)BLOCK_SIZE, in_buf, out_buf);
  }

  // Reduction stage
  params.learn_noise = 0;
  specbleach_load_parameters(handle, params);

  sf_count_t total_read = 0;
  sf_count_t read;
  while ((read = sf_readf_float(in_sf, in_buf, BLOCK_SIZE)) > 0) {
    specbleach_process(handle, (uint32_t)read, in_buf, out_buf);
    sf_readf_float(ref_sf, ref_buf, read);

    for (sf_count_t i = 0; i < read; i++) {
      // Use tolerance for floating point variations across platforms/compilers
      if (fabsf(out_buf[i] - ref_buf[i]) > 1e-3f) {
        fprintf(stderr, "Mismatch at sample %lld: %f != %f\n",
                (long long)(total_read + i), out_buf[i], ref_buf[i]);
        exit(1);
      }
    }
    total_read += read;
  }

  free(in_buf);
  free(out_buf);
  free(ref_buf);
  specbleach_free(handle);
  sf_close(in_sf);
  sf_close(ref_sf);

  printf("✓ Spectral denoiser file regression passed\n");
}

void test_adenoiser_file_regression(void) {
  printf("Testing adaptive denoiser file regression...\n");

  const char* input_path = TEST_DATA_DIR "Speech.wav";
  const char* reference_path = TEST_DATA_DIR "Speech_adaptive_denoised.wav";

  SF_INFO in_info, ref_info;
  SNDFILE* in_sf = sf_open(input_path, SFM_READ, &in_info);
  SNDFILE* ref_sf = sf_open(reference_path, SFM_READ, &ref_info);

  if (!ref_sf) {
    printf("  Note: Reference file %s not found. Skipping test.\n",
           reference_path);
    if (in_sf)
      sf_close(in_sf);
    return;
  }

  TEST_ASSERT(in_sf != NULL, "Failed to open input file");
  TEST_ASSERT(in_info.channels == 1, "Input file must be mono");

  SpectralBleachHandle handle = specbleach_initialize(
      (uint32_t)in_info.samplerate, FRAME_SIZE_ADAPTIVE_MS);
  TEST_ASSERT(handle != NULL, "Failed to initialize adaptive denoiser");

  SpectralBleachDenoiserParameters params = CANONICAL_ADENOISER_PARAMS;
  specbleach_load_parameters(handle, params);

  float* in_buf = malloc(BLOCK_SIZE * sizeof(float));
  float* out_buf = malloc(BLOCK_SIZE * sizeof(float));
  float* ref_buf = malloc(BLOCK_SIZE * sizeof(float));

  sf_count_t total_read = 0;
  sf_count_t read;
  while ((read = sf_readf_float(in_sf, in_buf, BLOCK_SIZE)) > 0) {
    specbleach_process(handle, (uint32_t)read, in_buf, out_buf);
    sf_readf_float(ref_sf, ref_buf, read);

    for (sf_count_t i = 0; i < read; i++) {
      if (fabsf(out_buf[i] - ref_buf[i]) > 1e-1f) {
        fprintf(stderr, "Mismatch at sample %lld: %f != %f\n",
                (long long)(total_read + i), out_buf[i], ref_buf[i]);
        exit(1);
      }
    }
    total_read += read;
  }

  free(in_buf);
  free(out_buf);
  free(ref_buf);
  specbleach_free(handle);
  sf_close(in_sf);
  sf_close(ref_sf);

  printf("✓ Adaptive denoiser file regression passed\n");
}

int main(void) {
  printf("Running audio file regression tests...\n\n");

  test_denoiser_file_regression();
  test_adenoiser_file_regression();

  printf("\n✅ All audio file regression tests passed!\n");
  return 0;
}
