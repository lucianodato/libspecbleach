/*
 * Unit tests for utility functions
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Include internal headers for testing
#include "../src/shared/utils/denoise_mixer.h"
#include "../src/shared/utils/general_utils.h"
#include "../src/shared/utils/spectral_features.h"
#include "../src/shared/utils/spectral_trailing_buffer.h"
#include "../src/shared/utils/spectral_utils.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_spectral_trailing_buffer(void) {
  printf("Testing spectral_trailing_buffer...\n");

  uint32_t spec_size = 512;
  uint32_t buf_size = 10;
  SpectralTrailingBuffer* stbuff =
      spectral_trailing_buffer_initialize(spec_size, buf_size);
  TEST_ASSERT(stbuff != NULL, "STBuff initialization should succeed");
  TEST_ASSERT(get_spectrum_size(stbuff) == spec_size, "Spec size check");
  TEST_ASSERT(get_spectrum_buffer_size(stbuff) == buf_size,
              "Buffer size check");

  float test_spec[512];
  for (int i = 0; i < 512; i++) {
    test_spec[i] = (float)i;
  }

  TEST_ASSERT(spectral_trailing_buffer_push_back(stbuff, test_spec),
              "Push should succeed");

  spectral_trailing_buffer_free(stbuff);
  printf("✓ spectral_trailing_buffer tests passed\n");
}

void test_denoise_mixer(void) {
  printf("Testing denoise_mixer...\n");

  uint32_t fft_size = 512;
  uint32_t sr = 44100;
  uint32_t hop = 128;

  DenoiseMixer* mixer = denoise_mixer_initialize(fft_size, sr, hop);
  TEST_ASSERT(mixer != NULL, "Mixer initialization should succeed");

  float* fft_spectrum = (float*)calloc(fft_size, sizeof(float));
  float* gain_spectrum = (float*)calloc(fft_size, sizeof(float));
  for (uint32_t i = 0; i < fft_size; i++) {
    fft_spectrum[i] = 1.0f;
    gain_spectrum[i] = 0.5f;
  }

  DenoiseMixerParameters params = {0};
  params.residual_listen = true;

  // Cover residual_listen branch
  TEST_ASSERT(denoise_mixer_run(mixer, fft_spectrum, gain_spectrum, params),
              "Mixer run should succeed");
  for (uint32_t i = 0; i < fft_size; i++) {
    TEST_FLOAT_CLOSE(fft_spectrum[i], 0.5f, 1e-6f); // 1.0 - (1.0 * 0.5)
  }

  // Cover NULL checks
  TEST_ASSERT(!denoise_mixer_run(NULL, fft_spectrum, gain_spectrum, params),
              "Mixer run should fail with NULL mixer");
  TEST_ASSERT(!denoise_mixer_run(mixer, NULL, gain_spectrum, params),
              "Mixer run should fail with NULL spectrum");
  TEST_ASSERT(!denoise_mixer_run(mixer, fft_spectrum, NULL, params),
              "Mixer run should fail with NULL gain");

  denoise_mixer_free(mixer);
  free(fft_spectrum);
  free(gain_spectrum);
  printf("✓ denoise_mixer tests passed\n");
}

void test_spectral_features(void) {
  printf("Testing spectral_features...\n");
  uint32_t fft_size = 512;
  uint32_t real_size = (fft_size / 2) + 1;
  SpectralFeatures* sf = spectral_features_initialize(real_size);
  TEST_ASSERT(sf != NULL, "SpectralFeatures initialization should succeed");

  float* fft_spectrum = (float*)calloc(fft_size, sizeof(float));
  for (uint32_t i = 0; i < fft_size; i++) {
    fft_spectrum[i] = (float)i / (float)fft_size;
  }

  // Cover Power Spectrum
  float* power =
      get_spectral_feature(sf, fft_spectrum, fft_size, POWER_SPECTRUM);
  TEST_ASSERT(power != NULL, "Power spectrum calculation should succeed");

  // Cover Magnitude Spectrum
  float* magnitude =
      get_spectral_feature(sf, fft_spectrum, fft_size, MAGNITUDE_SPECTRUM);
  TEST_ASSERT(magnitude != NULL,
              "Magnitude spectrum calculation should succeed");

  // Cover Phase Spectrum
  float* phase =
      get_spectral_feature(sf, fft_spectrum, fft_size, PHASE_SPECTRUM);
  TEST_ASSERT(phase != NULL, "Phase spectrum calculation should succeed");

  // Cover Invalid type
  float* invalid =
      get_spectral_feature(sf, fft_spectrum, fft_size, (SpectrumType)99);
  TEST_ASSERT(invalid == NULL, "Invalid type should return NULL");

  // Cover NULL checks
  TEST_ASSERT(get_spectral_feature(NULL, fft_spectrum, fft_size,
                                   POWER_SPECTRUM) == NULL,
              "Should fail with NULL handle");
  TEST_ASSERT(get_spectral_feature(sf, NULL, fft_size, POWER_SPECTRUM) == NULL,
              "Should fail with NULL spectrum");
  TEST_ASSERT(get_spectral_feature(sf, fft_spectrum, 0, POWER_SPECTRUM) == NULL,
              "Should fail with size 0");

  spectral_features_free(sf);
  free(fft_spectrum);
  printf("✓ spectral_features tests passed\n");
}

void test_sanitize_denormal(void) {
  printf("Testing sanitize_denormal...\n");

  // Test normal values
  TEST_FLOAT_CLOSE(sanitize_denormal(1.0f), 1.0f, 1e-6f);
  TEST_FLOAT_CLOSE(sanitize_denormal(-1.0f), -1.0f, 1e-6f);
  TEST_FLOAT_CLOSE(sanitize_denormal(0.5f), 0.5f, 1e-6f);

  // Test zero
  TEST_FLOAT_CLOSE(sanitize_denormal(0.0f), 0.0f, 1e-6f);

  // Test very small values (denormals)
  float denormal = 1e-40f; // This should be treated as denormal
  TEST_FLOAT_CLOSE(sanitize_denormal(denormal), 0.0f, 1e-6f);

  printf("✓ sanitize_denormal tests passed\n");
}

void test_from_db_to_coefficient(void) {
  printf("Testing from_db_to_coefficient...\n");

  // Test 0 dB = 1.0 coefficient
  TEST_FLOAT_CLOSE(from_db_to_coefficient(0.0f), 1.0f, 1e-6f);

  // Test -6 dB ≈ 0.501 coefficient
  TEST_FLOAT_CLOSE(from_db_to_coefficient(-6.0f), 0.501f, 0.01f);

  // Test -12 dB ≈ 0.251 coefficient
  TEST_FLOAT_CLOSE(from_db_to_coefficient(-12.0f), 0.251f, 0.01f);

  // Test positive dB
  TEST_ASSERT(from_db_to_coefficient(6.0f) > 1.5f,
              "Positive dB should give coefficient > 1.5");

  printf("✓ from_db_to_coefficient tests passed\n");
}

void test_get_next_divisible_two(void) {
  printf("Testing get_next_divisible_two...\n");

  TEST_ASSERT(get_next_divisible_two(1) == 2, "1 -> 2");
  TEST_ASSERT(get_next_divisible_two(2) == 2, "2 -> 2");
  TEST_ASSERT(get_next_divisible_two(3) == 4, "3 -> 4");
  TEST_ASSERT(get_next_divisible_two(4) == 4, "4 -> 4");
  TEST_ASSERT(get_next_divisible_two(5) == 6, "5 -> 6");
  TEST_ASSERT(get_next_divisible_two(7) == 8, "7 -> 8");
  TEST_ASSERT(get_next_divisible_two(8) == 8, "8 -> 8");
  TEST_ASSERT(get_next_divisible_two(9) == 10, "9 -> 10");

  printf("✓ get_next_divisible_two tests passed\n");
}

void test_get_next_power_two(void) {
  printf("Testing get_next_power_two...\n");

  TEST_ASSERT(get_next_power_two(1) == 1, "1 -> 1 (2^0)");
  TEST_ASSERT(get_next_power_two(2) == 2, "2 -> 2 (2^1)");
  TEST_ASSERT(get_next_power_two(3) == 4, "3 -> 4 (2^2)");
  TEST_ASSERT(get_next_power_two(4) == 4, "4 -> 4 (2^2)");
  TEST_ASSERT(get_next_power_two(5) == 8, "5 -> 8 (2^3)");
  TEST_ASSERT(get_next_power_two(7) == 8, "7 -> 8 (2^3)");
  TEST_ASSERT(get_next_power_two(8) == 8, "8 -> 8 (2^3)");
  TEST_ASSERT(get_next_power_two(9) == 16, "9 -> 16 (2^4)");
  TEST_ASSERT(get_next_power_two(15) == 16, "15 -> 16 (2^4)");
  TEST_ASSERT(get_next_power_two(16) == 16, "16 -> 16 (2^4)");
  TEST_ASSERT(get_next_power_two(17) == 32, "17 -> 32 (2^5)");

  printf("✓ get_next_power_two tests passed\n");
}

void test_remap_percentage_log_like_unity(void) {
  printf("Testing remap_percentage_log_like_unity...\n");

  // Test edge cases
  TEST_FLOAT_CLOSE(remap_percentage_log_like_unity(0.0f), 0.0f, 1e-6f);
  TEST_FLOAT_CLOSE(remap_percentage_log_like_unity(1.0f), 0.9502f, 0.001f);

  // Test monotonicity - function should be increasing
  float val1 = remap_percentage_log_like_unity(0.1f);
  float val2 = remap_percentage_log_like_unity(0.5f);
  float val3 = remap_percentage_log_like_unity(0.9f);
  TEST_ASSERT(val1 < val2 && val2 < val3,
              "Function should be monotonically increasing");

  printf("✓ remap_percentage_log_like_unity tests passed\n");
}

void test_get_fft_window(void) {
  printf("Testing get_fft_window...\n");

  uint32_t fft_size = 1024;
  float window[1024];

  // Test HANN
  TEST_ASSERT(get_fft_window(window, fft_size, HANN_WINDOW),
              "HANN should succeed");
  // Hann window at endpoints should be 0
  TEST_FLOAT_CLOSE(window[0], 0.0f, 1e-6f);
  // Hann window at middle should be around 1.0 (or slightly less depending on
  // implementation)
  TEST_ASSERT(window[512] > 0.9f, "Hann middle should be high");

  // Test HAMMING
  TEST_ASSERT(get_fft_window(window, fft_size, HAMMING_WINDOW),
              "HAMMING should succeed");
  // Hamming at endpoint is NOT zero (typically 0.08)
  TEST_ASSERT(window[0] > 0.07f && window[0] < 0.09f, "Hamming endpoint check");

  printf("✓ get_fft_window tests passed\n");
}

void test_initialize_spectrum_with_value(void) {
  printf("Testing initialize_spectrum_with_value...\n");

  uint32_t size = 100;
  float spectrum[100];

  TEST_ASSERT(initialize_spectrum_with_value(spectrum, size, 42.0f),
              "Init should succeed");
  for (uint32_t i = 0; i < size; i++) {
    TEST_FLOAT_CLOSE(spectrum[i], 42.0f, 1e-6f);
  }

  printf("✓ initialize_spectrum_with_value tests passed\n");
}

void test_spectral_min_max(void) {
  printf("Testing max_spectral_value and min_spectral_value...\n");

  float spectrum[] = {1.0f, 5.0f, 2.0f, 10.0f, -3.0f, 4.0f};
  uint32_t size = 6;

  TEST_FLOAT_CLOSE(max_spectral_value(spectrum, size), 10.0f, 1e-6f);
  TEST_FLOAT_CLOSE(min_spectral_value(spectrum, size), -3.0f, 1e-6f);

  printf("✓ spectral min/max tests passed\n");
}

void test_fft_frequency_conversions(void) {
  printf("Testing fft_bin_to_freq and freq_to_fft_bin...\n");

  uint32_t sr = 44100;
  uint32_t fft_size = 1024;

  // DC bin should be 0 Hz
  TEST_FLOAT_CLOSE(fft_bin_to_freq(0, sr, fft_size), 0.0f, 1e-6f);
  TEST_ASSERT(freq_to_fft_bin(0.0f, sr, fft_size) == 0, "0 Hz -> bin 0");

  // Nyquist bin (fft_size/2) should be sr/2
  float nyquist_freq = (float)sr / 2.0f;
  TEST_FLOAT_CLOSE(fft_bin_to_freq(512, sr, fft_size), nyquist_freq, 1e-6f);
  TEST_ASSERT(freq_to_fft_bin(nyquist_freq, sr, fft_size) == 512,
              "Nyquist -> bin 512");

  // Half way to Nyquist
  float mid_freq = nyquist_freq / 2.0f;
  TEST_FLOAT_CLOSE(fft_bin_to_freq(256, sr, fft_size), mid_freq, 1e-6f);
  TEST_ASSERT(freq_to_fft_bin(mid_freq, sr, fft_size) == 256,
              "Mid freq -> bin 256");

  printf("✓ fft frequency conversion tests passed\n");
}

void test_min_max_spectrum_float(void) {
  printf("Testing min_spectrum_float and max_spectrum_float...\n");

  uint32_t size = 4;
  float s1[] = {1.0f, 10.0f, 5.0f, 2.0f};
  float s2[] = {2.0f, 5.0f, 5.0f, 1.0f};

  float s1_min[] = {1.0f, 10.0f, 5.0f, 2.0f};
  min_spectrum_float(s1_min, s2, size);
  TEST_FLOAT_CLOSE(s1_min[0], 1.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_min[1], 5.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_min[2], 5.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_min[3], 1.0f, 1e-6f);

  float s1_max[] = {1.0f, 10.0f, 5.0f, 2.0f};
  max_spectrum_float(s1_max, s2, size);
  TEST_FLOAT_CLOSE(s1_max[0], 2.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_max[1], 10.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_max[2], 5.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_max[3], 2.0f, 1e-6f);

  printf("✓ min/max spectrum float tests passed\n");
}

void test_spectral_utils_null_and_edge_cases(void) {
  printf("Testing spectral_utils NULL pointer and edge cases...\n");

  float spectrum[10] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f,
                        6.0f, 7.0f, 8.0f, 9.0f, 10.0f};
  float spectrum2[10] = {10.0f, 9.0f, 8.0f, 7.0f, 6.0f,
                         5.0f,  4.0f, 3.0f, 2.0f, 1.0f};

  // Test NULL pointer handling
  TEST_ASSERT(!initialize_spectrum_with_value(NULL, 10, 1.0f),
              "initialize_spectrum_with_value should fail with NULL");
  TEST_ASSERT(!initialize_spectrum_with_value(spectrum, 0, 1.0f),
              "initialize_spectrum_with_value should fail with size 0");

  TEST_FLOAT_CLOSE(max_spectral_value(NULL, 10), 0.0f, 1e-6f);
  TEST_FLOAT_CLOSE(max_spectral_value(spectrum, 0), 0.0f, 1e-6f);

  TEST_FLOAT_CLOSE(min_spectral_value(NULL, 10), 0.0f, 1e-6f);
  TEST_FLOAT_CLOSE(min_spectral_value(spectrum, 0), 0.0f, 1e-6f);

  TEST_ASSERT(!min_spectrum_float(NULL, spectrum2, 10),
              "min_spectrum_float should fail with NULL spectrum_one");
  TEST_ASSERT(!min_spectrum_float(spectrum, NULL, 10),
              "min_spectrum_float should fail with NULL spectrum_two");
  TEST_ASSERT(!min_spectrum_float(spectrum, spectrum2, 0),
              "min_spectrum_float should fail with size 0");

  TEST_ASSERT(!max_spectrum_float(NULL, spectrum2, 10),
              "max_spectrum_float should fail with NULL spectrum_one");
  TEST_ASSERT(!max_spectrum_float(spectrum, NULL, 10),
              "max_spectrum_float should fail with NULL spectrum_two");
  TEST_ASSERT(!max_spectrum_float(spectrum, spectrum2, 0),
              "max_spectrum_float should fail with size 0");

  // Test double version edge cases
  double dspectrum[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  double dspectrum2[10] = {10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};

  TEST_ASSERT(!min_spectrum_double(NULL, dspectrum2, 10),
              "min_spectrum_double should fail with NULL");
  TEST_ASSERT(!min_spectrum_double(dspectrum, NULL, 10),
              "min_spectrum_double should fail with NULL");
  TEST_ASSERT(!min_spectrum_double(dspectrum, dspectrum2, 0),
              "min_spectrum_double should fail with size 0");

  TEST_ASSERT(!max_spectrum_double(NULL, dspectrum2, 10),
              "max_spectrum_double should fail with NULL");
  TEST_ASSERT(!max_spectrum_double(dspectrum, NULL, 10),
              "max_spectrum_double should fail with NULL");
  TEST_ASSERT(!max_spectrum_double(dspectrum, dspectrum2, 0),
              "max_spectrum_double should fail with size 0");

  // Test get_fft_window edge cases
  float window[10];
  TEST_ASSERT(!get_fft_window(NULL, 10, HANN_WINDOW),
              "get_fft_window should fail with NULL");
  TEST_ASSERT(!get_fft_window(window, 0, HANN_WINDOW),
              "get_fft_window should fail with size 0");

  // Test BLACKMAN and VORBIS windows
  float bwindow[64];
  TEST_ASSERT(get_fft_window(bwindow, 64, BLACKMAN_WINDOW),
              "BLACKMAN window should succeed");
  TEST_ASSERT(get_fft_window(bwindow, 64, VORBIS_WINDOW),
              "VORBIS window should succeed");

  printf("✓ spectral_utils NULL/edge case tests passed\n");
}

int main(void) {
  printf("Running utility function tests...\n\n");

  test_sanitize_denormal();
  test_from_db_to_coefficient();
  test_get_next_divisible_two();
  test_get_next_power_two();
  test_remap_percentage_log_like_unity();
  test_get_fft_window();
  test_initialize_spectrum_with_value();
  test_spectral_min_max();
  test_fft_frequency_conversions();
  test_min_max_spectrum_float();
  test_spectral_utils_null_and_edge_cases();
  test_denoise_mixer();
  test_spectral_trailing_buffer();
  test_spectral_features();

  printf("\n✅ All utility function tests passed!\n");
  return 0;
}
