/*
libspecbleach - A spectral processing library

Test suite for the frequency-adaptive tonal detector with parabolic
interpolation.
*/

#include "shared/utils/tonal_detector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_SAMPLE_RATE 48000
#define TEST_FFT_SIZE 2048
#define TEST_SPECTRUM_SIZE ((TEST_FFT_SIZE / 2) + 1)
#define BIN_HZ ((float)TEST_SAMPLE_RATE / (float)TEST_FFT_SIZE)

// Helper: create a flat noise profile with optional tonal peaks
static void create_test_profile(float* profile, uint32_t size,
                                float noise_floor, const int* peak_bins,
                                const float* peak_amplitudes, int n_peaks) {
  for (uint32_t k = 0; k < size; k++) {
    profile[k] = noise_floor;
  }
  for (int i = 0; i < n_peaks; i++) {
    if (peak_bins[i] >= 0 && peak_bins[i] < (int)size) {
      profile[peak_bins[i]] = peak_amplitudes[i];
    }
  }
}

// Helper: find the bin closest to a target frequency
static int freq_to_bin(float freq_hz) {
  return (int)roundf(freq_hz / BIN_HZ);
}

// --- Tests ---

void test_null_inputs(void) {
  printf("Testing NULL inputs...\n");
  float tonal_mask[TEST_SPECTRUM_SIZE];
  float profile[TEST_SPECTRUM_SIZE];
  float max_profile[TEST_SPECTRUM_SIZE];
  float median_profile[TEST_SPECTRUM_SIZE];

  memset(profile, 0, sizeof(profile));
  memset(max_profile, 0, sizeof(max_profile));
  memset(median_profile, 0, sizeof(median_profile));
  memset(tonal_mask, 0, sizeof(tonal_mask));

  // NULL profile should not crash
  detect_tonal_components(NULL, max_profile, median_profile, TEST_SPECTRUM_SIZE,
                          TEST_SAMPLE_RATE, TEST_FFT_SIZE, tonal_mask);

  // NULL tonal_mask should not crash
  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          NULL);

  // size < 5 should not crash
  detect_tonal_components(profile, max_profile, median_profile, 3,
                          TEST_SAMPLE_RATE, TEST_FFT_SIZE, tonal_mask);

  // Zero sample_rate should not crash
  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, 0, TEST_FFT_SIZE, tonal_mask);

  // Zero fft_size should not crash
  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, 0, tonal_mask);

  printf("✓ NULL/edge-case inputs handled safely\n");
}

void test_low_frequency_hum_detection(void) {
  printf("Testing low-frequency hum detection (60 Hz + 120 Hz)...\n");

  float* profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* max_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* median_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* tonal_mask = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));

  // 60 Hz fundamental + 120 Hz harmonic
  int bin_60 = freq_to_bin(60.0f);
  int bin_120 = freq_to_bin(120.0f);
  float noise_floor = 0.01f;

  int peak_bins[] = {bin_60, bin_120};
  float peak_amps[] = {0.1f, 0.05f}; // 60Hz stronger than 120Hz

  create_test_profile(profile, TEST_SPECTRUM_SIZE, noise_floor, peak_bins,
                      peak_amps, 2);
  // For stationarity: max ~ median means stable tone
  create_test_profile(max_profile, TEST_SPECTRUM_SIZE, noise_floor, peak_bins,
                      peak_amps, 2);
  create_test_profile(median_profile, TEST_SPECTRUM_SIZE, noise_floor,
                      peak_bins, peak_amps, 2);

  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  // Both bins should be detected
  if (tonal_mask[bin_60] <= 0.0f) {
    fprintf(stderr, "FAIL: 60 Hz bin %d not detected (mask = %f)\n", bin_60,
            tonal_mask[bin_60]);
    exit(1);
  }
  printf("  60 Hz (bin %d): mask = %.3f ✓\n", bin_60, tonal_mask[bin_60]);

  if (tonal_mask[bin_120] <= 0.0f) {
    fprintf(stderr, "FAIL: 120 Hz bin %d not detected (mask = %f)\n", bin_120,
            tonal_mask[bin_120]);
    exit(1);
  }
  printf("  120 Hz (bin %d): mask = %.3f ✓\n", bin_120, tonal_mask[bin_120]);

  // Test 30 Hz (Bin 1 at 48k/2048)
  int bin_30 = freq_to_bin(30.0f);
  create_test_profile(profile, TEST_SPECTRUM_SIZE, noise_floor, &bin_30,
                      peak_amps, 1);
  create_test_profile(max_profile, TEST_SPECTRUM_SIZE, noise_floor, &bin_30,
                      peak_amps, 1);
  create_test_profile(median_profile, TEST_SPECTRUM_SIZE, noise_floor, &bin_30,
                      peak_amps, 1);

  memset(tonal_mask, 0, TEST_SPECTRUM_SIZE * sizeof(float));
  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  if (tonal_mask[bin_30] <= 0.0f) {
    fprintf(stderr,
            "FAIL: 30 Hz (bin %d) not detected. Loop extension failed.\n",
            bin_30);
    exit(1);
  }
  printf("  30 Hz (bin %d): mask = %.3f ✓ (Loop extension to k=1 confirmed)\n",
         bin_30, tonal_mask[bin_30]);

  free(profile);
  free(max_profile);
  free(median_profile);
  free(tonal_mask);
  printf("✓ Low-frequency hum detection passed\n");
}

void test_high_frequency_detection(void) {
  printf("Testing high-frequency tone detection (5 kHz)...\n");

  float* profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* max_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* median_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* tonal_mask = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));

  int bin_5k = freq_to_bin(5000.0f);
  float noise_floor = 0.01f;

  int peak_bins[] = {bin_5k};
  float peak_amps[] = {0.1f};

  create_test_profile(profile, TEST_SPECTRUM_SIZE, noise_floor, peak_bins,
                      peak_amps, 1);
  create_test_profile(max_profile, TEST_SPECTRUM_SIZE, noise_floor, peak_bins,
                      peak_amps, 1);
  create_test_profile(median_profile, TEST_SPECTRUM_SIZE, noise_floor,
                      peak_bins, peak_amps, 1);

  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  if (tonal_mask[bin_5k] <= 0.0f) {
    fprintf(stderr, "FAIL: 5 kHz bin %d not detected (mask = %f)\n", bin_5k,
            tonal_mask[bin_5k]);
    exit(1);
  }
  printf("  5 kHz (bin %d): mask = %.3f ✓\n", bin_5k, tonal_mask[bin_5k]);

  free(profile);
  free(max_profile);
  free(median_profile);
  free(tonal_mask);
  printf("✓ High-frequency detection passed (no regression)\n");
}

void test_adaptive_width(void) {
  printf("Testing adaptive width detection via median filter...\n");

  float* profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* max_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* median_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* tonal_mask = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));

  float noise_floor = 0.01f;
  // Initialize backgrounds
  for (uint32_t i = 0; i < TEST_SPECTRUM_SIZE; i++) {
    profile[i] = noise_floor;
    max_profile[i] = noise_floor;
    median_profile[i] = noise_floor;
  }

  // Create a peak centered at bin 100
  // Peak: 0.1, Neighbors: 0.08, 0.04, 0.02 (still > PEAK_THRESHOLD * floor)
  int center = 100;
  profile[center] = 0.1f;
  max_profile[center] = 0.1f;
  median_profile[center] = 0.1f;

  profile[center - 1] = 0.08f;
  max_profile[center - 1] = 0.08f;
  median_profile[center - 1] = 0.08f;
  profile[center + 1] = 0.08f;
  max_profile[center + 1] = 0.08f;
  median_profile[center + 1] = 0.08f;

  profile[center - 2] = 0.04f;
  max_profile[center - 2] = 0.04f;
  median_profile[center - 2] = 0.04f;
  profile[center + 2] = 0.04f;
  max_profile[center + 2] = 0.04f;
  median_profile[center + 2] = 0.04f;

  profile[center - 3] = 0.02f;
  max_profile[center - 3] = 0.02f;
  median_profile[center - 3] = 0.02f;
  profile[center + 3] = 0.02f;
  max_profile[center + 3] = 0.02f;
  median_profile[center + 3] = 0.02f;

  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  // Verification:
  // Center should be detected
  if (tonal_mask[center] <= 0.0f) {
    fprintf(stderr, "FAIL: Center peak %d not detected\n", center);
    exit(1);
  }

  // +/- 3 should be detected because they are > PEAK_THRESHOLD * floor
  if (tonal_mask[center - 3] <= 0.0f || tonal_mask[center + 3] <= 0.0f) {
    fprintf(stderr,
            "FAIL: Peak edges (+/-3) not detected. Mask[-3]=%f, Mask[+3]=%f\n",
            tonal_mask[center - 3], tonal_mask[center + 3]);
    exit(1);
  }

  printf("  Center Mask: %.3f\n", tonal_mask[center]);
  printf("  Edge Mask (+/-3): %.3f / %.3f\n", tonal_mask[center - 3],
         tonal_mask[center + 3]);

  // +/- 4 should NOT be detected (floor)
  if (tonal_mask[center - 4] > 0.0f || tonal_mask[center + 4] > 0.0f) {
    fprintf(stderr, "FAIL: Noise floor bins masked! Mask[-4]=%f, Mask[+4]=%f\n",
            tonal_mask[center - 4], tonal_mask[center + 4]);
    exit(1);
  }

  free(profile);
  free(max_profile);
  free(median_profile);
  free(tonal_mask);
  printf("✓ Adaptive width detection passed\n");
}

void test_flat_noise_no_false_detection(void) {
  printf("Testing flat noise → no false tonal detections...\n");

  float* profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* max_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* median_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* tonal_mask = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));

  // Flat noise everywhere
  for (uint32_t k = 0; k < TEST_SPECTRUM_SIZE; k++) {
    profile[k] = 0.01f;
    max_profile[k] = 0.01f;
    median_profile[k] = 0.01f;
  }

  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  int detections = 0;
  for (uint32_t k = 0; k < TEST_SPECTRUM_SIZE; k++) {
    if (tonal_mask[k] > 0.0f) {
      detections++;
    }
  }

  if (detections > 0) {
    fprintf(stderr, "FAIL: Flat noise produced %d false detections\n",
            detections);
    exit(1);
  }

  free(profile);
  free(max_profile);
  free(median_profile);
  free(tonal_mask);
  printf("✓ No false detections on flat noise\n");
}

void test_off_center_peak_interpolation(void) {
  printf("Testing off-center peak detected via parabolic interpolation...\n");

  float* profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* max_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* median_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* tonal_mask = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));

  float noise_floor = 0.01f;

  // Simulate a tone between bins: peak at bin 100, but energy leaks
  // asymmetrically — more energy in bin 101 than bin 99
  for (uint32_t k = 0; k < TEST_SPECTRUM_SIZE; k++) {
    profile[k] = noise_floor;
    max_profile[k] = noise_floor;
    median_profile[k] = noise_floor;
  }

  // Create asymmetric peak (true frequency slightly above bin 100)
  profile[99] = 0.03f;
  profile[100] = 0.08f;
  profile[101] = 0.06f; // Asymmetric leakage
  max_profile[99] = 0.03f;
  max_profile[100] = 0.08f;
  max_profile[101] = 0.06f;
  median_profile[99] = 0.03f;
  median_profile[100] = 0.08f;
  median_profile[101] = 0.06f;

  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  // The peak bin should still be detected
  if (tonal_mask[100] <= 0.0f) {
    fprintf(stderr, "FAIL: Off-center peak at bin 100 not detected\n");
    exit(1);
  }

  // The asymmetric side (bin 101) should have higher mask than bin 99
  // due to parabolic bias
  printf("  Bin 99: %.3f, Bin 100: %.3f, Bin 101: %.3f\n", tonal_mask[99],
         tonal_mask[100], tonal_mask[101]);

  if (tonal_mask[101] <= 0.0f) {
    fprintf(stderr, "FAIL: Asymmetric sideband not detected at bin 101\n");
    exit(1);
  }

  free(profile);
  free(max_profile);
  free(median_profile);
  free(tonal_mask);
  printf("✓ Off-center peak interpolation passed\n");
}

void test_adaptive_mode_detection(void) {
  printf(
      "Testing adaptive mode detection (unlearned noise profile, fallback to "
      "profile)...\n");

  float* profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* max_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* median_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* tonal_mask = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));

  // max_profile and median_profile are kept at 0.0f (unlearned)
  // profile (representing the active noise spectrum) has a tone at 1 kHz
  int bin = freq_to_bin(1000.0f);
  float noise_floor = 0.01f;

  for (uint32_t k = 0; k < TEST_SPECTRUM_SIZE; k++) {
    profile[k] = noise_floor;
  }
  profile[bin] = 0.1f;
  profile[bin - 1] = 0.03f;
  profile[bin + 1] = 0.03f;

  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  if (tonal_mask[bin] <= 0.0f) {
    fprintf(stderr,
            "FAIL: Adaptive mode fallback failed to detect tone at bin %d\n",
            bin);
    exit(1);
  }
  printf("  Tone detected in adaptive mode at bin %d (mask=%.3f) ✓\n", bin,
         tonal_mask[bin]);

  free(profile);
  free(max_profile);
  free(median_profile);
  free(tonal_mask);
  printf("✓ Adaptive mode detection passed\n");
}

int main(void) {
  printf("=== Tonal Detector Tests ===\n\n");

  test_null_inputs();
  test_flat_noise_no_false_detection();
  test_low_frequency_hum_detection();
  test_high_frequency_detection();
  test_adaptive_width();
  test_off_center_peak_interpolation();
  test_adaptive_mode_detection();

  printf("\n=== All tonal detector tests passed ===\n");
  return 0;
}
