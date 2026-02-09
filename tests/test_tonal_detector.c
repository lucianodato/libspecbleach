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

void test_sideband_spread_scaling(void) {
  printf("Testing sideband spread scales with frequency...\n");

  float* profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* max_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* median_profile = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));
  float* tonal_mask = (float*)calloc(TEST_SPECTRUM_SIZE, sizeof(float));

  float noise_floor = 0.01f;

  // Test low-frequency tone (100 Hz) — should have wider spread
  int bin_low = freq_to_bin(100.0f);
  int peak_bins_low[] = {bin_low};
  float peak_amps_low[] = {0.1f};

  create_test_profile(profile, TEST_SPECTRUM_SIZE, noise_floor, peak_bins_low,
                      peak_amps_low, 1);
  create_test_profile(max_profile, TEST_SPECTRUM_SIZE, noise_floor,
                      peak_bins_low, peak_amps_low, 1);
  create_test_profile(median_profile, TEST_SPECTRUM_SIZE, noise_floor,
                      peak_bins_low, peak_amps_low, 1);

  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  // Count bins with non-zero mask around low-freq peak
  int low_spread = 0;
  for (int i = bin_low - 6; i <= bin_low + 6; i++) {
    if (i >= 0 && i < TEST_SPECTRUM_SIZE && tonal_mask[i] > 0.0f) {
      low_spread++;
    }
  }

  // Test high-frequency tone (8 kHz) — should have narrower spread
  int bin_high = freq_to_bin(8000.0f);
  int peak_bins_high[] = {bin_high};
  float peak_amps_high[] = {0.1f};

  memset(tonal_mask, 0, TEST_SPECTRUM_SIZE * sizeof(float));
  create_test_profile(profile, TEST_SPECTRUM_SIZE, noise_floor, peak_bins_high,
                      peak_amps_high, 1);
  create_test_profile(max_profile, TEST_SPECTRUM_SIZE, noise_floor,
                      peak_bins_high, peak_amps_high, 1);
  create_test_profile(median_profile, TEST_SPECTRUM_SIZE, noise_floor,
                      peak_bins_high, peak_amps_high, 1);

  detect_tonal_components(profile, max_profile, median_profile,
                          TEST_SPECTRUM_SIZE, TEST_SAMPLE_RATE, TEST_FFT_SIZE,
                          tonal_mask);

  int high_spread = 0;
  for (int i = bin_high - 6; i <= bin_high + 6; i++) {
    if (i >= 0 && i < TEST_SPECTRUM_SIZE && tonal_mask[i] > 0.0f) {
      high_spread++;
    }
  }

  printf("  Low-freq spread: %d bins, High-freq spread: %d bins\n", low_spread,
         high_spread);

  if (low_spread <= high_spread) {
    fprintf(stderr,
            "FAIL: Low-freq spread (%d) should be wider than high-freq (%d)\n",
            low_spread, high_spread);
    exit(1);
  }

  free(profile);
  free(max_profile);
  free(median_profile);
  free(tonal_mask);
  printf("✓ Sideband spread scaling passed\n");
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

int main(void) {
  printf("=== Tonal Detector Tests ===\n\n");

  test_null_inputs();
  test_flat_noise_no_false_detection();
  test_low_frequency_hum_detection();
  test_high_frequency_detection();
  test_sideband_spread_scaling();
  test_off_center_peak_interpolation();

  printf("\n=== All tonal detector tests passed ===\n");
  return 0;
}
