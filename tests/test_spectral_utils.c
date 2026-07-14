/*
libspecbleach - A spectral processing library

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "shared/utils/spectral_utils.h"
#include <math.h>
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

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_get_fft_window(void) {
  printf("Testing get_fft_window...\n");

  uint32_t fft_size = 64;
  float window[64];

  // Test NULL and size 0 checks
  TEST_ASSERT(!get_fft_window(NULL, fft_size, HANN_WINDOW),
              "NULL window check");
  TEST_ASSERT(!get_fft_window(window, 0, HANN_WINDOW), "Size 0 check");

  // Test Hann
  TEST_ASSERT(get_fft_window(window, fft_size, HANN_WINDOW), "HANN Window");
  TEST_FLOAT_CLOSE(window[0], 0.0f, 1e-6f);
  TEST_ASSERT(window[fft_size / 2] > 0.9f, "Hann middle check");

  // Test Hamming
  TEST_ASSERT(get_fft_window(window, fft_size, HAMMING_WINDOW),
              "HAMMING Window");
  TEST_ASSERT(window[0] > 0.07f && window[0] < 0.09f, "Hamming endpoint check");

  // Test Blackman
  TEST_ASSERT(get_fft_window(window, fft_size, BLACKMAN_WINDOW),
              "BLACKMAN Window");
  TEST_FLOAT_CLOSE(window[0], 0.0f, 1e-6f);

  // Test Vorbis
  TEST_ASSERT(get_fft_window(window, fft_size, VORBIS_WINDOW), "VORBIS Window");
  TEST_FLOAT_CLOSE(window[0], 0.0f, 1e-6f);

  // Test default window type (no-op loop but returns true)
  TEST_ASSERT(get_fft_window(window, fft_size, (WindowTypes)99),
              "Invalid window type check");

  printf("✓ get_fft_window tests passed\n");
}

void test_initialize_spectrum_with_value(void) {
  printf("Testing initialize_spectrum_with_value...\n");

  uint32_t size = 10;
  float spectrum[10];

  // Test NULL and size 0 checks
  TEST_ASSERT(!initialize_spectrum_with_value(NULL, size, 2.5f),
              "NULL spectrum check");
  TEST_ASSERT(!initialize_spectrum_with_value(spectrum, 0, 2.5f),
              "Size 0 check");

  // Test valid initialization
  TEST_ASSERT(initialize_spectrum_with_value(spectrum, size, 2.5f),
              "Init spectrum");
  for (uint32_t i = 0; i < size; i++) {
    TEST_FLOAT_CLOSE(spectrum[i], 2.5f, 1e-6f);
  }

  printf("✓ initialize_spectrum_with_value tests passed\n");
}

void test_min_max_spectrum_float(void) {
  printf("Testing min_spectrum_float and max_spectrum_float...\n");

  uint32_t size = 4;
  float s1[] = {1.0f, 10.0f, 5.0f, 2.0f};
  float s2[] = {2.0f, 5.0f, 5.0f, 1.0f};

  // Test NULL and size 0 checks for min_spectrum_float
  TEST_ASSERT(!min_spectrum_float(NULL, s2, size), "min_spectrum NULL 1");
  TEST_ASSERT(!min_spectrum_float(s1, NULL, size), "min_spectrum NULL 2");
  TEST_ASSERT(!min_spectrum_float(s1, s2, 0), "min_spectrum size 0");

  // Test NULL and size 0 checks for max_spectrum_float
  TEST_ASSERT(!max_spectrum_float(NULL, s2, size), "max_spectrum NULL 1");
  TEST_ASSERT(!max_spectrum_float(s1, NULL, size), "max_spectrum NULL 2");
  TEST_ASSERT(!max_spectrum_float(s1, s2, 0), "max_spectrum size 0");

  // Test operations
  float s1_min[] = {1.0f, 10.0f, 5.0f, 2.0f};
  TEST_ASSERT(min_spectrum_float(s1_min, s2, size), "min_spectrum execution");
  TEST_FLOAT_CLOSE(s1_min[0], 1.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_min[1], 5.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_min[2], 5.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_min[3], 1.0f, 1e-6f);

  float s1_max[] = {1.0f, 10.0f, 5.0f, 2.0f};
  TEST_ASSERT(max_spectrum_float(s1_max, s2, size), "max_spectrum execution");
  TEST_FLOAT_CLOSE(s1_max[0], 2.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_max[1], 10.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_max[2], 5.0f, 1e-6f);
  TEST_FLOAT_CLOSE(s1_max[3], 2.0f, 1e-6f);

  printf("✓ min/max spectrum float tests passed\n");
}

void test_fft_frequency_conversions(void) {
  printf("Testing fft_bin_to_freq and freq_to_fft_bin...\n");

  uint32_t sr = 44100;
  uint32_t fft_size = 1024;

  // DC bin
  TEST_FLOAT_CLOSE(fft_bin_to_freq(0, sr, fft_size), 0.0f, 1e-6f);
  TEST_ASSERT(freq_to_fft_bin(0.0f, sr, fft_size) == 0, "DC frequency check");

  // Nyquist bin
  float nyquist_freq = (float)sr / 2.0f;
  TEST_FLOAT_CLOSE(fft_bin_to_freq(512, sr, fft_size), nyquist_freq, 1e-6f);
  TEST_ASSERT(freq_to_fft_bin(nyquist_freq, sr, fft_size) == 512,
              "Nyquist check");

  printf("✓ fft frequency conversion tests passed\n");
}

void test_get_rolling_mean_spectrum(void) {
  printf("Testing get_rolling_mean_spectrum...\n");

  uint32_t size = 4;
  float current[] = {1.0f, 2.0f, 3.0f, 4.0f};
  float avg[4] = {0.0f};

  // Test NULL and size 0 checks
  TEST_ASSERT(!get_rolling_mean_spectrum(NULL, current, 1, size),
              "rolling_mean NULL avg");
  TEST_ASSERT(!get_rolling_mean_spectrum(avg, NULL, 1, size),
              "rolling_mean NULL current");
  TEST_ASSERT(!get_rolling_mean_spectrum(avg, current, 1, 0),
              "rolling_mean size 0");

  // Test blocks == 1 (should just copy)
  TEST_ASSERT(get_rolling_mean_spectrum(avg, current, 1, size),
              "rolling_mean blocks=1");
  for (uint32_t i = 0; i < size; i++) {
    TEST_FLOAT_CLOSE(avg[i], current[i], 1e-6f);
  }

  // Test blocks > 1
  float current2[] = {3.0f, 4.0f, 5.0f, 6.0f};
  TEST_ASSERT(get_rolling_mean_spectrum(avg, current2, 2, size),
              "rolling_mean blocks=2");
  // avg = avg + (current2 - avg) / 2 = 1 + (3-1)/2 = 2, etc.
  TEST_FLOAT_CLOSE(avg[0], 2.0f, 1e-6f);
  TEST_FLOAT_CLOSE(avg[1], 3.0f, 1e-6f);
  TEST_FLOAT_CLOSE(avg[2], 4.0f, 1e-6f);
  TEST_FLOAT_CLOSE(avg[3], 5.0f, 1e-6f);

  printf("✓ get_rolling_mean_spectrum tests passed\n");
}

void test_get_rolling_median_spectrum(void) {
  printf("Testing get_rolling_median_spectrum...\n");

  uint32_t size = 3;
  float frame1[] = {1.0f, 10.0f, 5.0f};
  float frame2[] = {3.0f, 2.0f, 8.0f};
  float frame3[] = {2.0f, 5.0f, 1.0f};
  const float* history3[] = {frame1, frame2, frame3};

  float median[3] = {0.0f};

  // Test NULL and size 0 checks
  TEST_ASSERT(!get_rolling_median_spectrum(NULL, history3, 3, size),
              "rolling_median NULL output");
  TEST_ASSERT(!get_rolling_median_spectrum(median, NULL, 3, size),
              "rolling_median NULL history");
  TEST_ASSERT(!get_rolling_median_spectrum(median, history3, 3, 0),
              "rolling_median size 0");

  // Odd number of blocks (3)
  // Sorted inputs for bin 0: 1.0, 2.0, 3.0 -> median is 2.0
  // Sorted inputs for bin 1: 2.0, 5.0, 10.0 -> median is 5.0
  // Sorted inputs for bin 2: 1.0, 5.0, 8.0 -> median is 5.0
  TEST_ASSERT(get_rolling_median_spectrum(median, history3, 3, size),
              "rolling_median blocks=3");
  TEST_FLOAT_CLOSE(median[0], 2.0f, 1e-6f);
  TEST_FLOAT_CLOSE(median[1], 5.0f, 1e-6f);
  TEST_FLOAT_CLOSE(median[2], 5.0f, 1e-6f);

  // Even number of blocks (2)
  const float* history2[] = {frame1, frame2};
  // Sorted inputs for bin 0: 1.0, 3.0 -> median is 2.0
  // Sorted inputs for bin 1: 2.0, 10.0 -> median is 6.0
  // Sorted inputs for bin 2: 5.0, 8.0 -> median is 6.5
  TEST_ASSERT(get_rolling_median_spectrum(median, history2, 2, size),
              "rolling_median blocks=2");
  TEST_FLOAT_CLOSE(median[0], 2.0f, 1e-6f);
  TEST_FLOAT_CLOSE(median[1], 6.0f, 1e-6f);
  TEST_FLOAT_CLOSE(median[2], 6.5f, 1e-6f);

  // Test NULL pointer in history fallback to 0.0f
  const float* history_null[] = {frame1, NULL, frame3};
  // Sorted inputs for bin 0: 0.0 (fallback), 1.0, 2.0 -> median is 1.0
  TEST_ASSERT(get_rolling_median_spectrum(median, history_null, 3, size),
              "rolling_median with NULL history element");
  TEST_FLOAT_CLOSE(median[0], 1.0f, 1e-6f);

  printf("✓ get_rolling_median_spectrum tests passed\n");
}

void test_smooth_spectrum(void) {
  printf("Testing smooth_spectrum...\n");

  float spec[] = {1.0f, 2.0f, 4.0f, 8.0f};
  uint32_t size = 4;

  // NULL check
  smooth_spectrum(NULL, size, 0.5f); // Should not crash

  // Size < 2 check
  float single_spec[] = {1.0f};
  smooth_spectrum(single_spec, 1, 0.5f);
  TEST_FLOAT_CLOSE(single_spec[0], 1.0f, 1e-6f);

  // Smoothing factor <= 0.0f check
  float spec_no_change[] = {1.0f, 2.0f, 4.0f, 8.0f};
  smooth_spectrum(spec_no_change, size, 0.0f);
  for (uint32_t i = 0; i < size; i++) {
    TEST_FLOAT_CLOSE(spec_no_change[i], spec[i], 1e-6f);
  }

  // Normal smoothing execution
  smooth_spectrum(spec, size, 0.5f);

  // Expected values calculation:
  // bin 0: avg(1.0, 2.0)*0.5 + 1.0*0.5 = 1.5*0.5 + 0.5 = 1.25
  // bin 1: avg(1.0, 2.0, 4.0)*0.5 + 2.0*0.5 = 2.3333*0.5 + 1.0 = 2.16667
  // bin 2: avg(2.0, 4.0, 8.0)*0.5 + 4.0*0.5 = 4.6667*0.5 + 2.0 = 4.33333
  // bin 3: avg(4.0, 8.0)*0.5 + 8.0*0.5 = 6.0*0.5 + 4.0 = 7.0
  TEST_FLOAT_CLOSE(spec[0], 1.25f, 1e-5f);
  TEST_FLOAT_CLOSE(spec[1], 2.16667f, 1e-5f);
  TEST_FLOAT_CLOSE(spec[2], 4.33333f, 1e-5f);
  TEST_FLOAT_CLOSE(spec[3], 7.00000f, 1e-5f);

  printf("✓ smooth_spectrum tests passed\n");
}

void test_interpolate_spectrum_gaps(void) {
  printf("Testing interpolate_spectrum_gaps...\n");

  float threshold = 0.1f;

  // NULL and size < 3 checks
  interpolate_spectrum_gaps(NULL, 10, threshold); // Should not crash
  float short_spec[] = {0.05f, 0.5f};
  interpolate_spectrum_gaps(short_spec, 2, threshold);
  TEST_FLOAT_CLOSE(short_spec[0], 0.05f, 1e-6f);

  // No gaps
  float no_gaps[] = {0.5f, 0.6f, 0.7f, 0.8f};
  interpolate_spectrum_gaps(no_gaps, 4, threshold);
  TEST_FLOAT_CLOSE(no_gaps[1], 0.6f, 1e-6f);

  // Gap in the middle (interpolated)
  float spec[] = {1.0f, 0.05f, 3.0f}; // Gap at index 1
  interpolate_spectrum_gaps(spec, 3, threshold);
  // Interpolated: start_val = spec[0] = 1.0, end_val = spec[2] = 3.0.
  // step = (3.0 - 1.0) / (2 - 0) = 1.0
  // spec[1] = 1.0 + 1.0 * (1 - 0) = 2.0
  TEST_FLOAT_CLOSE(spec[1], 2.0f, 1e-6f);

  // Gap extending to the end (should not be interpolated)
  float spec_end[] = {1.0f, 0.05f, 0.05f}; // Gaps at 1 and 2
  interpolate_spectrum_gaps(spec_end, 3, threshold);
  TEST_FLOAT_CLOSE(spec_end[1], 0.05f, 1e-6f);
  TEST_FLOAT_CLOSE(spec_end[2], 0.05f, 1e-6f);

  printf("✓ interpolate_spectrum_gaps tests passed\n");
}

void test_get_morphed_profile(void) {
  printf("Testing get_morphed_profile...\n");

  uint32_t size = 2;
  float mean[] = {2.0f, 4.0f};
  float median[] = {1.0f, 3.0f};
  float max[] = {4.0f, 6.0f};
  float min[] = {1.5f, 2.5f};

  float output[2] = {0.0f};

  // Test NULL and size 0 checks
  TEST_ASSERT(!get_morphed_profile(NULL, mean, median, max, min, size, 0.5f),
              "morphed NULL output");
  TEST_ASSERT(!get_morphed_profile(output, NULL, median, max, min, size, 0.5f),
              "morphed NULL mean");
  TEST_ASSERT(!get_morphed_profile(output, mean, median, max, min, 0, 0.5f),
              "morphed size 0");

  // Test positive aggressiveness (morph mean -> max)
  // aggressiveness = 0.5
  // output[0] = mean[0]*(1-0.5) + max[0]*0.5 = 2.0*0.5 + 4.0*0.5 = 3.0
  // Clamped with min[0] (1.5) -> 3.0
  // output[1] = mean[1]*(1-0.5) + max[1]*0.5 = 4.0*0.5 + 6.0*0.5 = 5.0
  // Clamped with min[1] (2.5) -> 5.0
  TEST_ASSERT(get_morphed_profile(output, mean, median, max, min, size, 0.5f),
              "morphed aggressiveness > 0");
  TEST_FLOAT_CLOSE(output[0], 3.0f, 1e-6f);
  TEST_FLOAT_CLOSE(output[1], 5.0f, 1e-6f);

  // Test negative aggressiveness (morph mean -> median)
  // aggressiveness = -0.5, t = 0.5
  // output[0] = mean[0]*(1-0.5) + median[0]*0.5 = 2.0*0.5 + 1.0*0.5 = 1.5
  // Clamped with min[0] (1.5) -> 1.5
  // output[1] = mean[1]*(1-0.5) + median[1]*0.5 = 4.0*0.5 + 3.0*0.5 = 3.5
  // Clamped with min[1] (2.5) -> 3.5
  TEST_ASSERT(get_morphed_profile(output, mean, median, max, min, size, -0.5f),
              "morphed aggressiveness < 0");
  TEST_FLOAT_CLOSE(output[0], 1.5f, 1e-6f);
  TEST_FLOAT_CLOSE(output[1], 3.5f, 1e-6f);

  // Test negative aggressiveness clamping with min profile
  // aggressiveness = -1.0, t = 1.0
  // output[0] = mean[0]*0.0 + median[0]*1.0 = 1.0
  // Clamped with min[0] (1.5) -> should clamp to 1.5
  TEST_ASSERT(get_morphed_profile(output, mean, median, max, min, size, -1.0f),
              "morphed aggressiveness clamping");
  TEST_FLOAT_CLOSE(output[0], 1.5f, 1e-6f);

  printf("✓ get_morphed_profile tests passed\n");
}

int main(void) {
  printf("Running spectral_utils tests...\n\n");

  test_get_fft_window();
  test_initialize_spectrum_with_value();
  test_min_max_spectrum_float();
  test_fft_frequency_conversions();
  test_get_rolling_mean_spectrum();
  test_get_rolling_median_spectrum();
  test_smooth_spectrum();
  test_interpolate_spectrum_gaps();
  test_get_morphed_profile();

  printf("\n✅ All spectral_utils tests passed!\n");
  return 0;
}
