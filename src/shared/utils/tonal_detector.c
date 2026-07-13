/*
libspecbleach - A spectral processing library
*/

#include "shared/utils/tonal_detector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"

static void insertion_sort(float* arr, int n) {
  for (int i = 1; i < n; i++) {
    float key = arr[i];
    int j = i - 1;
    while (j >= 0 && arr[j] > key) {
      arr[j + 1] = arr[j];
      j--;
    }
    arr[j + 1] = key;
  }
}

void detect_tonal_components(const float* profile, const float* max_profile,
                             const float* median_profile, uint32_t size,
                             uint32_t sample_rate, uint32_t fft_size,
                             float* tonal_mask) {
  if (!profile || !tonal_mask || size < 5 || !max_profile || !median_profile ||
      sample_rate == 0 || fft_size == 0) {
    return;
  }

  // Check if a pre-learned profile is available
  bool profile_learned = false;
  for (uint32_t k = 0U; k < size; k++) {
    if (median_profile[k] > 0.0f) {
      profile_learned = true;
      break;
    }
  }

  const float* detection_profile = profile_learned ? max_profile : profile;

  // 1. Perform frequency-domain median filtering to estimate the broadband
  // colored noise floor.
  //    To satisfy real-time safety constraints (no dynamic allocation), we use
  //    the output tonal_mask array as temporary storage for the broadband floor
  //    during the first pass.
  const int win_size = TONAL_MEDIAN_FILTER_WINDOW;
  const int half_win = win_size / 2;
  float win_buf[TONAL_MEDIAN_FILTER_WINDOW];

  for (uint32_t k = 0U; k < size; k++) {
    int n_count = 0;
    for (int i = -half_win; i <= half_win; i++) {
      int idx = (int)k + i;
      // Clamp boundaries
      if (idx < 0) {
        idx = 0;
      }
      if (idx >= (int)size) {
        idx = (int)size - 1;
      }
      win_buf[n_count++] = detection_profile[idx];
    }
    insertion_sort(win_buf, win_size);
    tonal_mask[k] = win_buf[half_win];
  }

  // 2. Compute the tonal mask based on the ratio of the spectrum to the
  // broadband floor.
  //    In this pass, we overwrite the tonal_mask in-place.
  for (uint32_t k = 0U; k < size; k++) {
    float floor_val = tonal_mask[k];
    float peak_val = detection_profile[k];

    // Prominence Guard: Peak must stand out in absolute terms
    if (peak_val - floor_val < MIN_PEAK_PROMINENCE) {
      tonal_mask[k] = 0.0f;
      continue;
    }

    float ratio = peak_val / (floor_val + 1e-20f);
    if (ratio > PEAK_THRESHOLD) {
      // In learned mode, the profile is captured from a noise-only segment, so
      // any peak is guaranteed to be a hum component. We do not need a
      // stationarity check.
      float stationarity_weight = 1.0f;

      // Calculate the fraction of energy belonging to the tone
      float tonal_energy = peak_val - floor_val;
      float tonal_fraction = tonal_energy / (peak_val + 1e-20f);

      // Final mask value is the product of stationarity and tonal fraction
      tonal_mask[k] = stationarity_weight * tonal_fraction;
    } else {
      tonal_mask[k] = 0.0f;
    }
  }
}
