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
  // colored noise floor using boundary-safe windowing (no DC padding duplication).
  const int win_size = TONAL_MEDIAN_FILTER_WINDOW;
  const int half_win = win_size / 2;
  float win_buf[TONAL_MEDIAN_FILTER_WINDOW];

  for (uint32_t k = 0U; k < size; k++) {
    int start_idx = (int)k - half_win;
    if (start_idx < 0) start_idx = 0;
    int end_idx = (int)k + half_win;
    if (end_idx >= (int)size) end_idx = (int)size - 1;
    int count = end_idx - start_idx + 1;

    for (int i = 0; i < count; i++) {
      win_buf[i] = detection_profile[start_idx + i];
    }
    insertion_sort(win_buf, count);
    tonal_mask[k] = win_buf[count / 2];
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

uint32_t tonal_detector_get_peaks(const float* tonal_mask, uint32_t size,
                                  uint32_t sample_rate, uint32_t fft_size,
                                  float* peak_freqs_hz, uint32_t max_peaks) {
  if (!tonal_mask || !peak_freqs_hz || size < 5 || sample_rate == 0 ||
      fft_size == 0 || max_peaks == 0) {
    return 0;
  }

  typedef struct {
    float freq_hz;
    float strength;
  } TonalPeakCandidate;

  TonalPeakCandidate candidates[MAX_TONAL_PEAKS_REPORTED];
  uint32_t candidate_count = 0;

  const float bin_width_hz = (float)sample_rate / (float)fft_size;

  for (uint32_t k = 1U; k < size - 1U; k++) {
    float mask_val = tonal_mask[k];

    // Peak center must be a local maximum above the significance threshold
    if (mask_val >= TONAL_PEAK_MIN_SIGNIFICANCE &&
        mask_val > tonal_mask[k - 1] && mask_val > tonal_mask[k + 1]) {

      // Parabolic interpolation for sub-bin frequency accuracy
      float left = tonal_mask[k - 1];
      float right = tonal_mask[k + 1];
      float denom = left - 2.0f * mask_val + right;
      float delta = 0.0f;
      if (fabsf(denom) > 1e-9f) {
        delta = 0.5f * (left - right) / denom;
        if (delta < -0.5f) delta = -0.5f;
        if (delta > 0.5f) delta = 0.5f;
      }

      float peak_bin = (float)k + delta;
      float freq_hz = peak_bin * bin_width_hz;

      if (freq_hz >= 20.0f && freq_hz <= (float)sample_rate * 0.48f) {
        if (candidate_count < MAX_TONAL_PEAKS_REPORTED) {
          candidates[candidate_count].freq_hz = freq_hz;
          candidates[candidate_count].strength = mask_val;
          candidate_count++;
        }
      }
    }
  }

  // Sort candidates by strength descending
  for (uint32_t i = 1; i < candidate_count; i++) {
    TonalPeakCandidate key = candidates[i];
    int j = (int)i - 1;
    while (j >= 0 && candidates[j].strength < key.strength) {
      candidates[j + 1] = candidates[j];
      j--;
    }
    candidates[j + 1] = key;
  }

  // Limit to max_peaks
  uint32_t return_count = candidate_count < max_peaks ? candidate_count : max_peaks;

  // Re-sort selected peaks by frequency ascending for ordered output
  for (uint32_t i = 1; i < return_count; i++) {
    TonalPeakCandidate key = candidates[i];
    int j = (int)i - 1;
    while (j >= 0 && candidates[j].freq_hz > key.freq_hz) {
      candidates[j + 1] = candidates[j];
      j--;
    }
    candidates[j + 1] = key;
  }

  for (uint32_t i = 0; i < return_count; i++) {
    peak_freqs_hz[i] = candidates[i].freq_hz;
  }

  return return_count;
}

uint32_t tonal_detector_get_peaks_from_profile(const float* profile, uint32_t size,
                                                uint32_t sample_rate, uint32_t fft_size,
                                                float* peak_freqs_hz, uint32_t max_peaks) {
  if (!profile || !peak_freqs_hz || size < 5 || sample_rate == 0 ||
      fft_size == 0 || max_peaks == 0) {
    return 0;
  }

  float* temp_mask = (float*)calloc(size, sizeof(float));
  if (!temp_mask) {
    return 0;
  }

  detect_tonal_components(profile, profile, profile, size, sample_rate, fft_size, temp_mask);

  uint32_t count = tonal_detector_get_peaks(temp_mask, size, sample_rate, fft_size, peak_freqs_hz, max_peaks);

  free(temp_mask);
  return count;
}
