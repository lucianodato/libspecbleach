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

#include "shared/utils/tonal_detector.h"
#include <math.h>
#include <stdlib.h>

#include "shared/configurations.h"
#include "shared/utils/simd_utils.h"

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

void detect_tonal_components(const float* profile, const float* median_profile,
                             uint32_t size, uint32_t sample_rate,
                             uint32_t fft_size, float* tonal_mask,
                             uint32_t* deque_workspace) {
  if (!profile || !tonal_mask || size < 5 || !median_profile ||
      sample_rate == 0 || fft_size == 0 ||
      (!deque_workspace && size > TONAL_DETECTOR_DEQUE_CAPACITY)) {
    return;
  }

  sb_simd_state_t old_simd_state = sb_simd_enable_ftz_daz();

  // Check if a pre-learned profile is available
  bool profile_learned = false;
  for (uint32_t k = 0U; k < size; k++) {
    if (median_profile[k] > 0.0f) {
      profile_learned = true;
      break;
    }
  }

  // The median profile is the most stable version of the noise profile.
  // We use it exclusively to find peaks, discarding max_profile completely.
  const float* detection_profile = profile_learned ? median_profile : profile;
  float base_threshold = PEAK_THRESHOLD;

  // 1. Calculate the absolute global maximum of the reference profile
  float global_max_val = 0.0f;
  for (uint32_t k = 0U; k < size; k++) {
    if (detection_profile[k] > global_max_val) {
      global_max_val = detection_profile[k];
    }
  }
  float global_max_power = global_max_val * global_max_val;

  // 2. Perform frequency-domain median filtering to estimate the broadband
  // colored noise floor using boundary-safe windowing (no DC padding
  // duplication).
  uint32_t half_win = TONAL_MEDIAN_FILTER_WINDOW / 2U;
  float win_buf[TONAL_MEDIAN_FILTER_WINDOW];

  for (uint32_t k = 0U; k < size; k++) {
    int start_idx = (int)k - (int)half_win;
    if (start_idx < 0) {
      start_idx = 0;
    }
    int end_idx = (int)k + (int)half_win;
    if (end_idx >= (int)size) {
      end_idx = (int)size - 1;
    }
    int count = end_idx - start_idx + 1;

    for (int i = 0; i < count; i++) {
      win_buf[i] = detection_profile[start_idx + i];
    }
    insertion_sort(win_buf, count);
    tonal_mask[k] = win_buf[count / 2];
  }

  // 2. Compute the tonal mask based on the ratio of the spectrum to the
  // broadband floor and octave-band relative dynamic range.
  //    In this pass, we overwrite the tonal_mask in-place.
  // 1-Octave sliding window max profile calculation via monotonic queue (O(N)
  // total)
  uint32_t deque_head = 0;
  uint32_t deque_tail = 0;
  uint32_t current_end = 0;

  uint32_t deque_stack[TONAL_DETECTOR_DEQUE_CAPACITY];
  uint32_t* deque = deque_workspace ? deque_workspace : deque_stack;

  for (uint32_t k = 0U; k < size; k++) {
    float floor_val = tonal_mask[k];
    float peak_val = detection_profile[k];

    // Prominence Guard: Peak must stand out in absolute terms
    if (peak_val - floor_val < MIN_PEAK_PROMINENCE) {
      tonal_mask[k] = 0.0f;
      continue;
    }

    // 1-Octave sliding window max profile calculation
    int start_octave = (int)((float)k * TONAL_OCTAVE_LOWER_RATIO);
    int end_octave = (int)((float)k * TONAL_OCTAVE_UPPER_RATIO);
    if (end_octave - start_octave < TONAL_OCTAVE_MIN_WIDTH_BINS) {
      start_octave = (int)k - TONAL_OCTAVE_FALLBACK_HALF_WIDTH;
      end_octave = (int)k + TONAL_OCTAVE_FALLBACK_HALF_WIDTH;
    }
    if (start_octave < 0) {
      start_octave = 0;
    }
    if (end_octave >= (int)size) {
      end_octave = (int)size - 1;
    }

    while (current_end <= (uint32_t)end_octave) {
      float val = detection_profile[current_end];
      while (deque_tail > deque_head &&
             detection_profile[deque[deque_tail - 1]] <= val) {
        deque_tail--;
      }
      deque[deque_tail++] = current_end++;
    }

    while (deque_head < deque_tail &&
           deque[deque_head] < (uint32_t)start_octave) {
      deque_head++;
    }

    float octave_max_val =
        (deque_head < deque_tail) ? detection_profile[deque[deque_head]] : 0.0f;

    float peak_power = peak_val * peak_val;
    float octave_max_power = octave_max_val * octave_max_val;

    // Reject candidates more than 20 dB below the max peak in their octave band
    if (peak_power < octave_max_power * TONAL_PEAK_MIN_OCTAVE_RELATIVE_POWER) {
      tonal_mask[k] = 0.0f;
      continue;
    }

    // Reject candidates that fall significantly below the global broadband peak
    // energy
    if (peak_power < global_max_power * TONAL_PEAK_MIN_GLOBAL_RELATIVE_POWER) {
      tonal_mask[k] = 0.0f;
      continue;
    }

    float ratio = peak_val / (floor_val + 1e-20f);

    // Since we are detecting directly on the median_profile (the most stable
    // version), we don't need a separate stationarity cross-check or dynamic
    // scaling.

    if (ratio > base_threshold) {
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

  sb_simd_restore_state(old_simd_state);
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

  for (uint32_t k = 2U; k < size - 2U; k++) {
    float mask_val = tonal_mask[k];

    float neighbor_avg = 0.5f * (tonal_mask[k - 1] + tonal_mask[k + 1]);
    float local_prominence = mask_val - neighbor_avg;

    float wider_neighbor_avg = 0.5f * (tonal_mask[k - 2] + tonal_mask[k + 2]);
    float wider_prominence = mask_val - wider_neighbor_avg;

    float approx_freq_hz = (float)k * bin_width_hz;

    // Scale wider prominence requirement for low frequencies where resolution
    // is poor
    float required_wider_prominence = TONAL_PEAK_MIN_WIDER_PROMINENCE;
    if (approx_freq_hz < 500.0f) {
      required_wider_prominence *= (approx_freq_hz / 500.0f);
    }

    // Peak center must be a local maximum above the significance threshold and
    // stand out from immediate neighbors
    if (mask_val >= TONAL_PEAK_MIN_SIGNIFICANCE &&
        mask_val > tonal_mask[k - 1] && mask_val >= tonal_mask[k + 1] &&
        local_prominence >= TONAL_PEAK_MIN_LOCAL_PROMINENCE &&
        wider_prominence >= required_wider_prominence) {

      // Parabolic interpolation for sub-bin frequency accuracy
      float left = tonal_mask[k - 1];
      float right = tonal_mask[k + 1];
      float denom = (left - (2.0f * mask_val)) + right;
      float delta = 0.0f;
      if (fabsf(denom) > 1e-9f) {
        delta = 0.5f * (left - right) / denom;
        if (delta < -0.5f) {
          delta = -0.5f;
        }
        if (delta > 0.5f) {
          delta = 0.5f;
        }
      }

      float peak_bin = (float)k + delta;
      float freq_hz = peak_bin * bin_width_hz;

      if (freq_hz >= TONAL_PEAK_MIN_FREQ_HZ &&
          freq_hz <= (float)sample_rate * TONAL_PEAK_NYQUIST_SAFETY_FACTOR) {
        if (candidate_count < MAX_TONAL_PEAKS_REPORTED) {
          candidates[candidate_count].freq_hz = freq_hz;
          candidates[candidate_count].strength = mask_val;
          candidate_count++;
        } else {
          uint32_t weakest_idx = 0;
          for (uint32_t c = 1; c < candidate_count; c++) {
            if (candidates[c].strength < candidates[weakest_idx].strength) {
              weakest_idx = c;
            }
          }
          if (mask_val > candidates[weakest_idx].strength) {
            candidates[weakest_idx].freq_hz = freq_hz;
            candidates[weakest_idx].strength = mask_val;
          }
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
  uint32_t return_count =
      candidate_count < max_peaks ? candidate_count : max_peaks;

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

uint32_t tonal_detector_get_peaks_from_profile(
    const float* profile, uint32_t size, uint32_t sample_rate,
    uint32_t fft_size, float* peak_freqs_hz, uint32_t max_peaks) {
  if (!profile || !peak_freqs_hz || size < 5 || sample_rate == 0 ||
      fft_size == 0 || max_peaks == 0) {
    return 0;
  }

  float* temp_mask = (float*)calloc(size, sizeof(float));
  if (!temp_mask) {
    return 0;
  }

  detect_tonal_components(profile, profile, size, sample_rate, fft_size,
                          temp_mask, NULL);

  uint32_t count = tonal_detector_get_peaks(temp_mask, size, sample_rate,
                                            fft_size, peak_freqs_hz, max_peaks);

  free(temp_mask);
  return count;
}
