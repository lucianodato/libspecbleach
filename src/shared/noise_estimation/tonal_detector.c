/*
libspecbleach - A spectral processing library
*/

#include "tonal_detector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PEAK_THRESHOLD 1.41f        // ~3dB above neighbor background
#define STATIONARITY_THRESHOLD 2.5f // Ratio of Max/Median spread

static int float_comparator(const void* a, const void* b) {
  float x = *(const float*)a;
  float y = *(const float*)b;
  return (x > y) - (x < y);
}

void detect_tonal_components(const float* profile, const float* max_profile,
                             const float* median_profile, uint32_t size,
                             float* tonal_mask) {
  if (!profile || !tonal_mask || size < 5 || !max_profile || !median_profile) {
    return;
  }

  memset(tonal_mask, 0, size * sizeof(float));

  // We use the median profile for detection as it is much more stable than
  // Mean/Instantaneous
  const float* detection_profile = median_profile;

  for (uint32_t k = 2; k < size - 2; k++) {
    // 1. Robust Local Maximum Check (5-bin window to capture broader peaks)
    if (detection_profile[k] >= detection_profile[k - 1] &&
        detection_profile[k] >= detection_profile[k + 1] &&
        detection_profile[k] > detection_profile[k - 2] &&
        detection_profile[k] > detection_profile[k + 2]) {

      // 2. Robust Local Background Estimation
      // We use a 15-bin median of the neighborhood (excluding the peak itself)
      // to avoid nearby harmonics biasing the background estimate.
      float neighbors[14];
      int n_count = 0;
      for (int i = -7; i <= 7; i++) {
        int idx = (int)k + i;
        if (i != 0 && idx >= 0 && idx < (int)size) {
          neighbors[n_count++] = detection_profile[idx];
        }
      }

      if (n_count > 0) {
        // Sort to find local background median
        qsort(neighbors, n_count, sizeof(float), float_comparator);
        float avg_background = neighbors[n_count / 2];

        // 3. Stationarity Index (Spread)
        float spread = (max_profile[k] + 1e-9f) / (median_profile[k] + 1e-9f);

        // We calculate a 'Stationarity Weight' [0.0 - 1.0]
        // A slightly more relaxed limit for stationarity
        float stationarity_weight = fmaxf(
            0.0f, 1.0f - ((spread - 1.0f) / (STATIONARITY_THRESHOLD - 1.0f)));

        // 4. Adaptive Peak Threshold
        // Relaxed multiplier (1.5x max) for better sensitivity to subtle tones
        float threshold =
            PEAK_THRESHOLD * (1.5f - (stationarity_weight * 0.5f));

        if (detection_profile[k] > avg_background * threshold) {
          // We found a tonal candidate.
          tonal_mask[k] = stationarity_weight;

          // Spread the influence to neighboring bins to cover sidebands
          // (leakage) Center: 1.0, ±1: 0.8, ±2: 0.4
          tonal_mask[k - 1] =
              fmaxf(tonal_mask[k - 1], stationarity_weight * 0.8f);
          tonal_mask[k + 1] =
              fmaxf(tonal_mask[k + 1], stationarity_weight * 0.8f);
          tonal_mask[k - 2] =
              fmaxf(tonal_mask[k - 2], stationarity_weight * 0.4f);
          tonal_mask[k + 2] =
              fmaxf(tonal_mask[k + 2], stationarity_weight * 0.4f);
        }
      }
    }
  }
}
