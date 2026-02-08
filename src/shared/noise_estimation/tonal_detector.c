/*
libspecbleach - A spectral processing library
*/

#include "tonal_detector.h"
#include <math.h>
#include <string.h>

#define PEAK_THRESHOLD 1.58f        // ~4dB above neighbor average
#define STATIONARITY_THRESHOLD 2.5f // Ratio of Max/Median spread

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

      // 2. Broad Neighborhood average (up to 15-bin window) to estimate
      // background Boundary-aware: clamp the window to [0, size-1] We exclude
      // the central 5 bins [k-2...k+2] to avoid peak-induced bias
      float sum_bg = 0.0f;
      int count_bg = 0;
      for (int i = -7; i <= 7; i++) {
        int idx = (int)k + i;
        if (idx >= 0 && idx < (int)size) {
          if (i < -2 || i > 2) {
            sum_bg += detection_profile[idx];
            count_bg++;
          }
        }
      }

      if (count_bg > 0) {
        float avg_background = sum_bg / (float)count_bg;

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
          if (k > 0) {
            tonal_mask[k - 1] =
                fmaxf(tonal_mask[k - 1], stationarity_weight * 0.8f);
          }
          if (k < size - 1) {
            tonal_mask[k + 1] =
                fmaxf(tonal_mask[k + 1], stationarity_weight * 0.8f);
          }
          if (k > 1) {
            tonal_mask[k - 2] =
                fmaxf(tonal_mask[k - 2], stationarity_weight * 0.4f);
          }
          if (k < size - 2) {
            tonal_mask[k + 2] =
                fmaxf(tonal_mask[k + 2], stationarity_weight * 0.4f);
          }
        }
      }
    }
  }
}
