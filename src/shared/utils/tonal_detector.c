/*
libspecbleach - A spectral processing library
*/

#include "shared/utils/tonal_detector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"

static int float_comparator(const void* a, const void* b) {
  float x = *(const float*)a;
  float y = *(const float*)b;
  return (x > y) - (x < y);
}

// Smooth linear interpolation between two breakpoints
static float lerp_by_freq(float freq_hz, float val_low, float val_high) {
  if (freq_hz <= LOW_FREQ_HZ) {
    return val_low;
  }
  if (freq_hz >= MID_FREQ_HZ) {
    return val_high;
  }
  float t = (freq_hz - LOW_FREQ_HZ) / (MID_FREQ_HZ - LOW_FREQ_HZ);
  return val_low + (t * (val_high - val_low));
}

void detect_tonal_components(const float* profile, const float* max_profile,
                             const float* median_profile, uint32_t size,
                             uint32_t sample_rate, uint32_t fft_size,
                             float* tonal_mask) {
  if (!profile || !tonal_mask || size < 5 || !max_profile || !median_profile ||
      sample_rate == 0 || fft_size == 0) {
    return;
  }

  memset(tonal_mask, 0, size * sizeof(float));

  // Hz per FFT bin
  const float bin_hz = (float)sample_rate / (float)fft_size;

  // We use the median profile for detection as it is much more stable than
  // Mean/Instantaneous
  const float* detection_profile = median_profile;

  for (uint32_t k = 2; k < size - 2; k++) {
    // 1. Robust Local Maximum Check (5-bin window to capture broader peaks)
    if (detection_profile[k] >= detection_profile[k - 1] &&
        detection_profile[k] >= detection_profile[k + 1] &&
        detection_profile[k] > detection_profile[k - 2] &&
        detection_profile[k] > detection_profile[k + 2]) {

      // Compute frequency for this bin
      float freq_hz = (float)k * bin_hz;

      // 2. Parabolic interpolation for sub-bin peak accuracy
      // Fit a parabola through log-magnitudes at (k-1, k, k+1)
      float log_prev = logf(detection_profile[k - 1] + 1e-20f);
      float log_peak = logf(detection_profile[k] + 1e-20f);
      float log_next = logf(detection_profile[k + 1] + 1e-20f);

      float denom_parab = log_prev - (2.0f * log_peak) + log_next;
      float delta = 0.0f;
      float true_amplitude = detection_profile[k];

      // Only interpolate if the parabola is concave (valid peak)
      if (denom_parab < -1e-10f) {
        delta = 0.5f * (log_prev - log_next) / denom_parab;
        // Clamp delta to [-0.5, 0.5]
        delta = fmaxf(-0.5f, fminf(0.5f, delta));
        true_amplitude =
            expf(log_peak - (0.25f * (log_prev - log_next) * delta));
      }

      // 3. Frequency-adaptive background estimation
      int bg_radius = (int)roundf(
          lerp_by_freq(freq_hz, (float)BG_RADIUS_LOW, (float)BG_RADIUS_HIGH));

      float neighbors[MAX_NEIGHBORS];
      int n_count = 0;
      for (int i = -bg_radius; i <= bg_radius; i++) {
        int idx = (int)k + i;
        if (i != 0 && idx >= 0 && idx < (int)size && n_count < MAX_NEIGHBORS) {
          neighbors[n_count++] = detection_profile[idx];
        }
      }

      if (n_count > 0) {
        // Sort to find local background median
        qsort(neighbors, n_count, sizeof(float), float_comparator);
        float avg_background = neighbors[n_count / 2];

        // 4. Stationarity Index (Spread)
        float spread = (max_profile[k] + 1e-9f) / (median_profile[k] + 1e-9f);

        // Stationarity Weight [0.0 - 1.0]
        float stationarity_weight = fmaxf(
            0.0f, 1.0f - ((spread - 1.0f) / (STATIONARITY_THRESHOLD - 1.0f)));

        // 5. Frequency-adaptive peak threshold
        float threshold_factor =
            lerp_by_freq(freq_hz, THRESHOLD_FACTOR_LOW, THRESHOLD_FACTOR_HIGH);

        float threshold = PEAK_THRESHOLD * threshold_factor *
                          (1.5f - (stationarity_weight * 0.5f));

        // Use parabolic-corrected amplitude for threshold comparison
        if (true_amplitude > avg_background * threshold) {
          // We found a tonal candidate
          tonal_mask[k] = stationarity_weight;

          // 6. Frequency-adaptive sideband spread
          int spread_bins = (int)roundf(
              lerp_by_freq(freq_hz, (float)SIDEBAND_LOW, (float)SIDEBAND_HIGH));

          // Apply asymmetric sideband spread based on parabolic offset
          for (int s = 1; s <= spread_bins; s++) {
            // Weight decreases with distance from peak
            float weight = stationarity_weight *
                           (1.0f - (float)s / ((float)spread_bins + 1.0f));

            // Bias the weight toward the side where the true peak sits
            float bias_minus = (delta < 0.0f)
                                   ? fminf(1.0f, 1.0f + (delta * 0.5f))
                                   : fmaxf(0.5f, 1.0f - (delta * 0.3f));
            float bias_plus = (delta > 0.0f)
                                  ? fminf(1.0f, (1.0f - (delta * 0.5f)) + 0.5f)
                                  : fmaxf(0.5f, 1.0f + (delta * 0.3f));

            int idx_minus = (int)k - s;
            int idx_plus = (int)k + s;

            if (idx_minus >= 0) {
              tonal_mask[idx_minus] =
                  fmaxf(tonal_mask[idx_minus], weight * bias_minus);
            }
            if (idx_plus < (int)size) {
              tonal_mask[idx_plus] =
                  fmaxf(tonal_mask[idx_plus], weight * bias_plus);
            }
          }
        }
      }
    }
  }
}
