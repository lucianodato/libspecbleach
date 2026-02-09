/*
libspecbleach - A spectral processing library
*/

#include "tonal_reducer.h"
#include <math.h>

void tonal_reducer_apply(float* alpha, const float* tonal_mask,
                         uint32_t real_spectrum_size,
                         float tonal_reduction_amount) {
  if (!alpha || !tonal_mask || real_spectrum_size == 0 ||
      tonal_reduction_amount >= 1.0f) {
    return;
  }

  // Alpha needed for the Wiener filter to produce the target gain:
  //   gain ≈ 1 - alpha  (when spectrum ≈ noise, i.e. noise-only)
  //   alpha = 1 - tonal_reduction_amount
  //
  // E.g. -48dB → tonal_reduction_amount = 0.004 → alpha_needed = 0.996
  // Each tonal bin gets the SAME relative dB reduction from its own level.
  float alpha_needed = 1.0f - tonal_reduction_amount;

  for (uint32_t k = 0U; k < real_spectrum_size; k++) {
    if (tonal_mask[k] <= 0.0f) {
      continue;
    }

    // Weight by mask strength (partial tonal bins get proportional boost)
    float target_alpha = alpha_needed * tonal_mask[k];

    // Only boost, never reduce (preserve existing suppression intent)
    alpha[k] = fmaxf(alpha[k], target_alpha);
  }
}
