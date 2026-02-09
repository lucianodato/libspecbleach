/*
libspecbleach - A spectral processing library
*/

#include "tonal_reducer.h"
#include "shared/utils/tonal_detector.h"
#include <math.h>
#include <stdlib.h>

struct TonalReducer {
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t fft_size;
  float* tonal_mask;
};

TonalReducer* tonal_reducer_initialize(uint32_t real_spectrum_size,
                                       uint32_t sample_rate,
                                       uint32_t fft_size) {
  TonalReducer* self = (TonalReducer*)calloc(1, sizeof(TonalReducer));
  if (!self) {
    return NULL;
  }

  self->real_spectrum_size = real_spectrum_size;
  self->sample_rate = sample_rate;
  self->fft_size = fft_size;

  self->tonal_mask = (float*)calloc(real_spectrum_size, sizeof(float));
  if (!self->tonal_mask) {
    free(self);
    return NULL;
  }

  return self;
}

void tonal_reducer_free(TonalReducer* self) {
  if (!self) {
    return;
  }
  free(self->tonal_mask);
  free(self);
}

void tonal_reducer_run(TonalReducer* self, const float* noise_spectrum,
                       const float* max_profile, const float* median_profile,
                       float* alpha, float tonal_reduction_amount) {
  if (!self || !noise_spectrum || !max_profile || !median_profile || !alpha) {
    return;
  }

  // 1. Detect tonal components in the noise profile
  detect_tonal_components(noise_spectrum, max_profile, median_profile,
                          self->real_spectrum_size, self->sample_rate,
                          self->fft_size, self->tonal_mask);

  // 2. Skip alpha boosting if reduction is at 0dB (no tonal suppression)
  if (tonal_reduction_amount >= 1.0f) {
    return;
  }

  // 3. Boost alpha at tonal bins
  //    gain ≈ 1 - alpha  (when spectrum ≈ noise)
  //    alpha = 1 - tonal_reduction_amount
  //    E.g. -48dB → 0.004 → alpha = 0.996
  float alpha_needed = 1.0f - tonal_reduction_amount;

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    if (self->tonal_mask[k] <= 0.0f) {
      continue;
    }

    // Weight by mask strength (partial tonal bins get proportional boost)
    float target_alpha = alpha_needed * self->tonal_mask[k];

    // Only boost, never reduce (preserve existing suppression intent)
    alpha[k] = fmaxf(alpha[k], target_alpha);
  }
}

const float* tonal_reducer_get_mask(const TonalReducer* self) {
  if (!self) {
    return NULL;
  }
  return self->tonal_mask;
}
