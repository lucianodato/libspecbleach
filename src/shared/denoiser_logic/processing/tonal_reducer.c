/*
libspecbleach - A spectral processing library
*/

#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/configurations.h"
#include "shared/utils/tonal_detector.h"
#include <math.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

struct TonalReducer {
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t fft_size;
  float* tonal_mask;
  float* tonal_mask_buffers[2];
  uint32_t* deque_workspace;
  atomic_int active_mask_idx;
  float last_profile_sum;
  const float* last_detection_profile;
  bool has_cached_mask;
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

  self->tonal_mask_buffers[0] =
      (float*)calloc(real_spectrum_size, sizeof(float));
  self->tonal_mask_buffers[1] =
      (float*)calloc(real_spectrum_size, sizeof(float));
  if (!self->tonal_mask_buffers[0] || !self->tonal_mask_buffers[1]) {
    free(self->tonal_mask_buffers[0]);
    free(self->tonal_mask_buffers[1]);
    free(self->tonal_mask);
    free(self);
    return NULL;
  }

  self->deque_workspace =
      (uint32_t*)calloc(real_spectrum_size, sizeof(uint32_t));
  if (!self->deque_workspace) {
    free(self->tonal_mask_buffers[0]);
    free(self->tonal_mask_buffers[1]);
    free(self->tonal_mask);
    free(self);
    return NULL;
  }

  atomic_init(&self->active_mask_idx, 0);
  self->last_profile_sum = -1.0f;
  self->last_detection_profile = NULL;
  self->has_cached_mask = false;

  return self;
}

void tonal_reducer_free(TonalReducer* self) {
  if (!self) {
    return;
  }
  free(self->deque_workspace);
  free(self->tonal_mask_buffers[0]);
  free(self->tonal_mask_buffers[1]);
  free(self->tonal_mask);
  free(self);
}

void tonal_reducer_run(TonalReducer* self, const float* noise_spectrum,
                       const float* max_profile, const float* median_profile,
                       float* alpha, float tonal_reduction_gain) {
  if (!self || !noise_spectrum || !max_profile || !median_profile || !alpha) {
    return;
  }

  // 1. Determine which profile is active for detection
  bool profile_learned = false;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    if (median_profile[k] > 0.0f) {
      profile_learned = true;
      break;
    }
  }
  const float* detection_profile =
      profile_learned ? max_profile : noise_spectrum;

  // Compute profile checksum
  float current_sum = 0.0f;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    current_sum += detection_profile[k];
  }

  // 2. Run detection only if the profile has changed
  if (!self->has_cached_mask || current_sum != self->last_profile_sum ||
      detection_profile != self->last_detection_profile) {
    detect_tonal_components(noise_spectrum, max_profile, median_profile,
                            self->real_spectrum_size, self->sample_rate,
                            self->fft_size, self->tonal_mask,
                            self->deque_workspace);
    // Write completed mask to inactive double-buffer and release-publish
    int published_idx =
        atomic_load_explicit(&self->active_mask_idx, memory_order_relaxed);
    int write_idx = 1 - published_idx;
    memcpy(self->tonal_mask_buffers[write_idx], self->tonal_mask,
           self->real_spectrum_size * sizeof(float));
    atomic_store_explicit(&self->active_mask_idx, write_idx,
                          memory_order_release);
    self->last_profile_sum = current_sum;
    self->last_detection_profile = detection_profile;
    self->has_cached_mask = true;
  }

  // 2. Skip alpha boosting if reduction is 0 (no tonal suppression)
  // Gain of 1.0 means no reduction.
  if (tonal_reduction_gain >= 0.999f) {
    return;
  }

  // 3. Boost alpha at tonal bins
  //    Internally we work with strength (0..1) where 1 is max reduction.
  float tonal_reduction_strength = 1.0f - tonal_reduction_gain;
  float alpha_needed =
      ALPHA_MIN + (tonal_reduction_strength * (ALPHA_MAX_TONAL - ALPHA_MIN));

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    if (self->tonal_mask[k] <= 0.0f) {
      continue;
    }

    // Weight by mask strength (partial tonal bins get proportional boost
    // between ALPHA_MIN and alpha_needed)
    float target_alpha =
        ALPHA_MIN + (self->tonal_mask[k] * (alpha_needed - ALPHA_MIN));

    // Only boost, never reduce (preserve existing suppression intent)
    alpha[k] = fmaxf(alpha[k], target_alpha);
  }
}

const float* tonal_reducer_get_mask(const TonalReducer* self) {
  if (!self) {
    return NULL;
  }
  int idx = atomic_load_explicit(&self->active_mask_idx, memory_order_acquire);
  return self->tonal_mask_buffers[idx];
}

uint32_t tonal_reducer_get_peaks(const TonalReducer* self, float* peak_freqs_hz,
                                 uint32_t max_peaks) {
  if (!self || !peak_freqs_hz || max_peaks == 0) {
    return 0;
  }
  int idx = atomic_load_explicit(&self->active_mask_idx, memory_order_acquire);
  const float* mask = self->tonal_mask_buffers[idx];
  if (!mask) {
    return 0;
  }
  return tonal_detector_get_peaks(mask, self->real_spectrum_size,
                                  self->sample_rate, self->fft_size,
                                  peak_freqs_hz, max_peaks);
}
