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

#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/configurations.h"
#include "shared/frame_rate_norm.h"
#include "shared/utils/simd_utils.h"
#include "shared/utils/spectral_utils.h"
#include "shared/utils/tonal_detector.h"
#include <float.h>
#include <math.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

static void tonal_reducer_compute_mask(TonalReducer* self,
                                       const float* noise_spectrum,
                                       const float* cv_mask_profile,
                                       bool cv_mask_available);

struct TonalReducer {
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t fft_size;
  float* tonal_mask;
  float* tonal_mask_buffers[2];
  uint32_t* deque_workspace;
  atomic_int active_mask_idx;

  // Dual-path decoupling state
  float* envelope;       // Opening(N): broadband baseline under the peaks
  float* alpha_tonal;    // Per-bin alpha for the tonal Wiener evaluation
  float* gain_memory[2]; // Tonal gain one-pole state per chain slot
  bool gain_seeded[2];
};

static void publish_mask(TonalReducer* self) {
  int published_idx =
      atomic_load_explicit(&self->active_mask_idx, memory_order_relaxed);
  int write_idx = 1 - published_idx;
  memcpy(self->tonal_mask_buffers[write_idx], self->tonal_mask,
         self->real_spectrum_size * sizeof(float));
  atomic_store_explicit(&self->active_mask_idx, write_idx,
                        memory_order_release);
}

void tonal_reducer_compute_mask(TonalReducer* self, const float* noise_spectrum,
                                const float* cv_mask_profile,
                                bool cv_mask_available) {
  if (!self || !noise_spectrum) {
    return;
  }

  if (cv_mask_available && cv_mask_profile) {
    // Manual Profile: copy CV mask directly
    memcpy(self->tonal_mask, cv_mask_profile,
           self->real_spectrum_size * sizeof(float));
  } else {
    // Idle (no profile, noise ~0) is the common steady state before learn:
    // detect_tonal_components does 1200x15 insertion sorts per frame and
    // dominates CPU while producing an all-zero mask. Early-out when noise is
    // negligible; adaptive case has non-zero noise and will not take this
    // path.
    float max_val = 0.0f;
    for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
      if (noise_spectrum[k] > max_val) {
        max_val = noise_spectrum[k];
      }
      if (max_val > TONAL_REDUCER_NEGLIGIBLE_NOISE_THRESHOLD) {
        break;
      }
    }
    if (max_val <= TONAL_REDUCER_NEGLIGIBLE_NOISE_THRESHOLD) {
      memset(self->tonal_mask, 0, self->real_spectrum_size * sizeof(float));
    } else {
      detect_tonal_components(noise_spectrum, noise_spectrum,
                              self->real_spectrum_size, self->sample_rate,
                              self->fft_size, self->tonal_mask,
                              self->deque_workspace);
    }
  }

  publish_mask(self);
}

void tonal_reducer_compute_split(TonalReducer* self,
                                 const float* noise_spectrum,
                                 const float* cv_mask_profile,
                                 bool cv_mask_available,
                                 float tonal_reduction_gain, float* noise_bb,
                                 float* noise_tonal) {
  if (!self || !noise_spectrum || !noise_bb || !noise_tonal) {
    return;
  }

  // The mask is always computed/published: the noise floor manager's
  // dual-path budget, the profile tonal offset and the peak reporter consume
  // it regardless of the requested reduction depth.
  tonal_reducer_compute_mask(self, noise_spectrum, cv_mask_profile,
                             cv_mask_available);

  // Inactive tonal path: identity split (bb = N, tonal = 0), so downstream
  // min() combining is a no-op and behavior matches the legacy single path.
  if (tonal_reduction_gain >= 0.999f) {
    memcpy(noise_bb, noise_spectrum, self->real_spectrum_size * sizeof(float));
    memset(noise_tonal, 0, self->real_spectrum_size * sizeof(float));
    return;
  }

  // Extract the broadband baseline under the tonal peaks and blend it in
  // gated by mask strength. Bins without tonal evidence keep the profile
  // exactly, so a missed detection stays on the legacy (safe) path.
  const float bin_hz = (self->fft_size > 0U)
                           ? ((float)self->sample_rate / (float)self->fft_size)
                           : 0.0F;
  const uint32_t window =
      sb_bins_for_hz(TONAL_DETONE_WINDOW_HZ, bin_hz, TONAL_DETONE_MIN_BINS,
                     TONAL_DETONE_MAX_BINS);
  if (!sb_spectral_envelope_opening(noise_spectrum, self->envelope, noise_bb,
                                    self->real_spectrum_size, window)) {
    // Envelope failure: keep the legacy profile everywhere
    memcpy(noise_bb, noise_spectrum, self->real_spectrum_size * sizeof(float));
    memset(noise_tonal, 0, self->real_spectrum_size * sizeof(float));
    return;
  }

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float w = fminf(self->tonal_mask[k], 1.0f);
    if (w <= 0.0f) {
      noise_bb[k] = noise_spectrum[k];
    } else {
      noise_bb[k] = (noise_spectrum[k] * (1.0f - w)) + (noise_bb[k] * w);
    }
    // Tonal residual handed to the parallel gain path (clamped >= 0)
    noise_tonal[k] = fmaxf(noise_spectrum[k] - noise_bb[k], 0.0f);
  }
}

void tonal_reducer_compute_tonal_gains(TonalReducer* self, uint32_t slot,
                                       const float* smoothed_magnitude,
                                       const float* noise_tonal,
                                       const float* tonal_mask,
                                       float tonal_reduction_gain,
                                       float* gain_tonal) {
  if (!self || !smoothed_magnitude || !noise_tonal || !gain_tonal ||
      slot >= 2U) {
    return;
  }

  // Per-bin sqrt/division on residual bins can decay into denormals; guard the
  // whole pass like the other per-frame DSP paths.
  sb_simd_state_t old_simd_state = sb_simd_enable_ftz_daz();

  // The mask must describe the same (possibly delayed) tile as noise_tonal.
  // Callers with a delayed residual pass the frame-aligned mask; a NULL mask
  // falls back to the current one for causal/low-latency callers.
  if (!tonal_mask) {
    tonal_mask = self->tonal_mask;
  }

  // Inactive tonal path (no residual anywhere): emit unity gains and keep the
  // one-pole state unused so reactivation starts from the fresh gains.
  bool any_tonal = false;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    if (noise_tonal[k] > FLT_EPSILON) {
      any_tonal = true;
      break;
    }
  }
  if (!any_tonal) {
    initialize_spectrum_with_value(gain_tonal, self->real_spectrum_size, 1.0f);
    self->gain_seeded[slot] = false;
    sb_simd_restore_state(old_simd_state);
    return;
  }

  // Alpha mapping identical to the legacy tonal alpha boost: mask-weighted
  // interpolation between ALPHA_MIN and the reduction-depth alpha.
  const float tonal_reduction_strength = 1.0f - tonal_reduction_gain;
  const float alpha_needed =
      ALPHA_MIN + (tonal_reduction_strength * (ALPHA_MAX_TONAL - ALPHA_MIN));
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float mask = fminf(tonal_mask[k], 1.0f);
    self->alpha_tonal[k] = ALPHA_MIN + (mask * (alpha_needed - ALPHA_MIN));
  }

  // Second Wiener evaluation on the tonal residual: plain subtraction curve
  // (no knee, no symmetry — the output is a half-spectrum whose upper bins
  // are ignored downstream). Notch is decisive at masked bins.
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float scaled_noise = noise_tonal[k] * self->alpha_tonal[k];
    const float mag = smoothed_magnitude[k];
    if (scaled_noise > FLT_MIN && mag > FLT_MIN) {
      const float diff = mag - scaled_noise;
      gain_tonal[k] = (diff > 0.0F) ? sqrtf(diff / mag) : 0.0F;
    } else {
      gain_tonal[k] = 1.0F;
    }
  }
  // Light one-pole stabilization: the notch is applied after the gain
  // smoothers, so mask flicker would otherwise enter the final gain raw.
  // Stationary tonal noise makes the small lag harmless. First-ever frame
  // per slot seeds the state with the raw gains (no ramp-in).
  const float alpha = expf(-1.0F / (float)TONAL_GAIN_STABILIZATION_HOPS);
  if (!self->gain_seeded[slot]) {
    memcpy(self->gain_memory[slot], gain_tonal,
           self->real_spectrum_size * sizeof(float));
    self->gain_seeded[slot] = true;
    sb_simd_restore_state(old_simd_state);
    return;
  }
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    gain_tonal[k] =
        (alpha * self->gain_memory[slot][k]) + ((1.0F - alpha) * gain_tonal[k]);
    self->gain_memory[slot][k] = gain_tonal[k];
  }

  sb_simd_restore_state(old_simd_state);
}

void tonal_reducer_promote_gain_slot(TonalReducer* self) {
  if (!self || !self->gain_memory[0] || !self->gain_memory[1]) {
    return;
  }
  memcpy(self->gain_memory[0], self->gain_memory[1],
         self->real_spectrum_size * sizeof(float));
  self->gain_seeded[0] = self->gain_seeded[1];
}

void tonal_reducer_swap_gain_slots(TonalReducer* self) {
  if (!self) {
    return;
  }
  float* gain_memory_tmp = self->gain_memory[0];
  self->gain_memory[0] = self->gain_memory[1];
  self->gain_memory[1] = gain_memory_tmp;
  const bool seeded_tmp = self->gain_seeded[0];
  self->gain_seeded[0] = self->gain_seeded[1];
  self->gain_seeded[1] = seeded_tmp;
}

void tonal_reducer_apply_alpha_boost(TonalReducer* self, float* alpha,
                                     float tonal_reduction_gain) {
  if (!self || !alpha || tonal_reduction_gain >= 0.999f) {
    return;
  }

  const float tonal_reduction_strength = 1.0f - tonal_reduction_gain;
  const float alpha_needed =
      ALPHA_MIN + (tonal_reduction_strength * (ALPHA_MAX_TONAL - ALPHA_MIN));
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    if (self->tonal_mask[k] <= 0.0f) {
      continue;
    }
    const float target_alpha = ALPHA_MIN + (fminf(self->tonal_mask[k], 1.0f) *
                                            (alpha_needed - ALPHA_MIN));
    alpha[k] = fmaxf(alpha[k], target_alpha);
  }
}

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
    tonal_reducer_free(self);
    return NULL;
  }

  self->deque_workspace =
      (uint32_t*)calloc(real_spectrum_size, sizeof(uint32_t));
  if (!self->deque_workspace) {
    tonal_reducer_free(self);
    return NULL;
  }

  self->envelope = (float*)calloc(real_spectrum_size, sizeof(float));
  self->alpha_tonal = (float*)calloc(real_spectrum_size, sizeof(float));
  self->gain_memory[0] = (float*)calloc(real_spectrum_size, sizeof(float));
  self->gain_memory[1] = (float*)calloc(real_spectrum_size, sizeof(float));
  if (!self->envelope || !self->alpha_tonal || !self->gain_memory[0] ||
      !self->gain_memory[1]) {
    tonal_reducer_free(self);
    return NULL;
  }

  atomic_init(&self->active_mask_idx, 0);

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
  free(self->envelope);
  free(self->alpha_tonal);
  free(self->gain_memory[0]);
  free(self->gain_memory[1]);
  free(self);
}

void tonal_reducer_reset(TonalReducer* self) {
  if (!self) {
    return;
  }
  if (self->tonal_mask) {
    memset(self->tonal_mask, 0, self->real_spectrum_size * sizeof(float));
  }
  if (self->tonal_mask_buffers[0]) {
    memset(self->tonal_mask_buffers[0], 0,
           self->real_spectrum_size * sizeof(float));
  }
  if (self->tonal_mask_buffers[1]) {
    memset(self->tonal_mask_buffers[1], 0,
           self->real_spectrum_size * sizeof(float));
  }
  if (self->gain_memory[0]) {
    memset(self->gain_memory[0], 0, self->real_spectrum_size * sizeof(float));
  }
  if (self->gain_memory[1]) {
    memset(self->gain_memory[1], 0, self->real_spectrum_size * sizeof(float));
  }
  self->gain_seeded[0] = false;
  self->gain_seeded[1] = false;
  atomic_store_explicit(&self->active_mask_idx, 0, memory_order_release);
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
