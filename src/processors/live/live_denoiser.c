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

/* Direct public implementation (no wrapper): unlike the spectral path,
 * whose wrapper owns STFT framing and the noise profile, this engine is
 * already sample-in/sample-out, so the public API maps onto it directly. */

#include "specbleach_live_denoiser.h"

#include "shared/configurations.h"
#include "shared/utils/simd_utils.h"
#include <math.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

#ifndef SB_PI_F
#define SB_PI_F (3.14159265358979323846F)
#endif

/* Signal-path synthesis: the DRY input runs through a SERIAL cascade
 * of dynamic peaking sections; see the struct comment. */

/* Adaptation speed presets for the noise floor tracker. */
static float live_speed_factor(SpecbleachLiveEstimationMethod method) {
  switch (method) {
    case SPECBLEACH_LIVE_SPP_MMSE:
      return 0.5F;
    case SPECBLEACH_LIVE_BRANDT:
      return 1.0F;
    case SPECBLEACH_LIVE_MARTIN:
    default:
      return 2.0F;
  }
}

struct specbleach_live_denoiser { // NOLINT(readability-identifier-naming)
  uint32_t sample_rate;
  SpecbleachLiveDenoiserParameters parameters;
  float gate_release_sec;

  /* Per-band bandpass (RBJ constant-peak biquad), analysis only, on
   * the dry input. The audio itself never runs through these: the
   * output path is the serial peaking cascade below. Single 2nd-order
   * stages: steeper cascades measured comb notches (+13..+25 dB sum
   * peaks) when used as a parallel split, so selectivity comes from
   * band count + masking reserve instead; analysis mismatch with the
   * synthesis curve is now harmless because reconstruction is exact
   * whenever gains are open. */
  float b0[LIVE_NUM_BANDS];
  float b1[LIVE_NUM_BANDS];
  float b2[LIVE_NUM_BANDS];
  float a1[LIVE_NUM_BANDS];
  float a2[LIVE_NUM_BANDS];
  float x1[LIVE_NUM_BANDS];
  float x2[LIVE_NUM_BANDS];
  float y1[LIVE_NUM_BANDS];
  float y2[LIVE_NUM_BANDS];
  uint8_t active[LIVE_NUM_BANDS];

  /* Per-band control state. env is a fast amplitude follower on the
   * raw band; noise is the learned *power* (mean env^2) of the noise
   * profile; r is the decision-directed smoothed power SNR
   * (env^2/floor) used for the gain law so noise-only periods yield a
   * steady reduction instead of the instantaneous envelope bobbing
   * above and below the floor. */
  float env[LIVE_NUM_BANDS];
  float noise[LIVE_NUM_BANDS];
  float r[LIVE_NUM_BANDS];
  float gain[LIVE_NUM_BANDS];
  /* Voice-gated HF residual follower: slow[] tracks HF band power
   * only while voiceband energy sits clearly above the frozen floor;
   * the HF gate then uses max(noise, slow) so hiss stays down under
   * speech while noise-only stretches (no voice -> follower parked
   * at/below the floor) are untouched. Calloc-zeroed at init,
   * cleared on Learn reset like noise[]. */
  float slow[LIVE_NUM_BANDS];
  /* Raw gate target for the current sample, written by pass 1 and read
   * back by the across-band smoothing step below (scratch, not state:
   * every active entry is rewritten each sample before it is read). */
  float gate_target[LIVE_NUM_BANDS];

  /* Design-time band edges in Hz (inactive bands store 0/0). */
  float band_lo_hz[LIVE_NUM_BANDS];
  float band_hi_hz[LIVE_NUM_BANDS];

  /* Floor re-learn request (Learn button): honored by process. */
  atomic_bool reset_floor;

  /* Delta monitoring: when true, process outputs removed noise. */
  atomic_bool delta_monitoring;

  /* Signal-path synthesis: the DRY input runs through a SERIAL
   * cascade of dynamic peaking sections (Bertom-style, no
   * split-and-sum bank). Section k peaking-filters at band k's center
   * with gain g[k]; at g = 1 every section is an exact identity, so
   * the all-open cascade is perfectly transparent (no comb, no phase
   * shift). Per-band asymmetric gains therefore cannot create
   * reconstruction ripple the way parallel splits do (Hoffman & King
   * 1978; Linkwitz-Riley allpass reconstruction breaks as soon as
   * bands differ). The peaking coefficients below are the
   * design-time (gain-independent) parts: cos(w) and sin(w)/(2Q). */
  float spk_cos[LIVE_NUM_BANDS];
  float spk_alpha[LIVE_NUM_BANDS];
  /* Cascade serial state (direct form I, one section per band; two
   * passes per band for depth). */
  float sx1[LIVE_NUM_BANDS];
  float sx2[LIVE_NUM_BANDS];
  float sy1[LIVE_NUM_BANDS];
  float sy2[LIVE_NUM_BANDS];

  /* Per-sample one-pole coefficients derived from time constants. */
  float env_attack;
  float env_release;
  float noise_up;
  float noise_down;
  float gate_attack;
  float gate_release;
  float snr_attack;
  float snr_attack_hf;
  float snr_release;

  float peak_a_floor;
  float follow_down;
  float follow_up;
  float follow_release;
  /* Voice/HF band split (design-time centers, ascending): voice is
   * the [0, n_voice) prefix, HF follow the [i_hf, NUM) suffix. */
  uint32_t n_voice;
  uint32_t i_hf;
  /* Per-band synthesis depth exponent: compensates the bank geometry
   * so the closed-gate static response is flat with frequency
   * (sparse edge regions otherwise cut shallower than dense mids).
   * Derived at design time from the section overlap, no tuning. */
  float peak_exp[LIVE_NUM_BANDS];
  /* Bark group index for the joint gate decision (see Pass 1b):
   * the bank tiles Bark uniformly, so consecutive bands form equal
   * Bark groups. Calloc-zeroed; filled at design time. */
  uint8_t group[LIVE_NUM_BANDS];
};

static float live_bark(float freq_hz) {
  const float f = freq_hz / 1000.0F;
  return 13.0F * atanf(0.76F * f) +
         3.5F * atanf((freq_hz / 7500.0F) * (freq_hz / 7500.0F));
}

static float live_bark_to_hz(float bark, float nyquist_hz) {
  float lo = 20.0F;
  float hi = nyquist_hz;
  for (int i = 0; i < 40; i++) {
    const float mid = 0.5F * (lo + hi);
    if (live_bark(mid) < bark) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5F * (lo + hi);
}

static float live_clamp01(float v) {
  if (v < 0.0F) {
    return 0.0F;
  }
  if (v > 1.0F) {
    return 1.0F;
  }
  return v;
}

static float live_clamp(float v, float lo, float hi) {
  if (v < lo) {
    return lo;
  }
  if (v > hi) {
    return hi;
  }
  return v;
}

/* Threshold dB offset to linear multiplier (matches the full denoiser
 * profile offset semantics). */
static float live_threshold_mult(float threshold_db) {
  return powf(10.0F, threshold_db / 20.0F);
}

static float live_time_to_coeff(float seconds, float sample_rate) {
  if (seconds <= 0.0F || sample_rate <= 0.0F) {
    return 1.0F;
  }
  return 1.0F - expf(-1.0F / (seconds * sample_rate));
}

static void live_design_bank(specbleach_live_denoiser* self);
static void live_flatten_synthesis(specbleach_live_denoiser* self);

static void live_design_bank(specbleach_live_denoiser* self) {
  const float sample_rate = (float)self->sample_rate;
  const float nyquist = 0.5F * sample_rate;
  /* Top edge margin: 3 % (0.97) rather than 5 % so the top-octave void
   the 5 % margin created (bands starving as edges clamp to max_edge) is
   pushed out of the audible sweep range. Leaves 23.28-24 kHz at 48 kHz
   with no covering band: measured -15.7 dB closed-gate depth there
   against -55 dB below 19 kHz.
   NOTE: max_edge sets bark_span, so changing it re-lays-out EVERY
   band, unlike LIVE_BAND_OVERLAP which preserves band edges. */
  const float max_edge = nyquist * 0.97F;

  /* Span the Bark range from 20 Hz (not 0): at 128 bands the bottom
   * bin would otherwise sit entirely below 20 Hz, invert its edges in
   * the Hz clamp into 0/0, and form a dead prefix that breaks the
   * documented active-prefix/inactive-suffix layout the display relies
   * on. Every band then satisfies hi > lo by construction. */
  /* Span the Bark range from 20 Hz to the top edge margin: a fixed
   * Bark cap (24 ~= 15.5 kHz) left everything above it uncovered, so
   * derive the span from max_edge and cover the full range at any
   * sample rate. */
  const float bark_floor = live_bark(20.0F);
  const float bark_span = live_bark(max_edge) - bark_floor;
  for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
    const float bark_lo =
        bark_floor + bark_span * (float)k / (float)LIVE_NUM_BANDS;
    const float bark_hi =
        bark_floor + bark_span * (float)(k + 1U) / (float)LIVE_NUM_BANDS;
    float freq_lo = live_bark_to_hz(bark_lo, nyquist);
    float freq_hi = live_bark_to_hz(bark_hi, nyquist);
    if (freq_lo < 20.0F) {
      freq_lo = 20.0F;
    }
    if (freq_hi > max_edge) {
      freq_hi = max_edge;
    }

    const float bark_center = 0.5F * (bark_lo + bark_hi);
    float center = live_bark_to_hz(bark_center, nyquist);
    if (center < 20.0F) {
      center = 20.0F;
    }

    float bandwidth = (freq_hi - freq_lo) * LIVE_BAND_OVERLAP;
    if (bandwidth < 20.0F) {
      bandwidth = 20.0F;
    }

    if (center >= max_edge || freq_hi <= freq_lo) {
      self->active[k] = 0U;
      self->b0[k] = 0.0F;
      self->b1[k] = 0.0F;
      self->b2[k] = 0.0F;
      self->a1[k] = 0.0F;
      self->a2[k] = 0.0F;
      self->spk_cos[k] = 0.0F;
      self->spk_alpha[k] = 0.0F;
      self->sx1[k] = 0.0F;
      self->sx2[k] = 0.0F;
      self->sy1[k] = 0.0F;
      self->sy2[k] = 0.0F;
      self->band_lo_hz[k] = 0.0F;
      self->band_hi_hz[k] = 0.0F;
      continue;
    }

    /* Clamp Q so LF bands stay stable and do not ring excessively. */
    float q = center / bandwidth;
    if (q < LIVE_BAND_Q_MIN) {
      q = LIVE_BAND_Q_MIN;
    }
    if (q > LIVE_BAND_Q_MAX) {
      q = LIVE_BAND_Q_MAX;
    }

    const float w = 2.0F * SB_PI_F * center / sample_rate;
    const float sin_w = sinf(w);
    const float cos_w = cosf(w);
    const float alpha = sin_w / (2.0F * q);
    const float a0 = 1.0F + alpha;
    const float a0_inv = 1.0F / a0;

    self->active[k] = 1U;
    self->b0[k] = alpha * a0_inv;
    self->b1[k] = 0.0F;
    self->b2[k] = -alpha * a0_inv;
    self->a1[k] = -2.0F * cos_w * a0_inv;
    self->a2[k] = (1.0F - alpha) * a0_inv;
    /* Same slot geometry for the synthesis peaking section. */
    self->spk_cos[k] = cos_w;
    /* Wider notch for synthesis (see LIVE_PEAK_Q_MULT). NOTE: easing
     * spk_alpha wider toward the top octave for closed-gate HF depth
     * was tried at x2.0 and x1.3 and reverted: every widened build
     * CRASHes the plugin host, every revert passes. HF depth stays
     * geometry-limited; fix by adding bands, not widening. */
    self->spk_alpha[k] = sin_w / (2.0F * q * LIVE_PEAK_Q_MULT);
    /* Bark group for the joint decision: uniform Bark tiling makes
     * consecutive bands equal Bark groups (see Pass 1b). */
    self->group[k] = (uint8_t)(((uint32_t)k * (uint32_t)LIVE_NUM_GROUPS) /
                               (uint32_t)LIVE_NUM_BANDS);
    self->band_lo_hz[k] = freq_lo;
    self->band_hi_hz[k] = freq_hi;
  }
  live_flatten_synthesis(self);
}

/* dB magnitude of an RBJ peaking section (a0-normalized coefficients)
 * at angular frequency w. Used once at design time to measure how much
 * cascade depth each band's region receives from its neighbors. */
static float live_peak_mag_db(float b0, float b1, float b2, float a1, float a2,
                              float cos_w, float sin_w, float cos_2w,
                              float sin_2w) {
  const float num_r = b0 + b1 * cos_w + b2 * cos_2w;
  const float num_i = -(b1 * sin_w + b2 * sin_2w);
  const float den_r = 1.0F + a1 * cos_w + a2 * cos_2w;
  const float den_i = -(a1 * sin_w + a2 * sin_2w);
  const float ratio =
      (num_r * num_r + num_i * num_i) / (den_r * den_r + den_i * den_i);
  return 10.0F * log10f(ratio + 1e-12F);
}

/* Per-band synthesis depth exponents: the closed-gate cascade depth at
 * band k is the dB sum of every section's skirt over k, so sparse
 * regions (edges, HF) cut shallower than dense mids at the same gate
 * gain. Scale each section's exponent inversely to its region overlap
 * (diagonal compensation, normalized to the 1 kHz region) so the
 * static response is flat. Works only where band density binds: above
 * ~19 kHz overlap is already shallow, so this pushes peak_exp toward
 * its 2.0 ceiling and the sections still cannot bridge the inter-band
 * gap at overlap 1.2 (-15 dB closed-gate depth at 23.5 kHz vs -65 dB
 * in the mids). Setup-only. */
static void live_flatten_synthesis(specbleach_live_denoiser* self) {
  /* Reference-depth peaking coefficients per section (A = -6 dB):
   * overlap ratios are depth-independent to first order. */
  float rb0[LIVE_NUM_BANDS];
  float rb1[LIVE_NUM_BANDS];
  float rb2[LIVE_NUM_BANDS];
  float ra1[LIVE_NUM_BANDS];
  float ra2[LIVE_NUM_BANDS];
  float center[LIVE_NUM_BANDS];
  for (uint32_t j = 0U; j < LIVE_NUM_BANDS; j++) {
    if (!self->active[j]) {
      rb0[j] = 1.0F;
      rb1[j] = 0.0F;
      rb2[j] = 0.0F;
      ra1[j] = 0.0F;
      ra2[j] = 0.0F;
      center[j] = 0.0F;
      continue;
    }
    const float a_ref = 0.5F;
    const float alpha = self->spk_alpha[j];
    const float cos_w = self->spk_cos[j];
    const float inv_a0 = 1.0F / (1.0F + alpha / a_ref);
    rb0[j] = (1.0F + alpha * a_ref) * inv_a0;
    rb1[j] = (-2.0F * cos_w) * inv_a0;
    rb2[j] = (1.0F - alpha * a_ref) * inv_a0;
    ra1[j] = rb1[j];
    ra2[j] = (1.0F - alpha / a_ref) * inv_a0;
    center[j] = self->band_lo_hz[j] > 0.0F
                    ? sqrtf(self->band_lo_hz[j] * self->band_hi_hz[j])
                    : 0.0F;
  }
  float overlap[LIVE_NUM_BANDS];
  for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
    if (!self->active[k] || center[k] <= 0.0F) {
      overlap[k] = 0.0F;
      continue;
    }
    const float w = 2.0F * SB_PI_F * center[k] / (float)self->sample_rate;
    const float cos_w = cosf(w);
    const float sin_w = sinf(w);
    const float cos_2w = cosf(2.0F * w);
    const float sin_2w = sinf(2.0F * w);
    float sum = 0.0F;
    for (uint32_t j = 0U; j < LIVE_NUM_BANDS; j++) {
      if (!self->active[j]) {
        continue;
      }
      sum += live_peak_mag_db(rb0[j], rb1[j], rb2[j], ra1[j], ra2[j], cos_w,
                              sin_w, cos_2w, sin_2w);
    }
    overlap[k] = sum;
  }
  uint32_t ref = 0U;
  float ref_dist = 1e30F;
  for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
    if (!self->active[k] || center[k] <= 0.0F) {
      continue;
    }
    const float d = fabsf(center[k] - 1000.0F);
    if (d < ref_dist) {
      ref_dist = d;
      ref = k;
    }
  }
  const float ref_overlap = overlap[ref];
  for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
    if (!self->active[k] || overlap[k] > -1e-6F) {
      self->peak_exp[k] = LIVE_PEAK_GAIN_EXP;
      continue;
    }
    float e = LIVE_PEAK_GAIN_EXP * ref_overlap / overlap[k];
    if (e < 0.05F) {
      e = 0.05F;
    } else if (e > 2.0F) {
      e = 2.0F;
    }
    self->peak_exp[k] = e;
  }
}

static void live_derive_control_coeffs(specbleach_live_denoiser* self) {
  const float sample_rate = (float)self->sample_rate;
  const float speed =
      live_speed_factor(self->parameters.noise_estimation_method);

  self->env_attack = live_time_to_coeff(LIVE_ENV_ATTACK_SEC, sample_rate);
  self->env_release = live_time_to_coeff(LIVE_ENV_RELEASE_SEC, sample_rate);
  /* Symmetric tracker time constant: an unbiased mean of the band
   * power (see LIVE_NOISE_TRACK_SEC). */
  self->noise_up =
      live_time_to_coeff(LIVE_NOISE_TRACK_SEC * speed, sample_rate);
  self->noise_down = self->noise_up;
  self->gate_attack =
      live_time_to_coeff(self->parameters.attack_time, sample_rate);
  self->gate_release_sec = self->parameters.release_time;
  self->gate_release = live_time_to_coeff(self->gate_release_sec, sample_rate);
  self->snr_attack = live_time_to_coeff(LIVE_SNR_ATTACK_SEC, sample_rate);
  /* NOTE: knob-tracked release (fast shallow, slow deep) was tried
   * and reverted: +0.03 R20 body for hotter onsets. Flat wins on
   * the complaint axis. */
  self->snr_release = live_time_to_coeff(LIVE_SNR_RELEASE_SEC, sample_rate);
  /* HF attack tracks the knob like anchor/tail: RX's HF decisions
   * quicken with Reduction (50 ms at R12 keeps the tuned match
   * exact, 8 ms at R20 opens harmonic peaks for binary contrast).
   * Derived from the slider gain floor, clamped both ends. */
  {
    const float gm = self->parameters.reduction_gain;
    const float rdb = -20.0F * log10f(gm > 1.0e-6F ? gm : 1.0e-6F);
    float asec =
        LIVE_SNR_SMOOTH_SEC - LIVE_HF_ATK_TRACK * (rdb - LIVE_ANCHOR_REF_DB);
    if (asec < LIVE_SNR_ATTACK_HF_MIN_SEC) {
      asec = LIVE_SNR_ATTACK_HF_MIN_SEC;
    } else if (asec > LIVE_SNR_SMOOTH_SEC) {
      asec = LIVE_SNR_SMOOTH_SEC;
    }
    self->snr_attack_hf = live_time_to_coeff(asec, sample_rate);
  }
  self->peak_a_floor = powf(10.0F, -LIVE_PEAK_CUT_MIN_DB / 20.0F);
  self->follow_down = live_time_to_coeff(LIVE_FOLLOW_DOWN_SEC, sample_rate);
  self->follow_up = live_time_to_coeff(LIVE_FOLLOW_UP_SEC, sample_rate);
  self->follow_release =
      live_time_to_coeff(LIVE_FOLLOW_RELEASE_SEC, sample_rate);
  /* Voice/HF split from design-time geometric centers (ascending). */
  uint32_t nv = 0U, ih = LIVE_NUM_BANDS;
  uint32_t k;
  for (k = 0U; k < LIVE_NUM_BANDS; k++) {
    if (!self->active[k] || self->band_lo_hz[k] <= 0.0F ||
        self->band_hi_hz[k] <= 0.0F) {
      continue;
    }
    const float c = sqrtf(self->band_lo_hz[k] * self->band_hi_hz[k]);
    if (c < LIVE_FOLLOW_VOICE_HZ) {
      nv = k + 1U;
    }
    if (c >= LIVE_FOLLOW_HF_HZ && ih == LIVE_NUM_BANDS) {
      ih = k;
    }
  }
  self->n_voice = nv;
  self->i_hf = ih;
}

specbleach_live_denoiser* specbleach_live_denoiser_initialize(
    uint32_t sample_rate) {
  if (sample_rate == 0U) {
    return NULL;
  }
  specbleach_live_denoiser* self =
      (specbleach_live_denoiser*)calloc(1U, sizeof(*self));
  if (self == NULL) {
    return NULL;
  }
  self->sample_rate = sample_rate;
  self->parameters.reduction_gain = powf(10.0F, LIVE_MIN_GAIN_DB / 20.0F);
  self->parameters.attack_time = LIVE_GATE_ATTACK_DEFAULT_SEC;
  self->parameters.release_time = LIVE_GATE_RELEASE_DEFAULT_SEC;
  self->parameters.adaptive_noise = true;
  self->parameters.noise_estimation_method = SPECBLEACH_LIVE_MARTIN;
  self->parameters.threshold_db = LIVE_THRESHOLD_DB_DEFAULT;
  self->parameters.knee_db = LIVE_GATE_KNEE_DEFAULT_DB;

  live_design_bank(self);
  live_derive_control_coeffs(self);
  atomic_init(&self->reset_floor, false);
  atomic_init(&self->delta_monitoring, false);

  /* Gates start open so the first samples pass through transparently;
   * they close only once a band envelope drops below its floor. */
  for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
    self->gain[k] = 1.0F;
  }
  return self;
}

void specbleach_live_denoiser_free(specbleach_live_denoiser* instance) {
  free(instance);
}

bool specbleach_live_denoiser_load_parameters(
    specbleach_live_denoiser* instance,
    const SpecbleachLiveDenoiserParameters* parameters,
    uint32_t parameters_size) {
  if (instance == NULL || parameters == NULL) {
    return false;
  }
  if (parameters_size != sizeof(SpecbleachLiveDenoiserParameters)) {
    return false;
  }
  instance->parameters.reduction_gain =
      live_clamp01(parameters->reduction_gain);
  instance->parameters.attack_time =
      live_clamp(parameters->attack_time, LIVE_GATE_ATTACK_MIN_SEC,
                 LIVE_GATE_ATTACK_MAX_SEC);
  instance->parameters.release_time =
      live_clamp(parameters->release_time, LIVE_GATE_RELEASE_MIN_SEC,
                 LIVE_GATE_RELEASE_MAX_SEC);
  instance->parameters.adaptive_noise = parameters->adaptive_noise;
  instance->parameters.noise_estimation_method =
      parameters->noise_estimation_method;
  instance->parameters.threshold_db = live_clamp(
      parameters->threshold_db, LIVE_THRESHOLD_DB_MIN, LIVE_THRESHOLD_DB_MAX);
  instance->parameters.knee_db = live_clamp(
      parameters->knee_db, LIVE_GATE_KNEE_MIN_DB, LIVE_GATE_KNEE_MAX_DB);
  live_derive_control_coeffs(instance);
  return true;
}

bool specbleach_live_denoiser_process(specbleach_live_denoiser* instance,
                                      uint32_t number_of_samples,
                                      const float* input, float* output) {
  if (instance == NULL || number_of_samples == 0U) {
    return false;
  }
  if (input == NULL || output == NULL) {
    return false;
  }

  const float g_min = instance->parameters.reduction_gain;
  const float threshold_mult =
      live_threshold_mult(instance->parameters.threshold_db);
  const float env_attack = instance->env_attack;
  const float env_release = instance->env_release;
  const float noise_up = instance->noise_up;
  const float noise_down = instance->noise_down;
  const float gate_attack = instance->gate_attack;
  const float gate_release = instance->gate_release;
  const float snr_attack = instance->snr_attack;
  const float snr_attack_hf = instance->snr_attack_hf;
  const float snr_release = instance->snr_release;
  const float follow_down = instance->follow_down;
  const float follow_up = instance->follow_up;
  const float follow_release = instance->follow_release;
  const uint32_t n_voice = instance->n_voice;
  const uint32_t i_hf = instance->i_hf;
  const bool adaptive = instance->parameters.adaptive_noise;
  const float knee_ratio = powf(10.0F, -instance->parameters.knee_db / 20.0F);
  /* Gain-law constants (block-level: every term depends only on the
   * slider set, so none of it belongs in the per-band/per-sample loop).
   *   reduction_db  the slider cut recovered from g_min, so the tail
   *                 below is expressed in dB against the knob itself.
   *   wener_anchor  the SNR (band power / learned floor) at which the
   *                 cut starts to relax: Threshold shifts it (raise
   *                 Threshold = believe more of the signal is noise =
   *                 cut harder), Knee widens it in the same direction
   *                 the old r_hi did. At the defaults (Threshold 0 dB,
   *                 Knee 0 dB) the factor is exactly 1.0, so a band
   *                 sitting at the learned floor is cut by exactly the
   *                 slider marking. */
  const float reduction_db = -20.0F * log10f(g_min > 1.0e-6F ? g_min : 1.0e-6F);
  const float span = 1.0F + (1.0F - knee_ratio) * (LIVE_GATE_OPEN_RATIO - 1.0F);
  const float wener_anchor =
      threshold_mult * threshold_mult * LIVE_GATE_RLO_MULT * span *
      powf(g_min / LIVE_ANCHOR_REF_LIN, LIVE_ANCHOR_TRACK_EXP);
  const float wener_inv_anchor = 1.0F / wener_anchor;
  /* Voice-gated HF anchor lift (power ratio): under voice the HF
   * tail stays deep on weak harmonics (RX holds -14.2 vs our -10.3
   * at R20) while fricatives, pauses and noise-only stretches never
   * engage the voice gate (fires 65% vowels / 0% fricatives /
   * never without lowband energy) so sibilance and floors are
   * untouched. R-tracked from the 12 dB reference, clamped at 0
   * below it, so the tuned 12 dB match is untouched. */
  float hf_lift_db = LIVE_HF_TRACK_DB * (reduction_db - LIVE_ANCHOR_REF_DB);
  if (hf_lift_db < 0.0F) {
    hf_lift_db = 0.0F;
  }
  const float hf_voice_lift = powf(10.0F, hf_lift_db / 10.0F);
  /* HF tail tracks the knob like the anchor: RX's expansion ratio
   * rises with Reduction (R12 HF tones release like a=0.6, R20
   * vowel peaks need a=0.9 for binary contrast). 0.6 at/under R12
   * so the tuned 12 dB match is untouched. */
  float hf_tail = LIVE_GAIN_TAIL_EXP +
                  LIVE_HF_TAIL_TRACK * (reduction_db - LIVE_ANCHOR_REF_DB);
  if (hf_tail < LIVE_GAIN_TAIL_EXP) {
    hf_tail = LIVE_GAIN_TAIL_EXP;
  }
  /* Lowband anchor lift: RX demands more SNR to open low gates (its
   * low node sits hotter), so buried voice lows stay shut instead of
   * leaking. Full lift at/below LOW_FULL_HZ, tapering to none at
   * LOW_TOP_HZ (log-frequency). Block-level like the anchor itself;
   * the floor (nu <= 1) is untouched. */
  float inv_anchor[LIVE_NUM_BANDS];
  {
    const float log_lo = logf(LIVE_LOW_LIFT_FULL_HZ);
    const float log_hi = logf(LIVE_LOW_LIFT_TOP_HZ);
    const float lift_track =
        LIVE_LIFT_TRACK_DB * (reduction_db - LIVE_ANCHOR_REF_DB);
    uint32_t k = 0U;
    for (k = 0U; k < LIVE_NUM_BANDS; k++) {
      float lift_db = 0.0F;
      if (instance->active[k] && instance->band_lo_hz[k] > 0.0F &&
          instance->band_hi_hz[k] > 0.0F) {
        const float c =
            sqrtf(instance->band_lo_hz[k] * instance->band_hi_hz[k]);
        if (c <= LIVE_LOW_LIFT_FULL_HZ) {
          lift_db = LIVE_LOW_LIFT_DB + lift_track;
        } else if (c < LIVE_LOW_LIFT_TOP_HZ) {
          lift_db = LIVE_LOW_LIFT_DB * (logf(c) - log_hi) / (log_lo - log_hi);
        }
      }
      if (lift_db < 0.0F) {
        lift_db = 0.0F;
      }
      inv_anchor[k] = wener_inv_anchor / powf(10.0F, lift_db / 10.0F);
    }
  }
  /* Learn: drop the floor so it re-converges on the current input. */
  if (atomic_exchange_explicit(&instance->reset_floor, false,
                               memory_order_acq_rel)) {
    memset(instance->noise, 0, sizeof(instance->noise));
    memset(instance->slow, 0, sizeof(instance->slow));
  }

  sb_simd_state_t simd_state = sb_simd_enable_ftz_daz();

  const bool delta_monitoring =
      atomic_load_explicit(&instance->delta_monitoring, memory_order_acquire);

  for (uint32_t n = 0U; n < number_of_samples; n++) {
    const float x = input[n];
    float v_env = 0.0F;
    float v_noise = 0.0F;
    bool voice = false;
    bool voice_lo = false;
    float m_env = 0.0F;
    float m_noise = 0.0F;

    /* Pass 1: per-band analysis + gain computation. */
    for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
      if (!instance->active[k]) {
        continue;
      }

      /* RBJ bandpass, direct form I on the dry input sample. */
      const float y = instance->b0[k] * x + instance->b1[k] * instance->x1[k] +
                      instance->b2[k] * instance->x2[k] -
                      instance->a1[k] * instance->y1[k] -
                      instance->a2[k] * instance->y2[k];
      instance->x2[k] = instance->x1[k];
      instance->x1[k] = x;
      instance->y2[k] = instance->y1[k];
      instance->y1[k] = y;

      /* Envelope follower with fast attack / slow release. */
      const float mag = fabsf(y);
      float env = instance->env[k];
      env += (mag > env ? env_attack : env_release) * (mag - env);
      instance->env[k] = env;

      /* Adaptive noise floor in the power domain: slow average of the
       * band power (the learned noise profile). Starts at zero so the
       * gate fails open until the floor learns. */
      const float power = env * env;
      if (adaptive) {
        float floor = instance->noise[k];
        floor += (power > floor ? noise_up : noise_down) * (power - floor);
        instance->noise[k] = floor;
      }

      /* Voice gate for the HF hold/follower: lowband power vs
       * frozen lowband floor (level-invariant ratio), AND midband
       * harmonic richness (speech-selective: steady sine gates are
       * loud down low but empty in the mids, so they must not hold
       * HF shut - vowels keep 60/64 coverage, sine-gates drop to
       * 0%). Low completes at the first non-voice band, mid at the
       * first HF band, both before any k >= i_hf use below. */
      if (k < n_voice) {
        v_env += power;
        v_noise += instance->noise[k];
      } else if (k == n_voice) {
        voice_lo = v_env > LIVE_FOLLOW_RATIO * (v_noise + SPECTRAL_EPSILON);
        if (i_hf <= n_voice) {
          voice = voice_lo;
        } else {
          m_env += power;
          m_noise += instance->noise[k];
        }
      } else if (k < i_hf) {
        m_env += power;
        m_noise += instance->noise[k];
      } else if (k == i_hf) {
        voice = voice_lo &&
                m_env > LIVE_FOLLOW_MIDRATIO * (m_noise + SPECTRAL_EPSILON);
      }

      /* Decision-directed smoothing of the instantaneous power SNR:
       * fast attack opens as soon as signal rises above the floor,
       * slow release collapses the noise's own fluctuation so the
       * gate target is steady (the time-domain analog of the full
       * denoiser's DD smoothing). */
      /* HF residual follower: voice-gated minima tracker. While
       * voice is present slow[] falls fast toward local minima (the
       * hiss under speech) and climbs back slowly, so it estimates
       * the HF noise floor even while sibilance fires above it.
       * Silent stretches relax it to the frozen floor; the gate uses
       * whichever is hotter. */
      float floor_eff = instance->noise[k];
      if (k >= i_hf) {
        float sl = instance->slow[k];
        if (voice) {
          const float c = (power < sl) ? follow_down : follow_up;
          sl += c * (power - sl);
        } else {
          sl += follow_release * (instance->noise[k] - sl);
        }
        instance->slow[k] = sl;
        if (sl > floor_eff) {
          floor_eff = sl;
        }
      }
      float lsnr = instance->r[k];
      const float inst_r = power / (floor_eff + SPECTRAL_EPSILON);
      /* HF attack runs faster (R-tracked): harmonic peaks open fully
       * for a frame or two while valleys stay shut - RX is
       * near-binary under voice. Release stays slow everywhere so
       * noise chatter still integrates out (a fast HF release was
       * tried: better vowel frame-mean but worse body at both
       * depths - it overfits the diagnostic). */
      const float sa = k >= i_hf ? snr_attack_hf : snr_attack;
      lsnr += (inst_r > lsnr ? sa : snr_release) * (inst_r - lsnr);
      instance->r[k] = lsnr;

      float target;
      /* RX-shaped gain law: the cut is expressed in dB against the
       * slider and relaxes as a POWER of the SNR above the anchor
       *
       *     cut_dB = Reduction_dB * nu^(-a),   nu = lsnr / anchor
       *     G      = 10^(-cut_dB/20),  clamped to [g_min, 1]
       *
       * Two properties were measured off RX 10 Voice De-noise by
       * sweeping a pure tone through its learned profile (knob = 12 dB,
       * 1 kHz, level relative to the point where the full cut ends):
       *   d (dB)   3     8    13    18    23    28    33    38
       *   RX cut  -10.7  -7.6  -5.3  -3.6  -2.5  -1.7  -1.1  -0.8
       * which is 12 * 10^(-d/32.4) to within 0.5 dB, i.e. a constant
       * -0.30 dB of cut per dB of SNR. The old Wiener rational
       * G = (nu/(1+nu))^p pinned the same floor but reached it far too
       * fast - it gave back half of RX's cut by d = 8 and all of it by
       * d = 20, which is precisely why the mixed part measured -0.4 dB
       * where RX held -1.3 dB and the background stayed audible under
       * speech. nu <= 1 keeps the exact floor, so the noise-only probe
       * and the slider meaning are unchanged. */
      float nu = lsnr * inv_anchor[k];
      /* Voice-gated HF hold: hotter anchor while voice is present so
       * weak harmonics cannot open the tail (fricatives never engage
       * the gate, so sibilance still opens). k >= i_hf always runs
       * after the gate completes (k == n_voice), so voice is valid. */
      if (voice && k >= i_hf) {
        nu /= hf_voice_lift;
      }
      if (nu <= 1.0F) {
        target = g_min;
      } else {
        /* Power-subtraction gain (Berouti-style, STFT-subtractor
         * family): G = sqrt(max(1 - alpha/nu, g_min^2)), alpha
         * SNR-adaptive (LIVE_SS_*). At/below the anchor the floor
         * keeps probe and slider meaning; just above it the
         * over-subtraction holds near-floor where the old nu^-tail
         * law already released to half-cut - that is the deeper
         * mixed-section cut (vowel frame-mean exact vs RX). In HF
         * alpha scales with the R-tracked tail ratio, preserving
         * hiss/harmonic contrast. Sibilance overcut is worked on
         * separately, not by cooling this law. */
        const float seg_snr = 10.0F * log10f(nu);
        float alpha;
        if (seg_snr <= LIVE_SS_SNR_LO_DB) {
          alpha = LIVE_SS_ALPHA_HI;
        } else if (seg_snr >= LIVE_SS_SNR_HI_DB) {
          alpha = 1.0F;
        } else {
          alpha = LIVE_SS_ALPHA_HI -
                  (LIVE_SS_ALPHA_HI - 1.0F) * (seg_snr - LIVE_SS_SNR_LO_DB) /
                      (LIVE_SS_SNR_HI_DB - LIVE_SS_SNR_LO_DB);
        }
        if (k >= i_hf) {
          alpha *= hf_tail / LIVE_GAIN_TAIL_EXP;
        }
        const float g_floor2 = g_min * g_min;
        const float g2 = 1.0F - alpha / nu;
        if (g2 <= g_floor2) {
          target = g_min;
        } else if (g2 >= 1.0F) {
          target = 1.0F;
        } else {
          target = sqrtf(g2);
        }
      }
      if (target < g_min) {
        target = g_min;
      } else if (target > 1.0F) {
        target = 1.0F;
      }

      instance->gate_target[k] = target;
    }

    /* Pass 1b: joint gate decision per Bark GROUP, then drive each
     * band's gain toward its group's value.
     *
     * 64 Bark groups over the 256 bands (RX's own "64 Bark gates"
     * architecture): a group opens only when a member is FULLY open
     * (a real harmonic/sibilant peak), never on mid-level chatter -
     * a plain max opened everything (R12 body 0.26 -> 0.38) because
     * smoothed noise excursions sit half-open. Peak-triggered groups
     * preserve lone harmonics while uniformly quiet groups shut
     * together: near-binary under voice. Smoothing the TARGET rather
     * than the gain is deliberate: gain is state, so re-blending it
     * every sample would diffuse all bands toward one another within
     * milliseconds. */
    float gmax[LIVE_NUM_GROUPS];
    for (uint32_t g = 0U; g < (uint32_t)LIVE_NUM_GROUPS; g++) {
      gmax[g] = 0.0F;
    }
    for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
      if (!instance->active[k]) {
        continue;
      }
      const float t = instance->gate_target[k];
      uint32_t g = instance->group[k];
      if (g >= (uint32_t)LIVE_NUM_GROUPS) {
        g = (uint32_t)LIVE_NUM_GROUPS - 1U;
      }
      if (t > gmax[g]) {
        gmax[g] = t;
      }
    }
    for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
      if (!instance->active[k]) {
        continue;
      }
      uint32_t g = instance->group[k];
      if (g >= (uint32_t)LIVE_NUM_GROUPS) {
        g = (uint32_t)LIVE_NUM_GROUPS - 1U;
      }
      const float own = instance->gate_target[k];
      const float shared = gmax[g];
      float desired = shared > LIVE_GROUP_PEAK_TRIG ? shared : own;
      /* NOTE: a voice-confidence mid trim lived here and was
       * removed: it fired on all voiced mids (corpus 0.93 -> 2.11)
       * while the hot moments it targeted never moved - the
       * controller cannot separate them. */
      /* NOTE: contrast lift was tried here and reverted (level sweep:
       * RX spares tones only above ~+12 dB SNR, voice lives below). */
      /* NOTE: a broadband transient bypass (25 ms full-open on flux)
       * lived here and was removed: the impulse probe showed it
       * pumping ~100 ms of noise after every attack while RX stays
       * shut. Each band now decides alone; attacks survive through
       * their own loud bins. Onset chop is the known cost. */

      /* De-clicked gain toward the target. NOTE: a tonal-absence
       * fast release was tried for pause residual and reverted: it
       * chops reverb tails (R20 pauses -19.9 vs RX -12.9) instead
       * of denoising - pauses were already at parity. */
      float gain = instance->gain[k];
      gain += (desired > gain ? gate_attack : gate_release) * (desired - gain);
      if (gain < g_min) {
        gain = g_min;
      } else if (gain > 1.0F) {
        gain = 1.0F;
      }
      instance->gain[k] = gain;
    }

    /* Pass 2: serial dynamic-EQ cascade on the dry signal. Each open
     * section is an exact identity (skipped); each cutting section is
     * an RBJ peaking biquad modulated to that band's gain. */
    float sig = x;
    for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
      if (!instance->active[k]) {
        continue;
      }
      const float gain = instance->gain[k];
      if (gain > LIVE_PEAK_BYPASS_GAIN) {
        continue;
      }
      /* RBJ peaking (boost A = linear cut gain): coefficients are
       * modulated per sample since the gain is smoothed per sample.
       * The command is exponent-mapped per band (see peak_exp: flat
       * closed-gate response) and the cut is floored
       * (LIVE_PEAK_CUT_MIN_DB) for numerical stability of the float
       * biquad state. */
      float a = powf(gain, instance->peak_exp[k]);
      if (a < instance->peak_a_floor) {
        a = instance->peak_a_floor;
      }
      const float alpha = instance->spk_alpha[k];
      const float cos_w = instance->spk_cos[k];
      const float inv_a0 = 1.0F / (1.0F + alpha / a);
      const float b0n = (1.0F + alpha * a) * inv_a0;
      const float b1n = (-2.0F * cos_w) * inv_a0;
      const float b2n = (1.0F - alpha * a) * inv_a0;
      const float a1n = b1n;
      const float a2n = (1.0F - alpha / a) * inv_a0;
      const float out = b0n * sig + b1n * instance->sx1[k] +
                        b2n * instance->sx2[k] - a1n * instance->sy1[k] -
                        a2n * instance->sy2[k];
      instance->sx2[k] = instance->sx1[k];
      instance->sx1[k] = sig;
      instance->sy2[k] = instance->sy1[k];
      instance->sy1[k] = out;
      sig = out;
    }
    if (delta_monitoring) {
      /* Removed-noise view: dry minus processed. */
      output[n] = x - sig;
    } else {
      output[n] = sig;
    }
    /* NOTE: an absence room-tone cap (fresh pauses capped at -13 dB
     * like RX) lived here and was removed: scoping it to long
     * silences never fires on the pause metric (short gaps), and
     * any wider scope heats onsets and opens vowels. R20 deep
     * pauses stay a known tradeoff. */
  }

  sb_simd_restore_state(simd_state);
  return true;
}

uint32_t specbleach_live_denoiser_get_latency(
    specbleach_live_denoiser* instance) {
  (void)instance;
  return 0U;
}

bool specbleach_live_denoiser_get_band_levels(
    specbleach_live_denoiser* instance, float* input, float* output,
    float* threshold) {
  if (instance == NULL) {
    return false;
  }
  if (input == NULL && output == NULL && threshold == NULL) {
    return false;
  }
  const float threshold_mult =
      live_threshold_mult(instance->parameters.threshold_db);
  for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
    if (!instance->active[k]) {
      if (input != NULL) {
        input[k] = 0.0F;
      }
      if (output != NULL) {
        output[k] = 0.0F;
      }
      if (threshold != NULL) {
        threshold[k] = 0.0F;
      }
      continue;
    }
    if (input != NULL) {
      input[k] = instance->env[k];
    }
    if (output != NULL) {
      output[k] = instance->gain[k] * instance->env[k];
    }
    if (threshold != NULL) {
      threshold[k] = sqrtf(instance->noise[k]) * threshold_mult;
    }
  }
  return true;
}

SPECBLEACH_API float specbleach_live_denoiser_debug_mean_gain(
    specbleach_live_denoiser* instance);
SPECBLEACH_API float specbleach_live_denoiser_debug_mean_gain(
    specbleach_live_denoiser* instance) {
  double sum = 0.0;
  uint32_t count = 0U;
  for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
    if (instance->active[k]) {
      sum += instance->gain[k];
      count++;
    }
  }
  return (float)(sum / (double)count);
}

bool specbleach_live_denoiser_reset_noise_floor(
    specbleach_live_denoiser* instance) {
  if (instance == NULL) {
    return false;
  }
  atomic_store_explicit(&instance->reset_floor, true, memory_order_release);
  return true;
}

bool specbleach_live_denoiser_get_band_edges(specbleach_live_denoiser* instance,
                                             float* lower_hz, float* upper_hz) {
  if (instance == NULL) {
    return false;
  }
  if (lower_hz == NULL && upper_hz == NULL) {
    return false;
  }
  for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
    if (lower_hz != NULL) {
      lower_hz[k] = instance->band_lo_hz[k];
    }
    if (upper_hz != NULL) {
      upper_hz[k] = instance->band_hi_hz[k];
    }
  }
  return true;
}

bool specbleach_live_denoiser_set_delta_monitoring(
    specbleach_live_denoiser* instance, bool enabled) {
  if (instance == NULL) {
    return false;
  }
  atomic_store_explicit(&instance->delta_monitoring, enabled,
                        memory_order_release);
  return true;
}
