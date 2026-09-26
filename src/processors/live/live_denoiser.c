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
  float snr_release;
  float peak_a_floor;
  /* Per-band synthesis depth exponent: compensates the bank geometry
   * so the closed-gate static response is flat with frequency
   * (sparse edge regions otherwise cut shallower than dense mids).
   * Derived at design time from the section overlap, no tuning. */
  float peak_exp[LIVE_NUM_BANDS];
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
    /* Wider notch for synthesis (see LIVE_PEAK_Q_MULT). */
    self->spk_alpha[k] = sin_w / (2.0F * q * LIVE_PEAK_Q_MULT);
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
  self->snr_attack =
      live_time_to_coeff(/* symmetric */ LIVE_SNR_SMOOTH_SEC, sample_rate);
  self->snr_release = self->snr_attack;
  self->peak_a_floor = powf(10.0F, -LIVE_PEAK_CUT_MIN_DB / 20.0F);
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
  const float snr_release = instance->snr_release;
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
      threshold_mult * threshold_mult * LIVE_GATE_RLO_MULT * span;
  const float wener_inv_anchor = 1.0F / wener_anchor;
  /* Learn: drop the floor so it re-converges on the current input. */
  if (atomic_exchange_explicit(&instance->reset_floor, false,
                               memory_order_acq_rel)) {
    memset(instance->noise, 0, sizeof(instance->noise));
  }

  sb_simd_state_t simd_state = sb_simd_enable_ftz_daz();

  const bool delta_monitoring =
      atomic_load_explicit(&instance->delta_monitoring, memory_order_acquire);

  for (uint32_t n = 0U; n < number_of_samples; n++) {
    const float x = input[n];

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

      /* Decision-directed smoothing of the instantaneous power SNR:
       * fast attack opens as soon as signal rises above the floor,
       * slow release collapses the noise's own fluctuation so the
       * gate target is steady (the time-domain analog of the full
       * denoiser's DD smoothing). */
      float lsnr = instance->r[k];
      const float inst_r = power / (instance->noise[k] + SPECTRAL_EPSILON);
      lsnr += (inst_r > lsnr ? snr_attack : snr_release) * (inst_r - lsnr);
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
      const float nu = lsnr * wener_inv_anchor;
      if (nu <= 1.0F) {
        target = g_min;
      } else {
        const float cut_db = reduction_db * powf(nu, -LIVE_GAIN_TAIL_EXP);
        target = powf(10.0F, -cut_db * 0.05F);
      }
      if (target < g_min) {
        target = g_min;
      } else if (target > 1.0F) {
        target = 1.0F;
      }

      instance->gate_target[k] = target;
    }

    /* Pass 1b: smooth the gain-law target ACROSS bands, then drive each
     * band's gain toward the smoothed value.
     *
     * The continuous law removed the hard 1.0/g_min cliff, but its slope
     * is still steepest where it matters most: right at the floor a 3 dB
     * SNR step between neighbours is a ~2.3 dB step in gain (12 -> 9.7 dB
     * of cut at 0 -> 3 dB past the anchor), and the noise profile decays
     * smoothly while the decision does not. Drawn and heard, that is a
     * visible band-to-band discontinuity. A binomial 5-tap [1 4 6 4 1]/16
     * cuts a step to 5/16 of its height, which puts the law's contribution
     * below the spectrum's own natural band-to-band steps (measured on
     * benchmark_suite_lead10s: gate-made step max 8.0 -> 2.78 dB against
     * an input-curve p95 of 2.11 dB).
     * Smoothing the TARGET rather than the gain is deliberate: gain is
     * state, so re-blending it every sample would diffuse all bands
     * toward one another within milliseconds. */
    for (uint32_t k = 0U; k < LIVE_NUM_BANDS; k++) {
      if (!instance->active[k]) {
        continue;
      }
      /* active[] is an active-prefix, so k-1/k-2 are always valid; the
       * upper taps clamp to k at the last active band. */
      const uint32_t km1 = (k > 0U) ? k - 1U : k;
      const uint32_t km2 = (km1 > 0U) ? km1 - 1U : km1;
      uint32_t kp1 = k + 1U;
      if (kp1 >= LIVE_NUM_BANDS || !instance->active[kp1]) {
        kp1 = k;
      }
      uint32_t kp2 = kp1 + 1U;
      if (kp2 >= LIVE_NUM_BANDS || !instance->active[kp2]) {
        kp2 = kp1;
      }
      const float* const t = instance->gate_target;
      const float desired =
          (t[km2] + 4.0F * t[km1] + 6.0F * t[k] + 4.0F * t[kp1] + t[kp2]) *
          (1.0F / 16.0F);

      /* De-clicked gain toward the target. */
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
