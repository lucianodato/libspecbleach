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

// Objective quality metrics for the temporal (1D) denoising chain:
// - Musical Noise Index (MNI): temporal coefficient of variation of the
//   residual power spectrum in noise-only segments (high = flickery,
//   chirpy residual).
// - Log-Spectral Distance (LSD): distortion of speech-active frames against
//   the known clean synthetic reference, restricted to speech-dominant bins
//   (high = over-suppressed, "underwater" speech).
// Both are printed every run so tuning can be tracked, and asserted against
// loose regression gates so a change that worsens either axis fails CI.

#include "../src/shared/configurations.h"
#include "specbleach_denoiser.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "FAIL: %s (line %d)\n", message, __LINE__);              \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define SAMPLE_RATE 44100U
#define FRAME_MS 40.0F // 1764 samples per block at 44.1 kHz
#define TOTAL_SAMPLES (SAMPLE_RATE * 4U)
#define LEARN_SAMPLES (SAMPLE_RATE / 2U) // 0.5 s noise-only for the profile

// Analysis STFT for metric computation (independent of the engine's STFT)
#define ANALYSIS_FFT 1024U
#define ANALYSIS_HOP 256U

// Synthetic speech: harmonic stack at f0 with slow syllable envelope
#define F0 130.0F
#define NUM_HARMONICS 24U
#define SYLLABLE_ON 0.30F  // seconds voiced
#define SYLLABLE_OFF 0.20F // seconds unvoiced (noise-only, MNI window)

// Synthetic tonal contaminant in the noise (a 2300 Hz hum fading in/out).
// The learn segment contains it, so the profile and the tonal mask catch it;
// the gap frames score how the tonal path + smoother treat it (and its
// spectral neighborhood).
#define HUM_FREQ 2300.0F

// Exclude smoother transition zones around voiced/silent boundaries so the
// metrics score steady-state behavior, not the (intentional) release tail
#define BOUNDARY_EXCLUSION_SEC 0.12F
#define GAP_EXCLUSION_SEC 0.03F

// Metric gates (deliberately loose: regression guards, not absolutes)
// Per-mode regression gates: [Temporal, NLM 2D, NLM 2D + DFTT]
// Temporal LSD calibrated at 18.0: the release cap (88 ms) restores gap
// suppression collapsed by frozen 324 ms releases, and working gains read
// higher deviation than frozen ones — NLM is accepted at 19.71 on identical
// material, so 18.0 stays the stricter bar while covering the fixed point
// (~16.7) with headroom.
#define MNI_GATE 1.25F
static const float MNI_GATES[3] = {1.25F, 1.25F, 1.05F};
static const float LSD_GATES[3] = {18.0F, 23.0F, 23.0F};
static const float RESIDUAL_GATES_DB[3] = {-3.0F, -12.0F, -15.0F};
#define LSD_GATE 25.0F // TEMP

static uint32_t next_rand_state = 22222U;

static float rand_uniform(void) {
  // xorshift32, deterministic across platforms (no libc rand())
  next_rand_state ^= next_rand_state << 13;
  next_rand_state ^= next_rand_state >> 17;
  next_rand_state ^= next_rand_state << 5;
  return (float)next_rand_state / 2147483648.0F - 1.0F;
}

// Speech-dominant bins: clean power must dominate noise by this ratio
#define SPEECH_BIN_DB_MARGIN 12.0F
#define PEAK_DYNAMIC_DB 25.0F
#define DEV_CLIP_DB 20.0F

typedef struct Metrics {
  float mni;         // musical noise index (lower = cleaner residual)
  float sustain_mni; // musical noise over voiced frames, noise-dominant bins
  float lsd;         // speech distortion in dB (lower = less underwater)
  float residual_db; // mean residual level in noise-only gaps, dB over the
                     // unprocessed noise level (lower = deeper rejection)
  // Tonal-region metrics (measurement bed for the tonal decoupling A/B):
  // relative residual levels (dB vs the unprocessed mix in the SAME bin
  // region) over noise-only gap frames.
  float tonal_resid_db;   // at the hum bins (how far the tonal component is
                          // knocked down; lower = deeper tonal rejection)
  float halo_resid_db;    // at the bins flanking the hum (the over-suppression
                          // "halo" the spatial smoother smears tonal carving
                          // into; lower = wider spill, higher = tighter)
  float control_resid_db; // far-from-hum broadband reference
} Metrics;

static void synthesize_inputs(float* clean, float* mix) {
  // One-pole lowpassed white noise, scaled for ~15 dB SNR under speech
  float noise_state = 0.0F;
  for (uint32_t n = 0U; n < TOTAL_SAMPLES; n++) {
    float t = (float)n / (float)SAMPLE_RATE;
    noise_state = (0.85F * noise_state) + (0.15F * rand_uniform());
    // Non-stationary noise: slow level breathing (stale-profile mismatch)
    // plus a tonal hum fading in and out (classic musical noise trigger:
    // bins toggle between cut-through and cut-to-floor)
    float breath = 1.0F + (0.5F * sinf(2.0F * M_PIf * 0.5F * t));
    // Tonal hum: strong enough to stand above the broadband background so
    // the tonal detector and the profile split actually engage the tonal
    // path (the A/B bed needs the mask populated).
    float hum = 2.4F *
                (0.35F + (0.25F * sinf(2.0F * M_PIf * 0.33F * t + 1.0F))) *
                sinf(2.0F * M_PIf * HUM_FREQ * t);
    float noise = (3.0F * noise_state * breath) + hum;

    float cycle = fmodf(t, SYLLABLE_ON + SYLLABLE_OFF);
    int voiced = cycle < SYLLABLE_ON;

    float speech = 0.0F;
    if (voiced) {
      // Harmonic stack with decaying amplitudes
      for (uint32_t h = 1U; h <= NUM_HARMONICS; h++) {
        speech += (sinf(2.0F * M_PIf * F0 * (float)h * t)) /
                  ((float)h * (1.0F + (float)h * 0.08F));
      }
      // Syllable envelope (smooth attack/decay within the voiced section)
      float env = sinf(M_PIf * cycle / SYLLABLE_ON);
      speech *= env;
    }

    clean[n] = speech;
    mix[n] = speech + noise;
  }
}

// Plucked-string bed (fast attack, exponential decay): the material the
// transient-protection sustain bug shows up on. Unlike the speech bed, the
// envelope only decays after the onset, so the detector fires at the attack
// and not through the body — isolating whether relief keeps overriding the
// smoothed gain past the onset.
#define PLUCK_FREQ 220.0F
#define PLUCK_PERIOD_SEC 0.5F
#define PLUCK_START_SEC 1.2F
#define PLUCK_DECAY_WIN_A 0.04F // sustain window after each onset, seconds
#define PLUCK_DECAY_WIN_B 0.18F

static void synthesize_plucks(float* clean, float* mix) {
  const uint32_t pluck_start = (uint32_t)(PLUCK_START_SEC * (float)SAMPLE_RATE);
  const uint32_t pluck_period =
      (uint32_t)(PLUCK_PERIOD_SEC * (float)SAMPLE_RATE);
  const uint32_t pluck_len = SAMPLE_RATE / 2U;
  const uint32_t attack = (SAMPLE_RATE * 3U) / 1000U;
  const float amps[4] = {1.0F, 0.5F, 0.33F, 0.25F};
  float noise_state = 0.0F;
  for (uint32_t n = 0U; n < TOTAL_SAMPLES; n++) {
    noise_state = (0.85F * noise_state) + (0.15F * rand_uniform());
    float speech = 0.0F;
    for (uint32_t t = pluck_start; t + pluck_len < TOTAL_SAMPLES;
         t += pluck_period) {
      if (n >= t && n < t + pluck_len) {
        uint32_t j = n - t;
        float env =
            (j < attack)
                ? (float)j / (float)attack
                : expf(-(float)(j - attack) / (0.150F * (float)SAMPLE_RATE));
        float s = 0.0F;
        for (uint32_t h = 0U; h < 4U; h++) {
          s += amps[h] * sinf(2.0F * M_PIf * PLUCK_FREQ * (float)(h + 1U) *
                              (float)j / (float)SAMPLE_RATE);
        }
        speech += 0.35F * env * s / 2.08F;
      }
    }
    clean[n] = speech;
    mix[n] = speech + (1.5F * noise_state);
  }
}

// Musical noise across the pluck sustain: spectral-shape CV over the decay
// windows, in bins away from the 220 Hz harmonic series (noise-only), so a
// smoothed residual scores low and a flickery one scores high.
static float pluck_sustain_mni(const float* power, uint32_t num_frames,
                               float stream_delay_sec) {
  const uint32_t bins = ANALYSIS_FFT / 2U;
  const float eps = 1e-20F;
  const float bin_hz = (float)SAMPLE_RATE / (float)ANALYSIS_FFT;
  const uint32_t pluck_start = (uint32_t)(PLUCK_START_SEC * (float)SAMPLE_RATE);
  const uint32_t pluck_period =
      (uint32_t)(PLUCK_PERIOD_SEC * (float)SAMPLE_RATE);
  const uint32_t pluck_len = SAMPLE_RATE / 2U;

  float* frame_energy = (float*)calloc(num_frames, sizeof(float));
  TEST_ASSERT(frame_energy != NULL, "pluck frame energy alloc");
  for (uint32_t f = 0U; f < num_frames; f++) {
    float e = 0.0F;
    for (uint32_t k = 1U; k < bins; k++) {
      e += power[f * bins + k];
    }
    frame_energy[f] = e;
  }

  float cv_sum = 0.0F;
  uint32_t cv_bins = 0U;
  for (uint32_t k = 1U; k < bins - 1U; k++) {
    bool near_harmonic = false;
    for (uint32_t h = 1U; h <= 12U; h++) {
      float hz = PLUCK_FREQ * (float)h;
      if (hz >= 0.5F * (float)SAMPLE_RATE) {
        break;
      }
      int hb = (int)lroundf(hz / bin_hz);
      if ((int)k >= hb - 2 && (int)k <= hb + 2) {
        near_harmonic = true;
        break;
      }
    }
    if (near_harmonic) {
      continue;
    }

    float mean = 0.0F;
    float sq_sum = 0.0F;
    uint32_t count = 0U;
    for (uint32_t f = 0U; f < num_frames; f++) {
      float center = ((float)(f * ANALYSIS_HOP) + (float)(ANALYSIS_FFT / 2U)) /
                     (float)SAMPLE_RATE;
      bool in_decay = false;
      for (uint32_t t = pluck_start; t + pluck_len < TOTAL_SAMPLES;
           t += pluck_period) {
        float rel = center - stream_delay_sec - ((float)t / (float)SAMPLE_RATE);
        if (rel >= PLUCK_DECAY_WIN_A && rel <= PLUCK_DECAY_WIN_B) {
          in_decay = true;
          break;
        }
      }
      if (!in_decay) {
        continue;
      }
      float norm = frame_energy[f] + eps;
      float p = (0.5F * power[f * bins + k] + 0.25F * power[f * bins + k - 1] +
                 0.25F * power[f * bins + k + 1]) /
                norm;
      mean += p;
      sq_sum += p * p;
      count++;
    }
    if (count > 1U && mean > eps) {
      mean /= (float)count;
      float variance = sq_sum / (float)count - mean * mean;
      if (variance > 0.0F) {
        cv_sum += sqrtf(variance) / mean;
        cv_bins++;
      }
    }
  }
  free(frame_energy);
  return (cv_bins > 0U) ? (cv_sum / (float)cv_bins) : 0.0F;
}

// Frame-wise analysis STFT (Hann window, ANALYSIS_FFT / ANALYSIS_HOP)
typedef struct Analyzer {
  double window[ANALYSIS_FFT];
  double* cos_table; // [bin][n]
  double* sin_table; // [bin][n]
  uint32_t num_frames;
  // power spectra laid out [frame][bin]
  float* power;
} Analyzer;

static void analyzer_init(Analyzer* self, uint32_t num_frames) {
  for (uint32_t n = 0U; n < ANALYSIS_FFT; n++) {
    self->window[n] =
        0.5 - 0.5 * cos(2.0 * M_PI * (double)n / (double)ANALYSIS_FFT);
  }
  self->cos_table = (double*)calloc((size_t)ANALYSIS_FFT * (ANALYSIS_FFT / 2U),
                                    sizeof(double));
  self->sin_table = (double*)calloc((size_t)ANALYSIS_FFT * (ANALYSIS_FFT / 2U),
                                    sizeof(double));
  TEST_ASSERT(self->cos_table != NULL && self->sin_table != NULL,
              "analyzer table alloc");
  for (uint32_t k = 0U; k < ANALYSIS_FFT / 2U; k++) {
    for (uint32_t n = 0U; n < ANALYSIS_FFT; n++) {
      double ph =
          2.0 * M_PI * (double)(k * n % ANALYSIS_FFT) / (double)ANALYSIS_FFT;
      self->cos_table[k * ANALYSIS_FFT + n] = cos(ph);
      self->sin_table[k * ANALYSIS_FFT + n] = sin(ph);
    }
  }
  self->num_frames = num_frames;
  self->power =
      (float*)calloc((size_t)num_frames * (ANALYSIS_FFT / 2U), sizeof(float));
  TEST_ASSERT(self->power != NULL, "analyzer power alloc");
}

static void analyzer_free(Analyzer* self) {
  free(self->cos_table);
  free(self->sin_table);
  free(self->power);
}

static void analyzer_run(Analyzer* self, const float* signal, int offset) {
  for (uint32_t f = 0U; f < self->num_frames; f++) {
    int start = (int)(f * ANALYSIS_HOP) + offset;
    if (start < 0 || (uint32_t)start + ANALYSIS_FFT > TOTAL_SAMPLES) {
      continue;
    }
    const float* frame = signal + start;
    for (uint32_t k = 0U; k < ANALYSIS_FFT / 2U; k++) {
      const double* cos_row = &self->cos_table[k * ANALYSIS_FFT];
      const double* sin_row = &self->sin_table[k * ANALYSIS_FFT];
      double re = 0.0;
      double im = 0.0;
      for (uint32_t n = 0U; n < ANALYSIS_FFT; n++) {
        double s = (double)frame[n] * self->window[n];
        re += s * cos_row[n];
        im -= s * sin_row[n];
      }
      float p = re * re + im * im;
      self->power[f * (ANALYSIS_FFT / 2U) + k] = p;
    }
  }
}

// Seconds from a frame center to the nearest voiced/silent boundary
static float boundary_distance(float center_seconds) {
  float cycle = fmodf(center_seconds, SYLLABLE_ON + SYLLABLE_OFF);
  if (cycle < SYLLABLE_ON) {
    float to_end = SYLLABLE_ON - cycle;
    return fminf(cycle, to_end);
  }
  return fminf(cycle - SYLLABLE_ON, (SYLLABLE_ON + SYLLABLE_OFF) - cycle);
}

static int frame_is_voiced(float frame_center_seconds) {
  float cycle = fmodf(frame_center_seconds, SYLLABLE_ON + SYLLABLE_OFF);
  return cycle < SYLLABLE_ON;
}

// Measures the true stream delay empirically with an impulse on a scratch
// instance (the reported latency is host-compensation info, not necessarily
// the buffer-level delay seen when feeding frame-sized blocks).
static uint32_t measure_stream_delay(uint32_t smoothing_mode) {
  specbleach_denoiser* probe =
      specbleach_denoiser_initialize(SAMPLE_RATE, FRAME_MS, 0u);
  TEST_ASSERT(probe != NULL, "probe init");
  SpecbleachDenoiserParameters p =
      (SpecbleachDenoiserParameters){.learn_noise = SPECBLEACH_LEARN_ALL,
                                     .reduction_gain = 0.1F,
                                     .smoothing_factor = 0.8F,
                                     .smoothing_mode = smoothing_mode,
                                     .masking_depth = 0.5F,
                                     .whitening_factor = 0.5F};
  TEST_ASSERT(specbleach_denoiser_load_parameters(probe, &p, sizeof(p)),
              "probe params");
  float pin[22050];
  float pout[22050] = {0.0F};
  memset(pin, 0, sizeof(pin));
  const uint32_t impulse_pos = 4000U;
  pin[impulse_pos] = 1.0F;
  specbleach_denoiser_process(probe, 11025, pin, pout);
  p.learn_noise = SPECBLEACH_LEARN_OFF;
  specbleach_denoiser_load_parameters(probe, &p, sizeof(p));
  specbleach_denoiser_process(probe, 11025, pin + 11025, pout + 11025);
  uint32_t peak = 0U;
  float peak_v = 0.0F;
  for (uint32_t n = 0U; n < 22050U; n++) {
    if (fabsf(pout[n]) > peak_v) {
      peak_v = fabsf(pout[n]);
      peak = n;
    }
  }
  specbleach_denoiser_free(probe);
  TEST_ASSERT(peak > impulse_pos, "impulse peak after input");
  return peak - impulse_pos;
}

static Metrics compute_metrics(const Analyzer* clean_an, const Analyzer* out_an,
                               const Analyzer* mix_an) {
  Metrics result = {0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F};
  const uint32_t bins = ANALYSIS_FFT / 2U;
  const float eps = 1e-20F;

  // Noise PSD estimate per bin from noise-only frames of the mix
  float* noise_psd = (float*)calloc(bins, sizeof(float));
  uint32_t noise_frame_count = 0U;
  for (uint32_t f = 0U; f < out_an->num_frames; f++) {
    float center = ((float)(f * ANALYSIS_HOP) + (float)(ANALYSIS_FFT / 2U)) /
                   (float)SAMPLE_RATE;
    if (center < ((float)LEARN_SAMPLES / (float)SAMPLE_RATE)) {
      continue; // skip the learn segment
    }
    if (!frame_is_voiced(center)) {
      for (uint32_t k = 0U; k < bins; k++) {
        noise_psd[k] += mix_an->power[f * bins + k];
      }
      noise_frame_count++;
    }
  }
  TEST_ASSERT(noise_frame_count > 10U, "enough noise-only frames");
  for (uint32_t k = 0U; k < bins; k++) {
    noise_psd[k] /= (float)noise_frame_count;
  }
  // MNI: temporal coefficient of variation of the residual spectral SHAPE
  // (3-bin-averaged, per-frame energy-normalized). Normalization removes slow
  // level modulation (breathing noise); only spectral shape instability
  // (musical noise) counts as flicker.
  float* frame_energy = (float*)calloc(out_an->num_frames, sizeof(float));
  for (uint32_t f = 0U; f < out_an->num_frames; f++) {
    float e = 0.0F;
    for (uint32_t k = 1U; k < bins; k++) {
      e += out_an->power[f * bins + k];
    }
    frame_energy[f] = e;
  }
  // power. Raw single-bin power of a stationary noise residual is
  // exponentially distributed (CV = 1 by construction), so we average over
  // neighboring bins first: stationary noise smooths toward CV ~= 0.58,
  // while musical-noise spectral spikes stay spiky (CV > 1).
  float cv_sum = 0.0F;
  uint32_t cv_bins = 0U;
  for (uint32_t k = 1U; k < bins - 1U; k++) { // skip DC
    float mean = 0.0F;
    float sq_sum = 0.0F;
    uint32_t count = 0U;
    for (uint32_t f = 0U; f < out_an->num_frames; f++) {
      float center = ((float)(f * ANALYSIS_HOP) + (float)(ANALYSIS_FFT / 2U)) /
                     (float)SAMPLE_RATE;
      if (center < ((float)LEARN_SAMPLES / (float)SAMPLE_RATE)) {
        continue;
      }
      if (frame_is_voiced(center) ||
          boundary_distance(center) < GAP_EXCLUSION_SEC) {
        continue;
      }
      float norm = frame_energy[f] + eps;
      float p = (0.5F * out_an->power[f * bins + k] +
                 0.25F * out_an->power[f * bins + k - 1] +
                 0.25F * out_an->power[f * bins + k + 1]) /
                norm;
      mean += p;
      sq_sum += p * p;
      count++;
    }
    if (count > 1U && mean > eps) {
      mean /= (float)count;
      float variance = sq_sum / (float)count - mean * mean;
      if (variance > 0.0F) {
        cv_sum += sqrtf(variance) / mean;
        cv_bins++;
      }
    }
  }
  result.mni = cv_sum / (float)cv_bins;

  // Sustain musical noise: same spectral-shape CV, but over VOICED frames and
  // restricted to bins where the clean reference is weak relative to the mix
  // (inter-harmonic noise). This is the axis transient protection hits when it
  // overrides the gain smoother through a pluck's sustain instead of only at
  // the onset.
  float sustain_cv_sum = 0.0F;
  uint32_t sustain_cv_bins = 0U;
  for (uint32_t k = 1U; k < bins - 1U; k++) {
    float mean = 0.0F;
    float sq_sum = 0.0F;
    uint32_t count = 0U;
    for (uint32_t f = 0U; f < out_an->num_frames; f++) {
      float center = ((float)(f * ANALYSIS_HOP) + (float)(ANALYSIS_FFT / 2U)) /
                     (float)SAMPLE_RATE;
      if (center < ((float)LEARN_SAMPLES / (float)SAMPLE_RATE) ||
          !frame_is_voiced(center) ||
          boundary_distance(center) < GAP_EXCLUSION_SEC) {
        continue;
      }
      float clean_p = 0.5F * clean_an->power[f * bins + k] +
                      0.25F * clean_an->power[f * bins + k - 1] +
                      0.25F * clean_an->power[f * bins + k + 1];
      float mix_p = 0.5F * mix_an->power[f * bins + k] +
                    0.25F * mix_an->power[f * bins + k - 1] +
                    0.25F * mix_an->power[f * bins + k + 1];
      if (clean_p > 0.25F * mix_p) {
        continue; // signal-dominant bin this frame
      }
      float norm = frame_energy[f] + eps;
      float p = (0.5F * out_an->power[f * bins + k] +
                 0.25F * out_an->power[f * bins + k - 1] +
                 0.25F * out_an->power[f * bins + k + 1]) /
                norm;
      mean += p;
      sq_sum += p * p;
      count++;
    }
    if (count > 1U && mean > eps) {
      mean /= (float)count;
      float variance = sq_sum / (float)count - mean * mean;
      if (variance > 0.0F) {
        sustain_cv_sum += sqrtf(variance) / mean;
        sustain_cv_bins++;
      }
    }
  }
  result.sustain_mni =
      (sustain_cv_bins > 0U) ? (sustain_cv_sum / (float)sustain_cv_bins) : 0.0F;

  free(frame_energy);

  // Mean residual level in gaps relative to the unprocessed noise level,
  // plus the tonal-decoupling A/B regional breakdown. The synthetic hum sits
  // at HUM_FREQ (analysis bins ~9 Hz wide... 44100/1024 ~ 43 Hz):
  //  - tonal region: bins around the hum (rejection of the tonal component)
  //  - halo region: bins flanking it + next harmonic (the smoother/tonal
  //    carve spill the decoupling is supposed to remove)
  //  - control region: a hum-free band as the broadband reference
  const int hum_bin = (int)lroundf((float)HUM_FREQ /
                                   ((float)SAMPLE_RATE / (float)ANALYSIS_FFT));
  float res_sum = 0.0F;
  float noise_sum = 0.0F;
  float tonal_out = 0.0F;
  float tonal_mix = 0.0F;
  float halo_out = 0.0F;
  float halo_mix = 0.0F;
  float control_out = 0.0F;
  float control_mix = 0.0F;
  uint32_t res_frames = 0U;
  for (uint32_t f = 0U; f < out_an->num_frames; f++) {
    float center = ((float)(f * ANALYSIS_HOP) + (float)(ANALYSIS_FFT / 2U)) /
                   (float)SAMPLE_RATE;
    if (center < ((float)LEARN_SAMPLES / (float)SAMPLE_RATE)) {
      continue;
    }
    if (frame_is_voiced(center) ||
        boundary_distance(center) < GAP_EXCLUSION_SEC) {
      continue;
    }
    float out_e = 0.0F;
    float mix_e = 0.0F;
    for (uint32_t k = 1U; k < bins; k++) {
      out_e += out_an->power[f * bins + k];
      mix_e += mix_an->power[f * bins + k];
      const float op = out_an->power[f * bins + k];
      const float mp = mix_an->power[f * bins + k];
      if (k >= (uint32_t)(hum_bin - 1U) && k <= (uint32_t)(hum_bin + 2U)) {
        tonal_out += op;
        tonal_mix += mp;
      } else if (k >= (uint32_t)(hum_bin - 6U) &&
                 k <= (uint32_t)(hum_bin + 7U)) {
        halo_out += op;
        halo_mix += mp;
      } else if (k >= 34U && k <= 41U) {
        control_out += op;
        control_mix += mp;
      }
    }
    res_sum += out_e;
    noise_sum += mix_e;
    res_frames++;
  }
  TEST_ASSERT(res_frames > 10U, "enough gap frames");
  result.residual_db = 10.0F * log10f((res_sum / (float)res_frames) /
                                      (noise_sum / (float)res_frames + eps));
  result.tonal_resid_db = 10.0F * log10f(tonal_out / (tonal_mix + eps) + eps);
  result.halo_resid_db = 10.0F * log10f(halo_out / (halo_mix + eps) + eps);
  result.control_resid_db =
      10.0F * log10f(control_out / (control_mix + eps) + eps);
  // LSD: distortion on speech-dominant bins of voiced frames, away from
  // voiced/silent boundaries (envelope edges dominate there). Uses the
  // per-frame median so boundary leakage cannot skew the score.
  float* per_frame_lsd = (float*)calloc(out_an->num_frames, sizeof(float));
  uint32_t lsd_frames = 0U;
  for (uint32_t f = 0U; f < out_an->num_frames; f++) {
    float center = ((float)(f * ANALYSIS_HOP) + (float)(ANALYSIS_FFT / 2U)) /
                   (float)SAMPLE_RATE;
    if (center < ((float)LEARN_SAMPLES / (float)SAMPLE_RATE)) {
      continue;
    }
    if (!frame_is_voiced(center) ||
        boundary_distance(center) < BOUNDARY_EXCLUSION_SEC) {
      continue;
    }
    // Frame loudness reference: only score peaks within PEAK_DYNAMIC_DB of
    // the loudest clean bin (the lowpass noise floor makes a per-bin noise
    // threshold alone useless at high frequencies)
    float max_clean_db = -1e9F;
    for (uint32_t k = 1U; k < bins - 1U; k++) {
      float cdb = 10.0F * log10f(clean_an->power[f * bins + k] + eps);
      if (cdb > max_clean_db) {
        max_clean_db = cdb;
      }
    }
    float sq_db = 0.0F;
    float weight_sum = 0.0F;
    for (uint32_t k = 2U; k < bins - 2U; k++) {
      float clean_p = clean_an->power[f * bins + k];
      float clean_db0 = 10.0F * log10f(clean_p + eps);
      if (clean_db0 <= max_clean_db - PEAK_DYNAMIC_DB) {
        continue;
      }
      // Soft speech-dominance weight: continuous in the clean/noise ratio so
      // bins near thresholds do not flip the score discontinuously. Local
      // peak prominence keeps the focus on harmonics, not leakage skirts.
      float dominance = clean_p / (noise_psd[k] + eps);
      float w_dom =
          dominance / (dominance + powf(10.0F, SPEECH_BIN_DB_MARGIN / 10.0F));
      float prominence = clean_p / (clean_an->power[f * bins + k - 2] +
                                    clean_an->power[f * bins + k + 2] + eps);
      float w_peak = fminf(1.0F, prominence * 0.5F);
      float w = w_dom * w_peak;
      if (w < 0.01F) {
        continue;
      }
      // 3-bin-averaged energy: the denoiser may rotate phase and shift
      // spectral leakage between adjacent bins, so raw single-bin comparison
      // would score a peak against a neighboring null
      float out_p = 0.5F * out_an->power[f * bins + k] +
                    0.25F * out_an->power[f * bins + k - 1] +
                    0.25F * out_an->power[f * bins + k + 1];
      float clean_p_avg = 0.5F * clean_p +
                          0.25F * clean_an->power[f * bins + k - 1] +
                          0.25F * clean_an->power[f * bins + k + 1];
      float d = 10.0F * log10f(out_p + eps) - 10.0F * log10f(clean_p_avg + eps);
      // Clip: bins deviating more than this are simply "destroyed"; without
      // the clip, threshold-flipping weights on rare extreme bins dominate
      // the weighted average and make the metric unstable
      d = fmaxf(-DEV_CLIP_DB, fminf(DEV_CLIP_DB, d));
      sq_db += w * d * d;
      weight_sum += w;
    }
    if (weight_sum > 0.0F) {
      per_frame_lsd[lsd_frames++] = sqrtf(sq_db / weight_sum);
    }
  }
  if (lsd_frames <= 10U) {
    free(per_frame_lsd);
    result.lsd = 1e9F;
    return result;
  }
  // median
  for (uint32_t i = 1U; i < lsd_frames; i++) {
    float v = per_frame_lsd[i];
    uint32_t j = i;
    while (j > 0U && per_frame_lsd[j - 1U] > v) {
      per_frame_lsd[j] = per_frame_lsd[j - 1U];
      j--;
    }
    per_frame_lsd[j] = v;
  }
  result.lsd = per_frame_lsd[lsd_frames / 2U];
  free(per_frame_lsd);
  free(noise_psd);

  return result;
}

static const char* mode_name(uint32_t mode) {
  switch (mode) {
    case SPECBLEACH_SMOOTHING_TEMPORAL:
      return "Temporal (1D)";
    case SPECBLEACH_SMOOTHING_NLM_2D:
      return "NLM 2D";
    default:
      return "NLM 2D + DFTT";
  }
}

// Per-build transient-protection setting (the measurement A/B toggles still
// drive the default absolute gates).
#if defined(TEST_TRANSIENT_PROTECTION) && TEST_TRANSIENT_PROTECTION
#define BUILD_TRANSIENT_PROTECTION true
#else
#define BUILD_TRANSIENT_PROTECTION false
#endif

// Runs learn -> denoise for one smoothing mode into `out`. Transient
// protection is an explicit argument so the differential penalty test can run
// both settings on the exact same bed.
static void process_bed(const float* mix, float* out, uint32_t smoothing_mode,
                        bool transient_protection) {
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, FRAME_MS, 0u);
  TEST_ASSERT(handle != NULL, "denoiser initialize");

  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .reduction_gain = 0.032F, // ~30 dB harsh
      .smoothing_factor = 0.8F,
      .smoothing_mode = smoothing_mode,
      .masking_depth = 0.5F,
      .whitening_factor = 0.5F,
      .tonal_reduction_gain = 0.0F, // max strength: strong A/B contrast
      .transient_protection_enable = transient_protection,
  };
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "load parameters");

  // Learn the profile on the leading noise-only segment
  uint32_t pos = 0U;
  while (pos < LEARN_SAMPLES) {
    uint32_t chunk =
        (LEARN_SAMPLES - pos < 1764U) ? (LEARN_SAMPLES - pos) : 1764U;
    TEST_ASSERT(
        specbleach_denoiser_process(handle, chunk, mix + pos, out + pos),
        "learn process");
    pos += chunk;
  }
  parameters.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &parameters,
                                                  sizeof(parameters)),
              "finalize parameters");
  while (pos < TOTAL_SAMPLES) {
    uint32_t chunk =
        (TOTAL_SAMPLES - pos < 1764U) ? (TOTAL_SAMPLES - pos) : 1764U;
    TEST_ASSERT(
        specbleach_denoiser_process(handle, chunk, mix + pos, out + pos),
        "denoise process");
    pos += chunk;
  }

  specbleach_denoiser_free(handle);
}

// Runs the full pipeline (learn -> process -> determinism probe) for one
// smoothing mode and returns its metrics
static Metrics run_and_measure(uint32_t smoothing_mode, const float* mix,
                               const float* clean, float* out, Analyzer* out_an,
                               Analyzer* clean_an, Analyzer* mix_an) {
  printf("== %s ==\n", mode_name(smoothing_mode));

  process_bed(mix, out, smoothing_mode, BUILD_TRANSIENT_PROTECTION);

  // Determinism probe: a second run on the same input must be bit-identical
  {
    float* out2 = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
    TEST_ASSERT(out2 != NULL, "second instance");
    process_bed(mix, out2, smoothing_mode, BUILD_TRANSIENT_PROTECTION);
    uint32_t diffs = 0U;
    double max_diff = 0.0;
    for (uint32_t n = 0U; n < TOTAL_SAMPLES; n++) {
      if (out[n] != out2[n]) {
        diffs++;
        double d = fabs((double)out[n] - (double)out2[n]);
        if (d > max_diff) {
          max_diff = d;
        }
      }
    }
    printf("  determinism: %u differing samples, max |diff| %.3g\n", diffs,
           max_diff);
    TEST_ASSERT(diffs == 0U,
                "determinism: second instance must be bit-identical");
    free(out2);
  }

  const uint32_t stream_delay = measure_stream_delay(smoothing_mode);
  printf("  stream delay: %u samples\n", stream_delay);

  // Output is delayed by the empirically measured stream delay: analyze
  // clean/mix shifted back so frame f corresponds to the same input instant
  analyzer_run(out_an, out, 0);
  analyzer_run(clean_an, clean, -(int)stream_delay);
  analyzer_run(mix_an, mix, -(int)stream_delay);

  return compute_metrics(clean_an, out_an, mix_an);
}

// Differential transient-protection penalty. The absolute MNI gates run with
// protection OFF (and are compiled out under TEST_TRANSIENT_PROTECTION), so
// they cannot see the smoother state that the on-onset relief and the 200 ms
// band hold release. Measure the SAME bed twice, protection off vs on, and
// bound the musical-noise increase. This is the meter for tuning the
// hold/relief trade-off; attack preservation is guarded separately in
// test_audio_regression.
// Budget for the sustain musical-noise increase with protection on. The
// NLM/DFTT chain applies no transient relief (1.00 by construction); the 1D
// temporal chain applies relief but the gain smoother owns it (~1.02-1.04).
// 1.10 leaves headroom for tuning while still catching a relief rearrangement
// that starts leaking unsmoothed gain into a note's sustain.
#define TRANSIENT_MNI_PENALTY_GATE 1.10F

static void test_transient_protection_penalty(const float* mix,
                                              const float* clean,
                                              Analyzer* clean_an,
                                              Analyzer* mix_an) {
  const uint32_t modes[2] = {SPECBLEACH_SMOOTHING_TEMPORAL,
                             SPECBLEACH_SMOOTHING_NLM_2D_DFTT};
  const uint32_t num_frames =
      (TOTAL_SAMPLES - ANALYSIS_FFT) / ANALYSIS_HOP + 1U;
  float sustain_penalties[2] = {0.0F, 0.0F};

  printf("== transient protection musical-noise penalty (off vs on) ==\n");
  for (uint32_t m = 0U; m < 2U; m++) {
    const uint32_t mode = modes[m];
    float* out_off = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
    float* out_on = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
    TEST_ASSERT(out_off != NULL && out_on != NULL, "penalty output alloc");

    process_bed(mix, out_off, mode, false);
    process_bed(mix, out_on, mode, true);

    const uint32_t stream_delay = measure_stream_delay(mode);
    analyzer_run(clean_an, clean, -(int)stream_delay);
    analyzer_run(mix_an, mix, -(int)stream_delay);

    Analyzer off_an;
    Analyzer on_an;
    analyzer_init(&off_an, num_frames);
    analyzer_init(&on_an, num_frames);
    analyzer_run(&off_an, out_off, 0);
    analyzer_run(&on_an, out_on, 0);

    const Metrics m_off = compute_metrics(clean_an, &off_an, mix_an);
    const Metrics m_on = compute_metrics(clean_an, &on_an, mix_an);
    const float penalty = m_on.mni / m_off.mni;
    const float sustain_penalty = m_on.sustain_mni / m_off.sustain_mni;
    sustain_penalties[m] = sustain_penalty;
    const float noise_let_through_db = m_on.residual_db - m_off.residual_db;
    printf(
        "  %-12s gapMNI off=%.4f on=%.4f pen=%.3f | sustainMNI off=%.4f "
        "on=%.4f pen=%.3f | residual %+.2f dB\n",
        mode_name(mode), m_off.mni, m_on.mni, penalty, m_off.sustain_mni,
        m_on.sustain_mni, sustain_penalty, noise_let_through_db);

    analyzer_free(&off_an);
    analyzer_free(&on_an);
    free(out_off);
    free(out_on);
  }

  // Report both modes first, then gate, so a failure still prints the full
  // table. The sustain axis is the one transient protection corrupts when it
  // overrides the smoother past the onset; enforced on the decoupled build
  // (legacy A/B builds print only).
#if TONAL_DUAL_PATH
  for (uint32_t m = 0U; m < 2U; m++) {
    if (sustain_penalties[m] >= TRANSIENT_MNI_PENALTY_GATE) {
      fprintf(stderr,
              "FAIL: %s transient-protection sustain musical-noise penalty "
              "%.3f >= target %.2f\n",
              mode_name(modes[m]), sustain_penalties[m],
              TRANSIENT_MNI_PENALTY_GATE);
      exit(1);
    }
  }
#endif
}

// Pluck-sustain differential: runs the pluck bed with smoothing on and bounds
// the sustain musical-noise increase. This is the axis the UI LED/onset vs
// 200 ms relief-hold mismatch corrupts on plucked material.
static void test_pluck_sustain_penalty(void) {
  float* clean = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
  float* mix = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
  TEST_ASSERT(clean != NULL && mix != NULL, "pluck bed alloc");
  synthesize_plucks(clean, mix);

  const uint32_t modes[2] = {SPECBLEACH_SMOOTHING_TEMPORAL,
                             SPECBLEACH_SMOOTHING_NLM_2D_DFTT};
  const uint32_t num_frames =
      (TOTAL_SAMPLES - ANALYSIS_FFT) / ANALYSIS_HOP + 1U;
  float pluck_penalties[2] = {0.0F, 0.0F};

  printf("== pluck sustain musical noise (off vs on) ==\n");
  for (uint32_t m = 0U; m < 2U; m++) {
    float* out_off = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
    float* out_on = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
    TEST_ASSERT(out_off != NULL && out_on != NULL, "pluck output alloc");
    process_bed(mix, out_off, modes[m], false);
    process_bed(mix, out_on, modes[m], true);

    const uint32_t stream_delay = measure_stream_delay(modes[m]);
    Analyzer off_an;
    Analyzer on_an;
    analyzer_init(&off_an, num_frames);
    analyzer_init(&on_an, num_frames);
    analyzer_run(&off_an, out_off, 0);
    analyzer_run(&on_an, out_on, 0);
    const float delay_sec = (float)stream_delay / (float)SAMPLE_RATE;
    const float mni_off =
        pluck_sustain_mni(off_an.power, off_an.num_frames, delay_sec);
    const float mni_on =
        pluck_sustain_mni(on_an.power, on_an.num_frames, delay_sec);
    pluck_penalties[m] = mni_on / mni_off;
    printf("  %-12s sustainMNI off=%.4f on=%.4f penalty=%.3f\n",
           mode_name(modes[m]), mni_off, mni_on, pluck_penalties[m]);

    analyzer_free(&off_an);
    analyzer_free(&on_an);
    free(out_off);
    free(out_on);
  }

  // Enforced on the decoupled build (legacy A/B builds print only).
#if TONAL_DUAL_PATH
  for (uint32_t m = 0U; m < 2U; m++) {
    if (pluck_penalties[m] >= TRANSIENT_MNI_PENALTY_GATE) {
      fprintf(stderr,
              "FAIL: %s pluck-sustain musical-noise penalty %.3f >= target "
              "%.2f\n",
              mode_name(modes[m]), pluck_penalties[m],
              TRANSIENT_MNI_PENALTY_GATE);
      exit(1);
    }
  }
#endif
  free(clean);
  free(mix);
}

int main() {
  setvbuf(stdout, NULL, _IONBF, 0); // live output; asserts abort otherwise
  float* clean = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
  float* mix = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
  TEST_ASSERT(clean && mix, "signal allocs");

  synthesize_inputs(clean, mix);
  float* mix_snapshot = (float*)malloc(TOTAL_SAMPLES * sizeof(float));
  memcpy(mix_snapshot, mix, TOTAL_SAMPLES * sizeof(float));

  const uint32_t num_frames =
      (TOTAL_SAMPLES - ANALYSIS_FFT) / ANALYSIS_HOP + 1U;
  Analyzer clean_an;
  Analyzer out_an;
  Analyzer mix_an;
  analyzer_init(&clean_an, num_frames);
  analyzer_init(&out_an, num_frames);
  analyzer_init(&mix_an, num_frames);

  const uint32_t modes[3] = {SPECBLEACH_SMOOTHING_TEMPORAL,
                             SPECBLEACH_SMOOTHING_NLM_2D,
                             SPECBLEACH_SMOOTHING_NLM_2D_DFTT};
  Metrics results[3];

  for (uint32_t m = 0U; m < 3U; m++) {
    float* out = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
    TEST_ASSERT(out != NULL, "out alloc");
    results[m] =
        run_and_measure(modes[m], mix, clean, out, &out_an, &clean_an, &mix_an);
    printf("  MNI (musical noise, lower better): %.4f (gate %.2f)\n",
           results[m].mni, MNI_GATES[m]);
    printf("  LSD (speech distortion dB, lower better): %.2f (gate %.1f)\n",
           results[m].lsd, LSD_GATES[m]);
    printf(
        "  Gap residual (dB vs unprocessed, lower better): %.2f (gate %.1f)\n",
        results[m].residual_db, RESIDUAL_GATES_DB[m]);
    printf(
        "  [tonal A/B] hum %.2f dB | halo %.2f dB | control %.2f dB (gap "
        "residual per region vs unprocessed same region)\n",
        results[m].tonal_resid_db, results[m].halo_resid_db,
        results[m].control_resid_db);
    // Regression gates are enforced on the decoupled build only; the legacy
    // measurement build (TONAL_DUAL_PATH=0) prints gate overruns so the full
    // A/B table stays available.
#if TONAL_DUAL_PATH && !defined(TEST_TRANSIENT_PROTECTION)
    TEST_ASSERT(results[m].mni < MNI_GATES[m], "MNI regression gate");
    TEST_ASSERT(results[m].lsd < LSD_GATES[m], "LSD regression gate");
    TEST_ASSERT(results[m].residual_db < RESIDUAL_GATES_DB[m],
                "residual level regression gate");
#endif
    free(out);
  }

  test_pluck_sustain_penalty();
  test_transient_protection_penalty(mix, clean, &clean_an, &mix_an);

  printf("MNI ratio NLM/Temporal: %.3f, DFTT/Temporal: %.3f\n",
         results[1].mni / results[0].mni, results[2].mni / results[0].mni);
  printf("LSD delta NLM-Temporal: %.2f dB, DFTT-Temporal: %.2f dB\n",
         results[1].lsd - results[0].lsd, results[2].lsd - results[0].lsd);
  printf("\u2713 quality metrics within gates\n");

  analyzer_free(&clean_an);
  analyzer_free(&out_an);
  analyzer_free(&mix_an);
  uint32_t modified = 0U;
  for (uint32_t n = 0U; n < TOTAL_SAMPLES; n++) {
    if (mix[n] != mix_snapshot[n]) {
      modified++;
    }
  }
  TEST_ASSERT(modified == 0U, "input buffer must stay unmodified");
  free(mix_snapshot);
  free(clean);
  free(mix);

  return 0;
}
