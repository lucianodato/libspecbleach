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

// Exclude smoother transition zones around voiced/silent boundaries so the
// metrics score steady-state behavior, not the (intentional) release tail
#define BOUNDARY_EXCLUSION_SEC 0.12F
#define GAP_EXCLUSION_SEC 0.03F

// Metric gates (deliberately loose: regression guards, not absolutes)
// Per-mode regression gates: [Temporal, NLM 2D, NLM 2D + DFTT]
#define MNI_GATE 1.25F
static const float MNI_GATES[3] = {1.25F, 1.25F, 1.05F};
static const float LSD_GATES[3] = {14.5F, 23.0F, 23.0F};
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
  float lsd;         // speech distortion in dB (lower = less underwater)
  float residual_db; // mean residual level in noise-only gaps, dB over the
                     // unprocessed noise level (lower = deeper rejection)
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
    float hum = 0.08F *
                (0.5F + (0.5F * sinf(2.0F * M_PIf * 0.33F * t + 1.0F))) *
                sinf(2.0F * M_PIf * 2300.0F * t);
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
      specbleach_denoiser_initialize(SAMPLE_RATE, FRAME_MS);
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
  Metrics result = {0.0F, 0.0F};
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
  free(frame_energy);

  // Mean residual level in gaps relative to the unprocessed noise level
  float res_sum = 0.0F;
  float noise_sum = 0.0F;
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
    }
    res_sum += out_e;
    noise_sum += mix_e;
    res_frames++;
  }
  TEST_ASSERT(res_frames > 10U, "enough gap frames");
  result.residual_db = 10.0F * log10f((res_sum / (float)res_frames) /
                                      (noise_sum / (float)res_frames + eps));
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

// Runs the full pipeline (learn -> process -> determinism probe) for one
// smoothing mode and returns its metrics
static Metrics run_and_measure(uint32_t smoothing_mode, const float* mix,
                               const float* clean, float* out, Analyzer* out_an,
                               Analyzer* clean_an, Analyzer* mix_an) {
  printf("== %s ==\n", mode_name(smoothing_mode));
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(SAMPLE_RATE, FRAME_MS);
  TEST_ASSERT(handle != NULL, "denoiser initialize");

  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .reduction_gain = 0.032F, // ~30 dB harsh
      .smoothing_factor = 0.8F,
      .smoothing_mode = smoothing_mode,
      .masking_depth = 0.5F,
      .whitening_factor = 0.5F,
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

  // Determinism probe: a second instance on the same input must produce
  // bit-identical output
  {
    float* out2 = (float*)calloc(TOTAL_SAMPLES, sizeof(float));
    specbleach_denoiser* handle2 =
        specbleach_denoiser_initialize(SAMPLE_RATE, FRAME_MS);
    TEST_ASSERT(handle2 != NULL, "second instance");
    SpecbleachDenoiserParameters p2 = parameters;
    p2.learn_noise = SPECBLEACH_LEARN_ALL; // parameters was already finalized
    TEST_ASSERT(specbleach_denoiser_load_parameters(handle2, &p2, sizeof(p2)),
                "load probe parameters");
    uint32_t pos2 = 0U;
    while (pos2 < LEARN_SAMPLES) {
      uint32_t chunk =
          (LEARN_SAMPLES - pos2 < 1764U) ? (LEARN_SAMPLES - pos2) : 1764U;
      specbleach_denoiser_process(handle2, chunk, mix + pos2, out2 + pos2);
      pos2 += chunk;
    }
    p2.learn_noise = SPECBLEACH_LEARN_OFF;
    specbleach_denoiser_load_parameters(handle2, &p2, sizeof(p2));
    while (pos2 < TOTAL_SAMPLES) {
      uint32_t chunk =
          (TOTAL_SAMPLES - pos2 < 1764U) ? (TOTAL_SAMPLES - pos2) : 1764U;
      specbleach_denoiser_process(handle2, chunk, mix + pos2, out2 + pos2);
      pos2 += chunk;
    }
    specbleach_denoiser_free(handle2);
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

  specbleach_denoiser_free(handle);

  const uint32_t stream_delay = measure_stream_delay(smoothing_mode);
  printf("  stream delay: %u samples\n", stream_delay);

  // Output is delayed by the empirically measured stream delay: analyze
  // clean/mix shifted back so frame f corresponds to the same input instant
  analyzer_run(out_an, out, 0);
  analyzer_run(clean_an, clean, -(int)stream_delay);
  analyzer_run(mix_an, mix, -(int)stream_delay);

  return compute_metrics(clean_an, out_an, mix_an);
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
    TEST_ASSERT(results[m].mni < MNI_GATES[m], "MNI regression gate");
    TEST_ASSERT(results[m].lsd < LSD_GATES[m], "LSD regression gate");
    TEST_ASSERT(results[m].residual_db < RESIDUAL_GATES_DB[m],
                "residual level regression gate");
    free(out);
  }

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
