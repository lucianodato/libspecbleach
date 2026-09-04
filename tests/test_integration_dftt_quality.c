/*
libspecbleach - A spectral processing library

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

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

/*
 * DFTT quality instrument (integration): quantifies the musical-noise
 * refinement of NLM+DFTT over plain NLM on a synthetic comb+noise stimulus.
 *
 * Metrics (computed on an analysis STFT of the plugin output):
 *  - Speckle flicker: mean per-bin temporal variance (dB^2) of the magnitude
 *    in noise-only bands. Musical noise = frame-to-frame flicker. DFTT must
 *    not increase it and (once tuned) must reduce it.
 *  - Comb preservation: mean magnitude (dB) at harmonic bins. DFTT must not
 *    muffle the tonal structure (within COMB_TOLERANCE_DB of plain NLM).
 *
 * Both engines see the identical deterministic stimulus; run in NLM_2D and
 * NLM_2D_DFTT modes and compare. Tuning knob: MIN_SPECKLE_IMPROVEMENT_PCT.
 */

#include "specbleach_denoiser.h"

#include "../src/shared/configurations.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define SAMPLE_RATE 44100U
#define ENGINE_FRAME_MS (46.0F)
#define LEARN_SECONDS (0.5F)
#define PROCESS_SECONDS (1.5F)
#define WARMUP_SECONDS (0.3F) // Skipped after the engine's latency alignment

// Analysis STFT (independent of the engine's internal one)
#define ANALYSIS_FFT 2048U
#define ANALYSIS_HOP 512U

// Comb at exact analysis bins: k = 10..60 (f0 = 10 * SR / ANALYSIS_FFT)
#define COMB_BIN_START 10U
#define COMB_BIN_STEP 10U
#define COMB_PARTIALS 6U
#define COMB_PARTIAL_AMP (0.08F)
#define NOISE_AMP (0.05F)

// Noise-only analysis bands: bins 150..350 (above the comb), every 7th bin
#define NOISE_BIN_START 150U
#define NOISE_BIN_END 350U
#define NOISE_BIN_STEP 7U

// Minimum relative speckle-flicker reduction DFTT must deliver vs plain NLM
// (measured 19-21% at K=32 on the reference stimulus; guard sits well below).
#define MIN_SPECKLE_IMPROVEMENT_PCT (10.0F)
#define COMB_TOLERANCE_DB (1.5F)

// Deterministic xorshift RNG so both runs get the identical noise
static uint32_t rng_state = 0x12345678U;
static float rng_uniform(void) {
  rng_state ^= rng_state << 13;
  rng_state ^= rng_state >> 17;
  rng_state ^= rng_state << 5;
  return (float)(rng_state >> 8) * (1.0F / 16777216.0F) - 0.5F;
}

typedef struct {
  float* data;
  size_t len;
} Signal;

static Signal make_signal(bool with_comb) {
  const size_t learn_len = (size_t)(LEARN_SECONDS * (float)SAMPLE_RATE);
  const size_t proc_len = (size_t)(PROCESS_SECONDS * (float)SAMPLE_RATE);
  Signal s = {.data = (float*)calloc(learn_len + proc_len, sizeof(float)),
              .len = learn_len + proc_len};
  TEST_ASSERT(s.data != NULL, "signal allocation failed");
  for (size_t i = 0; i < s.len; i++) {
    const float noise = NOISE_AMP * 2.0F * rng_uniform();
    float comb = 0.0F;
    if (with_comb && i >= learn_len) {
      const double t = (double)(i - learn_len) / (double)SAMPLE_RATE;
      for (uint32_t p = 0U; p < COMB_PARTIALS; p++) {
        const double f = (double)(COMB_BIN_START + p * COMB_BIN_STEP) *
                         (double)SAMPLE_RATE / (double)ANALYSIS_FFT;
        comb += COMB_PARTIAL_AMP * sin(6.283185307179586 * f * t);
      }
    }
    s.data[i] = noise + comb;
  }
  return s;
}

static void run_engine(SpecbleachSmoothingMode mode, const Signal* in,
                       float dftt_strength, float* out) {
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize((uint32_t)SAMPLE_RATE, ENGINE_FRAME_MS);
  TEST_ASSERT(handle != NULL, "engine init failed");

  SpecbleachDenoiserParameters params = {
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .reduction_gain = 0.05F,
      .smoothing_factor = 0.5F,
      .smoothing_mode = mode,
      .dftt_strength = dftt_strength,
      .whitening_factor = 0.0F,
      .residual_listen = false,
  };
  TEST_ASSERT(
      specbleach_denoiser_load_parameters(handle, &params, sizeof(params)),
      "load learn params failed");

  const size_t learn_len = (size_t)(LEARN_SECONDS * (float)SAMPLE_RATE);
  TEST_ASSERT(
      specbleach_denoiser_process(handle, (uint32_t)learn_len, in->data, out),
      "learn pass failed");

  params.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(
      specbleach_denoiser_load_parameters(handle, &params, sizeof(params)),
      "load reduce params failed");
  TEST_ASSERT(
      specbleach_denoiser_process(handle, (uint32_t)(in->len - learn_len),
                                  in->data + learn_len, out + learn_len),
      "reduce pass failed");

  specbleach_denoiser_free(handle);
}

// Mean per-bin temporal variance (dB^2) of the analysis-STFT magnitude
static float speckle_flicker(const float* x, size_t len, size_t offset) {
  const uint32_t nbins = (NOISE_BIN_END - NOISE_BIN_START) / NOISE_BIN_STEP;
  const size_t frames = (len - offset - ANALYSIS_FFT) / ANALYSIS_HOP;
  float* db = (float*)calloc((size_t)nbins * frames, sizeof(float));
  TEST_ASSERT(db != NULL, "flicker buffer failed");

  for (size_t t = 0U; t < frames; t++) {
    const float* frame = x + offset + t * ANALYSIS_HOP;
    for (uint32_t b = 0U; b < nbins; b++) {
      const uint32_t k = NOISE_BIN_START + b * NOISE_BIN_STEP;
      double re = 0.0;
      double im = 0.0;
      for (uint32_t n = 0U; n < ANALYSIS_FFT; n++) {
        const double w = 0.5 - 0.5 * cos(6.283185307179586 * (double)n /
                                         (double)ANALYSIS_FFT);
        const double s = (double)frame[n] * w;
        const double ph =
            6.283185307179586 * (double)(k * n) / (double)ANALYSIS_FFT;
        re += s * cos(ph);
        im -= s * sin(ph);
      }
      const float mag = (float)sqrt(re * re + im * im) + 1e-12F;
      db[(size_t)b * frames + t] = 20.0F * log10f(mag);
    }
  }

  double var_sum = 0.0;
  for (uint32_t b = 0U; b < nbins; b++) {
    double mean = 0.0;
    for (size_t t = 0U; t < frames; t++) {
      mean += db[(size_t)b * frames + t];
    }
    mean /= (double)frames;
    double var = 0.0;
    for (size_t t = 0U; t < frames; t++) {
      const double d = db[(size_t)b * frames + t] - mean;
      var += d * d;
    }
    var_sum += var / (double)frames;
  }
  free(db);
  return (float)(var_sum / (double)nbins);
}

// Mean magnitude (dB) at the comb bins across time
static float comb_level(const float* x, size_t len, size_t offset) {
  const size_t frames = (len - offset - ANALYSIS_FFT) / ANALYSIS_HOP;
  double sum = 0.0;
  uint32_t count = 0U;
  for (size_t t = 0U; t < frames; t++) {
    const float* frame = x + offset + t * ANALYSIS_HOP;
    for (uint32_t p = 0U; p < COMB_PARTIALS; p++) {
      const uint32_t k = COMB_BIN_START + p * COMB_BIN_STEP;
      double re = 0.0;
      double im = 0.0;
      for (uint32_t n = 0U; n < ANALYSIS_FFT; n++) {
        const double w = 0.5 - 0.5 * cos(6.283185307179586 * (double)n /
                                         (double)ANALYSIS_FFT);
        const double s = (double)frame[n] * w;
        const double ph =
            6.283185307179586 * (double)(k * n) / (double)ANALYSIS_FFT;
        re += s * cos(ph);
        im -= s * sin(ph);
      }
      sum += 20.0 * log10(sqrt(re * re + im * im) + 1e-12);
      count++;
    }
  }
  return (float)(sum / (double)count);
}

int main(void) {
  printf("Running DFTT quality instrument...\n\n");

  Signal in = make_signal(true);
  float* out_nlm = (float*)calloc(in.len, sizeof(float));
  float* out_dftt = (float*)calloc(in.len, sizeof(float));
  TEST_ASSERT(out_nlm && out_dftt, "output allocation failed");

  run_engine(SPECBLEACH_SMOOTHING_NLM_2D, &in, 1.0F, out_nlm);
  run_engine(SPECBLEACH_SMOOTHING_NLM_2D_DFTT, &in, 1.0F, out_dftt);

  const size_t offset =
      (size_t)(WARMUP_SECONDS * (float)SAMPLE_RATE) + ANALYSIS_FFT;

  float flicker_nlm = 0.0F;
  float flicker_dftt = 0.0F;
  float comb_nlm = 0.0F;
  float comb_dftt = 0.0F;
  for (size_t i = offset; i < in.len; i++) {
    TEST_ASSERT(isfinite(out_nlm[i]) && isfinite(out_dftt[i]),
                "outputs must be finite");
  }
  flicker_nlm = speckle_flicker(out_nlm, in.len, offset);
  flicker_dftt = speckle_flicker(out_dftt, in.len, offset);
  comb_nlm = comb_level(out_nlm, in.len, offset);
  comb_dftt = comb_level(out_dftt, in.len, offset);

  printf(
      "  Speckle flicker (dB^2, lower is better):  NLM %.3f  |  NLM+DFTT "
      "%.3f\n",
      (double)flicker_nlm, (double)flicker_dftt);
  printf(
      "  Comb level (dB):                           NLM %.3f  |  NLM+DFTT "
      "%.3f\n",
      (double)comb_nlm, (double)comb_dftt);

  TEST_ASSERT(flicker_nlm > 0.0F, "baseline flicker must be positive");
  TEST_ASSERT(flicker_dftt <= flicker_nlm, "DFTT must not increase flicker");

  const float improvement = 100.0F * (flicker_nlm - flicker_dftt) / flicker_nlm;
  printf("  DFTT flicker improvement: %.1f%% (min %.1f%%)\n",
         (double)improvement, (double)MIN_SPECKLE_IMPROVEMENT_PCT);
  TEST_ASSERT(improvement >= MIN_SPECKLE_IMPROVEMENT_PCT,
              "DFTT speckle improvement below required threshold");

  const float comb_delta = fabsf(comb_dftt - comb_nlm);
  TEST_ASSERT(comb_delta <= COMB_TOLERANCE_DB,
              "DFTT must not muffle the comb (comb level within tolerance)");

  /* Strength sweep (#155 plumbing): non-default dftt_strength values must
   * flow through load_parameters into the kill threshold. Contract, not
   * exact values: output stays finite, flicker never increases vs NLM, comb
   * stays preserved (same guard as the default strength). */
  for (uint32_t s = 0U; s < 2U; s++) {
    const float strength = (s == 0U) ? 0.5F : 2.0F;
    TEST_ASSERT(strength > 0.0F && strength < DFTT_STRENGTH_MAX,
                "sweep strength must be in (0, DFTT_STRENGTH_MAX)");
    float* out_s = (float*)calloc(in.len, sizeof(float));
    TEST_ASSERT(out_s != NULL, "strength sweep allocation failed");
    run_engine(SPECBLEACH_SMOOTHING_NLM_2D_DFTT, &in, strength, out_s);
    float ok = 1.0F;
    for (size_t i = offset; i < in.len; i++) {
      if (!isfinite(out_s[i])) {
        ok = 0.0F;
        break;
      }
    }
    TEST_ASSERT(ok > 0.5F, "strength-swept output must be finite");
    const float flicker_s = speckle_flicker(out_s, in.len, offset);
    const float comb_s = comb_level(out_s, in.len, offset);
    const float impr_s = 100.0F * (flicker_nlm - flicker_s) / flicker_nlm;
    printf("  strength %.1f: flicker %.3f (%.1f%%), comb %.3f\n",
           (double)strength, (double)flicker_s, (double)impr_s, (double)comb_s);
    TEST_ASSERT(flicker_s <= flicker_nlm,
                "strength sweep must not increase flicker vs NLM");
    TEST_ASSERT(fabsf(comb_s - comb_nlm) <= COMB_TOLERANCE_DB,
                "strength sweep must preserve comb");
    free(out_s);
  }

  free(in.data);
  free(out_nlm);
  free(out_dftt);

  printf("\n✅ DFTT quality instrument passed!\n");
  return 0;
}
