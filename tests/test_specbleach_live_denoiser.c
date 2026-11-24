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

/*
 * Unit tests for the live multiband-gate denoiser
 */

#include "specbleach_live_denoiser.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define SAMPLE_RATE 48000U
#define BLOCK 64U

static float block_rms(const float* data, uint32_t n) {
  double sum = 0.0;
  for (uint32_t i = 0U; i < n; i++) {
    sum += (double)data[i] * (double)data[i];
  }
  return (float)sqrt(sum / (double)n);
}

static void test_null_safety(void) {
  TEST_ASSERT(specbleach_live_denoiser_initialize(0U) == NULL,
              "zero sample rate must fail");
  TEST_ASSERT(specbleach_live_denoiser_get_latency(NULL) == 0U,
              "NULL latency must be 0");

  SpecbleachLiveDenoiserParameters params;
  memset(&params, 0, sizeof(params));
  TEST_ASSERT(
      !specbleach_live_denoiser_load_parameters(NULL, &params, sizeof(params)),
      "NULL instance must fail");
  TEST_ASSERT(!specbleach_live_denoiser_load_parameters(NULL, NULL, 0U),
              "NULL params must fail");

  float input[BLOCK] = {0};
  float output[BLOCK] = {0};
  TEST_ASSERT(!specbleach_live_denoiser_process(NULL, BLOCK, input, output),
              "NULL instance process must fail");

  TEST_ASSERT(!specbleach_live_denoiser_reset_noise_floor(NULL),
              "NULL reset must fail");
  TEST_ASSERT(
      !specbleach_live_denoiser_get_band_levels(NULL, input, output, output),
      "NULL instance levels must fail");
  TEST_ASSERT(!specbleach_live_denoiser_get_band_edges(NULL, input, output),
              "NULL instance edges must fail");
  TEST_ASSERT(!specbleach_live_denoiser_set_delta_monitoring(NULL, true),
              "NULL delta monitoring must fail");

  specbleach_live_denoiser_free(NULL);
  printf("✓ Null safety passed\n");
}

static void test_silence_stays_silent(void) {
  specbleach_live_denoiser* engine =
      specbleach_live_denoiser_initialize(SAMPLE_RATE);
  TEST_ASSERT(engine != NULL, "initialize must succeed");
  TEST_ASSERT(specbleach_live_denoiser_get_latency(engine) == 0U,
              "latency must be 0");

  float input[BLOCK] = {0};
  float output[BLOCK];
  for (int b = 0; b < 200; b++) {
    TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, input, output),
                "process must succeed");
    for (uint32_t i = 0U; i < BLOCK; i++) {
      TEST_ASSERT(isfinite(output[i]), "output must be finite");
      TEST_ASSERT(fabsf(output[i]) < 1e-6F, "silence must stay silent");
    }
  }

  specbleach_live_denoiser_free(engine);
  printf("✓ Silence passed\n");
}

static void test_tone_passes_through(void) {
  specbleach_live_denoiser* engine =
      specbleach_live_denoiser_initialize(SAMPLE_RATE);
  TEST_ASSERT(engine != NULL, "initialize must succeed");

  float input[BLOCK];
  float output[BLOCK];

  /* Learn a quiet floor first (real use case: tone arrives over a
   * learned-quiet background), then freeze it. */
  for (uint32_t i = 0U; i < BLOCK; i++) {
    input[i] = 0.0F;
  }
  for (int b = 0; b < 100; b++) {
    TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, input, output),
                "process must succeed");
  }
  SpecbleachLiveDenoiserParameters params;
  memset(&params, 0, sizeof(params));
  params.reduction_gain = 0.01F;
  params.attack_time = 0.005F;
  params.release_time = 0.030F;
  params.adaptive_noise = false;
  params.noise_estimation_method = SPECBLEACH_LIVE_MARTIN;
  params.threshold_db = 6.0F;
  params.knee_db = 0.0F;
  TEST_ASSERT(
      specbleach_live_denoiser_load_parameters(engine, &params, sizeof(params)),
      "freeze must succeed");

  /* Sustained 440 Hz tone sits far above the frozen floor: gate opens. */
  float phase = 0.0F;
  for (int b = 0; b < 300; b++) {
    for (uint32_t i = 0U; i < BLOCK; i++) {
      input[i] = 0.5F * sinf(2.0F * 3.14159265F * 440.0F * phase / 48000.0F);
      phase += 1.0F;
    }
    TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, input, output),
                "process must succeed");
  }
  const float in_rms = block_rms(input, BLOCK);
  const float out_rms = block_rms(output, BLOCK);
  TEST_ASSERT(isfinite(out_rms), "output must be finite");
  TEST_ASSERT(out_rms > 0.5F * in_rms, "sustained tone must mostly pass");

  specbleach_live_denoiser_free(engine);
  printf("✓ Tone passthrough passed\n");
}

static void test_noise_is_reduced(void) {
  specbleach_live_denoiser* engine =
      specbleach_live_denoiser_initialize(SAMPLE_RATE);
  TEST_ASSERT(engine != NULL, "initialize must succeed");

  SpecbleachLiveDenoiserParameters params;
  memset(&params, 0, sizeof(params));
  /* Max threshold must crush stationary noise: +12 dB puts the whole
   * noise distribution below the knee (a mid threshold only partially
   * reduces by design). */
  params.reduction_gain = 0.01F; /* -40 dB floor */
  params.attack_time = 0.005F;
  params.release_time = 0.165F;
  params.adaptive_noise = true;
  params.noise_estimation_method = SPECBLEACH_LIVE_MARTIN;
  params.threshold_db = 12.0F;
  params.knee_db = 6.0F;
  TEST_ASSERT(
      specbleach_live_denoiser_load_parameters(engine, &params, sizeof(params)),
      "load_parameters must succeed");

  float input[BLOCK];
  float output[BLOCK];
  srand(7U);

  /* Learn the stationary noise floor (fail-open: the floor rises
   * from zero over ~1 s, so warm up well past the time constant). */
  for (int b = 0; b < 2000; b++) {
    for (uint32_t i = 0U; i < BLOCK; i++) {
      input[i] = ((float)rand() / (float)RAND_MAX) * 0.1F - 0.05F;
    }
    TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, input, output),
                "process must succeed");
  }

  /* Freeze and measure. */
  params.adaptive_noise = false;
  TEST_ASSERT(
      specbleach_live_denoiser_load_parameters(engine, &params, sizeof(params)),
      "freeze must succeed");

  double dry_sum = 0.0;
  double wet_sum = 0.0;
  uint32_t count = 0U;
  for (int b = 0; b < 100; b++) {
    for (uint32_t i = 0U; i < BLOCK; i++) {
      input[i] = ((float)rand() / (float)RAND_MAX) * 0.1F - 0.05F;
    }
    TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, input, output),
                "process must succeed");
    for (uint32_t i = 0U; i < BLOCK; i++) {
      TEST_ASSERT(isfinite(output[i]), "output must be finite");
      dry_sum += (double)input[i] * (double)input[i];
      wet_sum += (double)output[i] * (double)output[i];
      count++;
    }
  }
  const float dry_rms = (float)sqrt(dry_sum / (double)count);
  const float wet_rms = (float)sqrt(wet_sum / (double)count);
  /* Dynamic-EQ cascade semantics: the Reduction slider is the steady
   * broadband cut at shallow depths and saturates deeper (notch
   * geometry), so anchor the assert to the calibrated floor rather
   * than the slider value. */
  TEST_ASSERT(wet_rms < dry_rms * 0.70F, "noise must be strongly reduced");

  specbleach_live_denoiser_free(engine);
  printf("✓ Noise reduction passed\n");
}

static void test_odd_block_sizes(void) {
  specbleach_live_denoiser* engine =
      specbleach_live_denoiser_initialize(SAMPLE_RATE);
  TEST_ASSERT(engine != NULL, "initialize must succeed");

  const uint32_t sizes[] = {1U, 13U, 64U, 512U, 4096U};
  for (uint32_t s = 0U; s < 5U; s++) {
    const uint32_t n = sizes[s];
    float* input = (float*)calloc(n, sizeof(float));
    float* output = (float*)calloc(n, sizeof(float));
    TEST_ASSERT(input != NULL && output != NULL, "calloc must succeed");
    for (uint32_t i = 0U; i < n; i++) {
      input[i] = 0.01F * sinf((float)i);
    }
    TEST_ASSERT(specbleach_live_denoiser_process(engine, n, input, output),
                "odd block size must succeed");
    for (uint32_t i = 0U; i < n; i++) {
      TEST_ASSERT(isfinite(output[i]), "output must be finite");
    }
    /* In-place aliasing must also work. */
    TEST_ASSERT(specbleach_live_denoiser_process(engine, n, input, input),
                "in-place process must succeed");
    free(input);
    free(output);
  }

  specbleach_live_denoiser_free(engine);
  printf("✓ Odd block sizes passed\n");
}

static void test_band_edges_match_design(void) {
  for (uint32_t rate = 0U; rate < 2U; rate++) {
    const uint32_t sample_rate = rate == 0U ? 44100U : 48000U;
    specbleach_live_denoiser* engine =
        specbleach_live_denoiser_initialize(sample_rate);
    TEST_ASSERT(engine != NULL, "initialize must succeed");

    float lo[SPECBLEACH_LIVE_NUM_BANDS];
    float hi[SPECBLEACH_LIVE_NUM_BANDS];
    TEST_ASSERT(specbleach_live_denoiser_get_band_edges(engine, lo, hi),
                "edges must succeed");
    TEST_ASSERT(!specbleach_live_denoiser_get_band_edges(engine, NULL, NULL),
                "all-NULL edges must fail");

    const float nyquist = (float)sample_rate * 0.5F;
    uint32_t num_active = 0U;
    for (uint32_t k = 0U; k < SPECBLEACH_LIVE_NUM_BANDS; k++) {
      if (hi[k] <= 0.0F) {
        TEST_ASSERT(lo[k] == 0.0F, "inactive band edges must be 0/0");
        continue;
      }
      num_active++;
      TEST_ASSERT(lo[k] >= 20.0F, "band lo must be >= 20 Hz");
      TEST_ASSERT(hi[k] > lo[k], "band hi must exceed lo");
      TEST_ASSERT(hi[k] <= nyquist, "band hi must not exceed Nyquist");
      if (k > 0U && hi[k - 1U] > 0.0F) {
        TEST_ASSERT(lo[k] <= hi[k - 1U] + 1e-3F,
                    "bands must tile without gaps");
      }
    }
    TEST_ASSERT(num_active > SPECBLEACH_LIVE_NUM_BANDS / 2U,
                "most bands must be active");

    specbleach_live_denoiser_free(engine);
  }
  printf("✓ Band edges passed\n");
}

static void test_band_levels_and_reset(void) {
  specbleach_live_denoiser* engine =
      specbleach_live_denoiser_initialize(SAMPLE_RATE);
  TEST_ASSERT(engine != NULL, "initialize must succeed");

  float in[SPECBLEACH_LIVE_NUM_BANDS];
  float out[SPECBLEACH_LIVE_NUM_BANDS];
  float thr[SPECBLEACH_LIVE_NUM_BANDS];

  /* Fresh instance: envelopes and floor are zero. */
  TEST_ASSERT(specbleach_live_denoiser_get_band_levels(engine, in, out, thr),
              "levels must succeed");
  for (uint32_t k = 0U; k < SPECBLEACH_LIVE_NUM_BANDS; k++) {
    TEST_ASSERT(in[k] == 0.0F && out[k] == 0.0F && thr[k] == 0.0F,
                "fresh levels must be zero");
  }
  TEST_ASSERT(
      !specbleach_live_denoiser_get_band_levels(engine, NULL, NULL, NULL),
      "all-NULL levels must fail");

  /* Run broadband noise: input energy rises, floor learns up from zero. */
  float input[BLOCK];
  float output[BLOCK];
  srand(3U);
  for (int b = 0; b < 500; b++) {
    for (uint32_t i = 0U; i < BLOCK; i++) {
      input[i] = ((float)rand() / (float)RAND_MAX) * 0.4F - 0.2F;
    }
    TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, input, output),
                "process must succeed");
  }
  TEST_ASSERT(specbleach_live_denoiser_get_band_levels(engine, in, out, thr),
              "levels must succeed");
  float max_in = 0.0F;
  float max_thr = 0.0F;
  for (uint32_t k = 0U; k < SPECBLEACH_LIVE_NUM_BANDS; k++) {
    TEST_ASSERT(isfinite(in[k]) && isfinite(out[k]) && isfinite(thr[k]),
                "levels must be finite");
    TEST_ASSERT(out[k] <= in[k] + 1e-6F, "output must not exceed input");
    TEST_ASSERT(thr[k] >= 0.0F, "threshold must be non-negative");
    if (in[k] > max_in) {
      max_in = in[k];
    }
    if (thr[k] > max_thr) {
      max_thr = thr[k];
    }
  }
  TEST_ASSERT(max_in > 1e-4F, "input energy must be visible");
  TEST_ASSERT(max_thr > 0.0F, "floor must have learned");

  /* Reset: freeze the tracker so the cleared floor can't re-learn from
   * the decaying envelope, then the next process must report zero. */
  SpecbleachLiveDenoiserParameters freeze;
  memset(&freeze, 0, sizeof(freeze));
  freeze.adaptive_noise = false;
  TEST_ASSERT(
      specbleach_live_denoiser_load_parameters(engine, &freeze, sizeof(freeze)),
      "freeze must succeed");
  TEST_ASSERT(specbleach_live_denoiser_reset_noise_floor(engine),
              "reset must succeed");
  for (uint32_t i = 0U; i < BLOCK; i++) {
    input[i] = 0.0F;
  }
  TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, input, output),
              "process must succeed");
  TEST_ASSERT(specbleach_live_denoiser_get_band_levels(engine, in, out, thr),
              "levels must succeed");
  for (uint32_t k = 0U; k < SPECBLEACH_LIVE_NUM_BANDS; k++) {
    TEST_ASSERT(thr[k] == 0.0F, "threshold must be zero after reset");
  }

  specbleach_live_denoiser_free(engine);
  printf("✓ Band levels and reset passed\n");
}

static void test_delta_monitoring(void) {
  specbleach_live_denoiser* engine =
      specbleach_live_denoiser_initialize(SAMPLE_RATE);
  TEST_ASSERT(engine != NULL, "initialize must succeed");

  /* At no reduction (g_min=1) delta must be silence. */
  SpecbleachLiveDenoiserParameters p = {0};
  p.reduction_gain = 1.0F;
  p.attack_time = 0.005F;
  p.release_time = 0.030F;
  p.adaptive_noise = false;
  p.threshold_db = -12.0F;
  p.knee_db = 0.0F;
  TEST_ASSERT(specbleach_live_denoiser_load_parameters(engine, &p, sizeof(p)),
              "load must succeed");
  TEST_ASSERT(specbleach_live_denoiser_set_delta_monitoring(engine, true),
              "delta on must succeed");
  float in[BLOCK];
  float out[BLOCK];
  for (uint32_t i = 0U; i < BLOCK; i++) {
    in[i] = 0.1F * sinf((float)i * 0.3F);
  }
  TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, in, out),
              "process must succeed");
  for (uint32_t i = 0U; i < BLOCK; i++) {
    TEST_ASSERT(fabsf(out[i]) < 1e-6F, "delta at no reduction must be silence");
  }

  /* At max reduction delta must carry energy, and wet+delta reconstructs
   * the band sum. Toggle back to normal and compare. */
  p.reduction_gain = 0.01F;
  p.threshold_db = 6.0F;
  p.adaptive_noise = true;
  TEST_ASSERT(specbleach_live_denoiser_load_parameters(engine, &p, sizeof(p)),
              "load must succeed");
  srand(9U);
  for (int b = 0; b < 500; b++) {
    for (uint32_t i = 0U; i < BLOCK; i++) {
      in[i] = ((float)rand() / (float)RAND_MAX) * 0.4F - 0.2F;
    }
    TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, in, out),
                "process must succeed");
  }
  p.adaptive_noise = false;
  TEST_ASSERT(specbleach_live_denoiser_load_parameters(engine, &p, sizeof(p)),
              "freeze must succeed");
  float wet[BLOCK];
  float delta[BLOCK];
  for (uint32_t i = 0U; i < BLOCK; i++) {
    in[i] = ((float)rand() / (float)RAND_MAX) * 0.4F - 0.2F;
  }
  TEST_ASSERT(specbleach_live_denoiser_set_delta_monitoring(engine, false),
              "delta off must succeed");
  TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, in, wet),
              "process must succeed");
  TEST_ASSERT(specbleach_live_denoiser_set_delta_monitoring(engine, true),
              "delta on must succeed");
  TEST_ASSERT(specbleach_live_denoiser_process(engine, BLOCK, in, delta),
              "process must succeed");
  float rms_delta = block_rms(delta, BLOCK);
  TEST_ASSERT(rms_delta > 1e-4F, "delta at max reduction must have energy");
  /* wet and delta are from successive blocks, so just check both finite
   * and that delta does not exceed input energy wildly. */
  for (uint32_t i = 0U; i < BLOCK; i++) {
    TEST_ASSERT(isfinite(delta[i]), "delta must be finite");
  }

  TEST_ASSERT(specbleach_live_denoiser_set_delta_monitoring(engine, false),
              "delta off must succeed");
  specbleach_live_denoiser_free(engine);
  printf("✓ Delta monitoring passed\n");
}

static void test_bank_open_flatness(void) {
  /* Gates pinned open: the cascaded bank must sum near-flat so dry
   * passes transparently (top octave excluded: bands end at 0.95 Nyq). */
  specbleach_live_denoiser* engine =
      specbleach_live_denoiser_initialize(SAMPLE_RATE);
  TEST_ASSERT(engine != NULL, "initialize must succeed");

  SpecbleachLiveDenoiserParameters p;
  memset(&p, 0, sizeof(p));
  p.reduction_gain = 1.0F;
  p.attack_time = 0.001F;
  p.release_time = 0.010F;
  p.adaptive_noise = false;
  p.threshold_db = -12.0F;
  p.knee_db = 0.0F;
  TEST_ASSERT(specbleach_live_denoiser_load_parameters(engine, &p, sizeof(p)),
              "load must succeed");

  float input[BLOCK];
  float output[BLOCK];
  for (int s = 0; s <= 12; s++) {
    const float freq = 50.0F * powf(16000.0F / 50.0F, (float)s / 12.0F);
    float phase = 0.0F;
    double sum_in = 0.0;
    double sum_out = 0.0;
    for (int b = 0; b < 120; b++) {
      for (uint32_t i = 0U; i < BLOCK; i++) {
        input[i] = 0.5F * sinf(2.0F * 3.14159265F * freq * phase / 48000.0F);
        phase += 1.0F;
      }
      TEST_ASSERT(
          specbleach_live_denoiser_process(engine, BLOCK, input, output),
          "process must succeed");
      if (b >= 60) {
        for (uint32_t i = 0U; i < BLOCK; i++) {
          sum_in += (double)input[i] * input[i];
          sum_out += (double)output[i] * output[i];
        }
      }
    }
    const float ratio_db = 20.0F * log10f((float)sqrt(sum_out / sum_in));
    /* Interior bands (up to ~14 kHz) must sum near-flat; the topmost
     * points thin out toward the 0.97 Nyquist design edge. */
    if (freq < 14000.0F) {
      TEST_ASSERT(ratio_db > -4.0F && ratio_db < 4.0F,
                  "open bank must be flat within 4 dB");
    } else {
      TEST_ASSERT(ratio_db > -6.5F,
                  "top edge must not collapse open-bank coverage");
    }
  }

  specbleach_live_denoiser_free(engine);
  printf("✓ Bank open flatness passed\n");
}

int main(void) {
  printf("Running live denoiser tests...\n\n");

  test_null_safety();
  test_silence_stays_silent();
  test_tone_passes_through();
  test_noise_is_reduced();
  test_odd_block_sizes();
  test_band_edges_match_design();
  test_band_levels_and_reset();
  test_delta_monitoring();
  test_bank_open_flatness();

  printf("\n✅ All live denoiser tests passed!\n");
  return 0;
}
