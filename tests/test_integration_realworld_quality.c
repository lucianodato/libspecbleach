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

// Objective quality metrics over real-world fixtures (see
// autoresearch/build_fixtures.py for fixture generation). Each case provides
// the clean signal and the noise separately; the mix is synthesized here at a
// designated SNR, so full references exist for every metric:
// - Att: noise attenuation (dB) over noise-only frames (lower bound on
//   usable rejection; higher = deeper).
// - SD: log-spectral distortion of signal-active frames against the clean
//   reference (dB, higher = more "underwater" / spectral damage).
// - MNI v2: musical-noise index of the output normalized by the unprocessed
//   noise's own fluctuation in the same frames (1.0 = no musical noise
//   added; > 1 = flickery residual introduced).
// - Decay deficit: post-transient decay shortfall in dB (higher = temporal
//   smearing / echoing after the clean signal drops).
// - SI-SDR: scale-invariant signal-to-distortion ratio of the output against
//   the clean reference (dB, higher = cleaner overall).
// - SAR: signal-to-artifacts ratio from a BSS_eval-style projection of the
//   output onto the span of {clean, noise} (dB, higher = less musical
//   noise/processing texture). SI-SDR/SAR are report-only (no CI gates yet).
// - sd/att: damage per dB rejected, from existing numbers (lower = cheaper
//   cleaning).
//
// Engine configuration mirrors the noise-repellent plugin defaults exactly
// (buildEngineParams at default slider values): that is the config the
// plugin actually runs and the one subjective listening tests judge. The
// sweep axis is the smoothing slider; reduction stays at the plugin default.
// CI runs the plugin-default smoothing (0.0) only; set
// SPECBLEACH_REALWORLD_FULL=1 to sweep all smoothing values locally.

#include "../src/shared/stft/fft_transform.h"
#include "specbleach_denoiser.h"

#include <math.h>
#include <sndfile.h>
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

#define TEST_DATA_DIR_PATTERN TEST_DATA_DIR "realworld/%s/%s.wav"

// Must match NOISE_LEAD_IN_SEC in autoresearch/build_fixtures.py
#define NOISE_LEAD_IN_SEC 1.5F
#define BODY_SEC 8.0F

#define FRAME_MS                                                               \
  46.0F // plugin default STFT frame size (index 2 of 23/32/46/64/93)

// Analysis STFT for metric computation (independent of the engine's STFT)
#define ANALYSIS_FFT 1024U
#define ANALYSIS_HOP 256U

// Frame gating against the clean reference
#define ACTIVE_DYNAMIC_DB 25.0F // frames within this of the clean peak = active
#define NOISE_DYNAMIC_DB                                                       \
  35.0F // frames below this of the clean peak = noise-only
#define BOUNDARY_EXCLUSION_FRAMES 2U // skip transition zones between the two
// The release tail after signal drops must not pollute noise-only metrics
#define RELEASE_TAIL_EXCLUSION_SEC 0.2F
// First reduce frames after finalization (smoother state warm-up)
#define INITIAL_REDUCE_EXCLUSION_SEC 0.1F
// Attenuation is measured on bins where the learned noise actually lives
#define NOISE_BIN_FLOOR_DB 20.0F
#define MIN_NOISE_BINS 4U
#define MIN_METRIC_FRAMES 10U

// SD metric tuning
#define PEAK_DYNAMIC_DB 25.0F // bin inclusion floor relative to frame peak
#define DEV_CLIP_DB 20.0F     // per-bin deviation clamp (outlier robustness)

// SD split diagnostic (report-only, no gates): where the damage lives.
// Bands in Hz; transient = clean frame-energy onset + short tail.
#define SD_BAND_MF_HZ 1000.0F
#define SD_BAND_HF_HZ 4000.0F
#define SD_TRANSIENT_ATTACK_DB 6.0F
#define SD_TRANSIENT_TAIL_FRAMES 1U

// Decay metric tuning
#define ENV_WINDOW 512U
#define ENV_HOP 256U
#define ECHO_DROP_DB 15.0F // clean decay that qualifies as a phrase/word end

// Plugin-default engine configuration (mirrors noise-repellent's
// buildEngineParams at default slider values). Reduction stays at the plugin
// default; the sweep axis is the smoothing slider.
#define PLUGIN_REDUCTION_DB 12.0F // Master Reduction default (linked tonal)
#define DFTT_STRENGTH_DOUBLING_DB 12.0F // kDfttStrengthDoublingReductionDb

// Synthetic white-noise case: deterministic LCG (high 24 bits) so the case is
// reproducible across platforms without shipping a fixture binary. Reproduces
// the low-SNR white-noise guitar scenario from the BM3D listening tests.
#define WHITE_NOISE_LCG_SEED 12345U
#define WHITE_NOISE_LCG_MULT 1664525U
#define WHITE_NOISE_LCG_INC 1013904223U
#define WHITE_NOISE_AMPLITUDE 0.2F

typedef struct SmoothingConfig {
  const char* name;
  float smoothing; // Smoothing slider 0..1 (plugin default 0.0)
} SmoothingConfig;

static const SmoothingConfig SMOOTHING_CONFIGS[3] = {
    {"smoothing-off", 0.0F},
    {"smoothing-30", 0.30F},
    {"smoothing-80", 0.80F},
};

#define NUM_CASES 5U

typedef struct Case {
  const char* name;
  float snr_db;
  const char* clean_source; // clean fixture dir (defaults to name)
  bool white_noise;         // synthetic deterministic noise (no fixture)
} Case;

// Fixture cases with their designated SNR (matrix: stationary/non-stationary
// noise x speech/music x moderate/hard SNR). case_guitar_white reuses the
// electric-guitar clean track against synthetic white noise at 3 dB: the
// low-SNR white-noise scenario where BM3D vs NLM differences actually show.
static const Case CASES[NUM_CASES] = {
    {"case_voice_fan", 12.0F, NULL, false},
    {"case_voice_city", 6.0F, NULL, false},
    {"case_guitar_electric", 18.0F, NULL, false},
    {"case_drums_fridge", 6.0F, NULL, false},
    {"case_guitar_white", 3.0F, "case_guitar_electric", true},
};

#define NUM_MODES 4U

// Per-mode regression gates for the CI (plugin-default smoothing) config,
// aggregated over cases: [Temporal, NLM 2D, NLM 2D + DFTT, BM3D]
// Calibrated on delay-compensated scores at 12 dB reduction (att ~10.9,
// sd ~3.0, mni ~1.22, decay ~0.0) with headroom for cross-platform variance.
static const float ATT_GATES_DB[4] = {8.0F, 8.0F, 8.0F,
                                      8.0F}; // minimum rejection
static const float SD_GATES_DB[4] = {5.0F, 5.0F, 5.0F, 5.0F}; // maximum damage
static const float MNI_GATES[4] = {1.4F, 1.4F, 1.4F, 1.4F}; // max musical noise
static const float DECAY_GATES_DB[4] = {2.0F, 2.0F, 2.0F, 2.0F}; // max smearing

typedef struct Metrics {
  float att_db;     // attenuation on learned-noise bins over noise-only frames
  float att_all_db; // attenuation over all bins (whitening-sensitive reference)
  float sd_db;      // log-spectral distortion on active frames
  float sd_lf_db;   // SD on LF bins (< 1 kHz), report-only
  float sd_mf_db;   // SD on MF bins (1-4 kHz), report-only
  float sd_hf_db;   // SD on HF bins (> 4 kHz), report-only
  float sd_tr_db;   // SD on transient active frames, report-only
  float sd_st_db;   // SD on steady active frames, report-only
  float mni_ratio;  // output shape-CV / unprocessed noise shape-CV
  float mni_out;    // raw output shape-CV (for reference when tuning)
  float decay_db;   // mean post-transient decay deficit
  float sisdr_db;   // SI-SDR of output vs clean (higher = cleaner)
  float sar_db;     // signal-to-artifacts ratio vs span{clean,noise}
} Metrics;

typedef struct Signals {
  const char* case_name;
  float* clean; // silence prefix + body
  float* noise; // scaled, covers the full duration
  float* mix;
  float* snapshot; // pristine mix copy for the immutability probe
  uint32_t total_samples;
  uint32_t lead_samples;
  uint32_t body_samples;
} Signals;

static char fixture_path[512];

static const char* fixture_file(const char* case_name, const char* file) {
  snprintf(fixture_path, sizeof(fixture_path), TEST_DATA_DIR_PATTERN, case_name,
           file);
  return fixture_path;
}

static float* load_wav(const char* path, uint32_t* sample_rate,
                       uint32_t* length) {
  SF_INFO info;
  memset(&info, 0, sizeof(info));
  SNDFILE* file = sf_open(path, SFM_READ, &info);
  TEST_ASSERT(file != NULL, "fixture open");
  TEST_ASSERT(info.channels == 1, "fixture must be mono");
  TEST_ASSERT(info.samplerate > 0, "fixture sample rate");
  float* data = (float*)malloc((size_t)info.frames * sizeof(float));
  TEST_ASSERT(data != NULL, "fixture alloc");
  TEST_ASSERT(sf_readf_float(file, data, info.frames) == info.frames,
              "fixture read");
  sf_close(file);
  *sample_rate = (uint32_t)info.samplerate;
  *length = (uint32_t)info.frames;
  return data;
}

static float* make_white_noise(uint32_t length) {
  float* data = (float*)malloc((size_t)length * sizeof(float));
  TEST_ASSERT(data != NULL, "white noise alloc");
  uint32_t state = WHITE_NOISE_LCG_SEED;
  for (uint32_t n = 0U; n < length; n++) {
    state = state * WHITE_NOISE_LCG_MULT + WHITE_NOISE_LCG_INC;
    data[n] = (((float)(state >> 8) / 16777216.0F) * 2.0F - 1.0F) *
              WHITE_NOISE_AMPLITUDE;
  }
  return data;
}

static double energy(const float* signal, uint32_t length) {
  double sum = 0.0;
  for (uint32_t n = 0U; n < length; n++) {
    sum += (double)signal[n] * (double)signal[n];
  }
  return sum;
}

// RMS-matched mixing at the designated SNR (mirrors the dataset generator,
// but the noise level is consistent across lead-in and body)
static void mix_case(const float* clean_body, const float* noise_raw,
                     uint32_t body_samples, uint32_t lead_samples, float snr_db,
                     Signals* out) {
  const double clean_rms =
      sqrt(energy(clean_body, body_samples) / (double)body_samples + 1e-12);
  const double noise_rms = sqrt(energy(noise_raw + lead_samples, body_samples) /
                                    (double)body_samples +
                                1e-12);
  const double factor =
      clean_rms / (noise_rms * pow(10.0, (double)snr_db / 20.0));

  const uint32_t total = lead_samples + body_samples;
  memset(out->clean, 0, total * sizeof(float));
  for (uint32_t n = 0U; n < total; n++) {
    out->noise[n] = (float)((double)noise_raw[n] * factor);
  }
  for (uint32_t n = 0U; n < body_samples; n++) {
    out->clean[lead_samples + n] = clean_body[n];
    out->mix[lead_samples + n] =
        clean_body[n] + (float)((double)noise_raw[lead_samples + n] * factor);
  }
  for (uint32_t n = 0U; n < lead_samples; n++) {
    out->mix[n] = out->noise[n];
  }

  float peak = 0.0F;
  for (uint32_t n = 0U; n < total; n++) {
    if (fabsf(out->mix[n]) > peak) {
      peak = fabsf(out->mix[n]);
    }
  }
  if (peak > 0.99F) {
    const float scale = 0.99F / peak;
    for (uint32_t n = 0U; n < total; n++) {
      out->clean[n] *= scale;
      out->noise[n] *= scale;
      out->mix[n] *= scale;
    }
  }
}

// Frame-wise analysis STFT (Hann window, ANALYSIS_FFT / ANALYSIS_HOP),
// computed with the library's own real FFT.
typedef struct Analyzer {
  FftTransform* fft;
  double window[ANALYSIS_FFT];
  uint32_t num_frames;
  float* power; // [frame][bin]
} Analyzer;

static void analyzer_init(Analyzer* self, uint32_t num_frames) {
  for (uint32_t n = 0U; n < ANALYSIS_FFT; n++) {
    self->window[n] =
        0.5 - 0.5 * cos(2.0 * M_PI * (double)n / (double)ANALYSIS_FFT);
  }
  self->fft = fft_transform_initialize_bins(ANALYSIS_FFT);
  TEST_ASSERT(self->fft != NULL, "analyzer fft alloc");
  self->num_frames = num_frames;
  self->power =
      (float*)calloc((size_t)num_frames * (ANALYSIS_FFT / 2U), sizeof(float));
  TEST_ASSERT(self->power != NULL, "analyzer power alloc");
}

static void analyzer_free(Analyzer* self) {
  fft_transform_free(self->fft);
  free(self->power);
}

static void analyzer_run(Analyzer* self, const float* signal, int offset,
                         uint32_t total_samples) {
  const uint32_t bins = ANALYSIS_FFT / 2U;
  float* input = get_fft_input_buffer(self->fft);
  for (uint32_t f = 0U; f < self->num_frames; f++) {
    int start = (int)(f * ANALYSIS_HOP) + offset;
    if (start < 0 || (uint32_t)start + ANALYSIS_FFT > total_samples) {
      continue;
    }
    const float* frame = signal + start;
    for (uint32_t n = 0U; n < ANALYSIS_FFT; n++) {
      input[n] = (float)((double)frame[n] * self->window[n]);
    }
    TEST_ASSERT(compute_forward_fft(self->fft), "analysis fft");
    // Halfcomplex layout: [re0, re1..re_{N/2-1}, re_{N/2}, im_{N/2-1}..im1]
    const float* spec = get_fft_output_buffer(self->fft);
    float* out = &self->power[f * bins];
    out[0] = spec[0] * spec[0];
    out[bins - 1U] = spec[bins - 1U] * spec[bins - 1U];
    for (uint32_t k = 1U; k + 1U < bins; k++) {
      const float re = spec[k];
      const float im = spec[ANALYSIS_FFT - k];
      out[k] = re * re + im * im;
    }
  }
}

// Exact mirror of the plugin's buildEngineParams at default slider values:
// suppression 1.0, aggressiveness 0.5, masking 100% (cubic -> 1.0), whitening
// 0, transient protection off, tonal reduction linked to the master reduction,
// DFTT strength coupled to reduction depth (Refinement mode).
// Diagnostic override: SPECBLEACH_REALWORLD_REDUCTION_DB (default 12).
static float plugin_reduction_db(void) {
  const char* env = getenv("SPECBLEACH_REALWORLD_REDUCTION_DB");
  if (env != NULL && env[0] != '\0') {
    float v = (float)atof(env);
    if (v > 0.0F && v <= 40.0F) {
      return v;
    }
  }
  return PLUGIN_REDUCTION_DB;
}

static SpecbleachDenoiserParameters make_parameters(float smoothing,
                                                    uint32_t mode) {
  const float reduction_db = plugin_reduction_db();
  const float reduction_gain = powf(10.0F, -reduction_db / 20.0F);
  return (SpecbleachDenoiserParameters){
      .residual_listen = false,
      .aggressiveness = 0.5F,
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .reduction_gain = reduction_gain,
      .smoothing_factor = smoothing,
      .smoothing_mode = (SpecbleachSmoothingMode)mode,
      .whitening_factor = 0.0F,
      .adaptive_noise = false,
      .noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_MARTIN,
      .masking_depth = 1.0F,
      .suppression_strength = 1.0F,
      .tonal_reduction_gain = reduction_gain,
      .dftt_strength =
          mode == 2U ? 1.0F + reduction_db / DFTT_STRENGTH_DOUBLING_DB : 1.0F,
      .transient_protection_enable = false};
}

// Measures the true stream delay empirically: noise burst on a scratch
// instance, then cross-correlate output against input (the reported latency
// is host-compensation info). Peak-picking an impulse response lands on the
// first OLA lobe (one block); correlation finds the energy centroid the dry
// path must align to (two blocks here).
static uint32_t measure_stream_delay(uint32_t sample_rate,
                                     uint32_t smoothing_mode) {
  specbleach_denoiser* probe =
      specbleach_denoiser_initialize(sample_rate, FRAME_MS, 0u);
  TEST_ASSERT(probe != NULL, "probe init");
  SpecbleachDenoiserParameters p = make_parameters(0.0F, smoothing_mode);
  TEST_ASSERT(specbleach_denoiser_load_parameters(probe, &p, sizeof(p)),
              "probe params");
  const uint32_t half = sample_rate / 2U;
  float* pin = (float*)calloc((size_t)2U * half, sizeof(float));
  float* pout = (float*)calloc((size_t)2U * half, sizeof(float));
  TEST_ASSERT(pin != NULL && pout != NULL, "probe buffers");
  uint32_t rng = 0x12345678U; // deterministic burst (no rand())
  for (uint32_t n = 0U; n < 2U * half; n++) {
    rng ^= rng << 13U;
    rng ^= rng >> 17U;
    rng ^= rng << 5U;
    pin[n] = (float)(int32_t)rng / 2147483648.0F * 0.5F;
  }
  specbleach_denoiser_process(probe, half, pin, pout);
  p.learn_noise = SPECBLEACH_LEARN_OFF;
  specbleach_denoiser_load_parameters(probe, &p, sizeof(p));
  specbleach_denoiser_process(probe, half, pin + half, pout + half);
  const uint32_t block =
      (uint32_t)((float)sample_rate * FRAME_MS / 1000.0F + 0.5F);
  const uint32_t max_lag = 4U * block;
  // Skip the learn->off switch transient at `half` when correlating.
  const uint32_t corr_first = half + 2U * block;
  double best = 0.0;
  uint32_t delay = 0U;
  for (uint32_t lag = 0U; lag <= max_lag; lag++) {
    double acc = 0.0;
    for (uint32_t n = corr_first; n < 2U * half; n++) {
      acc += (double)pout[n] * (double)pin[n - lag];
    }
    if (acc > best) {
      best = acc;
      delay = lag;
    }
  }
  free(pin);
  free(pout);
  specbleach_denoiser_free(probe);
  TEST_ASSERT(best > 0.0 && delay > 0U && delay < max_lag,
              "correlation peak strictly inside lag range");
  const char* dbg = getenv("SPECBLEACH_REALWORLD_DEBUG");
  if (dbg != NULL && dbg[0] != '0' && dbg[0] != '\0') {
    printf("  [probe] correlation delay: %u samples (block=%u)\n", delay,
           block);
  }
  return delay;
}

// Denoises the mix: learns on the lead-in, finalizes, then reduces.
// Returns the allocated output; the caller frees it.
static float* process_mix(const Signals* signals, uint32_t sample_rate,
                          const SmoothingConfig* config, uint32_t mode) {
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(sample_rate, FRAME_MS, 0u);
  TEST_ASSERT(handle != NULL, "denoiser initialize");
  SpecbleachDenoiserParameters p = make_parameters(config->smoothing, mode);
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &p, sizeof(p)),
              "parameters load");

  // Diagnostic override: cap the learn segment (e.g. to mimic short learns)
  uint32_t learn_samples = signals->lead_samples;
  const char* learn_env = getenv("SPECBLEACH_REALWORLD_LEARN_SAMPLES");
  if (learn_env != NULL) {
    uint32_t capped = (uint32_t)strtoul(learn_env, NULL, 10);
    if (capped > 0U && capped < learn_samples) {
      learn_samples = capped;
    }
  }

  const uint32_t total = signals->total_samples;
  uint32_t block = (uint32_t)((float)sample_rate * FRAME_MS / 1000.0F + 0.5F);
  // Diagnostic override: block size (e.g. 512 to mimic the demo CLI)
  const char* block_env = getenv("SPECBLEACH_REALWORLD_BLOCK");
  if (block_env != NULL) {
    uint32_t custom = (uint32_t)strtoul(block_env, NULL, 10);
    if (custom > 0U) {
      block = custom;
    }
  }
  float* out = (float*)calloc(total, sizeof(float));
  TEST_ASSERT(out != NULL, "output alloc");

  uint32_t pos = 0U;
  while (pos < learn_samples) {
    uint32_t chunk = block;
    if (pos + chunk > learn_samples) {
      chunk = learn_samples - pos;
    }
    TEST_ASSERT(specbleach_denoiser_process(handle, chunk, signals->mix + pos,
                                            out + pos),
                "learn process");
    pos += chunk;
  }
  // Finalize the profile (capture modes only finalize on learn -> off)
  TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(handle, 1),
              "noise profile was not learned");
  p.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(specbleach_denoiser_load_parameters(handle, &p, sizeof(p)),
              "parameters finalize");

  while (pos < total) {
    uint32_t chunk = block;
    if (pos + chunk > total) {
      chunk = total - pos;
    }
    TEST_ASSERT(specbleach_denoiser_process(handle, chunk, signals->mix + pos,
                                            out + pos),
                "reduce process");
    pos += chunk;
  }

  // Capture modes 3-4 (STD_DEV/CV_MASK) finalize on the learn -> off
  // transition; the tonal reducer is load-bearing on CV_MASK availability.
  TEST_ASSERT(specbleach_denoiser_noise_profile_available_for_mode(handle, 4),
              "CV_MASK profile was not finalized");

  specbleach_denoiser_free(handle);
  return out;
}

// Classifies body frames as active (signal dominates) or noise-only from the
// clean reference power. Excludes transition zones and the release tail that
// follows active stretches (the denoiser keeps decaying signal energy there,
// which must not count as steady-state noise behavior).
static void frame_masks(const Analyzer* clean_an, float active_threshold,
                        float noise_threshold, uint32_t body_start_frame,
                        uint32_t sample_rate, uint8_t* active,
                        uint8_t* noise_only) {
  const uint32_t bins = ANALYSIS_FFT / 2U;
  const uint32_t tail_frames =
      (uint32_t)(RELEASE_TAIL_EXCLUSION_SEC * (float)sample_rate /
                 (float)ANALYSIS_HOP);
  const uint32_t warmup_frames =
      (uint32_t)(INITIAL_REDUCE_EXCLUSION_SEC * (float)sample_rate /
                 (float)ANALYSIS_HOP);
  const uint32_t noise_start =
      body_start_frame + warmup_frames > body_start_frame
          ? body_start_frame + warmup_frames
          : body_start_frame;
  uint32_t last_active = 0U;
  uint8_t has_active = 0U;
  for (uint32_t f = 0U; f < clean_an->num_frames; f++) {
    float p = 0.0F;
    for (uint32_t k = 1U; k < bins; k++) {
      p += clean_an->power[f * bins + k];
    }
    if (p > active_threshold) {
      active[f] = 1U;
      noise_only[f] = 0U;
      last_active = f;
      has_active = 1U;
      continue;
    }
    active[f] = 0U;
    // Noise-only: quiet enough, in the body, not inside the release tail of
    // a previous active stretch, and past the post-finalization warm-up.
    noise_only[f] = (p < noise_threshold && f >= noise_start &&
                     (!has_active || f - last_active > tail_frames));
  }
}

static void remove_boundary_frames(uint8_t* mask, const uint8_t* other,
                                   uint32_t num_frames) {
  for (uint32_t f = 0U; f < num_frames; f++) {
    if (!mask[f]) {
      continue;
    }
    for (uint32_t e = 1U; e <= BOUNDARY_EXCLUSION_FRAMES; e++) {
      if ((f >= e && other[f - e]) || (f + e < num_frames && other[f + e])) {
        mask[f] = 0U;
        break;
      }
    }
  }
}

// MNI v2: mean temporal CV of 3-bin-averaged, per-frame energy-normalized
// spectral shape over noise-only frames.
static float shape_cv(const Analyzer* an, const uint8_t* noise_only) {
  const uint32_t bins = ANALYSIS_FFT / 2U;
  const uint32_t group_size = 3U;
  const uint32_t num_groups = (bins - 1U) / group_size;

  double* mean = (double*)calloc(num_groups, sizeof(double));
  double* m2 = (double*)calloc(num_groups, sizeof(double));
  TEST_ASSERT(mean != NULL && m2 != NULL, "mni alloc");

  uint32_t count = 0U;
  for (uint32_t f = 0U; f < an->num_frames; f++) {
    if (!noise_only[f]) {
      continue;
    }
    float total = 0.0F;
    for (uint32_t g = 0U; g < num_groups; g++) {
      float group_power = 0.0F;
      for (uint32_t j = 0U; j < group_size; j++) {
        group_power += an->power[f * bins + 1U + g * group_size + j];
      }
      total += group_power;
    }
    if (total <= 0.0F) {
      continue;
    }
    count++;
    for (uint32_t g = 0U; g < num_groups; g++) {
      float group_power = 0.0F;
      for (uint32_t j = 0U; j < group_size; j++) {
        group_power += an->power[f * bins + 1U + g * group_size + j];
      }
      const double shape = (double)group_power / (double)total;
      const double delta = shape - mean[g];
      mean[g] += delta / (double)count;
      m2[g] += delta * (shape - mean[g]);
    }
  }

  double cv_sum = 0.0;
  uint32_t cv_count = 0U;
  for (uint32_t g = 0U; g < num_groups; g++) {
    if (mean[g] > 0.0) {
      cv_sum += sqrt(m2[g] / (double)count) / mean[g];
      cv_count++;
    }
  }
  free(mean);
  free(m2);
  TEST_ASSERT(count >= MIN_METRIC_FRAMES && cv_count > 0U,
              "enough noise-only frames for MNI");
  return (float)(cv_sum / (double)cv_count);
}

// SD split diagnostic: overall + LF/MF/HF band medians + transient/steady
// medians over signal-active frames. Band medians reuse the frame's global
// peak floor (signal-dominant bins only); transient = clean energy onset.
// Subsets below MIN_METRIC_FRAMES fall back to the overall median.
typedef struct SdSplit {
  float lf_db;
  float mf_db;
  float hf_db;
  float tr_db;
  float st_db;
} SdSplit;

static float median_in_place(float* values, uint32_t count) {
  // median via insertion sort (count is small)
  for (uint32_t i = 1U; i < count; i++) {
    float key = values[i];
    int j = (int)i - 1;
    while (j >= 0 && values[j] > key) {
      values[j + 1] = values[j];
      j--;
    }
    values[j + 1] = key;
  }
  return values[count / 2U];
}

// SD: median per-frame log-spectral distance over signal-active frames,
// restricted to bins within PEAK_DYNAMIC_DB of the frame's clean peak.
static float spectral_distortion(const Analyzer* clean_an,
                                 const Analyzer* out_an, const uint8_t* active,
                                 uint32_t sample_rate, SdSplit* split) {
  const uint32_t bins = ANALYSIS_FFT / 2U;
  const float bin_floor_ratio = powf(10.0F, -PEAK_DYNAMIC_DB / 10.0F);
  const float clip = DEV_CLIP_DB * DEV_CLIP_DB;
  const float bin_hz = (float)sample_rate / (float)ANALYSIS_FFT;
  const uint32_t mf_bin = (uint32_t)(SD_BAND_MF_HZ / bin_hz);
  const uint32_t hf_bin = (uint32_t)(SD_BAND_HF_HZ / bin_hz);
  const uint32_t nf = clean_an->num_frames;

  // Transient mask: clean frame-energy onsets (+ tail), intersected with
  // active at accumulation time.
  uint8_t* transient = (uint8_t*)calloc(nf, sizeof(uint8_t));
  TEST_ASSERT(transient != NULL, "transient alloc");
  {
    float prev_p = 0.0F;
    for (uint32_t f = 0U; f < nf; f++) {
      float p = 0.0F;
      for (uint32_t k = 1U; k < bins; k++) {
        p += clean_an->power[f * bins + k];
      }
      const float onset_db = 10.0F * log10f((p + 1e-12F) / (prev_p + 1e-12F));
      prev_p = p;
      if (onset_db > SD_TRANSIENT_ATTACK_DB) {
        const uint32_t last = f + SD_TRANSIENT_TAIL_FRAMES < nf
                                  ? f + SD_TRANSIENT_TAIL_FRAMES
                                  : nf - 1U;
        for (uint32_t t = f; t <= last; t++) {
          transient[t] = 1U;
        }
      }
    }
  }

  float* frame_sd = (float*)calloc(nf, sizeof(float));
  float* frame_lf = (float*)calloc(nf, sizeof(float));
  float* frame_mf = (float*)calloc(nf, sizeof(float));
  float* frame_hf = (float*)calloc(nf, sizeof(float));
  float* frame_tr = (float*)calloc(nf, sizeof(float));
  float* frame_st = (float*)calloc(nf, sizeof(float));
  TEST_ASSERT(frame_sd != NULL && frame_lf != NULL && frame_mf != NULL &&
                  frame_hf != NULL && frame_tr != NULL && frame_st != NULL,
              "sd alloc");
  uint32_t count = 0U;
  uint32_t count_lf = 0U;
  uint32_t count_mf = 0U;
  uint32_t count_hf = 0U;
  uint32_t count_tr = 0U;
  uint32_t count_st = 0U;
  for (uint32_t f = 0U; f < nf; f++) {
    if (!active[f]) {
      continue;
    }
    float peak = 0.0F;
    for (uint32_t k = 1U; k < bins; k++) {
      if (clean_an->power[f * bins + k] > peak) {
        peak = clean_an->power[f * bins + k];
      }
    }
    if (peak <= 0.0F) {
      continue;
    }
    const float floor_power = peak * bin_floor_ratio;
    double sum = 0.0;
    double sum_lf = 0.0;
    double sum_mf = 0.0;
    double sum_hf = 0.0;
    uint32_t used = 0U;
    uint32_t used_lf = 0U;
    uint32_t used_mf = 0U;
    uint32_t used_hf = 0U;
    for (uint32_t k = 1U; k < bins; k++) {
      float cp = clean_an->power[f * bins + k];
      float op = out_an->power[f * bins + k];
      if (cp < floor_power) {
        continue;
      }
      double dev = 10.0 * log10(((double)op + 1e-20) / ((double)cp + 1e-20));
      if (dev > DEV_CLIP_DB) {
        dev = DEV_CLIP_DB;
      }
      if (dev < -DEV_CLIP_DB) {
        dev = -DEV_CLIP_DB;
      }
      sum += dev * dev;
      used++;
      if (k < mf_bin) {
        sum_lf += dev * dev;
        used_lf++;
      } else if (k < hf_bin) {
        sum_mf += dev * dev;
        used_mf++;
      } else {
        sum_hf += dev * dev;
        used_hf++;
      }
    }
    if (used > 0U) {
      const float frame_value = (float)sqrt(sum / (double)used);
      frame_sd[count++] = frame_value;
      if (transient[f]) {
        frame_tr[count_tr++] = frame_value;
      } else {
        frame_st[count_st++] = frame_value;
      }
    }
    if (used_lf > 0U) {
      frame_lf[count_lf++] = (float)sqrt(sum_lf / (double)used_lf);
    }
    if (used_mf > 0U) {
      frame_mf[count_mf++] = (float)sqrt(sum_mf / (double)used_mf);
    }
    if (used_hf > 0U) {
      frame_hf[count_hf++] = (float)sqrt(sum_hf / (double)used_hf);
    }
  }
  TEST_ASSERT(count >= MIN_METRIC_FRAMES, "enough active frames for SD");
  const float overall = median_in_place(frame_sd, count);
  split->lf_db = count_lf >= MIN_METRIC_FRAMES
                     ? median_in_place(frame_lf, count_lf)
                     : overall;
  split->mf_db = count_mf >= MIN_METRIC_FRAMES
                     ? median_in_place(frame_mf, count_mf)
                     : overall;
  split->hf_db = count_hf >= MIN_METRIC_FRAMES
                     ? median_in_place(frame_hf, count_hf)
                     : overall;
  split->tr_db = count_tr >= MIN_METRIC_FRAMES
                     ? median_in_place(frame_tr, count_tr)
                     : overall;
  split->st_db = count_st >= MIN_METRIC_FRAMES
                     ? median_in_place(frame_st, count_st)
                     : overall;
  free(transient);
  free(frame_sd);
  free(frame_lf);
  free(frame_mf);
  free(frame_hf);
  free(frame_tr);
  free(frame_st);
  return overall;
}

// Decay deficit: at clean-envelope drop points (phrase/word endings well
// above the floor), how much energy the processed signal carries above the
// simultaneous unprocessed noise level (dB, mean over drop points).
// Self-normalizing: 0 = decays at least as fast as the noise floor; > 0 =
// temporal smear echoing the dropped signal.
static float decay_deficit(const float* clean, const float* out,
                           const float* mix, uint32_t total_samples,
                           uint32_t lead_samples, uint32_t delay_samples) {
  const uint32_t num_env =
      (total_samples - delay_samples - ENV_WINDOW) / ENV_HOP + 1U;
  float* clean_db = (float*)malloc((size_t)num_env * sizeof(float));
  float* out_db = (float*)malloc((size_t)num_env * sizeof(float));
  float* mix_db = (float*)malloc((size_t)num_env * sizeof(float));
  TEST_ASSERT(clean_db != NULL && out_db != NULL && mix_db != NULL,
              "envelope alloc");

  for (uint32_t i = 0U; i < num_env; i++) {
    const uint32_t start = i * ENV_HOP;
    const uint32_t out_start = start + delay_samples;
    double cs = 0.0;
    double os = 0.0;
    double ms = 0.0;
    for (uint32_t n = 0U; n < ENV_WINDOW; n++) {
      const double c = (double)clean[start + n];
      cs += c * c;
      const double o = (double)out[out_start + n];
      os += o * o;
      const double m = (double)mix[start + n];
      ms += m * m;
    }
    clean_db[i] = (float)(10.0 * log10(cs + 1e-20));
    out_db[i] = (float)(10.0 * log10(os + 1e-20));
    mix_db[i] = (float)(10.0 * log10(ms + 1e-20));
  }

  float mix_level = -1000.0F;
  for (uint32_t i = 0U; i < num_env; i++) {
    if (mix_db[i] > mix_level) {
      mix_level = mix_db[i];
    }
  }
  const float floor_db = mix_level - 25.0F;

  double error_sum = 0.0;
  uint32_t count = 0U;
  for (uint32_t i = lead_samples / ENV_HOP; i + 2U < num_env; i++) {
    if (clean_db[i] < floor_db) {
      continue; // drop point must start above the noise floor
    }
    const float clean_drop = clean_db[i] - clean_db[i + 2U];
    if (clean_drop < ECHO_DROP_DB) {
      continue;
    }
    const float excess = out_db[i + 2U] - mix_db[i + 2U];
    if (excess > 0.0F) {
      error_sum += (double)excess;
      count++;
    }
  }
  free(clean_db);
  free(out_db);
  free(mix_db);
  if (count == 0U) {
    return 0.0F;
  }
  return (float)(error_sum / (double)count);
}

static const char* mode_name(uint32_t mode) {
  static const char* names[4] = {"Temporal (1D)", "NLM 2D", "NLM 2D + DFTT",
                                 "BM3D"};
  return names[mode];
}

// Diagnostic: dump mix/out wavs when SPECBLEACH_REALWORLD_DUMP_DIR is set
static void dump_wav(const char* dir, const char* name, const float* data,
                     uint32_t length, uint32_t sample_rate) {
  char path[512];
  snprintf(path, sizeof(path), "%s/%s.wav", dir, name);
  SF_INFO info;
  memset(&info, 0, sizeof(info));
  info.samplerate = (int)sample_rate;
  info.channels = 1;
  info.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
  SNDFILE* file = sf_open(path, SFM_WRITE, &info);
  TEST_ASSERT(file != NULL, "dump open");
  TEST_ASSERT(
      sf_writef_float(file, data, (sf_count_t)length) == (sf_count_t)length,
      "dump write");
  sf_close(file);
}

// Full pipeline for one (case, smoothing config, mode): mix -> learn ->
// reduce -> analyze -> probes -> metrics.
// SI-SDR (Le Roux et al.): score est[i] against ref[i - delay] over
// [lo, hi), both zero-meaned. Higher = cleaner overall.
static float si_sdr_db(const float* est, const float* ref, uint32_t lo,
                       uint32_t hi, uint32_t delay) {
  double sum_e = 0.0, sum_r = 0.0;
  uint32_t n = 0U;
  for (uint32_t i = lo; i < hi; i++) {
    sum_e += (double)est[i];
    sum_r += (double)ref[i - delay];
    n++;
  }
  if (n == 0U) {
    return -60.0F;
  }
  const double mean_e = sum_e / (double)n;
  const double mean_r = sum_r / (double)n;
  double cc = 0.0, ec = 0.0, ee = 0.0;
  for (uint32_t i = lo; i < hi; i++) {
    const double e = (double)est[i] - mean_e;
    const double c = (double)ref[i - delay] - mean_r;
    cc += c * c;
    ec += e * c;
    ee += e * e;
  }
  if (!(cc > 0.0)) {
    return -60.0F;
  }
  const double denom = cc * ee - ec * ec;
  if (!(denom > 0.0)) {
    return 100.0F;
  }
  double v = 10.0 * log10(ec * ec / denom);
  if (v > 100.0) {
    v = 100.0;
  } else if (v < -60.0) {
    v = -60.0;
  }
  return (float)v;
}

// SAR (BSS_eval style): project est[i] onto the span of {clean, noise} at
// [i - delay]; the residual is processing artifact (musical noise). Higher =
// less artifact. Valid because mix == clean + noise exactly here.
static float sar_db(const float* est, const float* clean, const float* noise,
                    uint32_t lo, uint32_t hi, uint32_t delay) {
  double cc = 0.0, nn = 0.0, cn = 0.0, ec = 0.0, en = 0.0, ee = 0.0;
  for (uint32_t i = lo; i < hi; i++) {
    const double e = (double)est[i];
    const double c = (double)clean[i - delay];
    const double v = (double)noise[i - delay];
    cc += c * c;
    nn += v * v;
    cn += c * v;
    ec += e * c;
    en += e * v;
    ee += e * e;
  }
  const double det = cc * nn - cn * cn;
  if (!(det > 1e-12 * cc * nn)) {
    return -60.0F; // collinear or silent references
  }
  const double wc = (ec * nn - en * cn) / det;
  const double wn = (en * cc - ec * cn) / det;
  const double allowed = wc * ec + wn * en; // ||projection||^2
  const double artif = ee - allowed;
  if (!(artif > 1e-12 * ee)) {
    return 100.0F; // (near) bit-exact projection, no artifacts
  }
  if (!(allowed > 0.0)) {
    return -60.0F;
  }
  double val = 10.0 * log10(allowed / artif);
  if (val > 100.0) {
    val = 100.0;
  } else if (val < -60.0) {
    val = -60.0;
  }
  return (float)val;
}

static Metrics run_case_mode(const Signals* signals, uint32_t sample_rate,
                             const SmoothingConfig* config, uint32_t mode) {
  // Determinism + immutability probes
  float* out = process_mix(signals, sample_rate, config, mode);
  float* out2 = process_mix(signals, sample_rate, config, mode);
  const uint32_t total = signals->total_samples;
  uint32_t diffs = 0U;
  double max_diff = 0.0;
  for (uint32_t n = 0U; n < total; n++) {
    if (out[n] != out2[n]) {
      diffs++;
      double d = fabs((double)out[n] - (double)out2[n]);
      if (d > max_diff) {
        max_diff = d;
      }
    }
  }
  free(out2);
  printf("  determinism: %u differing samples, max |diff| %.3g\n", diffs,
         max_diff);
  TEST_ASSERT(diffs == 0U,
              "determinism: second instance must be bit-identical");

  uint32_t modified = 0U;
  for (uint32_t n = 0U; n < total; n++) {
    if (signals->mix[n] != signals->snapshot[n]) {
      modified++;
    }
  }
  TEST_ASSERT(modified == 0U, "input immutability");

  const char* dump_dir = getenv("SPECBLEACH_REALWORLD_DUMP_DIR");
  if (dump_dir != NULL && dump_dir[0] != '\0') {
    char name[128];
    snprintf(name, sizeof(name), "%s_%s_mode%u_out", config->name,
             signals->case_name, mode);
    dump_wav(dump_dir, name, out, total, sample_rate);
    snprintf(name, sizeof(name), "%s_mix", signals->case_name);
    dump_wav(dump_dir, name, signals->mix, total, sample_rate);
    snprintf(name, sizeof(name), "%s_clean", signals->case_name);
    dump_wav(dump_dir, name, signals->clean, total, sample_rate);
  }

  const uint32_t stream_delay = measure_stream_delay(sample_rate, mode);

  const uint32_t num_frames = (total - ANALYSIS_FFT) / ANALYSIS_HOP + 1U;
  Analyzer clean_an;
  Analyzer out_an;
  Analyzer mix_an;
  analyzer_init(&clean_an, num_frames);
  analyzer_init(&out_an, num_frames);
  analyzer_init(&mix_an, num_frames);
  analyzer_run(&out_an, out, 0, total);
  analyzer_run(&clean_an, signals->clean, -(int)stream_delay, total);
  analyzer_run(&mix_an, signals->mix, -(int)stream_delay, total);

  // Frame classification against the clean reference
  const uint32_t bins = ANALYSIS_FFT / 2U;
  float clean_peak = 0.0F;
  const uint32_t body_start_frame =
      ((uint32_t)(signals->lead_samples +
                  (uint32_t)(0.05F * (float)sample_rate)) -
       ANALYSIS_FFT / 2U) /
      ANALYSIS_HOP;
  for (uint32_t f = body_start_frame; f < num_frames; f++) {
    float p = 0.0F;
    for (uint32_t k = 1U; k < bins; k++) {
      p += clean_an.power[f * bins + k];
    }
    if (p > clean_peak) {
      clean_peak = p;
    }
  }
  const float active_threshold =
      clean_peak * powf(10.0F, -ACTIVE_DYNAMIC_DB / 10.0F);
  const float noise_threshold =
      clean_peak * powf(10.0F, -NOISE_DYNAMIC_DB / 10.0F);
  uint8_t* active = (uint8_t*)calloc(num_frames, 1);
  uint8_t* noise_only = (uint8_t*)calloc(num_frames, 1);
  TEST_ASSERT(active != NULL && noise_only != NULL, "mask alloc");
  frame_masks(&clean_an, active_threshold, noise_threshold, body_start_frame,
              sample_rate, active, noise_only);
  remove_boundary_frames(active, noise_only, num_frames);
  remove_boundary_frames(noise_only, active, num_frames);

  Metrics m = {0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
               0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F};

  // Learned-noise bin selection: per-bin PSD from lead-in mix frames (the
  // same region the engine learned from). Energy-ratio attenuation over all
  // bins is blind to whitening on colored noise; restricting to the bins
  // where the noise actually lives measures true rejection.
  uint32_t lead_first =
      (uint32_t)((0.2F * (float)sample_rate + (float)stream_delay) /
                 (float)ANALYSIS_HOP);
  uint32_t lead_last =
      (uint32_t)((1.35F * (float)sample_rate + (float)stream_delay -
                  (float)(ANALYSIS_FFT / 2U)) /
                 (float)ANALYSIS_HOP);
  if (lead_last >= num_frames) {
    lead_last = num_frames - 1U;
  }
  TEST_ASSERT(lead_last > lead_first, "lead-in analysis frames exist");
  double* noise_psd = (double*)calloc(bins, sizeof(double));
  TEST_ASSERT(noise_psd != NULL, "noise psd alloc");
  uint32_t lead_frames = 0U;
  for (uint32_t f = lead_first; f <= lead_last; f++) {
    lead_frames++;
    for (uint32_t k = 0U; k < bins; k++) {
      noise_psd[k] += (double)mix_an.power[f * bins + k];
    }
  }
  for (uint32_t k = 0U; k < bins; k++) {
    noise_psd[k] /= (double)lead_frames;
  }
  double psd_peak = 0.0;
  for (uint32_t k = 1U; k < bins; k++) {
    if (noise_psd[k] > psd_peak) {
      psd_peak = noise_psd[k];
    }
  }
  const double psd_floor =
      psd_peak * pow(10.0, -(double)NOISE_BIN_FLOOR_DB / 10.0);
  uint8_t* noise_bins = (uint8_t*)calloc(bins, 1);
  TEST_ASSERT(noise_bins != NULL, "noise bin mask alloc");
  uint32_t noise_bin_count = 0U;
  for (uint32_t k = 1U; k < bins; k++) {
    if (noise_psd[k] > psd_floor) {
      noise_bins[k] = 1U;
      noise_bin_count++;
    }
  }
  TEST_ASSERT(noise_bin_count >= MIN_NOISE_BINS, "learned noise bins found");

  // Att: attenuation over noise-only frames
  double mix_e = 0.0;
  double out_e = 0.0;
  double mix_e_all = 0.0;
  double out_e_all = 0.0;
  uint32_t att_frames = 0U;
  for (uint32_t f = 0U; f < num_frames; f++) {
    if (!noise_only[f]) {
      continue;
    }
    att_frames++;
    for (uint32_t k = 1U; k < bins; k++) {
      mix_e_all += (double)mix_an.power[f * bins + k];
      out_e_all += (double)out_an.power[f * bins + k];
      if (noise_bins[k]) {
        mix_e += (double)mix_an.power[f * bins + k];
        out_e += (double)out_an.power[f * bins + k];
      }
    }
  }
  TEST_ASSERT(att_frames >= MIN_METRIC_FRAMES, "enough noise-only frames");
  m.att_db = (float)(10.0 * log10(mix_e / out_e));
  m.att_all_db = (float)(10.0 * log10(mix_e_all / out_e_all));
  free(noise_psd);
  free(noise_bins);

  // SD over active frames (+ band/transient split diagnostic)
  SdSplit sd_split = {0.0F, 0.0F, 0.0F, 0.0F, 0.0F};
  m.sd_db =
      spectral_distortion(&clean_an, &out_an, active, sample_rate, &sd_split);
  m.sd_lf_db = sd_split.lf_db;
  m.sd_mf_db = sd_split.mf_db;
  m.sd_hf_db = sd_split.hf_db;
  m.sd_tr_db = sd_split.tr_db;
  m.sd_st_db = sd_split.st_db;

  // MNI v2
  m.mni_out = shape_cv(&out_an, noise_only);
  const float mni_mix = shape_cv(&mix_an, noise_only);
  m.mni_ratio = m.mni_out / mni_mix;

  // Decay deficit
  m.decay_db = decay_deficit(signals->clean, out, signals->mix, total,
                             signals->lead_samples, stream_delay);

  // SI-SDR / SAR over the delay-compensated body (report-only, no gates).
  // The 1.5 s lead-in settles all filter fill, so the body window is clean
  // for every mode; out[i] aligns with clean/noise at [i - stream_delay].
  {
    const uint32_t lo = signals->lead_samples + stream_delay;
    TEST_ASSERT(lo < total, "scored body window exists");
    m.sisdr_db = si_sdr_db(out, signals->clean, lo, total, stream_delay);
    m.sar_db =
        sar_db(out, signals->clean, signals->noise, lo, total, stream_delay);
  }

  analyzer_free(&clean_an);
  analyzer_free(&out_an);
  analyzer_free(&mix_an);
  free(active);
  free(noise_only);
  free(out);
  return m;
}

int main(void) {
  setvbuf(stdout, NULL, _IONBF, 0); // live output; asserts abort otherwise

  int full_mode = 0;
  const char* full_env = getenv("SPECBLEACH_REALWORLD_FULL");
  if (full_env != NULL && full_env[0] != '\0' && full_env[0] != '0') {
    full_mode = 1;
  }
  int dense_mode = 0;
  const char* dense_env = getenv("SPECBLEACH_REALWORLD_DENSE");
  if (dense_env != NULL && dense_env[0] != '\0' && dense_env[0] != '0') {
    dense_mode = 1;
  }
  // CI default: plugin-default smoothing (0.0); FULL sweeps 3 values; DENSE
  // sweeps 0.0..1.0 in 0.1 steps for knob-response curves (audit only)
#define SMOOTHING_DENSE_STEPS 11U
  const uint32_t num_configs =
      dense_mode ? SMOOTHING_DENSE_STEPS : (full_mode ? 3U : 1U);

  uint32_t sample_rate = 0U;
  uint32_t body_samples = 0U;

  // Load fixtures (all cases must share the sample rate)
  float* clean_body[NUM_CASES];
  float* noise_raw[NUM_CASES];
  Signals signals[NUM_CASES];
  for (uint32_t c = 0U; c < NUM_CASES; c++) {
    uint32_t sr = 0U;
    uint32_t len = 0U;
    clean_body[c] = load_wav(
        fixture_file(CASES[c].clean_source != NULL ? CASES[c].clean_source
                                                   : CASES[c].name,
                     "clean"),
        &sr, &len);
    if (c == 0U) {
      sample_rate = sr;
      body_samples = len;
    }
    TEST_ASSERT(sr == sample_rate && len == body_samples,
                "fixture sample rates and lengths must match");
    if (CASES[c].white_noise) {
      noise_raw[c] =
          make_white_noise((uint32_t)(NOISE_LEAD_IN_SEC * (float)sr) + len);
      continue;
    }
    uint32_t noise_len = 0U;
    noise_raw[c] =
        load_wav(fixture_file(CASES[c].name, "noise"), &sr, &noise_len);
    TEST_ASSERT(sr == sample_rate, "noise fixture sample rate");
    TEST_ASSERT(
        noise_len >=
            (uint32_t)(NOISE_LEAD_IN_SEC * (float)sample_rate) + body_samples,
        "noise fixture must cover lead-in + body");
  }
  const uint32_t lead_samples =
      (uint32_t)(NOISE_LEAD_IN_SEC * (float)sample_rate);
  const uint32_t total = lead_samples + body_samples;

  printf(
      "== real-world quality metrics (sample rate %u, %u cases, plugin-default "
      "config at %.0f dB reduction, %s smoothing config%s) ==\n",
      sample_rate, NUM_CASES, (double)plugin_reduction_db(),
      full_mode ? "all" : "default", full_mode ? "s" : "");

  for (uint32_t ci = 0U; ci < num_configs; ci++) {
    char dense_name[32];
    SmoothingConfig dense_config;
    const SmoothingConfig* config;
    if (dense_mode) {
      snprintf(dense_name, sizeof(dense_name), "smoothing-%u", ci * 10U);
      dense_config.name = dense_name;
      dense_config.smoothing = (float)ci / 10.0F;
      config = &dense_config;
    } else {
      config = &SMOOTHING_CONFIGS[ci];
    }
    for (uint32_t mode = 0U; mode < NUM_MODES; mode++) {
      printf("== %s | %s ==\n", config->name, mode_name(mode));
      Metrics agg = {0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
                     0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F};
      for (uint32_t c = 0U; c < NUM_CASES; c++) {
        Signals* sig = &signals[c];
        sig->case_name = CASES[c].name;
        sig->total_samples = total;
        sig->lead_samples = lead_samples;
        sig->body_samples = body_samples;
        sig->clean = (float*)calloc(total, sizeof(float));
        sig->noise = (float*)calloc(total, sizeof(float));
        sig->mix = (float*)calloc(total, sizeof(float));
        sig->snapshot = (float*)malloc(total * sizeof(float));
        TEST_ASSERT(sig->clean != NULL && sig->noise != NULL &&
                        sig->mix != NULL && sig->snapshot != NULL,
                    "signal allocs");
        mix_case(clean_body[c], noise_raw[c], body_samples, lead_samples,
                 CASES[c].snr_db, sig);
        memcpy(sig->snapshot, sig->mix, total * sizeof(float));

        Metrics m = run_case_mode(sig, sample_rate, config, mode);
        agg.att_db += m.att_db;
        agg.att_all_db += m.att_all_db;
        agg.sd_db += m.sd_db;
        agg.sd_lf_db += m.sd_lf_db;
        agg.sd_mf_db += m.sd_mf_db;
        agg.sd_hf_db += m.sd_hf_db;
        agg.sd_tr_db += m.sd_tr_db;
        agg.sd_st_db += m.sd_st_db;
        agg.mni_ratio += m.mni_ratio;
        agg.mni_out += m.mni_out;
        agg.decay_db += m.decay_db;
        agg.sisdr_db += m.sisdr_db;
        agg.sar_db += m.sar_db;
        printf(
            "  %s: att %.1f dB (all %.1f), sd %.1f dB, mni %.2f (raw "
            "%.2f), decay %.1f dB, sisdr %.1f dB, sar %.1f dB\n",
            CASES[c].name, m.att_db, m.att_all_db, m.sd_db, m.mni_ratio,
            m.mni_out, m.decay_db, m.sisdr_db, m.sar_db);
        printf(
            "    sd split: LF %.1f dB, MF %.1f dB, HF %.1f dB | tr %.1f dB, "
            "st %.1f dB\n",
            m.sd_lf_db, m.sd_mf_db, m.sd_hf_db, m.sd_tr_db, m.sd_st_db);
        free(sig->clean);
        free(sig->noise);
        free(sig->mix);
        free(sig->snapshot);
      }
      agg.att_db /= (float)NUM_CASES;
      agg.att_all_db /= (float)NUM_CASES;
      agg.sd_db /= (float)NUM_CASES;
      agg.sd_lf_db /= (float)NUM_CASES;
      agg.sd_mf_db /= (float)NUM_CASES;
      agg.sd_hf_db /= (float)NUM_CASES;
      agg.sd_tr_db /= (float)NUM_CASES;
      agg.sd_st_db /= (float)NUM_CASES;
      agg.mni_ratio /= (float)NUM_CASES;
      agg.mni_out /= (float)NUM_CASES;
      agg.decay_db /= (float)NUM_CASES;
      agg.sisdr_db /= (float)NUM_CASES;
      agg.sar_db /= (float)NUM_CASES;
      const float sd_per_att =
          (agg.att_db > 0.01F) ? agg.sd_db / agg.att_db : 0.0F;
      printf(
          "  AGGREGATE: att %.1f dB (all %.1f), sd %.1f dB, mni %.2f, decay "
          "%.1f dB, sisdr %.1f dB, sar %.1f dB, sd/att %.2f\n",
          agg.att_db, agg.att_all_db, agg.sd_db, agg.mni_ratio, agg.decay_db,
          agg.sisdr_db, agg.sar_db, sd_per_att);
      printf(
          "    sd split: LF %.1f dB, MF %.1f dB, HF %.1f dB | tr %.1f dB, "
          "st %.1f dB\n",
          agg.sd_lf_db, agg.sd_mf_db, agg.sd_hf_db, agg.sd_tr_db, agg.sd_st_db);

      if (!full_mode && !dense_mode) {
        TEST_ASSERT(agg.att_db > ATT_GATES_DB[mode], "attenuation gate");
        TEST_ASSERT(agg.sd_db < SD_GATES_DB[mode], "spectral distortion gate");
        TEST_ASSERT(agg.mni_ratio < MNI_GATES[mode], "musical noise gate");
        TEST_ASSERT(agg.decay_db < DECAY_GATES_DB[mode], "decay smearing gate");
      }
    }
  }

  for (uint32_t c = 0U; c < NUM_CASES; c++) {
    free(clean_body[c]);
    free(noise_raw[c]);
  }
  printf("\n[OK] real-world quality metrics passed\n");
  return 0;
}
