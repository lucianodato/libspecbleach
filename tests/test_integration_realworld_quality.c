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
// CI runs the moderate preset by default; set SPECBLEACH_REALWORLD_FULL=1 to
// sweep light/moderate/heavy locally.

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

#define FRAME_MS 40.0F // 1920 samples per block at 48 kHz

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

// Decay metric tuning
#define ENV_WINDOW 512U
#define ENV_HOP 256U
#define ECHO_DROP_DB 15.0F // clean decay that qualifies as a phrase/word end

// Presets (mirroring the historical autoresearch evaluator)
typedef struct Preset {
  const char* name;
  float reduction_db;
  float whitening; // 0..1
  float smoothing; // 0..1
} Preset;

static const Preset PRESETS[3] = {
    {"light", 12.0F, 0.10F, 0.05F},
    {"moderate", 18.0F, 0.20F, 0.15F},
    {"heavy", 24.0F, 0.30F, 0.25F},
};

typedef struct Case {
  const char* name;
  float snr_db;
} Case;

// Fixture cases with their designated SNR (matrix: stationary/non-stationary
// noise x speech/music x moderate/hard SNR)
static const Case CASES[4] = {
    {"case_voice_fan", 12.0F},
    {"case_voice_city", 6.0F},
    {"case_guitar_electric", 18.0F},
    {"case_drums_fridge", 6.0F},
};

#define NUM_MODES 3U

// Per-mode regression gates for the CI (moderate) preset, aggregated over
// cases: [Temporal, NLM 2D, NLM 2D + DFTT]
static const float ATT_GATES_DB[3] = {3.0F, 3.0F, 3.0F};   // minimum rejection
static const float SD_GATES_DB[3] = {25.0F, 25.0F, 25.0F}; // maximum damage
static const float MNI_GATES[3] = {2.5F, 2.5F, 2.5F};      // max musical noise
static const float DECAY_GATES_DB[3] = {8.0F, 8.0F, 8.0F}; // max smearing

typedef struct Metrics {
  float att_db;     // attenuation on learned-noise bins over noise-only frames
  float att_all_db; // attenuation over all bins (whitening-sensitive reference)
  float sd_db;      // log-spectral distortion on active frames
  float mni_ratio;  // output shape-CV / unprocessed noise shape-CV
  float mni_out;    // raw output shape-CV (for reference when tuning)
  float decay_db;   // mean post-transient decay deficit
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

// Measures the true stream delay empirically with an impulse on a scratch
// instance (the reported latency is host-compensation info).
static uint32_t measure_stream_delay(uint32_t sample_rate,
                                     uint32_t smoothing_mode) {
  specbleach_denoiser* probe =
      specbleach_denoiser_initialize(sample_rate, FRAME_MS);
  TEST_ASSERT(probe != NULL, "probe init");
  SpecbleachDenoiserParameters p =
      (SpecbleachDenoiserParameters){.aggressiveness = 1.0F,
                                     .learn_noise = SPECBLEACH_LEARN_ALL,
                                     .reduction_gain = 0.1F,
                                     .smoothing_factor = 0.15F,
                                     .smoothing_mode = smoothing_mode,
                                     .masking_depth = 0.5F,
                                     .whitening_factor = 0.2F,
                                     .tonal_reduction_gain = 1.0F,
                                     .hpss_enable = true};
  TEST_ASSERT(specbleach_denoiser_load_parameters(probe, &p, sizeof(p)),
              "probe params");
  const uint32_t half = sample_rate / 2U;
  float* pin = (float*)calloc((size_t)2U * half, sizeof(float));
  float* pout = (float*)calloc((size_t)2U * half, sizeof(float));
  TEST_ASSERT(pin != NULL && pout != NULL, "probe buffers");
  const uint32_t impulse_pos = 4000U;
  pin[impulse_pos] = 1.0F;
  specbleach_denoiser_process(probe, half, pin, pout);
  p.learn_noise = SPECBLEACH_LEARN_OFF;
  specbleach_denoiser_load_parameters(probe, &p, sizeof(p));
  specbleach_denoiser_process(probe, half, pin + half, pout + half);
  uint32_t peak = 0U;
  float peak_v = 0.0F;
  for (uint32_t n = 0U; n < 2U * half; n++) {
    if (fabsf(pout[n]) > peak_v) {
      peak_v = fabsf(pout[n]);
      peak = n;
    }
  }
  free(pin);
  free(pout);
  specbleach_denoiser_free(probe);
  TEST_ASSERT(peak > impulse_pos, "impulse peak after input");
  const char* dbg = getenv("SPECBLEACH_REALWORLD_DEBUG");
  if (dbg != NULL && dbg[0] != '0' && dbg[0] != '\0') {
    printf("  [probe] impulse delay: %u samples (block=%u)\n",
           peak - impulse_pos,
           (uint32_t)((float)sample_rate * FRAME_MS / 1000.0F + 0.5F));
  }
  return peak - impulse_pos;
}

static SpecbleachDenoiserParameters make_parameters(const Preset* preset,
                                                    uint32_t mode) {
  return (SpecbleachDenoiserParameters){
      .residual_listen = false,
      .aggressiveness = 1.0F,
      .learn_noise = SPECBLEACH_LEARN_ALL,
      .reduction_gain = powf(10.0F, -preset->reduction_db / 20.0F),
      .smoothing_factor = preset->smoothing,
      .smoothing_mode = (SpecbleachSmoothingMode)mode,
      .whitening_factor = preset->whitening,
      .masking_depth = 0.5F,
      .tonal_reduction_gain = 1.0F,
      .hpss_enable = true};
}

// Denoises the mix: learns on the lead-in, finalizes, then reduces.
// Returns the allocated output; the caller frees it.
static float* process_mix(const Signals* signals, uint32_t sample_rate,
                          const Preset* preset, uint32_t mode) {
  specbleach_denoiser* handle =
      specbleach_denoiser_initialize(sample_rate, FRAME_MS);
  TEST_ASSERT(handle != NULL, "denoiser initialize");
  SpecbleachDenoiserParameters p = make_parameters(preset, mode);
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

// SD: median per-frame log-spectral distance over signal-active frames,
// restricted to bins within PEAK_DYNAMIC_DB of the frame's clean peak.
static float spectral_distortion(const Analyzer* clean_an,
                                 const Analyzer* out_an,
                                 const uint8_t* active) {
  const uint32_t bins = ANALYSIS_FFT / 2U;
  const float bin_floor_ratio = powf(10.0F, -PEAK_DYNAMIC_DB / 10.0F);
  const float clip = DEV_CLIP_DB * DEV_CLIP_DB;

  float* frame_sd = (float*)calloc(clean_an->num_frames, sizeof(float));
  TEST_ASSERT(frame_sd != NULL, "sd alloc");
  uint32_t count = 0U;
  for (uint32_t f = 0U; f < clean_an->num_frames; f++) {
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
    uint32_t used = 0U;
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
    }
    if (used > 0U) {
      frame_sd[count++] = (float)sqrt(sum / (double)used);
    }
  }
  TEST_ASSERT(count >= MIN_METRIC_FRAMES, "enough active frames for SD");
  // median via insertion sort (count is small)
  for (uint32_t i = 1U; i < count; i++) {
    float key = frame_sd[i];
    int j = (int)i - 1;
    while (j >= 0 && frame_sd[j] > key) {
      frame_sd[j + 1] = frame_sd[j];
      j--;
    }
    frame_sd[j + 1] = key;
  }
  float median = frame_sd[count / 2U];
  free(frame_sd);
  return median;
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
  static const char* names[3] = {"Temporal (1D)", "NLM 2D", "NLM 2D + DFTT"};
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

// Full pipeline for one (case, preset, mode): mix -> learn -> reduce ->
// analyze -> probes -> metrics.
static Metrics run_case_mode(const Signals* signals, uint32_t sample_rate,
                             const Preset* preset, uint32_t mode) {
  // Determinism + immutability probes
  float* out = process_mix(signals, sample_rate, preset, mode);
  float* out2 = process_mix(signals, sample_rate, preset, mode);
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
    snprintf(name, sizeof(name), "%s_%s_mode%u_out", preset->name,
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

  Metrics m = {0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F};

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

  // SD over active frames
  m.sd_db = spectral_distortion(&clean_an, &out_an, active);

  // MNI v2
  m.mni_out = shape_cv(&out_an, noise_only);
  const float mni_mix = shape_cv(&mix_an, noise_only);
  m.mni_ratio = m.mni_out / mni_mix;

  // Decay deficit
  m.decay_db = decay_deficit(signals->clean, out, signals->mix, total,
                             signals->lead_samples, stream_delay);

  analyzer_free(&clean_an);
  analyzer_free(&out_an);
  analyzer_free(&mix_an);
  free(active);
  free(noise_only);
  free(out);
  return m;
}

int main() {
  setvbuf(stdout, NULL, _IONBF, 0); // live output; asserts abort otherwise

  int full_mode = 0;
  const char* full_env = getenv("SPECBLEACH_REALWORLD_FULL");
  if (full_env != NULL && full_env[0] != '\0' && full_env[0] != '0') {
    full_mode = 1;
  }
  const uint32_t num_presets = full_mode ? 3U : 1U; // CI default: moderate
  const uint32_t preset_offset = full_mode ? 0U : 1U;

  uint32_t sample_rate = 0U;
  uint32_t body_samples = 0U;

  // Load fixtures (all cases must share the sample rate)
  float* clean_body[4];
  float* noise_raw[4];
  Signals signals[4];
  for (uint32_t c = 0U; c < 4U; c++) {
    uint32_t sr = 0U;
    uint32_t len = 0U;
    clean_body[c] = load_wav(fixture_file(CASES[c].name, "clean"), &sr, &len);
    if (c == 0U) {
      sample_rate = sr;
      body_samples = len;
    }
    TEST_ASSERT(sr == sample_rate && len == body_samples,
                "fixture sample rates and lengths must match");
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
      "== real-world quality metrics (sample rate %u, %u cases, %s "
      "preset%s) ==\n",
      sample_rate, 4U, full_mode ? "all" : "moderate", full_mode ? "s" : "");

  for (uint32_t pi = 0U; pi < num_presets; pi++) {
    const Preset* preset = &PRESETS[preset_offset + pi];
    for (uint32_t mode = 0U; mode < NUM_MODES; mode++) {
      printf("== %s | %s ==\n", preset->name, mode_name(mode));
      Metrics agg = {0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F};
      for (uint32_t c = 0U; c < 4U; c++) {
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

        Metrics m = run_case_mode(sig, sample_rate, preset, mode);
        agg.att_db += m.att_db;
        agg.att_all_db += m.att_all_db;
        agg.sd_db += m.sd_db;
        agg.mni_ratio += m.mni_ratio;
        agg.mni_out += m.mni_out;
        agg.decay_db += m.decay_db;
        printf(
            "  %s: att %.1f dB (all %.1f), sd %.1f dB, mni %.2f (raw "
            "%.2f), decay %.1f dB\n",
            CASES[c].name, m.att_db, m.att_all_db, m.sd_db, m.mni_ratio,
            m.mni_out, m.decay_db);
        free(sig->clean);
        free(sig->noise);
        free(sig->mix);
        free(sig->snapshot);
      }
      agg.att_db /= 4.0F;
      agg.att_all_db /= 4.0F;
      agg.sd_db /= 4.0F;
      agg.mni_ratio /= 4.0F;
      agg.mni_out /= 4.0F;
      agg.decay_db /= 4.0F;
      printf(
          "  AGGREGATE: att %.1f dB (all %.1f), sd %.1f dB, mni %.2f, decay "
          "%.1f dB\n",
          agg.att_db, agg.att_all_db, agg.sd_db, agg.mni_ratio, agg.decay_db);

      if (!full_mode) {
        TEST_ASSERT(agg.att_db > ATT_GATES_DB[mode], "attenuation gate");
        TEST_ASSERT(agg.sd_db < SD_GATES_DB[mode], "spectral distortion gate");
        TEST_ASSERT(agg.mni_ratio < MNI_GATES[mode], "musical noise gate");
        TEST_ASSERT(agg.decay_db < DECAY_GATES_DB[mode], "decay smearing gate");
      }
    }
  }

  for (uint32_t c = 0U; c < 4U; c++) {
    free(clean_body[c]);
    free(noise_raw[c]);
  }
  printf("\n[OK] real-world quality metrics passed\n");
  return 0;
}
