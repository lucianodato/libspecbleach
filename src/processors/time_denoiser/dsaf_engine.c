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

#include "dsaf_engine.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/gain_calculator.h"
#include "shared/stft/fft_transform.h"
#include "shared/time_domain/fir_history_convolver.h"
#include "shared/time_domain/minimum_phase_fir.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/spectral_utils.h"
#include <math.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(_WIN32)
#include <windows.h>
#define SB_SLEEP_US(us) Sleep((us) / 1000U)
#else
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#define SB_SLEEP_US(us)                                                        \
  do {                                                                         \
    const struct timespec sb_ts = {(us) / 1000000, ((us) % 1000000) * 1000};   \
    nanosleep(&sb_ts, NULL);                                                   \
  } while (0)
#endif

#if defined(_WIN32)
#define SB_THREAD_RET DWORD
#define SB_THREAD_CALL WINAPI
typedef HANDLE sb_thread_t;
static sb_thread_t sb_thread_spawn(DWORD(WINAPI* fn)(LPVOID), void* arg) {
  return CreateThread(NULL, 0, fn, arg, 0, NULL);
}
static void sb_thread_join(sb_thread_t t) {
  WaitForSingleObject(t, INFINITE);
}
static sb_thread_t SB_THREAD_INVALID = NULL;
#define SB_THREAD_MAIN_DECL DWORD SB_THREAD_CALL
#else
#define SB_THREAD_RET void*
#define SB_THREAD_CALL
typedef pthread_t sb_thread_t;
static sb_thread_t sb_thread_invalid(void) {
  sb_thread_t t;
  memset(&t, 0, sizeof(t));
  return t;
}
static bool sb_thread_valid(sb_thread_t t) {
  sb_thread_t zero = sb_thread_invalid();
  return memcmp(&t, &zero, sizeof(t)) != 0;
}
static sb_thread_t sb_thread_spawn(void* (*fn)(void*), void* arg) {
  sb_thread_t t;
  if (pthread_create(&t, NULL, fn, arg) != 0) {
    return sb_thread_invalid();
  }
  return t;
}
static void sb_thread_join(sb_thread_t t) {
  pthread_join(t, NULL);
}
#define SB_THREAD_MAIN_DECL void*
#endif

struct SbDsafEngine {
  uint32_t sample_rate;
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t taps;

  // Real-time path (audio thread only)
  SbFirHistoryConvolver* convolver;
  int rt_coeff_idx;

  // Lock-free FIR coefficient publication (worker -> RT)
  float* coeffs_buffers[2];
  atomic_int active_coeff_idx;      // 0 or 1, published by worker
  atomic_bool coeff_buffer_free[2]; // RT may free a buffer it stopped using

  // SPSC hop ring (RT producer -> worker consumer)
  float* ring;
  uint32_t ring_capacity; // power of two
  uint32_t ring_mask;
  // Only producer writes write_idx, only consumer writes read_idx
  atomic_uint ring_write_idx;
  atomic_uint ring_read_idx;

  // Worker thread
  sb_thread_t worker_thread;
  bool worker_started;
  atomic_bool worker_running;

  // Worker-owned analysis state
  FftTransform* analysis_fft;
  FftTransform* cepstral_fft;
  AdaptiveNoiseEstimator* noise_estimator;
  CriticalBands* critical_bands;
  float* analysis_window;
  float* frame_buffer;
  float* power_spectrum;
  float* noise_power;
  float* frozen_noise_scratch;
  float* bin_gains;
  float* prev_clean_power;
  float* band_gains_raw;
  float* band_gains_smooth;
  float* taper_window;

  DsafParameters parameters;
};

static SB_THREAD_RET SB_THREAD_CALL dsaf_worker_main(void* arg);
static void dsaf_analyze_frame(SbDsafEngine* self);
static bool dsaf_ring_push(SbDsafEngine* self, uint32_t n, const float* input);
static uint32_t dsaf_ring_pop(SbDsafEngine* self, uint32_t n, float* output);

SbDsafEngine* dsaf_engine_initialize(const uint32_t sample_rate) {
  if (sample_rate == 0U) {
    return NULL;
  }

  SbDsafEngine* self = (SbDsafEngine*)calloc(1U, sizeof(SbDsafEngine));
  if (!self) {
    return NULL;
  }

  self->sample_rate = sample_rate;
  self->fft_size = DSAF_FFT_SIZE;
  self->real_spectrum_size = (self->fft_size / 2U) + 1U;
  self->taps = DSAF_FIR_TAPS;
  self->parameters.reduction_gain =
      powf(10.0F, DSAF_MIN_GAIN_DB / 20.0F); // ponytail: silent G_min default
  self->parameters.smoothing_factor = 0.5F;
  self->parameters.adaptive_noise = true;
  self->parameters.noise_estimation_method = SPP_MMSE_METHOD;
  self->parameters.suppression_strength = 0.5F;

  self->convolver = fir_history_convolver_initialize(self->taps);
  self->analysis_fft = fft_transform_initialize_bins(self->fft_size);
  self->cepstral_fft = fft_transform_initialize_bins(self->fft_size);
  self->noise_estimator = adaptive_estimator_initialize(
      self->real_spectrum_size, sample_rate, self->fft_size,
      self->parameters.noise_estimation_method);
  self->critical_bands =
      critical_bands_initialize(sample_rate, self->fft_size, ERB_SCALE);

  const uint32_t band_count =
      get_number_of_critical_bands(self->critical_bands);

  self->coeffs_buffers[0] = (float*)calloc(self->taps, sizeof(float));
  self->coeffs_buffers[1] = (float*)calloc(self->taps, sizeof(float));
  self->ring = (float*)calloc(DSAF_RING_CAPACITY, sizeof(float));
  self->analysis_window = (float*)calloc(self->fft_size, sizeof(float));
  self->frame_buffer = (float*)calloc(self->fft_size, sizeof(float));
  self->power_spectrum =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->noise_power = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->frozen_noise_scratch =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->bin_gains = (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->prev_clean_power =
      (float*)calloc(self->real_spectrum_size, sizeof(float));
  self->band_gains_raw = (float*)calloc(band_count, sizeof(float));
  self->band_gains_smooth = (float*)calloc(band_count, sizeof(float));
  self->taper_window = (float*)calloc(self->taps, sizeof(float));

  if (!self->convolver || !self->analysis_fft || !self->cepstral_fft ||
      !self->noise_estimator || !self->critical_bands ||
      !self->coeffs_buffers[0] || !self->coeffs_buffers[1] || !self->ring ||
      !self->analysis_window || !self->frame_buffer || !self->power_spectrum ||
      !self->noise_power || !self->frozen_noise_scratch || !self->bin_gains ||
      !self->prev_clean_power || !self->band_gains_raw ||
      !self->band_gains_smooth || !self->taper_window) {
    dsaf_engine_free(self);
    return NULL;
  }

  get_fft_window(self->analysis_window, self->fft_size, HANN_WINDOW);
  minimum_phase_fir_taper_window(self->taps, self->taper_window);

  // Delta pass-through filter until the first coefficient update lands
  self->coeffs_buffers[0][0] = 1.0F;
  self->coeffs_buffers[1][0] = 1.0F;

  self->ring_capacity = DSAF_RING_CAPACITY;
  self->ring_mask = self->ring_capacity - 1U;
  atomic_init(&self->ring_write_idx, 0U);
  atomic_init(&self->ring_read_idx, 0U);
  atomic_init(&self->active_coeff_idx, 0);
  atomic_init(&self->coeff_buffer_free[0], false);
  atomic_init(&self->coeff_buffer_free[1], true);
  atomic_init(&self->worker_running, true);
  self->rt_coeff_idx = 0;
  self->worker_started = false;

  self->worker_thread = sb_thread_spawn(dsaf_worker_main, self);
  if (!sb_thread_valid(self->worker_thread)) {
    dsaf_engine_free(self);
    return NULL;
  }
  self->worker_started = true;

  return self;
}

void dsaf_engine_free(SbDsafEngine* self) {
  if (!self) {
    return;
  }

  atomic_store(&self->worker_running, false);
  if (self->worker_started) {
    sb_thread_join(self->worker_thread);
  }

  fir_history_convolver_free(self->convolver);
  fft_transform_free(self->analysis_fft);
  fft_transform_free(self->cepstral_fft);
  adaptive_estimator_free(self->noise_estimator);
  critical_bands_free(self->critical_bands);
  free(self->coeffs_buffers[0]);
  free(self->coeffs_buffers[1]);
  free(self->ring);
  free(self->analysis_window);
  free(self->frame_buffer);
  free(self->power_spectrum);
  free(self->noise_power);
  free(self->frozen_noise_scratch);
  free(self->bin_gains);
  free(self->prev_clean_power);
  free(self->band_gains_raw);
  free(self->band_gains_smooth);
  free(self->taper_window);
  free(self);
}

void dsaf_engine_load_parameters(SbDsafEngine* self,
                                 const DsafParameters* parameters) {
  if (!self || !parameters) {
    return;
  }
  // Load/process never run concurrently on the same instance
  DsafParameters clamped = *parameters;
  clamped.reduction_gain = fmaxf(fminf(clamped.reduction_gain, 1.0F), 0.0F);
  clamped.smoothing_factor = fmaxf(fminf(clamped.smoothing_factor, 1.0F), 0.0F);
  clamped.suppression_strength =
      fmaxf(fminf(clamped.suppression_strength, 1.0F), 0.0F);

  // Estimator internals are method-specific: rebuild on change. Setup-only
  // context, allocation is allowed here.
  if (clamped.noise_estimation_method !=
      self->parameters.noise_estimation_method) {
    adaptive_estimator_free(self->noise_estimator);
    self->noise_estimator = adaptive_estimator_initialize(
        self->real_spectrum_size, self->sample_rate, self->fft_size,
        clamped.noise_estimation_method);
  }

  self->parameters = clamped;
}

uint32_t dsaf_engine_get_latency(SbDsafEngine* self) {
  (void)self;
  return 0U; // Delayless by construction
}

bool dsaf_engine_process(SbDsafEngine* self, const uint32_t number_of_samples,
                         const float* input, float* output) {
  if (!self || number_of_samples == 0U || !input || !output) {
    return false;
  }

  // Feed the analysis worker (lock-free SPSC ring). Chunked so host block
  // sizes larger than the ring capacity still reach the worker; overflow
  // drops only the slices that don't fit (worker temporarily behind).
  const float* input_slice = input;
  uint32_t remaining = number_of_samples;
  while (remaining > 0U) {
    const uint32_t chunk =
        remaining > self->fft_size ? self->fft_size : remaining;
    dsaf_ring_push(self, chunk, input_slice);
    input_slice += chunk;
    remaining -= chunk;
  }

  // Pick up newly published FIR coefficients (lock-free double buffer)
  const int idx = atomic_load(&self->active_coeff_idx);
  if (idx != self->rt_coeff_idx) {
    atomic_store(&self->coeff_buffer_free[self->rt_coeff_idx], true);
    self->rt_coeff_idx = idx;
  }

  fir_history_convolver_process(self->convolver, number_of_samples, input,
                                output,
                                self->coeffs_buffers[self->rt_coeff_idx]);

  return true;
}

/* --------------------------------------------------------------------- */
/* --------------------------- SPSC hop ring --------------------------- */
/* --------------------------------------------------------------------- */

static bool dsaf_ring_push(SbDsafEngine* self, const uint32_t n,
                           const float* input) {
  const uint32_t write =
      atomic_load_explicit(&self->ring_write_idx, memory_order_relaxed);
  const uint32_t read =
      atomic_load_explicit(&self->ring_read_idx, memory_order_acquire);
  const uint32_t free_space = self->ring_capacity - (write - read);

  if (n > free_space) {
    return false; // Worker is behind: skip this analysis feed
  }

  for (uint32_t i = 0U; i < n; i++) {
    self->ring[(write + i) & self->ring_mask] = input[i];
  }
  atomic_store_explicit(&self->ring_write_idx, write + n, memory_order_release);
  return true;
}

static uint32_t dsaf_ring_pop(SbDsafEngine* self, const uint32_t n,
                              float* output) {
  const uint32_t read =
      atomic_load_explicit(&self->ring_read_idx, memory_order_relaxed);
  const uint32_t write =
      atomic_load_explicit(&self->ring_write_idx, memory_order_acquire);
  const uint32_t available = write - read;

  if (available < n) {
    return 0U;
  }

  for (uint32_t i = 0U; i < n; i++) {
    output[i] = self->ring[(read + i) & self->ring_mask];
  }
  atomic_store_explicit(&self->ring_read_idx, read + n, memory_order_release);
  return n;
}

/* --------------------------------------------------------------------- */
/* ------------------------ Worker thread body ------------------------- */
/* --------------------------------------------------------------------- */

static SB_THREAD_RET SB_THREAD_CALL dsaf_worker_main(void* arg) {
  SbDsafEngine* self = (SbDsafEngine*)arg;

  while (atomic_load(&self->worker_running)) {
    const int target = 1 - atomic_load(&self->active_coeff_idx);

    if (!atomic_load(&self->coeff_buffer_free[target])) {
      // RT thread has not released the inactive buffer yet: drop this update
      SB_SLEEP_US(DSAF_WORKER_POLL_US);
      continue;
    }

    if (dsaf_ring_pop(self, self->fft_size, self->frame_buffer) == 0U) {
      SB_SLEEP_US(DSAF_WORKER_POLL_US);
      continue;
    }

    // Reserve the buffer before synthesizing into it
    atomic_store(&self->coeff_buffer_free[target], false);
    dsaf_analyze_frame(self);

    if (minimum_phase_fir_synthesize(
            self->fft_size, self->taps, self->bin_gains, self->cepstral_fft,
            self->taper_window, self->coeffs_buffers[target])) {
      atomic_store_explicit(&self->active_coeff_idx, target,
                            memory_order_release);
    } else {
      atomic_store(&self->coeff_buffer_free[target], true);
    }
  }

  return NULL;
}

static void dsaf_analyze_frame(SbDsafEngine* self) {
  float* fft_in = get_fft_input_buffer(self->analysis_fft);

  // Windowed analysis frame
  for (uint32_t i = 0U; i < self->fft_size; i++) {
    fft_in[i] = self->frame_buffer[i] * self->analysis_window[i];
  }

  if (!compute_forward_fft(self->analysis_fft)) {
    return;
  }

  const float* fft_out = get_fft_output_buffer(self->analysis_fft);
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const uint32_t im_idx = (k == 0U || k == self->real_spectrum_size - 1U)
                                ? k
                                : self->fft_size - k;
    const float re = fft_out[k];
    const float im = fft_out[im_idx];
    self->power_spectrum[k] = re * re + im * im;
  }

  // Per-bin noise floor tracking (reuses the shared statistical estimators)
  float aggressiveness = 0.F;
  const float param_aggressiveness = self->parameters.suppression_strength;
  adaptive_estimator_run(self->noise_estimator, self->power_spectrum,
                         self->frozen_noise_scratch, &aggressiveness,
                         param_aggressiveness);
  if (self->parameters.adaptive_noise) {
    memcpy(self->noise_power, self->frozen_noise_scratch,
           self->real_spectrum_size * sizeof(float));
  }
  // When adaptive_noise is false, the last tracked floor is held in
  // noise_power and only refreshed into a scratch buffer.

  // Decision-directed a priori SNR + Wiener gain per bin
  const float g_min = self->parameters.reduction_gain;
  const float alpha = DSAF_A_PRIORI_SNR_ALPHA;
  const float one_minus_alpha = 1.0F - alpha;
  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    const float px = self->power_spectrum[k];
    const float pn = fmaxf(
        self->noise_power[k] * (1.0F + self->parameters.suppression_strength),
        SPECTRAL_EPSILON);
    const float gamma = px / pn;
    const float a_posteriori = fmaxf(gamma - 1.0F, 0.0F);
    const float a_priori = alpha * (self->prev_clean_power[k] / pn) +
                           one_minus_alpha * a_posteriori;
    const float wiener = a_priori / (1.0F + a_priori);
    self->bin_gains[k] = fmaxf(wiener, g_min);

    // Clean-speech power estimate for the next decision-directed step
    const float clean = self->bin_gains[k] * self->bin_gains[k] * px;
    self->prev_clean_power[k] =
        fmaxf(self->prev_clean_power[k] * DSAF_CLEAN_POWER_DECAY, clean);
  }

  // Pool bin gains into ERB bands (mean gain per band)
  compute_critical_bands_spectrum(self->critical_bands, self->bin_gains,
                                  self->band_gains_raw);
  const uint32_t band_count =
      get_number_of_critical_bands(self->critical_bands);
  for (uint32_t b = 0U; b < band_count; b++) {
    const CriticalBandIndexes idx = get_band_indexes(self->critical_bands, b);
    const float bins = (float)(idx.end_position - idx.start_position);
    self->band_gains_raw[b] /= fmaxf(bins, 1.0F);

    // Temporal smoothing per band
    const float s = self->parameters.smoothing_factor;
    self->band_gains_smooth[b] =
        s * self->band_gains_smooth[b] + (1.0F - s) * self->band_gains_raw[b];
  }

  // Expand band gains to bins (piecewise constant) and smooth boundaries
  for (uint32_t b = 0U; b < band_count; b++) {
    const CriticalBandIndexes idx = get_band_indexes(self->critical_bands, b);
    for (uint32_t k = idx.start_position; k < idx.end_position; k++) {
      self->bin_gains[k] = self->band_gains_smooth[b];
    }
  }
  smooth_spectrum(self->bin_gains, self->real_spectrum_size,
                  DSAF_ERB_GAIN_SMOOTHING);

  // Gentle high-frequency taper toward unity to avoid edge artifacts
  const uint32_t taper_start =
      self->real_spectrum_size -
      (uint32_t)((float)self->real_spectrum_size * DSAF_HF_TAPER_RATIO);
  for (uint32_t k = taper_start; k < self->real_spectrum_size; k++) {
    const float t = (float)(k - taper_start) /
                    (float)(self->real_spectrum_size - taper_start);
    self->bin_gains[k] = self->bin_gains[k] + (1.0F - self->bin_gains[k]) * t;
  }
}
