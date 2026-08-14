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

#include "hpss_filter.h"
#include "shared/configurations.h"
#include "shared/utils/spectral_circular_buffer.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct HpssFilter {
  HpssConfig config;
  uint32_t latency_frames;
  SbSpectralCircularBuffer* circular_buffer;
  uint32_t mag_layer_id;
  uint32_t rect_mag_layer_id;

  float* rectified_mag;
  float* prev_mag;
  float* m_h;
  float* m_p;
  float* time_sort_buffer;
  float* onset_boost_buffer;
  float* onset_ratio_buffer;

  uint32_t write_pos;
  float delayed_onset_ratio;
  bool is_initialized_flux;
};

static float fast_median(float* arr, uint32_t n) {
  for (uint32_t i = 1U; i < n; ++i) {
    float key = arr[i];
    int32_t j = (int32_t)i - 1;
    while (j >= 0 && arr[j] > key) {
      arr[j + 1] = arr[j];
      j--;
    }
    arr[j + 1] = key;
  }
  if (n % 2U == 1U) {
    return arr[n / 2U];
  }
  return 0.5f * (arr[(n / 2U) - 1U] + arr[n / 2U]);
}

static void compute_adaptive_freq_median(const float* src, float* dst,
                                         int num_bins, int max_win_size) {
  int max_half = max_win_size / 2;
  if (max_half > 31) {
    max_half = 31;
  }
  float scratch[64];

  for (int k = 0; k < num_bins; k++) {
    // Scale half-window size: 1 bin at bass, scaling up to max_half at high frequencies
    int half_win = k / 12;
    if (half_win < 1) {
      half_win = 1; // 3-bin window for bass (k-1, k, k+1)
    }
    if (half_win > max_half) {
      half_win = max_half; // Cap at max_half (e.g. 8 for 17-bin)
    }

    int start_idx = k - half_win;
    if (start_idx < 0) {
      start_idx = 0;
    }

    int end_idx = k + half_win;
    if (end_idx >= num_bins) {
      end_idx = num_bins - 1;
    }

    int count = end_idx - start_idx + 1;

    for (int i = 0; i < count; i++) {
      scratch[i] = src[start_idx + i];
    }

    // Fast insertion sort for small array
    for (int i = 1; i < count; i++) {
      float key = scratch[i];
      int j = i - 1;
      while (j >= 0 && scratch[j] > key) {
        scratch[j + 1] = scratch[j];
        j--;
      }
      scratch[j + 1] = key;
    }

    dst[k] = scratch[count / 2];
  }
}

HpssFilter* hpss_filter_initialize(HpssConfig config) {
  if (config.real_spectrum_size == 0U || config.time_window_size == 0U ||
      config.freq_window_size == 0U) {
    return NULL;
  }

  // Enforce odd window sizes
  if (config.time_window_size % 2U == 0U) {
    config.time_window_size += 1U;
  }
  if (config.freq_window_size % 2U == 0U) {
    config.freq_window_size += 1U;
  }

  if (config.noise_oversubtraction <= 0.0f) {
    config.noise_oversubtraction = HPSS_NOISE_OVERSUBTRACTION_DEFAULT;
  }

  HpssFilter* self = (HpssFilter*)calloc(1U, sizeof(HpssFilter));
  if (!self) {
    return NULL;
  }

  self->config = config;
  self->latency_frames = (config.time_window_size - 1U) / 2U;
  self->write_pos = 0U;
  self->delayed_onset_ratio = 0.0f;
  self->is_initialized_flux = false;

  // Pre-allocate circular buffer to maximum window size
  self->circular_buffer =
      spectral_circular_buffer_create(HPSS_TIME_WINDOW_MAX);
  if (!self->circular_buffer) {
    hpss_filter_free(self);
    return NULL;
  }

  self->mag_layer_id = spectral_circular_buffer_add_layer(
      self->circular_buffer, config.real_spectrum_size);
  self->rect_mag_layer_id = spectral_circular_buffer_add_layer(
      self->circular_buffer, config.real_spectrum_size);

  if (self->mag_layer_id == 0xFFFFFFFFU ||
      self->rect_mag_layer_id == 0xFFFFFFFFU) {
    hpss_filter_free(self);
    return NULL;
  }

  self->rectified_mag = (float*)calloc(config.real_spectrum_size, sizeof(float));
  self->prev_mag = (float*)calloc(config.real_spectrum_size, sizeof(float));
  self->m_h = (float*)calloc(config.real_spectrum_size, sizeof(float));
  self->m_p = (float*)calloc(config.real_spectrum_size, sizeof(float));
  self->time_sort_buffer =
      (float*)calloc(HPSS_TIME_WINDOW_MAX, sizeof(float));
  self->onset_boost_buffer =
      (float*)calloc(HPSS_TIME_WINDOW_MAX, sizeof(float));
  self->onset_ratio_buffer =
      (float*)calloc(HPSS_TIME_WINDOW_MAX, sizeof(float));

  if (!self->rectified_mag || !self->prev_mag || !self->m_h || !self->m_p ||
      !self->time_sort_buffer || !self->onset_boost_buffer ||
      !self->onset_ratio_buffer) {
    hpss_filter_free(self);
    return NULL;
  }

  for (uint32_t i = 0U; i < HPSS_TIME_WINDOW_MAX; ++i) {
    self->onset_boost_buffer[i] = 0.0f;
    self->onset_ratio_buffer[i] = 0.0f;
  }

  return self;
}

void hpss_filter_free(HpssFilter* self) {
  if (!self) {
    return;
  }

  if (self->circular_buffer) {
    spectral_circular_buffer_free(self->circular_buffer);
  }
  if (self->rectified_mag) {
    free(self->rectified_mag);
  }
  if (self->prev_mag) {
    free(self->prev_mag);
  }
  if (self->m_h) {
    free(self->m_h);
  }
  if (self->m_p) {
    free(self->m_p);
  }
  if (self->time_sort_buffer) {
    free(self->time_sort_buffer);
  }
  if (self->onset_boost_buffer) {
    free(self->onset_boost_buffer);
  }
  if (self->onset_ratio_buffer) {
    free(self->onset_ratio_buffer);
  }
  free(self);
}

void hpss_filter_set_quality_mode(HpssFilter* self, HpssQualityMode mode) {
  if (!self) {
    return;
  }

  uint32_t new_time_win = HPSS_TIME_WINDOW_MEDIUM;
  uint32_t new_freq_win = HPSS_FREQ_WINDOW_MEDIUM;
  switch (mode) {
    case HPSS_QUALITY_LOW:
      new_time_win = HPSS_TIME_WINDOW_LOW;
      new_freq_win = HPSS_FREQ_WINDOW_LOW;
      break;
    case HPSS_QUALITY_MEDIUM:
      new_time_win = HPSS_TIME_WINDOW_MEDIUM;
      new_freq_win = HPSS_FREQ_WINDOW_MEDIUM;
      break;
    case HPSS_QUALITY_HIGH:
      new_time_win = HPSS_TIME_WINDOW_HIGH;
      new_freq_win = HPSS_FREQ_WINDOW_HIGH;
      break;
    default:
      new_time_win = HPSS_TIME_WINDOW_MEDIUM;
      new_freq_win = HPSS_FREQ_WINDOW_MEDIUM;
      break;
  }

  if (self->config.time_window_size == new_time_win &&
      self->config.freq_window_size == new_freq_win) {
    return;
  }

  self->config.time_window_size = new_time_win;
  self->config.freq_window_size = new_freq_win;
  self->latency_frames = (new_time_win - 1U) / 2U;

  if (self->circular_buffer) {
    spectral_circular_buffer_clear(self->circular_buffer);
  }
  if (self->onset_boost_buffer) {
    memset(self->onset_boost_buffer, 0,
           HPSS_TIME_WINDOW_MAX * sizeof(float));
  }
  if (self->onset_ratio_buffer) {
    memset(self->onset_ratio_buffer, 0,
           HPSS_TIME_WINDOW_MAX * sizeof(float));
  }
  self->write_pos = 0U;
  self->delayed_onset_ratio = 0.0f;
  self->is_initialized_flux = false;
}

uint32_t hpss_filter_get_latency_frames(const HpssFilter* self) {
  return self ? self->latency_frames : 0U;
}

float hpss_filter_get_onset_ratio(const HpssFilter* self) {
  return self ? self->delayed_onset_ratio : 0.0f;
}

bool hpss_filter_process(HpssFilter* self, const float* current_magnitude,
                         const float* noise_profile,
                         float* delayed_magnitude_out, float* mask_harmonic_out,
                         float* mask_percussive_out) {
  if (!self || !current_magnitude) {
    return false;
  }

  const uint32_t K = self->config.real_spectrum_size;
  const uint32_t time_win = self->config.time_window_size;
  const uint32_t freq_win = self->config.freq_window_size;
  const float gamma = self->config.noise_oversubtraction;

  // 1. Noise Rectification on incoming frame
  if (noise_profile) {
    for (uint32_t k = 0U; k < K; ++k) {
      float rect = current_magnitude[k] - gamma * noise_profile[k];
      self->rectified_mag[k] = (rect > 0.0f) ? rect : 0.0f;
    }
  } else {
    memcpy(self->rectified_mag, current_magnitude, K * sizeof(float));
  }

  // 2. Frame-level Onset / Spectral Flux Detection
  float onset_ratio = 0.0f;
  float boost_add = 0.0f;
  if (self->is_initialized_flux) {
    float flux = 0.0f;
    float prev_sum = 0.0f;
    for (uint32_t k = 0U; k < K; ++k) {
      float diff = current_magnitude[k] - self->prev_mag[k];
      if (diff > 0.0f) {
        flux += diff;
      }
      prev_sum += self->prev_mag[k];
    }
    onset_ratio = flux / (prev_sum + 1e-6f);
    boost_add = 2.0f * onset_ratio;
    if (boost_add > 3.0f) {
      boost_add = 3.0f;
    } else if (boost_add < 0.0f) {
      boost_add = 0.0f;
    }
  } else {
    self->is_initialized_flux = true;
  }
  memcpy(self->prev_mag, current_magnitude, K * sizeof(float));
  self->onset_boost_buffer[self->write_pos] = boost_add;
  self->onset_ratio_buffer[self->write_pos] = onset_ratio;

  // 3. Push unrectified and rectified frames into circular buffer
  spectral_circular_buffer_push(self->circular_buffer, self->mag_layer_id,
                                current_magnitude);
  spectral_circular_buffer_push(self->circular_buffer, self->rect_mag_layer_id,
                                self->rectified_mag);

  // 4. Retrieve delayed frames
  const float* delayed_mag = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->mag_layer_id, self->latency_frames);
  const float* delayed_rect_mag = spectral_circular_buffer_retrieve(
      self->circular_buffer, self->rect_mag_layer_id, self->latency_frames);

  if (!delayed_mag || !delayed_rect_mag) {
    return false;
  }

  if (delayed_magnitude_out) {
    memcpy(delayed_magnitude_out, delayed_mag, K * sizeof(float));
  }

  uint32_t delayed_idx =
      (self->write_pos + HPSS_TIME_WINDOW_MAX - self->latency_frames) %
      HPSS_TIME_WINDOW_MAX;
  float frame_onset_boost = self->onset_boost_buffer[delayed_idx];
  self->delayed_onset_ratio = self->onset_ratio_buffer[delayed_idx];

  // 5. 1D Median Filtering on Rectified Spectrum
  // 5a. Harmonic median (M_H) along time axis for each frequency bin
  for (uint32_t k = 0U; k < K; ++k) {
    for (uint32_t t = 0U; t < time_win; ++t) {
      const float* frame_t = spectral_circular_buffer_retrieve(
          self->circular_buffer, self->rect_mag_layer_id, t);
      self->time_sort_buffer[t] = frame_t ? frame_t[k] : 0.0f;
    }
    self->m_h[k] = fast_median(self->time_sort_buffer, time_win);
  }

  // 5b. Percussive median (M_P) along frequency axis for delayed frame
  compute_adaptive_freq_median(delayed_rect_mag, self->m_p, (int)K, (int)freq_win);

  // 6. Onset-Boosted Soft Masking with Frequency Exponential Decay
  for (uint32_t k = 0U; k < K; ++k) {
    float m_h = self->m_h[k];
    float bass_decay = expf(-(float)k / HPSS_BASS_CUTOFF_BINS);
    float m_p_boosted = self->m_p[k] * (1.0f + frame_onset_boost * bass_decay);

    float m_h_sq = m_h * m_h;
    float m_p_sq = m_p_boosted * m_p_boosted;
    float sum_sq = m_h_sq + m_p_sq;

    float w_p = 0.0f;
    float w_h = 1.0f;
    if (sum_sq > SPECTRAL_EPSILON) {
      w_h = m_h_sq / sum_sq;
      w_p = m_p_sq / sum_sq;
    }

    if (mask_harmonic_out) {
      mask_harmonic_out[k] = w_h;
    }
    if (mask_percussive_out) {
      mask_percussive_out[k] = w_p;
    }
  }

  // 7. Advance circular buffer and write position
  self->write_pos = (self->write_pos + 1U) % HPSS_TIME_WINDOW_MAX;
  spectral_circular_buffer_advance(self->circular_buffer);

  return true;
}
