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

#include "shared/utils/spectral_neighbors.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void spectral_neighbors_compute_power_row(float* power_row,
                                          const float* fft_spectrum,
                                          uint32_t fft_size) {
  if (!power_row || !fft_spectrum || fft_size < 2U) {
    return;
  }

  const uint32_t n = fft_size;
  const uint32_t n2 = n / 2U;
  const bool is_even = (n % 2U == 0);

  // DC bin (purely real)
  power_row[0] = fft_spectrum[0] * fft_spectrum[0];

  // Complex bins (imag mirrored at n - k)
  for (uint32_t k = 1U; k < n2; k++) {
    const float real = fft_spectrum[k];
    const float imag = fft_spectrum[n - k];
    power_row[k] = (real * real) + (imag * imag);
  }

  // Nyquist bin (purely real, even transform only)
  if (is_even) {
    power_row[n2] = fft_spectrum[n2] * fft_spectrum[n2];
  }
}

// Clamp a signed neighbor tap bin to the valid spectrum range
// (edge-extension at DC and Nyquist)
static inline uint32_t spectral_neighbors_clamp_bin(int32_t bin,
                                                    uint32_t spectrum_size) {
  if (bin < 0) {
    return 0U;
  }
  if ((uint32_t)bin >= spectrum_size) {
    return spectrum_size - 1U;
  }
  return (uint32_t)bin;
}

// Frequency axis with per-bin Erb-scaled variable-radius Hann. Edge taps
// clamp to DC/Nyquist and the used weights renormalize so the output stays
// a true average there too.
static void spectral_neighbors_average_variable(
    const SpectralNeighborFreqKernel* freq, const float* acc_scratch,
    uint32_t spectrum_size, float* out) {
  for (uint32_t k = 0U; k < spectrum_size; k++) {
    const uint32_t halo = freq->halo[k];
    const uint32_t first = freq->start[k];
    float sum = 0.0F;
    float used = 0.0F;
    for (uint32_t i = 0U; i < 2U * halo + 1U; i++) {
      const uint32_t bin = spectral_neighbors_clamp_bin(
          (int32_t)k + ((int32_t)i - (int32_t)halo), spectrum_size);
      sum += freq->weights[first + i] * acc_scratch[bin];
      used += freq->weights[first + i];
    }
    out[k] = sum * (1.0F / used);
  }
}

// Frequency axis: symmetric FIR, edges clamped (edge-extension)
static void spectral_neighbors_average_fixed(const float* kernel,
                                             const float* acc_scratch,
                                             uint32_t spectrum_size,
                                             float* out) {
  for (uint32_t k = 0U; k < spectrum_size; k++) {
    float sum = 0.0F;
    for (uint32_t j = 0U; j < NEIGHBOR_KERNEL_TAPS; j++) {
      uint32_t idx = k + j - NEIGHBOR_KERNEL_HALO;
      if (idx >= spectrum_size) { // unsigned clamp: negatives wrap high
        idx = (j < NEIGHBOR_KERNEL_HALO) ? 0U : (spectrum_size - 1U);
      }
      sum += kernel[j] * acc_scratch[idx];
    }
    out[k] = sum;
  }
}

bool spectral_neighbor_average(const SpectralNeighborParameters* parameters,
                               float* acc_scratch, float* out) {
  if (!parameters || !parameters->rows || !acc_scratch || !out ||
      parameters->spectrum_size == 0U || parameters->rows_count == 0U ||
      (parameters->rows_count % 2U) == 0U ||
      parameters->rows_count > NEIGHBOR_KERNEL_TAPS) {
    return false;
  }

  const float kernel[NEIGHBOR_KERNEL_TAPS] = NEIGHBOR_KERNEL_W;
  const uint32_t spectrum_size = parameters->spectrum_size;

  // Time axis: every row participates fully or not at all, so the presence
  // weighting is a single scalar and the accumulation stays per-bin linear.
  float time_weight = 0.0F;
  memset(acc_scratch, 0, spectrum_size * sizeof(float));
  for (uint32_t i = 0U; i < parameters->rows_count; i++) {
    const float* row = parameters->rows[i];
    if (!row) {
      continue;
    }
    const float weight = kernel[i];
    for (uint32_t k = 0U; k < spectrum_size; k++) {
      acc_scratch[k] += weight * row[k];
    }
    time_weight += weight;
  }
  if (time_weight <= 0.0F) {
    return false;
  }

  // Reciprocal multiply instead of division inside the bin loops
  const float inv_weight = 1.0F / time_weight;
  for (uint32_t k = 0U; k < spectrum_size; k++) {
    acc_scratch[k] *= inv_weight;
  }

  if (parameters->freq_kernel != NULL) {
    spectral_neighbors_average_variable(parameters->freq_kernel, acc_scratch,
                                        spectrum_size, out);
    return true;
  }

  spectral_neighbors_average_fixed(kernel, acc_scratch, spectrum_size, out);

  return true;
}

SpectralNeighborFreqKernel* spectral_neighbors_freq_kernel_create(
    float bin_hz, uint32_t spectrum_size) {
  if (!(bin_hz > 0.0F) || spectrum_size == 0U) {
    return NULL;
  }

  SpectralNeighborFreqKernel* freq = (SpectralNeighborFreqKernel*)calloc(
      1U, sizeof(SpectralNeighborFreqKernel));
  if (!freq) {
    return NULL;
  }

  const size_t bins = (size_t)spectrum_size;
  freq->spectrum_size = spectrum_size;
  freq->halo = (uint32_t*)calloc(bins, sizeof(uint32_t));
  freq->start = (uint32_t*)calloc(bins, sizeof(uint32_t));

  uint32_t total = 0U;
  for (uint32_t k = 0U; k < spectrum_size; k++) {
    // Erb half-width grows with center frequency; converts to bins
    const float center_khz = ((float)k * bin_hz) / 1000.0F;
    const float erb_hz =
        NEIGHBOR_ERB_ROOT_HZ * (NEIGHBOR_ERB_PER_KHZ * center_khz + 1.0F);
    const float halo_bins = (NEIGHBOR_ERB_HALO_SHARE * erb_hz) / bin_hz;
    uint32_t halo = (uint32_t)(halo_bins + 0.5F);
    if (halo < NEIGHBOR_FREQ_HALO_MIN) {
      halo = NEIGHBOR_FREQ_HALO_MIN;
    }
    if (halo > NEIGHBOR_FREQ_HALO_MAX) {
      halo = NEIGHBOR_FREQ_HALO_MAX;
    }
    freq->halo[k] = halo;
    freq->start[k] = total;
    total += 2U * halo + 1U;
  }
  freq->total_taps = total;

  freq->weights = (float*)calloc((size_t)total, sizeof(float));
  if (!freq->halo || !freq->start || !freq->weights) {
    spectral_neighbors_freq_kernel_free(freq);
    return NULL;
  }

  // Flat per-bin Hann taper, each bin's taps summing to 1
  for (uint32_t k = 0U; k < spectrum_size; k++) {
    const uint32_t halo = freq->halo[k];
    const uint32_t start = freq->start[k];
    float sum = 0.0F;
    for (uint32_t i = 0U; i < 2U * halo + 1U; i++) {
      const float d = (float)i - (float)halo; // signed offset from center
      const float w = 0.5F * (1.0F + cosf((M_PIf * d) / (float)halo));
      freq->weights[start + i] = w;
      sum += w;
    }
    const float inv = 1.0F / sum;
    for (uint32_t i = 0U; i < 2U * halo + 1U; i++) {
      freq->weights[start + i] *= inv;
    }
  }

  return freq;
}

void spectral_neighbors_freq_kernel_free(SpectralNeighborFreqKernel* self) {
  if (!self) {
    return;
  }
  free(self->halo);
  free(self->start);
  free(self->weights);
  free(self);
}
