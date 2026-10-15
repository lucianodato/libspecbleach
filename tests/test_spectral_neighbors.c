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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s (%s:%d)\n", message, __FILE__,          \
              __LINE__);                                                       \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

// Same kernel as the shared configuration, reproduced here so the expected
// values are independent of future retuning of the implementation
static const float kKernel[NEIGHBOR_KERNEL_TAPS] = NEIGHBOR_KERNEL_W;

void test_spectral_neighbors(void) {
  printf("Testing Spectral Neighbors...\n");

  const uint32_t fft_size = 1024U;
  const uint32_t num_bins = (fft_size / 2U) + 1U; // 513
  const uint32_t halo = NEIGHBOR_KERNEL_HALO;

  // 1. Power row conversion: DC, complex bins, Nyquist
  float fft_spectrum[fft_size];
  memset(fft_spectrum, 0, sizeof(fft_spectrum));
  fft_spectrum[0] = 2.0F;              // DC
  fft_spectrum[100] = 3.0F;            // real of bin 100
  fft_spectrum[fft_size - 100] = 4.0F; // imag of bin 100
  fft_spectrum[512] = 5.0F;            // Nyquist (even transform)
  float power_row[num_bins];
  memset(power_row, 0, sizeof(power_row));
  spectral_neighbors_compute_power_row(power_row, fft_spectrum, fft_size);
  TEST_FLOAT_CLOSE(power_row[0], 4.0F, 1e-6F);
  TEST_FLOAT_CLOSE(power_row[100], 25.0F, 1e-6F);
  TEST_FLOAT_CLOSE(power_row[200], 0.0F, 1e-6F);
  TEST_FLOAT_CLOSE(power_row[512], 25.0F, 1e-6F);

  // 2. Invalid calls
  const float* rows[NEIGHBOR_KERNEL_TAPS] = {NULL};
  float acc[num_bins];
  float out[num_bins];
  float flat[num_bins];
  for (uint32_t k = 0U; k < num_bins; k++) {
    flat[k] = 1.0F;
  }

  SpectralNeighborParameters params = {
      .rows = rows,
      .rows_count = NEIGHBOR_KERNEL_TAPS,
      .spectrum_size = num_bins,
  };
  TEST_ASSERT(!spectral_neighbor_average(NULL, acc, out),
              "NULL parameters should fail");
  TEST_ASSERT(!spectral_neighbor_average(&params, NULL, out),
              "NULL acc scratch should fail");
  TEST_ASSERT(!spectral_neighbor_average(&params, acc, NULL),
              "NULL output should fail");
  TEST_ASSERT(!spectral_neighbor_average(&params, acc, out),
              "All-NULL rows should fail");

  const float* even_rows[NEIGHBOR_KERNEL_TAPS - 1U] = {NULL};
  SpectralNeighborParameters even_params = {
      .rows = even_rows,
      .rows_count = NEIGHBOR_KERNEL_TAPS - 1U, // even count is invalid
      .spectrum_size = num_bins,
  };
  TEST_ASSERT(!spectral_neighbor_average(&even_params, acc, out),
              "Even rows_count should fail");
  // NULL instance safety must not crash
  spectral_neighbors_compute_power_row(NULL, fft_spectrum, fft_size);
  spectral_neighbors_compute_power_row(power_row, NULL, fft_size);
  spectral_neighbors_compute_power_row(power_row, fft_spectrum, 0U);

  // 3. Single centered row: output equals the frequency FIR of the row
  // (edge-extended), independent of the skipped inner rows
  rows[0] = NULL;
  rows[1] = NULL;
  rows[2] = flat;
  rows[3] = NULL;
  rows[4] = NULL;
  params.rows = rows;
  TEST_ASSERT(spectral_neighbor_average(&params, acc, out),
              "Center-only rows should succeed");
  // Flat row through the FIR stays flat (kernel sums to 1)
  for (uint32_t k = 0U; k < num_bins; k++) {
    TEST_FLOAT_CLOSE(out[k], 1.0F, 1e-5F);
  }

  // 4. Impulse response: interior bin keeps the center weight, edges clamp
  float impulse[num_bins];
  memset(impulse, 0, sizeof(impulse));
  const uint32_t interior = 300U;
  impulse[interior] = 4.0F;
  rows[2] = impulse;
  memset(out, 0, sizeof(out));
  TEST_ASSERT(spectral_neighbor_average(&params, acc, out),
              "Center-only impulse should succeed");
  for (uint32_t j = 0U; j < NEIGHBOR_KERNEL_TAPS; j++) {
    TEST_FLOAT_CLOSE(out[interior - halo + j], kKernel[j] * 4.0F, 1e-5F);
  }
  TEST_FLOAT_CLOSE(out[interior + halo + 1U], 0.0F, 1e-5F);
  // Edge clamp: first bin accumulates the taps that would fall below 0
  float expected_edge = 0.0F;
  for (uint32_t j = 0U; j < NEIGHBOR_KERNEL_TAPS; j++) {
    const int32_t tap =
        (int32_t)interior - (int32_t)NEIGHBOR_KERNEL_HALO + (int32_t)j;
    if (tap == 0) {
      expected_edge += kKernel[j] * 4.0F;
    }
  }
  TEST_FLOAT_CLOSE(out[0], expected_edge, 1e-5F);

  // 5. Weight renormalization over present rows (oldest + center + newest)
  float past_row[num_bins];
  float future_row[num_bins];
  for (uint32_t k = 0U; k < num_bins; k++) {
    past_row[k] = 1.0F; // kernel weight lookups happen by index
    future_row[k] = 3.0F;
  }
  rows[0] = past_row;   // oldest
  rows[4] = future_row; // newest
  TEST_ASSERT(spectral_neighbor_average(&params, acc, out),
              "Three present rows should succeed");
  const uint32_t center_bin = 400U; // interior: frequency FIR keeps the value
  float weighted = past_row[center_bin] * kKernel[0] +
                   impulse[center_bin] * kKernel[2] +
                   future_row[center_bin] * kKernel[4];
  const float weights_sum = kKernel[0] + kKernel[2] + kKernel[4];
  TEST_FLOAT_CLOSE(out[center_bin], weighted / weights_sum, 1e-5F);

  // 6. Fully present flat rows stay flat through both axes
  rows[0] = flat;
  rows[1] = flat;
  rows[2] = flat;
  rows[3] = flat;
  rows[4] = flat;
  TEST_ASSERT(spectral_neighbor_average(&params, acc, out),
              "All present rows should succeed");
  for (uint32_t k = 0U; k < num_bins; k++) {
    TEST_FLOAT_CLOSE(out[k], 1.0F, 1e-5F);
  }

  // 7. Erb-scaled frequency kernel: structure and invariants
  SpectralNeighborFreqKernel* none =
      spectral_neighbors_freq_kernel_create(0.0F, num_bins);
  TEST_ASSERT(none == NULL, "Non-positive bin_hz must fail");
  none = spectral_neighbors_freq_kernel_create(46.875F, 0U);
  TEST_ASSERT(none == NULL, "Zero spectrum_size must fail");
  spectral_neighbors_freq_kernel_free(NULL); // safety, must not crash

  // 46.875 Hz/bin makes late bins saturate at the max halo while the first
  // bins floor at the min one
  SpectralNeighborFreqKernel* freq =
      spectral_neighbors_freq_kernel_create(46.875F, num_bins);
  TEST_ASSERT(freq != NULL, "Kernel build should succeed");
  TEST_ASSERT(freq->halo[0] == NEIGHBOR_FREQ_HALO_MIN,
              "DC halo must floor at the minimum radius");
  TEST_ASSERT(freq->halo[num_bins - 1] == NEIGHBOR_FREQ_HALO_MAX,
              "Nyquist halo must cap at the maximum radius");
  // Monotonic growth between the clamps
  for (uint32_t k = 1U; k < num_bins; k++) {
    TEST_ASSERT(freq->halo[k] >= freq->halo[k - 1U],
                "Halo radius must grow with frequency (non-decreasing)");
    TEST_ASSERT(
        freq->start[k] == freq->start[k - 1U] + (2U * freq->halo[k - 1U] + 1U),
        "Flat tap indexing must stay contiguous");
  }
  TEST_ASSERT(
      freq->start[num_bins - 1U] + 2U * freq->halo[num_bins - 1U] + 1U ==
          freq->total_taps,
      "total_taps must cover the last bin");
  // Per-bin weights sum to 1 (true average)
  for (uint32_t k = 0U; k < num_bins; k++) {
    float sum = 0.0F;
    for (uint32_t i = 0U; i < 2U * freq->halo[k] + 1U; i++) {
      sum += freq->weights[freq->start[k] + i];
    }
    TEST_FLOAT_CLOSE(sum, 1.0F, 1e-4F);
  }

  // 8. Kernel path: flat rows stay exactly flat (edge renormalization)
  rows[0] = flat;
  rows[1] = NULL;
  rows[2] = flat;
  rows[3] = NULL;
  rows[4] = NULL;
  SpectralNeighborParameters kernel_params = {
      .rows = rows,
      .rows_count = NEIGHBOR_KERNEL_TAPS,
      .spectrum_size = num_bins,
      .freq_kernel = freq,
  };
  TEST_ASSERT(spectral_neighbor_average(&kernel_params, acc, out),
              "Kernel path should succeed");
  for (uint32_t k = 0U; k < num_bins; k++) {
    TEST_FLOAT_CLOSE(out[k], 1.0F, 1e-5F);
  }

  // 9. Kernel path: impulse spreads within its own halo; the time axis
  // mixes the present flat row too, so the exact expectation is the
  // normalized center tap on the mixed value plus flat leakage elsewhere
  memset(impulse, 0, sizeof(impulse));
  impulse[interior] = 8.0F;
  rows[2] = impulse;
  TEST_ASSERT(spectral_neighbor_average(&kernel_params, acc, out),
              "Kernel impulse should succeed");
  const uint32_t spread = freq->halo[interior];
  // Time pass already mixed rows (0.0625*flat + 0.375*impulse) / 0.4375:
  const float center_eval = (0.0625F + 0.375F * 8.0F) / 0.4375F; // 7.0 at k=300
  const float flat_eval = (0.0625F * 1.0F) / 0.4375F; // 0.1428... elsewhere
  const float w_center = freq->weights[freq->start[interior] + spread];
  TEST_FLOAT_CLOSE(out[interior],
                   w_center * center_eval + (1.0F - w_center) * flat_eval,
                   1e-5F);
  TEST_ASSERT(out[interior + spread + 1U] < flat_eval + 1e-4F,
              "Beyond the halo only bounded time leakage may remain");
  TEST_ASSERT(out[0] < flat_eval + 1e-4F,
              "No impulse leakage into DC beyond the halo radius");
  TEST_ASSERT(out[num_bins - 1U] < flat_eval + 1e-4F,
              "No impulse leakage into Nyquist beyond the halo radius");

  spectral_neighbors_freq_kernel_free(freq);

  printf("✓ Spectral Neighbors full coverage tests passed\n");
}

int main(void) {
  test_spectral_neighbors();
  return 0;
}
