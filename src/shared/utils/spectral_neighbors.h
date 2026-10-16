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

#ifndef SHARED_UTILS_SPECTRAL_NEIGHBORS_H
#define SHARED_UTILS_SPECTRAL_NEIGHBORS_H

#include "shared/configurations.h"
#include <math.h>
#include <stdbool.h>
#include <stdint.h>

/**
 * @brief Convert a complex FFT spectrum to a power row (one bin per band).
 *
 * Same domain convention as the shared spectral features (DC and Nyquist are
 * purely real); fft_spectrum uses the interleaved layout of the library FFT
 * (imag mirrored at fft_spectrum[n - k]).
 */
void spectral_neighbors_compute_power_row(float* power_row,
                                          const float* fft_spectrum,
                                          uint32_t fft_size);

/**
 * @brief Erb-scaled frequency-axis kernel (build once at init, RT-safe use).
 *
 * For every bin the frequency halo radius follows the critical band width:
 * radius = clamp(round(NEIGHBOR_ERB_HALO_SHARE * ERB(center) / bin_hz),
 * NEIGHBOR_FREQ_HALO_MIN, NEIGHBOR_FREQ_HALO_MAX) bins around the center.
 * Taps per bin form a Hann taper normalized so their sum is 1, and edge
 * clamping renormalizes over the taps actually used: the output stays a
 * scale-preserving average at every radius, including DC and Nyquist.
 * Setup-only; the built kernel is read-only afterwards.
 */
typedef struct SpectralNeighborFreqKernel {
  uint32_t spectrum_size;
  uint32_t total_taps; /**< Total stored tap values (weights[] length) */
  uint32_t* halo;      /**< Per-bin frequency halo radius (bins) */
  uint32_t* start;     /**< Per-bin index of its first tap in weights[] */
  float* weights;      /**< Flat per-bin Hann tapers, each summing to 1 */
} SpectralNeighborFreqKernel;

/**
 * @brief Build the Erb-scaled frequency kernel. Setup-only: allocates.
 *
 * @param bin_hz Frequency of one bin in Hz (> 0)
 * @param spectrum_size Number of bins per row (> 0)
 * @return Kernel instance or NULL on invalid parameters / allocation failure
 */
SpectralNeighborFreqKernel* spectral_neighbors_freq_kernel_create(
    float bin_hz, uint32_t spectrum_size);

/** Frees a kernel built by spectral_neighbors_freq_kernel_create (NULL-safe).
 */
void spectral_neighbors_freq_kernel_free(SpectralNeighborFreqKernel* self);

/**
 * @brief Neighborhood average parameters.
 *
 * rows are time ordered (oldest first) and must hold an odd count; the
 * center entry describes the frame being processed. NULL entries (e.g.
 * look-ahead rows not yet filled) are skipped and the kernel weights are
 * renormalized over the present rows.
 */
typedef struct SpectralNeighborParameters {
  const float** rows;     /**< Time ordered power rows, NULL entries skipped */
  uint32_t rows_count;    /**< Odd number of entries, <= NEIGHBOR_KERNEL_TAPS */
  uint32_t spectrum_size; /**< Number of bins per row */
  /**
   * Optional Erb-scaled frequency kernel. NULL renders the fixed 5-tap
   * Hann neighborhood with frequency radius NEIGHBOR_KERNEL_HALO bins.
   */
  const SpectralNeighborFreqKernel* freq_kernel;
} SpectralNeighborParameters;

/**
 * @brief Separable time-frequency neighborhood average of power rows.
 *
 * Time axis: weighted accumulation over the present rows followed by
 * reciprocal normalization. Frequency axis: when parameters->freq_kernel is
 * set, per-bin Erb-scaled variable-radius averaging; otherwise a symmetric
 * fixed 5-tap FIR. Edge taps are clamped; the kernel path renormalizes the
 * clamped taps. All buffers are caller owned; acc_scratch holds the
 * time-averaged row on return and out receives the averaged row.
 * Real-time safe: no allocation.
 *
 * @param parameters Row bundle (see SpectralNeighborParameters)
 * @param acc_scratch Intermediate buffer, spectrum_size floats
 * @param out Output buffer, spectrum_size floats
 * @return false on invalid parameters or when no row is present
 */
bool spectral_neighbor_average(const SpectralNeighborParameters* parameters,
                               float* acc_scratch, float* out);

#endif // SHARED_UTILS_SPECTRAL_NEIGHBORS_H
