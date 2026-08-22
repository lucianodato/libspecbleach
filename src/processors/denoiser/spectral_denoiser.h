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

#ifndef SPECTRAL_DENOISER_H
#define SPECTRAL_DENOISER_H

#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/spectral_processor.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct DenoiserParameters {
  float reduction_amount;
  bool residual_listen;
  int learn_noise;
  float smoothing_factor;
  float whitening_factor;
  int adaptive_noise;
  int noise_estimation_method; /**< 0=SPP-MMSE, 1=Brandt, 2=Martin MS */
  float masking_depth;
  float suppression_strength;
  float aggressiveness;  /**< -1.0 (Median/Min) to 1.0 (Max), 0.0 (Mean) */
  float tonal_reduction; /**< 0.0 to 1.0 (Phase 3) */
  int hpss_enable;       /**< 0=disabled, 1=enabled */
  float noise_profile_offset_linear; /**< Linear scalar for noise profile */
  const float* reduction_curve_bias; /**< Per-bin dB bias, NULL = disabled */
  bool reduction_curve_enabled;
} DenoiserParameters;

SpectralProcessorHandle spectral_denoiser_initialize(
    uint32_t sample_rate, uint32_t fft_size, uint32_t overlap_factor,
    NoiseProfile* noise_profile);
void spectral_denoiser_free(SpectralProcessorHandle instance);
bool load_reduction_parameters(SpectralProcessorHandle instance,
                               DenoiserParameters parameters);
bool spectral_denoiser_run(SpectralProcessorHandle instance,
                           float* fft_spectrum);
const float* spectral_denoiser_get_tonal_mask(SpectralProcessorHandle instance);
uint32_t spectral_denoiser_get_peaks(SpectralProcessorHandle instance,
                                     float* peak_freqs_hz, uint32_t max_peaks);
const float* spectral_denoiser_get_active_noise_profile(
    SpectralProcessorHandle instance);
void spectral_denoiser_reset_noise_profile(SpectralProcessorHandle instance);
uint32_t spectral_denoiser_get_latency_frames(SpectralProcessorHandle instance);
bool spectral_denoiser_is_transient_detected(SpectralProcessorHandle instance);

#endif
