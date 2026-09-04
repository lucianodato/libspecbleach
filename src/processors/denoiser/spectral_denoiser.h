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
#include "specbleach_denoiser.h"
#include <stdbool.h>
#include <stdint.h>

/**
 * Parameters for the unified spectral denoiser. The smoothing strategy is
 * selected with smoothing_mode and can be changed at runtime (allocation-free,
 * see load_reduction_parameters).
 */
typedef struct DenoiserParameters {
  int learn_noise;        /**< Learning mode: 0=disabled, 1=learn all modes */
  bool residual_listen;   /**< Output residue instead of denoised signal */
  float reduction_amount; /**< Gain floor / reduction amount (linear) */
  int smoothing_mode;     /**< 0=temporal (1D), 1=NLM 2D, 2=NLM 2D + DFTT */
  float smoothing_factor; /**< Temporal smoothing strength / NLM h parameter */
  float whitening_factor; /**< Whitening factor (0.0 to 1.0) */
  int adaptive_noise;     /**< Adaptive noise mode: 0=disabled, 1=enabled */
  int noise_estimation_method; /**< 0=SPP-MMSE, 1=Brandt, 2=Martin MS */
  float masking_depth;        /**< Masking veto/protection depth (0.0 to 1.0) */
  float suppression_strength; /**< Suppression aggressiveness (0.0 to 1.0) */
  float aggressiveness;       /**< -1.0 (Median/Min) to 1.0 (Max), 0.0 (Mean) */
  float tonal_reduction;      /**< 0.0 to 1.0 */
  int hpss_enable;            /**< 0=disabled, 1=enabled */
  float noise_profile_offset_linear; /**< Linear scalar for noise profile */
  float tonal_noise_profile_offset_linear; /**< Linear scalar at tonal bins */
  const float* reduction_curve_bias;       /**< Per-bin dB bias, NULL = off */
  bool reduction_curve_enabled;
} DenoiserParameters;

SpectralProcessorHandle spectral_denoiser_initialize(
    uint32_t sample_rate, uint32_t fft_size, uint32_t overlap_factor,
    NoiseProfile* noise_profile);
/* Preferred: pass the true STFT hop (frame/overlap) so time constants stay
 * frame-rate independent even when the FFT is padded above the frame. */
SpectralProcessorHandle spectral_denoiser_initialize_with_hop(
    uint32_t sample_rate, uint32_t fft_size, uint32_t overlap_factor,
    uint32_t hop_samples, NoiseProfile* noise_profile);
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
float spectral_denoiser_get_transient_intensity(
    SpectralProcessorHandle instance);

#endif
