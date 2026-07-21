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

#include <math.h>

#include "shared/denoiser_logic/core/denoiser_post_process.h"
#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/processing/tonal_reducer.h"
#include "shared/utils/spectral_utils.h"

static float get_bin_magnitude(const float* fft_spectrum, uint32_t k, uint32_t fft_size) {
  if (k == 0U) {
    return fabsf(fft_spectrum[0]);
  } else if (k == fft_size / 2U && (fft_size % 2U == 0)) {
    return fabsf(fft_spectrum[k]);
  } else {
    float real = fft_spectrum[k];
    float imag = fft_spectrum[fft_size - k];
    return sqrtf(real * real + imag * imag);
  }
}

void denoiser_post_process_apply(DenoiserPostProcessParams params) {
  // 1. Apply noise floor management
  const float* tonal_mask = tonal_reducer_get_mask(params.tonal_reducer);
  if (params.noise_floor_manager && params.noise_spectrum) {
    noise_floor_manager_apply(
        params.noise_floor_manager, params.real_spectrum_size, params.fft_size,
        params.gain_spectrum, params.noise_spectrum, params.reduction_amount,
        params.tonal_reduction, tonal_mask, params.whitening_factor);

    // 1.1. Spectral Gain Patching (Tonal Gain Interpolation)
    //      For bins where tonal hum noise was aggressively reduced, we interpolate the target
    //      gains from adjacent non-tonal bins to avoid notches/holes in the spectrum.
    //      - In silence, the gain at the hum bin is matched to the active surrounding broadband
    //        noise magnitude to preserve the dynamically moving, grainy texture of the remaining noise.
    //      - During signal playback, the gain is restored to the neighbor gain level to protect
    //        the wanted audio.
    if (tonal_mask) {
      for (uint32_t k = 0U; k < params.real_spectrum_size; k++) {
        if (tonal_mask[k] > 0.0f) {
          // Find nearest non-tonal neighbor on the left
          int left_idx = -1;
          for (int idx = (int)k - 1; idx >= 0; idx--) {
            if (tonal_mask[idx] == 0.0f) {
              left_idx = idx;
              break;
            }
          }

          // Find nearest non-tonal neighbor on the right
          int right_idx = -1;
          for (int idx = (int)k + 1; idx < (int)params.real_spectrum_size; idx++) {
            if (tonal_mask[idx] == 0.0f) {
              right_idx = idx;
              break;
            }
          }

          // Calculate the interpolated total gain and dynamic broadband noise magnitude from neighbors
          float g_floor = params.gain_spectrum[k];
          float s_in_broadband = 0.0f;
          bool has_neighbor = false;

          float s_left = (left_idx >= 0) ? get_bin_magnitude(params.fft_spectrum, left_idx, params.fft_size) : 0.0f;
          float s_right = (right_idx >= 0) ? get_bin_magnitude(params.fft_spectrum, right_idx, params.fft_size) : 0.0f;

          if (left_idx >= 0 && right_idx >= 0) {
            float t = (float)(k - left_idx) / (float)(right_idx - left_idx);
            g_floor = params.gain_spectrum[left_idx] +
                      t * (params.gain_spectrum[right_idx] - params.gain_spectrum[left_idx]);
            s_in_broadband = s_left + t * (s_right - s_left);
            has_neighbor = true;
          } else if (left_idx >= 0) {
            g_floor = params.gain_spectrum[left_idx];
            s_in_broadband = s_left;
            has_neighbor = true;
          } else if (right_idx >= 0) {
            g_floor = params.gain_spectrum[right_idx];
            s_in_broadband = s_right;
            has_neighbor = true;
          }

          // Compute signal magnitude at bin k
          float s_in = get_bin_magnitude(params.fft_spectrum, k, params.fft_size);

          // Compute signal-presence weight (w_sig) based on the neighbor gains
          float g_base = params.reduction_amount;
          float g_max = 0.0f;
          if (left_idx >= 0) {
            g_max = fmaxf(g_max, params.gain_spectrum[left_idx]);
          }
          if (right_idx >= 0) {
            g_max = fmaxf(g_max, params.gain_spectrum[right_idx]);
          }

          float w_sig = 0.0f;
          if (g_max > g_base) {
            w_sig = (g_max - g_base) / 0.2f;
            if (w_sig > 1.0f) {
              w_sig = 1.0f;
            }
          }

          // Calculate dynamic noise target ratio
          float ratio = 1.0f;
          if (s_in > 1e-12f) {
            ratio = s_in_broadband / s_in;
            if (ratio > 1.0f) {
              ratio = 1.0f;
            }
          }

          if (has_neighbor) {
            float target_factor = (1.0f - w_sig) * ratio + w_sig;
            params.gain_spectrum[k] = g_floor * target_factor;
          }
        }
      }

      // Apply symmetric copy of the patched gains to the upper half of the
      // spectrum (gain_spectrum size is fft_size)
      for (uint32_t k = 1U; k < params.fft_size - k; k++) {
        params.gain_spectrum[params.fft_size - k] = params.gain_spectrum[k];
      }
    }
  }

  // 2. Mixing Logic (formerly denoise_mixer)
  // We perform mixing in-place or with stack variables avoiding extra heap
  // allocation if possible However, to calculate residual we need the denoised
  // version and original. params.fft_spectrum holds the original.

  // We can optimize this by calculating per-bin.

  if (params.residual_listen) {
    // Residual Listen Mode: Output = Original - Denoised
    // Original - (Original * Gain) = Original * (1 - Gain)
    for (uint32_t k = 0U; k < params.fft_size; k++) {
      float original = params.fft_spectrum[k];
      float gain = params.gain_spectrum[k];
      float denoised = original * gain;
      params.fft_spectrum[k] = original - denoised;
    }
  } else {
    // Normal Mode: Output = Original * Gain
    for (uint32_t k = 0U; k < params.fft_size; k++) {
      params.fft_spectrum[k] *= params.gain_spectrum[k];
    }
  }
}
