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

#include "minimum_phase_fir.h"
#include "shared/configurations.h"
#include <math.h>
#include <string.h>

bool minimum_phase_fir_taper_window(const uint32_t taps, float* window) {
  if (taps == 0U || !window) {
    return false;
  }

  const uint32_t taper_start =
      (uint32_t)((float)taps * DSAF_ASYMMETRIC_TAPER_START);

  for (uint32_t i = 0U; i < taps; i++) {
    if (i < taper_start || taper_start + 1U >= taps) {
      window[i] = 1.0F;
    } else {
      // Raised-cosine decay across the tail to smooth truncation
      const float t = (float)(i - taper_start) / (float)(taps - taper_start);
      window[i] = 0.5F * (1.0F + cosf(M_PIf * t));
    }
  }

  return true;
}

bool minimum_phase_fir_synthesize(const uint32_t fft_size, const uint32_t taps,
                                  const float* gain_bins, FftTransform* fft,
                                  const float* taper_window, float* out_ir) {
  if (fft_size == 0U || taps == 0U || taps > fft_size || !gain_bins || !fft ||
      !taper_window || !out_ir) {
    return false;
  }

  const uint32_t n = fft_size;
  const uint32_t half = n / 2U;

  // Step 1: log-magnitude spectrum in halfcomplex layout. The target gain
  // spectrum is real, so imaginary parts (out[n - k]) are zero.
  float* out = get_fft_output_buffer(fft);
  out[0] = logf(fmaxf(gain_bins[0], DSAF_CEPSTRAL_EPSILON));
  out[half] = logf(fmaxf(gain_bins[half], DSAF_CEPSTRAL_EPSILON));
  for (uint32_t k = 1U; k < half; k++) {
    out[k] = logf(fmaxf(gain_bins[k], DSAF_CEPSTRAL_EPSILON));
    out[n - k] = 0.0F;
  }

  // Step 2: real cepstrum (backward real FFT is unnormalized -> scale 1/N)
  if (!compute_backward_fft(fft)) {
    return false;
  }
  float* in = get_fft_input_buffer(fft);
  const float scale = 1.0F / (float)n;
  for (uint32_t i = 0U; i < n; i++) {
    in[i] *= scale;
  }

  // Step 3: minimum-phase cepstral fold (zero negative quefrency)
  for (uint32_t i = 1U; i < half; i++) {
    in[i] = 2.0F * in[i];
    in[n - i] = 0.0F;
  }

  // Step 4: forward FFT -> complex logarithm (halfcomplex)
  if (!compute_forward_fft(fft)) {
    return false;
  }

  // Step 5: exp in place -> minimum-phase complex spectrum W[k]
  out[0] = expf(out[0]);
  out[half] = expf(out[half]);
  for (uint32_t k = 1U; k < half; k++) {
    const float cre = out[k];
    const float cim = out[n - k];
    const float mag = expf(cre);
    out[k] = mag * cosf(cim);
    out[n - k] = mag * sinf(cim);
  }

  // Step 6: inverse FFT -> time-domain impulse response (scale 1/N)
  if (!compute_backward_fft(fft)) {
    return false;
  }

  // Step 7: truncate to M taps with the precomputed asymmetric taper
  // (backward real FFT is unnormalized -> scale 1/N)
  const float* ir = get_fft_input_buffer(fft);
  for (uint32_t i = 0U; i < taps; i++) {
    out_ir[i] = ir[i] * scale * taper_window[i];
  }

  return true;
}
