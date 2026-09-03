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

#ifndef SPECBLEACH_FRAME_RATE_NORM_H
#define SPECBLEACH_FRAME_RATE_NORM_H

#include "shared/configurations.h"
#include "shared/utils/general_utils.h"
#include <math.h>
#include <stdint.h>

/* Frame-rate normalization helpers (init-time only, never RT audio thread).
 * Converts legacy frame-count / per-frame-alpha tuning into seconds so
 * behavior is invariant across STFT frame sizes (23-93ms). */

static inline SB_UNUSED float sb_hop_sec(uint32_t hop_samples,
                                         uint32_t sample_rate) {
  if (hop_samples == 0U || sample_rate == 0U) {
    return 0.0F;
  }
  return (float)hop_samples / (float)sample_rate;
}

static inline SB_UNUSED float sb_bin_hz(uint32_t sample_rate,
                                        uint32_t fft_size) {
  if (fft_size == 0U) {
    return 0.0F;
  }
  return (float)sample_rate / (float)fft_size;
}

/* tau <-> alpha for a one-pole per-hop smoother: alpha = exp(-hop/tau). */
static inline SB_UNUSED float sb_alpha_from_tau(float tau_sec, float hop_sec) {
  if (!(tau_sec > 0.0F) || !(hop_sec > 0.0F)) {
    return 0.0F;
  }
  return expf(-hop_sec / tau_sec);
}

static inline SB_UNUSED float sb_tau_from_alpha(float alpha_ref,
                                                float hop_ref_sec) {
  if (!(alpha_ref > 0.0F) || alpha_ref >= 1.0F || !(hop_ref_sec > 0.0F)) {
    return 0.0F;
  }
  return -hop_ref_sec / logf(alpha_ref);
}

/* alpha_ref tuned at hop_ref -> equivalent alpha at hop_sec. */
static inline SB_UNUSED float sb_alpha_retuned(float alpha_ref,
                                               float hop_ref_sec,
                                               float hop_sec) {
  const float tau = sb_tau_from_alpha(alpha_ref, hop_ref_sec);
  if (!(tau > 0.0F)) {
    return alpha_ref;
  }
  return sb_alpha_from_tau(tau, hop_sec);
}

static inline SB_UNUSED uint32_t sb_frames_for_ms(float target_ms,
                                                  float hop_sec, uint32_t lo,
                                                  uint32_t hi) {
  if (!(target_ms > 0.0F) || !(hop_sec > 0.0F)) {
    return lo;
  }
  const float hop_ms = hop_sec * 1000.0F;
  uint32_t n = (uint32_t)(target_ms / hop_ms + 0.5F);
  if (n < lo) {
    n = lo;
  }
  if (n > hi) {
    n = hi;
  }
  return n;
}

static inline SB_UNUSED uint32_t sb_bins_for_hz(float target_hz, float bin_hz,
                                                uint32_t lo, uint32_t hi) {
  if (!(target_hz > 0.0F) || !(bin_hz > 0.0F)) {
    return lo;
  }
  uint32_t n = (uint32_t)(target_hz / bin_hz + 0.5F);
  if (n < lo) {
    n = lo;
  }
  if (n > hi) {
    n = hi;
  }
  return n;
}

typedef struct SbNlmGeometry {
  uint32_t patch;
  uint32_t past;
  uint32_t future;
  uint32_t search_freq;
  uint32_t paste;
} SbNlmGeometry;

/* Listen-validated per-option table for the 5 plugin frame sizes (time axis).
 * Values == round(ms_target / hop) with symmetric 128ms search context;
 * frozen so init is deterministic and reviewable. Frequency axis is always
 * recomputed from Hz (SR/FFT dependent) and is not tabled. */
static inline SB_UNUSED SbNlmGeometry
sb_nlm_geometry_for_frame_ms(float frame_ms, float hop_sec, float bin_hz) {
  SbNlmGeometry g = {NLM_PATCH_SIZE, NLM_SEARCH_RANGE_TIME_PAST,
                     NLM_SEARCH_RANGE_TIME_FUTURE, NLM_SEARCH_RANGE_FREQ,
                     NLM_PASTE_BLOCK_SIZE};
  /* Known plugin options: match within 0.6ms tolerance. */
  if (frame_ms > 22.4F && frame_ms < 23.6F) {
    g.patch = 16U;
    g.past = 22U;
    g.future = 22U;
  } else if (frame_ms > 31.4F && frame_ms < 32.6F) {
    g.patch = 12U;
    g.past = 16U;
    g.future = 16U;
  } else if (frame_ms > 45.4F && frame_ms < 46.6F) {
    g.patch = 8U;
    g.past = 11U;
    g.future = 11U;
  } else if (frame_ms > 63.4F && frame_ms < 64.6F) {
    g.patch = 6U;
    g.past = 8U;
    g.future = 8U;
  } else if (frame_ms > 92.4F && frame_ms < 93.6F) {
    g.patch = 4U;
    g.past = 6U;
    g.future = 6U;
  } else if (hop_sec > 0.0F) {
    g.patch = sb_frames_for_ms(NLM_PATCH_TIME_MS, hop_sec, 4U, 16U);
    g.past = sb_frames_for_ms(NLM_SEARCH_PAST_MS, hop_sec, 8U, 32U);
    g.future = sb_frames_for_ms(NLM_SEARCH_FUTURE_MS, hop_sec, 8U, 32U);
  }
  if (bin_hz > 0.0F) {
    g.search_freq = sb_bins_for_hz(NLM_SEARCH_FREQ_HZ, bin_hz, 2U, 32U);
    g.paste = sb_bins_for_hz(NLM_PASTE_FREQ_HZ, bin_hz, 2U, 16U);
  }
  /* Paste must stay <= patch for the block algorithm; keep SIMD sizes. */
  if (g.paste > g.patch) {
    g.paste = g.patch;
  }
  return g;
}

#endif /* SPECBLEACH_FRAME_RATE_NORM_H */
