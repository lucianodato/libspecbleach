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

#include "shared/denoiser_logic/processing/dftt_filter.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"
#include "shared/utils/simd_utils.h"

#ifndef DFTT_TWO_PI
#define DFTT_TWO_PI (6.283185307179586F)
#endif

/* ponytail: DFT cost profile. Tile sizes are <=32x16 and every output frame
 * recomputes ~spec/hop tiles, so the naive O(N^2) DFT dominated the stage
 * (~4x NLM cost measured in-host). Now: radix-2 FFT dispatch per dimension
 * when the size is a power of two (test/host configs hit bf=32, bt=16), the
 * naive SIMD DFT stays as fallback for odd Hz/ms-derived sizes. The inverse
 * collapses the time axis straight to the single emitted row index (only the
 * newest row is ever reconstructed), cutting it to bt*bf MACs + one 1D
 * transform. Gather/emit run through the sb_vec8 kit on the clamp-free
 * interior. Float-level results differ from the naive path only through
 * summation reassociation (~1e-7, same class as the existing SIMD DFT path);
 * the old-vs-new A/B harness pins the exact delta. */
#define DFTT_MAX_DIM (32U)

struct DfttFilter {
  uint32_t spectrum_size;
  uint32_t time_span;
  uint32_t block_freq;
  uint32_t block_hop;
  float kill_k;
  float** noisy_ring;
  float** smooth_ring;
  uint32_t head;
  uint32_t filled;
  bool fresh;
  float* win_freq;
  float* win_time;
  float* cos_freq;
  float* sin_freq;
  float* cos_time;
  float* sin_time;
  uint32_t* br_freq;
  uint32_t* br_time;
  float* mod_cos;
  float* mod_sin;
  bool pow2_freq;
  bool pow2_time;
  float* wsum;
  float* tile_re;
  float* tile_im;
  float* ref_re;
  float* ref_im;
};

static uint32_t dftt_clamp_bin(int32_t bin, uint32_t size) {
  if (bin < 0) {
    return 0U;
  }
  if ((uint32_t)bin >= size) {
    return size - 1U;
  }
  return (uint32_t)bin;
}

/* Periodic Hann (STFT-style) floored so overlap-add weights stay strictly
 * positive everywhere, including the DC bin. Exact reconstruction under
 * unity gain is preserved by win^2 normalization. */
#define DFTT_WIN_FLOOR (1e-3F)

static void dftt_build_hann(float* win, uint32_t n) {
  if (n == 0U) {
    return;
  }
  if (n == 1U) {
    win[0] = 1.0F;
    return;
  }
  for (uint32_t i = 0U; i < n; i++) {
    const float w = 0.5F - (0.5F * cosf((DFTT_TWO_PI * (float)i) / (float)n));
    win[i] = w < DFTT_WIN_FLOOR ? DFTT_WIN_FLOOR : w;
  }
}

/* Rising half-Hann for the causal time axis: 0 at the oldest row, 1.0 at the
 * newest (emitted) row. A symmetric window would taper the emitted row to
 * near zero, and threshold redistribution from the tile interior would then
 * be amplified by win^2 normalization. Recency weighting is also the right
 * prior for a past-only tile. */
static void dftt_build_time_window(float* win, uint32_t n) {
  if (n == 0U) {
    return;
  }
  if (n == 1U) {
    win[0] = 1.0F;
    return;
  }
  for (uint32_t i = 0U; i < n; i++) {
    const float w = 0.5F - (0.5F * cosf((M_PIf * (float)i) / (float)(n - 1U)));
    win[i] = w < DFTT_WIN_FLOOR ? DFTT_WIN_FLOOR : w;
  }
}

static void dftt_build_twiddles(float* cos_tab, float* sin_tab, uint32_t n) {
  for (uint32_t k = 0U; k < n; k++) {
    for (uint32_t m = 0U; m < n; m++) {
      const float ph = ((DFTT_TWO_PI * (float)k) * (float)m) / (float)n;
      cos_tab[(k * n) + m] = cosf(ph);
      sin_tab[(k * n) + m] = sinf(ph);
    }
  }
}

static bool dftt_is_pow2(uint32_t n) {
  return (n != 0U) && ((n & (n - 1U)) == 0U);
}

static uint32_t dftt_bit_reverse(uint32_t x, uint32_t bits) {
  uint32_t r = 0U;
  for (uint32_t i = 0U; i < bits; i++) {
    r = (r << 1U) | (x & 1U);
    x >>= 1U;
  }
  return r;
}

/* In-place 1D complex DFT over strided data. inverse flips the sine sign and
 * scales by 1/n. The stride-1 path (frequency rows, ~2/3 of the work) runs
 * 8-wide MACs through the shared sb_vec8 kit; strided time columns stay
 * scalar. Bit-level results differ from scalar reassociation (~1e-7). */
static void dftt_dft_1d(float* re, float* im, uint32_t stride, uint32_t n,
                        const float* cos_tab, const float* sin_tab,
                        bool inverse) {
  float tmp_re[DFTT_MAX_DIM];
  float tmp_im[DFTT_MAX_DIM];
  for (uint32_t k = 0U; k < n; k++) {
    float sr = 0.0F;
    float si = 0.0F;
    uint32_t m = 0U;
    if (stride == 1U) {
      const sb_vec8_t vzero = sb_set8(0.0F);
      sb_vec8_t acc_rc = vzero;
      sb_vec8_t acc_is = vzero;
      sb_vec8_t acc_ic = vzero;
      sb_vec8_t acc_rs = vzero;
      for (; m + 8U <= n; m += 8U) {
        sb_vec8_t xr = sb_load8(re + m);
        sb_vec8_t xi = sb_load8(im + m);
        sb_vec8_t c = sb_load8(cos_tab + ((size_t)k * n) + m);
        sb_vec8_t s = sb_load8(sin_tab + ((size_t)k * n) + m);
        if (inverse) {
          s = sb_sub8(vzero, s);
        }
        acc_rc = sb_add8(acc_rc, sb_mul8(xr, c));
        acc_is = sb_add8(acc_is, sb_mul8(xi, s));
        acc_ic = sb_add8(acc_ic, sb_mul8(xi, c));
        acc_rs = sb_add8(acc_rs, sb_mul8(xr, s));
      }
      sr = sb_vec8_hsum(acc_rc) + sb_vec8_hsum(acc_is);
      si = sb_vec8_hsum(acc_ic) - sb_vec8_hsum(acc_rs);
    }
    for (; m < n; m++) {
      const float c = cos_tab[((size_t)k * n) + m];
      float s = sin_tab[((size_t)k * n) + m];
      if (inverse) {
        s = -s;
      }
      const float xr = re[((size_t)m * stride)];
      const float xi = im[((size_t)m * stride)];
      sr += (xr * c) + (xi * s);
      si += (xi * c) - (xr * s);
    }
    tmp_re[k] = sr;
    tmp_im[k] = si;
  }
  const float scale = inverse ? (1.0F / (float)n) : 1.0F;
  for (uint32_t k = 0U; k < n; k++) {
    re[((size_t)k * stride)] = tmp_re[k] * scale;
    im[((size_t)k * stride)] = tmp_im[k] * scale;
  }
}

/* Iterative radix-2 Cooley-Tukey (DIT) over strided data. Twiddles reuse row
 * 1 of the shared table (e^{-2pi i j / n}); inverse flips the sine sign and
 * scales by 1/n, matching the naive DFT conventions. */
static void dftt_fft_1d(float* re, float* im, uint32_t stride, uint32_t n,
                        const uint32_t* br, const float* cos_row,
                        const float* sin_row, bool inverse) {
  for (uint32_t i = 0U; i < n; i++) {
    const uint32_t j = br[i];
    if (j > i) {
      const size_t a = (size_t)i * stride;
      const size_t b = (size_t)j * stride;
      float t = re[a];
      re[a] = re[b];
      re[b] = t;
      t = im[a];
      im[a] = im[b];
      im[b] = t;
    }
  }
  for (uint32_t len = 2U; len <= n; len <<= 1U) {
    const uint32_t half = len >> 1U;
    const uint32_t step = n / len;
    for (uint32_t base = 0U; base < n; base += len) {
      for (uint32_t j = 0U; j < half; j++) {
        const size_t tw = (size_t)j * step;
        const float c = cos_row[tw];
        const float s = inverse ? sin_row[tw] : -sin_row[tw];
        const size_t p = (size_t)(base + j) * stride;
        const size_t q = (size_t)(base + j + half) * stride;
        const float tr = (re[q] * c) - (im[q] * s);
        const float ti = (re[q] * s) + (im[q] * c);
        re[q] = re[p] - tr;
        im[q] = im[p] - ti;
        re[p] += tr;
        im[p] += ti;
      }
    }
  }
  if (inverse) {
    const float scale = 1.0F / (float)n;
    for (uint32_t i = 0U; i < n; i++) {
      const size_t a = (size_t)i * stride;
      re[a] *= scale;
      im[a] *= scale;
    }
  }
}

static void dftt_fwd_rows(DfttFilter* f, float* re, float* im) {
  const uint32_t bf = f->block_freq;
  for (uint32_t r = 0U; r < f->time_span; r++) {
    float* row_r = re + ((size_t)r * bf);
    float* row_i = im + ((size_t)r * bf);
    if (f->pow2_freq) {
      dftt_fft_1d(row_r, row_i, 1U, bf, f->br_freq, f->cos_freq + bf,
                  f->sin_freq + bf, false);
    } else {
      dftt_dft_1d(row_r, row_i, 1U, bf, f->cos_freq, f->sin_freq, false);
    }
  }
}

static void dftt_fwd_cols(DfttFilter* f, float* re, float* im) {
  const uint32_t bf = f->block_freq;
  const uint32_t bt = f->time_span;
  for (uint32_t c = 0U; c < bf; c++) {
    if (f->pow2_time) {
      dftt_fft_1d(re + c, im + c, bf, bt, f->br_time, f->cos_time + bt,
                  f->sin_time + bt, false);
    } else {
      dftt_dft_1d(re + c, im + c, bf, bt, f->cos_time, f->sin_time, false);
    }
  }
}

/* Inverse restricted to the emitted (newest) time row: collapse the time axis
 * directly at output index bt-1 (one sum per quefrency column —
 * e^{+2pi i (bt-1) r / bt} reduces to e^{-2pi i r / bt}), then one inverse
 * transform along quefrency frequency. Same separable math as the full 2D
 * inverse; costs bt*bf MACs + T(bf) instead of bt*T(bf) + bf*T(bt). */
static void dftt_inv_last_row(DfttFilter* f, const float* re, const float* im,
                              float* out_re, float* out_im) {
  const uint32_t bf = f->block_freq;
  const uint32_t bt = f->time_span;
  const float inv_bt = 1.0F / (float)bt;
  for (uint32_t k = 0U; k < bf; k++) {
    float sr = 0.0F;
    float si = 0.0F;
    for (uint32_t r = 0U; r < bt; r++) {
      const float xr = re[((size_t)r * bf) + k];
      const float xi = im[((size_t)r * bf) + k];
      const float mc = f->mod_cos[r];
      const float ms = f->mod_sin[r];
      sr += (xr * mc) + (xi * ms);
      si += (xi * mc) - (xr * ms);
    }
    out_re[k] = sr * inv_bt;
    out_im[k] = si * inv_bt;
  }
  if (f->pow2_freq) {
    dftt_fft_1d(out_re, out_im, 1U, bf, f->br_freq, f->cos_freq + bf,
                f->sin_freq + bf, true);
  } else {
    dftt_dft_1d(out_re, out_im, 1U, bf, f->cos_freq, f->sin_freq, true);
  }
}

DfttFilter* dftt_filter_initialize(uint32_t spectrum_size,
                                   uint32_t time_span_frames,
                                   uint32_t block_freq_bins) {
  if (spectrum_size == 0U || time_span_frames == 0U || block_freq_bins == 0U) {
    return NULL;
  }
  if (time_span_frames > DFTT_MAX_TIME_FRAMES) {
    time_span_frames = DFTT_MAX_TIME_FRAMES;
  }
  if (block_freq_bins > DFTT_MAX_DIM) {
    block_freq_bins = DFTT_MAX_DIM;
  }
  /* Fixed-size stack temporaries cap both DFT dimensions. */
  if (time_span_frames > DFTT_MAX_DIM) {
    return NULL;
  }

  DfttFilter* f = (DfttFilter*)calloc(1U, sizeof(DfttFilter));
  if (!f) {
    return NULL;
  }
  f->spectrum_size = spectrum_size;
  f->time_span = time_span_frames;
  f->block_freq = block_freq_bins;
  f->block_hop = block_freq_bins / DFTT_FREQ_OVERLAP;
  f->kill_k = DFTT_KILL_K;
  f->pow2_freq = dftt_is_pow2(block_freq_bins);
  f->pow2_time = dftt_is_pow2(time_span_frames);
  if (f->block_hop == 0U) {
    f->block_hop = 1U;
  }

  const uint32_t bt = f->time_span;
  const uint32_t bf = f->block_freq;
  const size_t tile = (size_t)bt * bf;

  f->noisy_ring = (float**)calloc(bt, sizeof(float*));
  f->smooth_ring = (float**)calloc(bt, sizeof(float*));
  f->win_freq = (float*)calloc(bf, sizeof(float));
  f->win_time = (float*)calloc(bt, sizeof(float));
  f->cos_freq = (float*)calloc((size_t)bf * bf, sizeof(float));
  f->sin_freq = (float*)calloc((size_t)bf * bf, sizeof(float));
  f->cos_time = (float*)calloc((size_t)bt * bt, sizeof(float));
  f->sin_time = (float*)calloc((size_t)bt * bt, sizeof(float));
  f->mod_cos = (float*)calloc(bt, sizeof(float));
  f->mod_sin = (float*)calloc(bt, sizeof(float));
  f->wsum = (float*)calloc(spectrum_size, sizeof(float));
  f->tile_re = (float*)calloc(tile, sizeof(float));
  f->tile_im = (float*)calloc(tile, sizeof(float));
  f->ref_re = (float*)calloc(tile, sizeof(float));
  f->ref_im = (float*)calloc(tile, sizeof(float));
  if (f->pow2_freq) {
    f->br_freq = (uint32_t*)calloc(bf, sizeof(uint32_t));
  }
  if (f->pow2_time) {
    f->br_time = (uint32_t*)calloc(bt, sizeof(uint32_t));
  }
  if (!f->noisy_ring || !f->smooth_ring || !f->win_freq || !f->win_time ||
      !f->cos_freq || !f->sin_freq || !f->cos_time || !f->sin_time ||
      !f->mod_cos || !f->mod_sin || !f->wsum || !f->tile_re || !f->tile_im ||
      !f->ref_re || !f->ref_im || (f->pow2_freq && !f->br_freq) ||
      (f->pow2_time && !f->br_time)) {
    dftt_filter_free(f);
    return NULL;
  }
  for (uint32_t r = 0U; r < bt; r++) {
    f->noisy_ring[r] = (float*)calloc(spectrum_size, sizeof(float));
    f->smooth_ring[r] = (float*)calloc(spectrum_size, sizeof(float));
    if (!f->noisy_ring[r] || !f->smooth_ring[r]) {
      dftt_filter_free(f);
      return NULL;
    }
  }

  dftt_build_hann(f->win_freq, bf);
  dftt_build_time_window(f->win_time, bt);
  dftt_build_twiddles(f->cos_freq, f->sin_freq, bf);
  dftt_build_twiddles(f->cos_time, f->sin_time, bt);
  for (uint32_t r = 0U; r < bt; r++) {
    const float ph = (DFTT_TWO_PI * (float)r) / (float)bt;
    f->mod_cos[r] = cosf(ph);
    f->mod_sin[r] = sinf(ph);
  }
  if (f->pow2_freq) {
    uint32_t bits = 0U;
    while ((1U << bits) < bf) {
      bits++;
    }
    for (uint32_t i = 0U; i < bf; i++) {
      f->br_freq[i] = dftt_bit_reverse(i, bits);
    }
  }
  if (f->pow2_time) {
    uint32_t bits = 0U;
    while ((1U << bits) < bt) {
      bits++;
    }
    for (uint32_t i = 0U; i < bt; i++) {
      f->br_time[i] = dftt_bit_reverse(i, bits);
    }
  }

  /* Overlap-add normalization: only the newest time row is emitted, so each
   * tile contributes syn*ana = win^2 weights. Tiles are centered (first tile
   * starts at -(BF-HF)) so edge bins are covered by a full-weight tile
   * instead of a near-zero window skirt that normalization would blow up. */
  const float wt = f->win_time[bt - 1U];
  const int32_t tile_start = -((int32_t)bf - (int32_t)f->block_hop);
  for (int32_t fs = tile_start; fs < (int32_t)spectrum_size;
       fs += (int32_t)f->block_hop) {
    for (uint32_t i = 0U; i < bf; i++) {
      const uint32_t bin = dftt_clamp_bin(fs + (int32_t)i, spectrum_size);
      const float w = wt * f->win_freq[i];
      f->wsum[bin] += w * w;
    }
  }

  return f;
}

void dftt_filter_free(DfttFilter* f) {
  if (!f) {
    return;
  }
  if (f->noisy_ring) {
    for (uint32_t r = 0U; r < f->time_span; r++) {
      free(f->noisy_ring[r]);
    }
    free((void*)f->noisy_ring);
  }
  if (f->smooth_ring) {
    for (uint32_t r = 0U; r < f->time_span; r++) {
      free(f->smooth_ring[r]);
    }
    free((void*)f->smooth_ring);
  }
  free(f->win_freq);
  free(f->win_time);
  free(f->cos_freq);
  free(f->sin_freq);
  free(f->cos_time);
  free(f->sin_time);
  free(f->br_freq);
  free(f->br_time);
  free(f->mod_cos);
  free(f->mod_sin);
  free(f->wsum);
  free(f->tile_re);
  free(f->tile_im);
  free(f->ref_re);
  free(f->ref_im);
  free(f);
}

void dftt_filter_set_strength(DfttFilter* f, float strength) {
  if (!f || strength <= 0.0F) {
    return;
  }
  f->kill_k = DFTT_KILL_K * fminf(strength, DFTT_STRENGTH_MAX);
}

void dftt_filter_push(DfttFilter* f, const float* noisy_snr,
                      const float* smoothed_snr) {
  if (!f || !noisy_snr || !smoothed_snr) {
    return;
  }
  memcpy(f->noisy_ring[f->head], noisy_snr, f->spectrum_size * sizeof(float));
  memcpy(f->smooth_ring[f->head], smoothed_snr,
         f->spectrum_size * sizeof(float));
  f->head = (f->head + 1U) % f->time_span;
  f->filled++;
  f->fresh = true;
}

bool dftt_filter_is_ready(DfttFilter* f) {
  return f && f->filled >= f->time_span;
}

void dftt_filter_reset(DfttFilter* f) {
  if (!f) {
    return;
  }
  for (uint32_t r = 0U; r < f->time_span; r++) {
    memset(f->noisy_ring[r], 0, f->spectrum_size * sizeof(float));
    memset(f->smooth_ring[r], 0, f->spectrum_size * sizeof(float));
  }
  f->head = 0U;
  f->filled = 0U;
  f->fresh = false;
}

bool dftt_filter_process(DfttFilter* f, float* refined_snr) {
  if (!f || !refined_snr) {
    return false;
  }
  if (!f->fresh || f->filled < f->time_span) {
    return false;
  }
  f->fresh = false;

  sb_simd_state_t old_simd_state = sb_simd_enable_ftz_daz();

  const uint32_t spec = f->spectrum_size;
  const uint32_t bt = f->time_span;
  const uint32_t bf = f->block_freq;
  const sb_vec8_t vzero = sb_set8(0.0F);
  float crow_re[DFTT_MAX_DIM];
  float crow_im[DFTT_MAX_DIM];

  memset(refined_snr, 0, spec * sizeof(float));
  const float* newest_smooth = f->smooth_ring[((f->head + bt) - 1U) % bt];

  const int32_t tile_start = -((int32_t)bf - (int32_t)f->block_hop);
  for (int32_t fs = tile_start; fs < (int32_t)spec;
       fs += (int32_t)f->block_hop) {
    /* Gather past-only tile, oldest row first. Interior tiles need no bin
     * clamping and run through the 8-wide kit; edge tiles stay scalar. */
    const bool clamped = (fs < 0) || ((fs + (int32_t)bf) > (int32_t)spec);
    float esum = 0.0F;
    float wsum_r2 = 0.0F;
    for (uint32_t r = 0U; r < bt; r++) {
      const float* row_n = f->noisy_ring[(f->head + r) % bt];
      const float* row_s = f->smooth_ring[(f->head + r) % bt];
      const float wt = f->win_time[r];
      const size_t row_off = (size_t)r * bf;
      uint32_t i = 0U;
      if (!clamped) {
        const sb_vec8_t vw = sb_set8(wt);
        for (; i + 8U <= bf; i += 8U) {
          const sb_vec8_t vn = sb_load8(row_n + fs + i);
          const sb_vec8_t vs = sb_load8(row_s + fs + i);
          esum += sb_vec8_hsum(vn);
          const sb_vec8_t w = sb_mul8(sb_load8(f->win_freq + i), vw);
          const sb_vec8_t rd = sb_sub8(vn, vs);
          wsum_r2 += sb_vec8_hsum(sb_mul8(sb_mul8(w, w), sb_mul8(rd, rd)));
          sb_store8(f->tile_re + row_off + i, sb_mul8(vn, w));
          sb_store8(f->tile_im + row_off + i, vzero);
          sb_store8(f->ref_re + row_off + i, sb_mul8(vs, w));
          sb_store8(f->ref_im + row_off + i, vzero);
        }
      }
      for (; i < bf; i++) {
        const uint32_t bin = dftt_clamp_bin(fs + (int32_t)i, spec);
        esum += row_n[bin];
        const float w = wt * f->win_freq[i];
        const float resid = row_n[bin] - row_s[bin];
        wsum_r2 += (w * w) * (resid * resid);
        const size_t at = row_off + i;
        f->tile_re[at] = row_n[bin] * w;
        f->tile_im[at] = 0.0F;
        f->ref_re[at] = row_s[bin] * w;
        f->ref_im[at] = 0.0F;
      }
    }

    const float wt_last = f->win_time[bt - 1U];
    if (esum < DFTT_SILENCE_EPS) {
      /* Silent tile: carry the smoothed row through, OLA-weighted. */
      for (uint32_t i = 0U; i < bf; i++) {
        const uint32_t bin = dftt_clamp_bin(fs + (int32_t)i, spec);
        const float w = wt_last * f->win_freq[i];
        refined_snr[bin] += (w * w) * newest_smooth[bin];
      }
      continue;
    }

    dftt_fwd_rows(f, f->tile_re, f->tile_im);
    dftt_fwd_cols(f, f->tile_re, f->tile_im);
    dftt_fwd_rows(f, f->ref_re, f->ref_im);
    dftt_fwd_cols(f, f->ref_re, f->ref_im);
    /* Per-coefficient quefrency-domain rule (paper S4.2): the NLM-smoothed
     * tile provides the per-coefficient SNR estimate of the suppression
     * rule. Speckle is white in the tile-DFT domain, so its per-coefficient
     * power follows directly from Parseval on the spatial residual
     * (noisy - NLM): sigma2 = sum(w^2 * r^2) over the tile. Gain: Wiener
     * against the structure prior — coefficients where the NLM tile shows
     * structure (pr >> sigma2) pass the noisy (sharp) value, coefficients
     * where it shows none (pr ~ 0) die. The tile's flat level and slow
     * envelopes live at huge pr, so they pass without exemptions. */
    const float sigma2 = wsum_r2;
    const float speckle_power = f->kill_k * sigma2;
    for (size_t k = 0U; k < (size_t)bt * bf; k++) {
      const float pr =
          (f->ref_re[k] * f->ref_re[k]) + (f->ref_im[k] * f->ref_im[k]);
      const float g = pr / (pr + speckle_power);
      f->tile_re[k] *= g;
      f->tile_im[k] *= g;
    }

    /* Emit the newest time row only; past rows were emitted before. The
     * gathered tile already carries the analysis weight w, so emitting with
     * w gives syn*ana = w^2 per tile, matching the wsum normalization below
     * (unity gain reconstructs exactly). */
    dftt_inv_last_row(f, f->tile_re, f->tile_im, crow_re, crow_im);
    uint32_t i = 0U;
    if (!clamped) {
      const sb_vec8_t vw = sb_set8(wt_last);
      for (; i + 8U <= bf; i += 8U) {
        const sb_vec8_t w = sb_mul8(sb_load8(f->win_freq + i), vw);
        const sb_vec8_t v = sb_mul8(sb_load8(crow_re + i), w);
        sb_store8(refined_snr + fs + i,
                  sb_add8(sb_load8(refined_snr + fs + i), v));
      }
    }
    for (; i < bf; i++) {
      const uint32_t bin = dftt_clamp_bin(fs + (int32_t)i, spec);
      refined_snr[bin] += (wt_last * f->win_freq[i]) * crow_re[i];
    }
  }

  for (uint32_t k = 0U; k < spec; k++) {
    refined_snr[k] = f->wsum[k] > 0.0F ? refined_snr[k] / f->wsum[k] : 0.0F;
  }

  /* Monotone safety: the refined map may only remove energy relative to
   * NLM's verdict, never add it — windowing cross-talk from neighbouring
   * tiles can never creep energy back in (no hum creep). With the
   * structure-prior rule this is a no-op except in edge cases: refined
   * already sits below NLM on the diffuse floor and matches it on
   * structure. */
  for (uint32_t k = 0U; k < spec; k++) {
    if (refined_snr[k] > newest_smooth[k]) {
      refined_snr[k] = newest_smooth[k];
    }
  }

  sb_simd_restore_state(old_simd_state);
  return true;
}
