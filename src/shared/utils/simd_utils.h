/*
libspecbleach - A spectral processing library

Copyright 2024 Luciano Dato <lucianodato@gmail.com>

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

#ifndef SHARED_UTILS_SIMD_UTILS_H
#define SHARED_UTILS_SIMD_UTILS_H

#include <stdint.h>
#include <string.h>

#ifdef __ARM_NEON
#include <arm_neon.h>
#elif defined(__SSE__)
#include <xmmintrin.h>
#ifdef __AVX__
#include <immintrin.h>
#endif
#ifdef __SSE4_1__
#include <smmintrin.h>
#endif
#endif

#ifdef __GNUC__
#define SB_SIMD_INLINE static inline __attribute__((unused))
#else
#define SB_SIMD_INLINE static inline
#endif

/* ------------------------------------------------------------------------- */
/* DENORMAL HANDLING (FTZ/DAZ)                                               */
/* ------------------------------------------------------------------------- */

typedef uint32_t sb_simd_state_t;

/**
 * Enables "Flush-to-Zero" and "Denormals-are-Zero" modes.
 * returns the previous state to be restored later.
 */
SB_SIMD_INLINE sb_simd_state_t sb_simd_enable_ftz_daz(void) {
  sb_simd_state_t old_state = 0;
#ifdef __SSE__
  old_state = _mm_getcsr();
  _mm_setcsr(old_state | 0x8040); // MXCSR: bits 15 (FTZ) and 6 (DAZ)
#elif defined(__ARM_NEON) && defined(__aarch64__)
  // On ARM64, we manipulate the FPCR (Floating-point Control Register)
  // Bit 24 is FZ (Flush-to-zero)
  __asm__ __volatile__("mrs %x0, fpcr" : "=r"(old_state));
  sb_simd_state_t new_state = old_state | (1U << 24);
  __asm__ __volatile__("msr fpcr, %x0" : : "r"(new_state));
#endif
  return old_state;
}

/**
 * Restores a previously saved SIMD state.
 */
SB_SIMD_INLINE void sb_simd_restore_state(sb_simd_state_t state) {
#ifdef __SSE__
  _mm_setcsr(state);
#elif defined(__ARM_NEON) && defined(__aarch64__)
  __asm__ __volatile__("msr fpcr, %x0" : : "r"(state));
#else
  (void)state;
#endif
}

/* ------------------------------------------------------------------------- */
/* 8-WIDE VECTOR ABSTRACTION (Ideal for 8x8 patches)                         */
/* ------------------------------------------------------------------------- */

#ifdef __AVX__
typedef __m256 sb_vec8_t;
#elif defined(__SSE__)
typedef struct {
  __m128 v1, v2;
} sb_vec8_t;
#elif defined(__ARM_NEON)
typedef struct {
  float32x4_t v1, v2;
} sb_vec8_t;
#else
typedef struct {
  float v[8];
} sb_vec8_t;
#endif

/**
 * Loads 8 unaligned floats into an 8-wide vector.
 */
SB_SIMD_INLINE sb_vec8_t sb_load8(const float* p) {
#ifdef __AVX__
  return _mm256_loadu_ps(p);
#elif defined(__SSE__)
  sb_vec8_t r;
  r.v1 = _mm_loadu_ps(p);
  r.v2 = _mm_loadu_ps(p + 4);
  return r;
#elif defined(__ARM_NEON)
  sb_vec8_t r;
  r.v1 = vld1q_f32(p);
  r.v2 = vld1q_f32(p + 4);
  return r;
#else
  sb_vec8_t r;
  memcpy(r.v, p, 8 * sizeof(float));
  return r;
#endif
}

/**
 * Computes Sum of Squared Differences (SSD) between two 8-wide vectors.
 */
SB_SIMD_INLINE float sb_vec8_ssd(sb_vec8_t a, sb_vec8_t b) {
#ifdef __AVX__
  __m256 d = _mm256_sub_ps(a, b);
  d = _mm256_mul_ps(d, d);

  // Horizontal Sum
  __m128 low = _mm256_castps256_ps128(d);
  __m128 high = _mm256_extractf128_ps(d, 1);
  __m128 sum128 = _mm_add_ps(low, high);
  __m128 shuf = _mm_shuffle_ps(sum128, sum128, _MM_SHUFFLE(2, 3, 0, 1));
  __m128 sums = _mm_add_ps(sum128, shuf);
  shuf = _mm_movehl_ps(shuf, sums);
  sums = _mm_add_ss(sums, shuf);
  float f;
  _mm_store_ss(&f, sums);
  return f;
#elif defined(__SSE__)
  __m128 d1 = _mm_sub_ps(a.v1, b.v1);
  __m128 d2 = _mm_sub_ps(a.v2, b.v2);
  d1 = _mm_mul_ps(d1, d1);
  d2 = _mm_mul_ps(d2, d2);
  __m128 sum = _mm_add_ps(d1, d2);

  // Horizontal Sum
  __m128 shuf = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(2, 3, 0, 1));
  __m128 sums = _mm_add_ps(sum, shuf);
  shuf = _mm_movehl_ps(shuf, sums);
  sums = _mm_add_ss(sums, shuf);
  float f;
  _mm_store_ss(&f, sums);
  return f;
#elif defined(__ARM_NEON)
  float32x4_t d1 = vsubq_f32(a.v1, b.v1);
  float32x4_t d2 = vsubq_f32(a.v2, b.v2);
  d1 = vmulq_f32(d1, d1);
  d2 = vmulq_f32(d2, d2);
  float32x4_t sum = vaddq_f32(d1, d2);
  return vgetq_lane_f32(sum, 0) + vgetq_lane_f32(sum, 1) +
         vgetq_lane_f32(sum, 2) + vgetq_lane_f32(sum, 3);
#else
  float ssd = 0.0f;
  for (int i = 0; i < 8; i++) {
    float d = a.v[i] - b.v[i];
    ssd += d * d;
  }
  return ssd;
#endif
}

/**
 * Accumulates Sum of Squared Differences into an accumulator vector.
 * Useful for loops like register blocking.
 */
#ifdef __AVX__
typedef __m256 sb_acc8_t;
#elif defined(__SSE__)
typedef __m128 sb_acc8_t;
#elif defined(__ARM_NEON)
typedef float32x4_t sb_acc8_t;
#else
typedef float sb_acc8_t;
#endif

SB_SIMD_INLINE sb_acc8_t sb_acc8_zero(void) {
#ifdef __AVX__
  return _mm256_setzero_ps();
#elif defined(__SSE__)
  return _mm_setzero_ps();
#elif defined(__ARM_NEON)
  return vdupq_n_f32(0.0f);
#else
  return 0.0f;
#endif
}

SB_SIMD_INLINE sb_acc8_t sb_acc8_add_ssd(sb_acc8_t acc, sb_vec8_t a,
                                         sb_vec8_t b) {
#ifdef __AVX__
  __m256 d = _mm256_sub_ps(a, b);
  return _mm256_add_ps(acc, _mm256_mul_ps(d, d));
#elif defined(__SSE__)
  __m128 d1 = _mm_sub_ps(a.v1, b.v1);
  __m128 d2 = _mm_sub_ps(a.v2, b.v2);
  return _mm_add_ps(acc, _mm_add_ps(_mm_mul_ps(d1, d1), _mm_mul_ps(d2, d2)));
#elif defined(__ARM_NEON)
  float32x4_t d1 = vsubq_f32(a.v1, b.v1);
  float32x4_t d2 = vsubq_f32(a.v2, b.v2);
  return vaddq_f32(acc, vaddq_f32(vmulq_f32(d1, d1), vmulq_f32(d2, d2)));
#else
  float ssd = 0.0f;
  for (int i = 0; i < 8; i++) {
    float d = a.v[i] - b.v[i];
    ssd += d * d;
  }
  return acc + ssd;
#endif
}

SB_SIMD_INLINE float sb_acc8_hsum(sb_acc8_t acc) {
#ifdef __AVX__
  __m128 low = _mm256_castps256_ps128(acc);
  __m128 high = _mm256_extractf128_ps(acc, 1);
  __m128 sum128 = _mm_add_ps(low, high);
  __m128 shuf = _mm_shuffle_ps(sum128, sum128, _MM_SHUFFLE(2, 3, 0, 1));
  __m128 sums = _mm_add_ps(sum128, shuf);
  shuf = _mm_movehl_ps(shuf, sums);
  sums = _mm_add_ss(sums, shuf);
  float f;
  _mm_store_ss(&f, sums);
  return f;
#elif defined(__SSE__)
  __m128 shuf = _mm_shuffle_ps(acc, acc, _MM_SHUFFLE(2, 3, 0, 1));
  __m128 sums = _mm_add_ps(acc, shuf);
  shuf = _mm_movehl_ps(shuf, sums);
  sums = _mm_add_ss(sums, shuf);
  float f;
  _mm_store_ss(&f, sums);
  return f;
#elif defined(__ARM_NEON)
  return vgetq_lane_f32(acc, 0) + vgetq_lane_f32(acc, 1) +
         vgetq_lane_f32(acc, 2) + vgetq_lane_f32(acc, 3);
#else
  return acc;
#endif
}

/**
 * Basic Arithmetic for 8-wide vectors
 */
SB_SIMD_INLINE sb_vec8_t sb_add8(sb_vec8_t a, sb_vec8_t b) {
#ifdef __AVX__
  return _mm256_add_ps(a, b);
#elif defined(__SSE__)
  sb_vec8_t r;
  r.v1 = _mm_add_ps(a.v1, b.v1);
  r.v2 = _mm_add_ps(a.v2, b.v2);
  return r;
#elif defined(__ARM_NEON)
  sb_vec8_t r;
  r.v1 = vaddq_f32(a.v1, b.v1);
  r.v2 = vaddq_f32(a.v2, b.v2);
  return r;
#else
  sb_vec8_t r;
  for (int i = 0; i < 8; i++)
    r.v[i] = a.v[i] + b.v[i];
  return r;
#endif
}

SB_SIMD_INLINE sb_vec8_t sb_sub8(sb_vec8_t a, sb_vec8_t b) {
#ifdef __AVX__
  return _mm256_sub_ps(a, b);
#elif defined(__SSE__)
  sb_vec8_t r;
  r.v1 = _mm_sub_ps(a.v1, b.v1);
  r.v2 = _mm_sub_ps(a.v2, b.v2);
  return r;
#elif defined(__ARM_NEON)
  sb_vec8_t r;
  r.v1 = vsubq_f32(a.v1, b.v1);
  r.v2 = vsubq_f32(a.v2, b.v2);
  return r;
#else
  sb_vec8_t r;
  for (int i = 0; i < 8; i++)
    r.v[i] = a.v[i] - b.v[i];
  return r;
#endif
}

SB_SIMD_INLINE sb_vec8_t sb_mul8(sb_vec8_t a, sb_vec8_t b) {
#ifdef __AVX__
  return _mm256_mul_ps(a, b);
#elif defined(__SSE__)
  sb_vec8_t r;
  r.v1 = _mm_mul_ps(a.v1, b.v1);
  r.v2 = _mm_mul_ps(a.v2, b.v2);
  return r;
#elif defined(__ARM_NEON)
  sb_vec8_t r;
  r.v1 = vmulq_f32(a.v1, b.v1);
  r.v2 = vmulq_f32(a.v2, b.v2);
  return r;
#else
  sb_vec8_t r;
  for (int i = 0; i < 8; i++)
    r.v[i] = a.v[i] * b.v[i];
  return r;
#endif
}

SB_SIMD_INLINE sb_vec8_t sb_div8(sb_vec8_t a, sb_vec8_t b) {
#ifdef __AVX__
  return _mm256_div_ps(a, b);
#elif defined(__SSE__)
  sb_vec8_t r;
  r.v1 = _mm_div_ps(a.v1, b.v1);
  r.v2 = _mm_div_ps(a.v2, b.v2);
  return r;
#elif defined(__ARM_NEON)
  // NEON doesn't have a direct div instruction, use reciprocal or scalar
  float ta1[4];
  float tb1[4];
  float tr1[4];
  float ta2[4];
  float tb2[4];
  float tr2[4];
  vst1q_f32(ta1, a.v1);
  vst1q_f32(tb1, b.v1);
  vst1q_f32(ta2, a.v2);
  vst1q_f32(tb2, b.v2);
  for (int i = 0; i < 4; i++) {
    tr1[i] = ta1[i] / tb1[i];
    tr2[i] = ta2[i] / tb2[i];
  }
  sb_vec8_t r;
  r.v1 = vld1q_f32(tr1);
  r.v2 = vld1q_f32(tr2);
  return r;
#else
  sb_vec8_t r;
  for (int i = 0; i < 8; i++)
    r.v[i] = a.v[i] / b.v[i];
  return r;
#endif
}

SB_SIMD_INLINE sb_vec8_t sb_max8(sb_vec8_t a, sb_vec8_t b) {
#ifdef __AVX__
  return _mm256_max_ps(a, b);
#elif defined(__SSE__)
  sb_vec8_t r;
  r.v1 = _mm_max_ps(a.v1, b.v1);
  r.v2 = _mm_max_ps(a.v2, b.v2);
  return r;
#elif defined(__ARM_NEON)
  sb_vec8_t r;
  r.v1 = vmaxq_f32(a.v1, b.v1);
  r.v2 = vmaxq_f32(a.v2, b.v2);
  return r;
#else
  sb_vec8_t r;
  for (int i = 0; i < 8; i++)
    r.v[i] = (a.v[i] > b.v[i]) ? a.v[i] : b.v[i];
  return r;
#endif
}

SB_SIMD_INLINE void sb_store8(float* p, sb_vec8_t v) {
#ifdef __AVX__
  _mm256_storeu_ps(p, v);
#elif defined(__SSE__)
  _mm_storeu_ps(p, v.v1);
  _mm_storeu_ps(p + 4, v.v2);
#elif defined(__ARM_NEON)
  vst1q_f32(p, v.v1);
  vst1q_f32(p + 4, v.v2);
#else
  memcpy(p, v.v, 8 * sizeof(float));
#endif
}

SB_SIMD_INLINE sb_vec8_t sb_set8(float f) {
#ifdef __AVX__
  return _mm256_set1_ps(f);
#elif defined(__SSE__)
  sb_vec8_t r;
  r.v1 = _mm_set1_ps(f);
  r.v2 = r.v1;
  return r;
#elif defined(__ARM_NEON)
  sb_vec8_t r;
  r.v1 = vdupq_n_f32(f);
  r.v2 = r.v1;
  return r;
#else
  sb_vec8_t r;
  for (int i = 0; i < 8; i++)
    r.v[i] = f;
  return r;
#endif
}

/* ------------------------------------------------------------------------- */
/* 4-WIDE VECTOR ABSTRACTION                                                 */
/* ------------------------------------------------------------------------- */

#if defined(__SSE__) || defined(__AVX__)
typedef __m128 sb_vec4_t;
#elif defined(__ARM_NEON)
typedef float32x4_t sb_vec4_t;
#else
typedef struct {
  float v[4];
} sb_vec4_t;
#endif

SB_SIMD_INLINE sb_vec4_t sb_load4(const float* p) {
#if defined(__SSE__) || defined(__AVX__)
  return _mm_loadu_ps(p);
#elif defined(__ARM_NEON)
  return vld1q_f32(p);
#else
  sb_vec4_t r;
  memcpy(r.v, p, 4 * sizeof(float));
  return r;
#endif
}

SB_SIMD_INLINE float sb_vec4_ssd(sb_vec4_t a, sb_vec4_t b) {
#if defined(__SSE__) || defined(__AVX__)
  __m128 d = _mm_sub_ps(a, b);
  d = _mm_mul_ps(d, d);
  __m128 shuf = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2, 3, 0, 1));
  __m128 sums = _mm_add_ps(d, shuf);
  shuf = _mm_movehl_ps(shuf, sums);
  sums = _mm_add_ss(sums, shuf);
  float f;
  _mm_store_ss(&f, sums);
  return f;
#elif defined(__ARM_NEON)
  float32x4_t d = vsubq_f32(a, b);
  d = vmulq_f32(d, d);
  return vgetq_lane_f32(d, 0) + vgetq_lane_f32(d, 1) + vgetq_lane_f32(d, 2) +
         vgetq_lane_f32(d, 3);
#else
  float ssd = 0.0f;
  for (int i = 0; i < 4; i++) {
    float d = a.v[i] - b.v[i];
    ssd += d * d;
  }
  return ssd;
#endif
}

/**
 * 4-Wide Arithmetic
 */
SB_SIMD_INLINE sb_vec4_t sb_add4(sb_vec4_t a, sb_vec4_t b) {
#if defined(__SSE__) || defined(__AVX__)
  return _mm_add_ps(a, b);
#elif defined(__ARM_NEON)
  return vaddq_f32(a, b);
#else
  sb_vec4_t r;
  for (int i = 0; i < 4; i++)
    r.v[i] = a.v[i] + b.v[i];
  return r;
#endif
}

SB_SIMD_INLINE sb_vec4_t sb_sub4(sb_vec4_t a, sb_vec4_t b) {
#if defined(__SSE__) || defined(__AVX__)
  return _mm_sub_ps(a, b);
#elif defined(__ARM_NEON)
  return vsubq_f32(a, b);
#else
  sb_vec4_t r;
  for (int i = 0; i < 4; i++)
    r.v[i] = a.v[i] - b.v[i];
  return r;
#endif
}

SB_SIMD_INLINE sb_vec4_t sb_mul4(sb_vec4_t a, sb_vec4_t b) {
#if defined(__SSE__) || defined(__AVX__)
  return _mm_mul_ps(a, b);
#elif defined(__ARM_NEON)
  return vmulq_f32(a, b);
#else
  sb_vec4_t r;
  for (int i = 0; i < 4; i++)
    r.v[i] = a.v[i] * b.v[i];
  return r;
#endif
}

SB_SIMD_INLINE sb_vec4_t sb_set4(float f) {
#if defined(__SSE__) || defined(__AVX__)
  return _mm_set1_ps(f);
#elif defined(__ARM_NEON)
  return vdupq_n_f32(f);
#else
  sb_vec4_t r;
  for (int i = 0; i < 4; i++)
    r.v[i] = f;
  return r;
#endif
}

SB_SIMD_INLINE void sb_store4(float* p, sb_vec4_t v) {
#if defined(__SSE__) || defined(__AVX__)
  _mm_storeu_ps(p, v);
#elif defined(__ARM_NEON)
  vst1q_f32(p, v);
#else
  memcpy(p, v.v, 4 * sizeof(float));
#endif
}

/**
 * 8-Wide Bitwise/Logic (for masks)
 */
SB_SIMD_INLINE sb_vec8_t sb_gt8(sb_vec8_t a, sb_vec8_t b) {
#ifdef __AVX__
  return _mm256_cmp_ps(a, b, _CMP_GT_OQ);
#elif defined(__SSE__)
  sb_vec8_t r;
  r.v1 = _mm_cmpgt_ps(a.v1, b.v1);
  r.v2 = _mm_cmpgt_ps(a.v2, b.v2);
  return r;
#elif defined(__ARM_NEON)
  sb_vec8_t r;
  r.v1 = (float32x4_t)vcgtq_f32(a.v1, b.v1);
  r.v2 = (float32x4_t)vcgtq_f32(a.v2, b.v2);
  return r;
#else
  sb_vec8_t r;
  for (int i = 0; i < 8; i++)
    r.v[i] = (a.v[i] > b.v[i]) ? 1.0f : 0.0f;
  return r;
#endif
}

SB_SIMD_INLINE sb_vec8_t sb_sel8(sb_vec8_t mask, sb_vec8_t a, sb_vec8_t b) {
#ifdef __AVX__
  return _mm256_blendv_ps(b, a, mask);
#elif defined(__SSE__)
// SSE4.1 blendv_ps or manual bitwise
#ifdef __SSE4_1__
  return (sb_vec8_t){_mm_blendv_ps(b.v1, a.v1, mask.v1),
                     _mm_blendv_ps(b.v2, a.v2, mask.v2)};
#else
  __m128 res1 =
      _mm_or_ps(_mm_and_ps(mask.v1, a.v1), _mm_andnot_ps(mask.v1, b.v1));
  __m128 res2 =
      _mm_or_ps(_mm_and_ps(mask.v2, a.v2), _mm_andnot_ps(mask.v2, b.v2));
  return (sb_vec8_t){res1, res2};
#endif
#elif defined(__ARM_NEON)
  return (sb_vec8_t){vbslq_f32((uint32x4_t)mask.v1, a.v1, b.v1),
                     vbslq_f32((uint32x4_t)mask.v2, a.v2, b.v2)};
#else
  sb_vec8_t r;
  for (int i = 0; i < 8; i++)
    r.v[i] = (mask.v[i] != 0.0f) ? a.v[i] : b.v[i];
  return r;
#endif
}

#endif // SHARED_UTILS_SIMD_UTILS_H
