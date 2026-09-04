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

/* Unit tests for frame-rate normalization helpers (#152). Fast and exact:
 * pins the per-option NLM table, the ms/Hz fallback formulas, clamps, and
 * alpha retuning (alpha_new == alpha_ref^(hop_new/hop_ref)). */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "shared/configurations.h"
#include "shared/frame_rate_norm.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

static void test_hop_and_bin_helpers(void) {
  printf("Testing hop/bin helpers...\n");
  TEST_ASSERT(
      sb_hop_sec(552U, 48000U) > 0.0114F && sb_hop_sec(552U, 48000U) < 0.0116F,
      "hop_sec(552 @48k) should be ~11.5ms");
  TEST_ASSERT(sb_hop_sec(0U, 48000U) == 0.0F, "zero hop must give 0");
  TEST_ASSERT(sb_hop_sec(552U, 0U) == 0.0F, "zero rate must give 0");
  TEST_ASSERT(
      sb_bin_hz(48000U, 2048U) > 23.4F && sb_bin_hz(48000U, 2048U) < 23.5F,
      "bin_hz(48k/2048) should be ~23.4Hz");
  TEST_ASSERT(sb_bin_hz(48000U, 0U) == 0.0F, "zero fft must give 0");
  printf("✓ Hop/bin helpers passed\n");
}

static void test_alpha_retuning(void) {
  printf("Testing alpha retuning...\n");
  /* Doubling the hop squares the alpha: 0.8^(2) == 0.64. */
  float a =
      sb_alpha_retuned(0.8F, LEGACY_REF_HOP_SEC, 2.0F * LEGACY_REF_HOP_SEC);
  TEST_ASSERT(fabsf(a - 0.64F) < 1e-6F, "doubled hop must square alpha");
  /* Same hop is identity. */
  float b = sb_alpha_retuned(0.8F, LEGACY_REF_HOP_SEC, LEGACY_REF_HOP_SEC);
  TEST_ASSERT(fabsf(b - 0.8F) < 1e-6F, "same hop must be identity");
  /* Degenerate inputs fall back to the reference alpha. */
  TEST_ASSERT(sb_alpha_retuned(0.0F, LEGACY_REF_HOP_SEC, 0.02F) == 0.0F,
              "zero alpha must pass through");
  TEST_ASSERT(sb_alpha_retuned(1.5F, LEGACY_REF_HOP_SEC, 0.02F) == 1.5F,
              "alpha >= 1 must pass through");
  TEST_ASSERT(sb_alpha_retuned(0.8F, 0.0F, 0.02F) == 0.8F,
              "zero ref hop must pass through");
  TEST_ASSERT(sb_alpha_from_tau(0.0F, 0.01F) == 0.0F, "zero tau must give 0");
  TEST_ASSERT(sb_tau_from_alpha(1.0F, 0.01F) == 0.0F,
              "unit alpha must give 0 tau");
  printf("✓ Alpha retuning passed\n");
}

static void test_frames_and_bins(void) {
  printf("Testing frames/bins conversion...\n");
  TEST_ASSERT(sb_frames_for_ms(128.0F, 0.0115F, 8U, 32U) == 11U,
              "128ms @11.5ms hop must be 11 frames");
  TEST_ASSERT(sb_frames_for_ms(0.0F, 0.0115F, 8U, 32U) == 8U,
              "zero ms must clamp to lo");
  TEST_ASSERT(sb_frames_for_ms(10000.0F, 0.0115F, 8U, 32U) == 32U,
              "huge ms must clamp to hi");
  TEST_ASSERT(sb_frames_for_ms(128.0F, 0.0F, 8U, 32U) == 8U,
              "zero hop must give lo");
  TEST_ASSERT(sb_bins_for_hz(170.0F, 21.5F, 2U, 32U) == 8U,
              "170Hz @21.5Hz bins must be 8");
  TEST_ASSERT(sb_bins_for_hz(0.0F, 21.5F, 2U, 32U) == 2U,
              "zero Hz must clamp to lo");
  TEST_ASSERT(sb_bins_for_hz(170.0F, 0.0F, 2U, 32U) == 2U,
              "zero bin width must give lo");
  TEST_ASSERT(sb_bins_for_hz(170.0F, 1.0F, 2U, 32U) == 32U,
              "tiny bins must clamp to hi");
  printf("✓ Frames/bins conversion passed\n");
}

static void check_table(float frame_ms, uint32_t patch, uint32_t past,
                        uint32_t future) {
  const float hop_sec = frame_ms / 4000.0F; /* hop = frame/4 */
  SbNlmGeometry g = sb_nlm_geometry_for_frame_ms(frame_ms, hop_sec, 0.0F);
  if (g.patch != patch || g.past != past || g.future != future) {
    fprintf(stderr,
            "TEST FAILED: %.0fms geometry is %u/%u/%u, expected %u/%u/%u\n",
            frame_ms, g.patch, g.past, g.future, patch, past, future);
    exit(1);
  }
}

static void test_nlm_table(void) {
  printf("Testing NLM per-option table...\n");
  check_table(23.0F, 16U, 32U, 8U);
  check_table(32.0F, 12U, 23U, 6U);
  check_table(46.0F, 8U, 16U, 4U);
  check_table(64.0F, 6U, 12U, 3U);
  check_table(93.0F, 4U, 8U, 2U);

  /* Fallback formula for non-table sizes: 50ms -> 92/184/46ms targets. */
  SbNlmGeometry g50 = sb_nlm_geometry_for_frame_ms(50.0F, 0.0125F, 21.5F);
  TEST_ASSERT(g50.patch == 7U && g50.past == 15U && g50.future == 4U,
              "50ms fallback must be 7/15/4");
  TEST_ASSERT(g50.search_freq == 8U && g50.paste == 4U,
              "50ms @21.5Hz bins must be 8 search / 4 paste");

  /* Paste never exceeds patch: tiny bins + 93ms frame. */
  SbNlmGeometry gclamp =
      sb_nlm_geometry_for_frame_ms(93.0F, 93.0F / 4000.0F, 1.95F);
  TEST_ASSERT(gclamp.search_freq == 32U, "tiny bins must clamp search to 32");
  TEST_ASSERT(gclamp.paste == gclamp.patch,
              "paste must clamp down to patch size");

  /* Degenerate hop keeps legacy defaults. */
  SbNlmGeometry gdeg = sb_nlm_geometry_for_frame_ms(50.0F, 0.0F, 0.0F);
  TEST_ASSERT(gdeg.patch == NLM_PATCH_SIZE &&
                  gdeg.past == NLM_SEARCH_RANGE_TIME_PAST &&
                  gdeg.future == NLM_SEARCH_RANGE_TIME_FUTURE,
              "zero hop must keep legacy defaults");
  printf("✓ NLM per-option table passed\n");
}

int main(void) {
  printf("Running frame-rate normalization tests...\n\n");
  test_hop_and_bin_helpers();
  test_alpha_retuning();
  test_frames_and_bins();
  test_nlm_table();
  printf("\n✅ All frame-rate normalization tests passed!\n");
  return 0;
}
