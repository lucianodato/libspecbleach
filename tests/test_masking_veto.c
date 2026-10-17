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

/*
 * Unit tests for the Masking Veto (psychoacoustic audibility floor)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/configurations.h"
#include "shared/denoiser_logic/processing/masking_veto.h"
#include "shared/utils/spectral_features.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

static void fill_uniform(float* buffer, uint32_t size, float value) {
  for (uint32_t i = 0U; i < size; i++) {
    buffer[i] = value;
  }
}

/* Run the veto for N frames with a caller-computed (reset) base gain, as in
 * the engine. The floor attack settles over ~one hop*20, so callers must
 * provide enough frames for asymptotic checks. */
static void run_frames(MaskingVeto* mv, const float* spectrum,
                       const float* noise, const float* future, float* gain,
                       uint32_t real_size, float base, float depth,
                       int frames) {
  for (int f = 0; f < frames; f++) {
    fill_uniform(gain, real_size, base);
    masking_veto_apply_floor(mv, spectrum, noise, future, gain, depth);
  }
}

void test_masking_veto_logic(void) {
  printf("Testing MaskingVeto floor logic...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE,
                                            POWER_SPECTRUM, true, true);
  TEST_ASSERT(mv != NULL, "Init failed");

  float spectrum[513];
  float noise[513];
  float gain[513];
  fill_uniform(spectrum, real_size, 1000.0F);
  fill_uniform(noise, real_size, 1.0F);

  // 1. Normal run: strong signal, tiny noise -> relief converges to the
  // ceiling (legacy alpha = 1.0 anchor); attack settles over ~a hop*20
  run_frames(mv, spectrum, noise, NULL, gain, real_size, 0.2F, 1.0F, 40);
  TEST_FLOAT_CLOSE(gain[100], MASKING_VETO_RELIEF_CEILING, 0.02F);

  // 2. Depth 0 must be an exact no-op
  fill_uniform(gain, real_size, 0.2F);
  masking_veto_apply_floor(mv, spectrum, noise, NULL, gain, 0.0F);
  TEST_FLOAT_CLOSE(gain[100], 0.2F, 1e-6F);

  // 3. Edge/Invalid Depth (negative treated as no-op)
  fill_uniform(gain, real_size, 0.5F);
  masking_veto_apply_floor(mv, spectrum, noise, NULL, gain, -1.0F);
  TEST_FLOAT_CLOSE(gain[100], 0.5F, 1e-6F);

  // 4. NULL guards
  masking_veto_apply_floor(NULL, spectrum, noise, NULL, gain, 1.0F);
  masking_veto_apply_floor(mv, NULL, noise, NULL, gain, 1.0F);
  masking_veto_apply_floor(mv, spectrum, NULL, NULL, gain, 1.0F);
  masking_veto_apply_floor(mv, spectrum, noise, NULL, NULL, 1.0F);

  // 5. Frame-rate normalization guards
  masking_veto_set_hop_sec(NULL, 0.0115F); // must not crash
  masking_veto_set_hop_sec(mv, 0.0115F);
  masking_veto_set_hop_sec(mv, 0.0F); // edge: ignored

  masking_veto_free(mv);
  printf("OK Masking Veto floor logic tests passed\n");
}

void test_audibility_floor_knobs(void) {
  printf("Testing audibility floor bounds...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;

  // Scenario: noise-only (no masker present anywhere) -> no relief
  {
    MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE,
                                              POWER_SPECTRUM, false, true);
    TEST_ASSERT(mv != NULL, "Init failed (noise-only scenario)");
    float* spectrum = (float*)calloc(real_size, sizeof(float));
    float* noise = (float*)calloc(real_size, sizeof(float));
    float* gain = (float*)calloc(real_size, sizeof(float));
    fill_uniform(spectrum, real_size, 0.0F); // Near silence
    fill_uniform(noise, real_size, 100.0F);  // Loud noise
    run_frames(mv, spectrum, noise, NULL, gain, real_size, 0.1F, 1.0F, 40);
    TEST_FLOAT_CLOSE(gain[100], 0.1F, 1e-3F);
    TEST_ASSERT(gain[42] >= 0.1F, "Relief must be upward-only");
    free(spectrum);
    free(noise);
    free(gain);
    masking_veto_free(mv);
  }

  // Scenario: masked noise -> relief, and depth lerps the amount
  {
    MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE,
                                              POWER_SPECTRUM, false, true);
    TEST_ASSERT(mv != NULL, "Init failed (masked scenario)");
    float* spectrum = (float*)calloc(real_size, sizeof(float));
    float* noise = (float*)calloc(real_size, sizeof(float));
    float* gain_full = (float*)calloc(real_size, sizeof(float));
    float* gain_half = (float*)calloc(real_size, sizeof(float));
    fill_uniform(spectrum, real_size, 1000.0F); // Strong signal
    fill_uniform(noise, real_size, 1.0F);       // Masked noise
    run_frames(mv, spectrum, noise, NULL, gain_half, real_size, 0.1F, 0.5F, 40);
    run_frames(mv, spectrum, noise, NULL, gain_full, real_size, 0.1F, 1.0F, 40);
    TEST_ASSERT(gain_full[100] > gain_half[100],
                "Full depth must relieve more than half depth");
    TEST_ASSERT(gain_half[100] > 0.1F,
                "Half depth must still provide some relief");
    TEST_FLOAT_CLOSE(gain_full[100], MASKING_VETO_RELIEF_CEILING, 0.02F);
    TEST_ASSERT(gain_full[100] <= MASKING_VETO_RELIEF_CEILING + 1e-3F,
                "Relief must never exceed the ceiling");
    free(spectrum);
    free(noise);
    free(gain_full);
    free(gain_half);
    masking_veto_free(mv);
  }

  printf("OK Audibility floor tests passed\n");
}

void test_temporal_stability(void) {
  printf("Testing temporal stability (attack/release/no crash)...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise = (float*)calloc(real_size, sizeof(float));
  float* gain = (float*)calloc(real_size, sizeof(float));

  MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE,
                                            POWER_SPECTRUM, false, true);
  TEST_ASSERT(mv != NULL, "Init failed");
  masking_veto_set_hop_sec(mv, 0.0116F);

  // A: sustained signal -> relief converges to the ceiling (never overshoots)
  fill_uniform(spectrum, real_size, 1000.0F);
  fill_uniform(noise, real_size, 1.0F);
  run_frames(mv, spectrum, noise, NULL, gain, real_size, 0.1F, 1.0F, 60);
  TEST_ASSERT(gain[100] > 0.20F, "Protected band should get relief");
  TEST_ASSERT(gain[100] <= MASKING_VETO_RELIEF_CEILING + 1e-3F,
              "Relief must never overshoot the ceiling");

  // B: signal vanishes. The floor must decay smoothly (release), never
  // crash back to zero on the first weak frame (the alpha veto's
  // musical-noise generator at note decay). Persistence-only relief
  // preserves low-level detail, capped at the stricter persistence ceiling
  // — the masking model shapes reduction, it never overrides it.
  fill_uniform(spectrum, real_size, 0.0F);
  fill_uniform(noise, real_size, 100.0F);
  for (int frame = 0; frame < 3; frame++) {
    fill_uniform(gain, real_size, 0.1F);
    masking_veto_apply_floor(mv, spectrum, noise, NULL, gain, 1.0F);
    TEST_ASSERT(gain[100] > 0.10F,
                "Floor must hold through note-off (release)");
  }
  // ...until the release settles (the floor decays smoothly through the
  // estimator's forward-masking decay, no crash back to zero)
  for (int frame = 0; frame < 150; frame++) {
    fill_uniform(gain, real_size, 0.1F);
    masking_veto_apply_floor(mv, spectrum, noise, NULL, gain, 1.0F);
  }
  TEST_ASSERT(gain[100] < 0.20F, "Floor must decay to no relief eventually");

  // C: signal returns -> attack re-claims relief gradually (no flicker
  // between carved and relieved states)
  fill_uniform(spectrum, real_size, 1000.0F);
  fill_uniform(noise, real_size, 1.0F);
  run_frames(mv, spectrum, noise, NULL, gain, real_size, 0.1F, 1.0F, 30);
  TEST_ASSERT(gain[100] > 0.20F, "Attack must re-claim relief gradually");

  free(spectrum);
  free(noise);
  free(gain);
  masking_veto_free(mv);

  printf("OK Temporal stability tests passed\n");
}

void test_backward_masking_future_frame(void) {
  printf("Testing backward masking from the future frame...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;
  float* spectrum = (float*)calloc(real_size, sizeof(float));
  float* noise = (float*)calloc(real_size, sizeof(float));
  float* future = (float*)calloc(real_size, sizeof(float));
  float* gain = (float*)calloc(real_size, sizeof(float));

  // Noise-only bands (clean gated to zero): no masker anywhere, so
  // suppression must run at full depth (no relief).
  fill_uniform(noise, real_size, 1.0F);

  MaskingVeto* mv = masking_veto_initialize(fft_size, 44100, OPUS_SCALE,
                                            POWER_SPECTRUM, false, true);
  TEST_ASSERT(mv != NULL, "Init failed");
  masking_veto_set_hop_sec(mv, 0.0116F);

  run_frames(mv, spectrum, noise, NULL, gain, real_size, 0.1F, 1.0F, 40);
  TEST_FLOAT_CLOSE(gain[100], 0.1F, 0.001F);

  // Loud impulse in the next frame: pre-masking (backward masking) must
  // relieve suppression in the current frame ahead of the hit.
  fill_uniform(future, real_size, 1000.0F);
  run_frames(mv, spectrum, noise, future, gain, real_size, 0.1F, 1.0F, 60);
  TEST_ASSERT(gain[100] > 0.1F + 0.05F,
              "Future masker must relieve suppression via backward masking");

  free(spectrum);
  free(noise);
  free(future);
  free(gain);
  masking_veto_free(mv);

  printf("OK Backward masking tests passed\n");
}

int main(void) {
  printf("Running Masking Veto tests...\n");
  test_masking_veto_logic();
  test_audibility_floor_knobs();
  test_temporal_stability();
  test_backward_masking_future_frame();
  printf("OK All Masking Veto tests passed!\n");
  return 0;
}
