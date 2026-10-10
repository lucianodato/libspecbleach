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

#include "shared/denoiser_logic/processing/patch_filter_context.h"
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

static void test_history_lifecycle_and_indexing(void) {
  PatchFilterContext context = {0};
  TEST_ASSERT(patch_filter_context_initialize(&context, 2U, 5U, 2U, 1U, 1U),
              "Context initialization should succeed");
  TEST_ASSERT(!patch_filter_context_is_ready(&context),
              "Empty history should not be ready");
  TEST_ASSERT(patch_filter_context_get_latency_frames(&context) == 1U,
              "Latency should equal the future search range");

  for (uint32_t frame_index = 0U; frame_index < 5U; frame_index++) {
    float frame[2] = {(float)frame_index, (float)(frame_index + 100U)};
    patch_filter_context_push_frame(&context, frame);
  }
  TEST_ASSERT(patch_filter_context_is_ready(&context),
              "Full history should be ready");
  TEST_ASSERT(context.buffer_head == 0U && context.frames_filled == 5U,
              "Full history should wrap and saturate the frame count");
  TEST_ASSERT(patch_filter_context_get_frame(&context, 0)[0] == 3.0F,
              "Target index should retain the future-frame offset");
  TEST_ASSERT(patch_filter_context_get_frame(&context, -2)[0] == 1.0F,
              "Past frame lookup should wrap correctly");
  TEST_ASSERT(patch_filter_context_get_frame(&context, 1)[0] == 4.0F,
              "Future frame lookup should wrap correctly");

  patch_filter_context_populate_frame_ptrs(&context);
  for (int32_t dt = -10; dt <= 9; dt++) {
    TEST_ASSERT(patch_filter_context_cached_get_frame(&context, dt) ==
                    patch_filter_context_get_frame(&context, dt),
                "Cached frame pointer should match ring lookup");
  }

  float wrapped_frame[2] = {5.0F, 105.0F};
  patch_filter_context_push_frame(&context, wrapped_frame);
  TEST_ASSERT(context.buffer_head == 1U && context.frames_filled == 5U,
              "Additional pushes should wrap without growing readiness");
  TEST_ASSERT(patch_filter_context_get_frame(&context, 0)[0] == 4.0F &&
                  patch_filter_context_get_frame(&context, 1)[0] == 5.0F,
              "Wrapped write should retain the configured future offset");

  patch_filter_context_reset(&context);
  TEST_ASSERT(!patch_filter_context_is_ready(&context) &&
                  context.buffer_head == 0U && context.frames_filled == 0U,
              "Reset should clear history state");
  for (uint32_t i = 0U; i < context.time_buffer_size; i++) {
    TEST_ASSERT(context.frame_buffer[i][0] == 0.0F &&
                    context.frame_buffer[i][1] == 0.0F,
                "Reset should zero every history frame");
  }

  patch_filter_context_free(&context);
  patch_filter_context_free(NULL);
}

static void test_clamped_patch_row_ssd(void) {
  const float target[3] = {1.0F, 2.0F, 3.0F};
  const float candidate[3] = {0.0F, 2.0F, 4.0F};
  const float distance = patch_filter_accumulate_patch_row_ssd(
      7.0F, (PatchRowArgs){.target_row = target,
                           .candidate_row = candidate,
                           .target_freq = 0U,
                           .candidate_freq = 2U,
                           .patch_size = 3U,
                           .half_patch_size = 1U,
                           .spectrum_size = 3U});
  TEST_ASSERT(distance == 21.0F,
              "Clamped SSD should retain the initial sum and bin order");
}

int main(void) {
  test_history_lifecycle_and_indexing();
  test_clamped_patch_row_ssd();
  printf("Patch filter context tests passed.\n");
  return 0;
}
