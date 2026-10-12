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
 * Unit tests for BM3D-lite Filter (1-pass collaborative)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "shared/denoiser_logic/processing/bm3d_filter.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

static void test_bm3d_lifecycle(void) {
  Bm3dFilterConfig config = {
      .spectrum_size = 64,
      .time_buffer_size = 21,
      .patch_size = 8,
      .search_range_freq = 8,
      .search_range_time_past = 16,
      .search_range_time_future = 4,
      .h_parameter = 1.0f,
  };
  Bm3dFilter* filter = bm3d_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "BM3D init should succeed");
  TEST_ASSERT(bm3d_filter_get_latency_frames(filter) == 4,
              "BM3D latency must equal future search range (NLM parity)");
  bm3d_filter_free(filter);

  Bm3dFilterConfig invalid = {.spectrum_size = 0};
  TEST_ASSERT(bm3d_filter_initialize(invalid) == NULL,
              "BM3D init with zero spectrum_size should fail");
}

static void test_bm3d_geometry_validation(void) {
  // Nonzero time_buffer_size smaller than past+future+1 would alias the
  // frame ring and must be rejected.
  Bm3dFilterConfig short_buffer = {
      .spectrum_size = 64,
      .time_buffer_size = 20, // past(16) + future(4) + 1 == 21 required
      .patch_size = 8,
      .search_range_freq = 8,
      .search_range_time_past = 16,
      .search_range_time_future = 4,
      .h_parameter = 1.0f,
  };
  TEST_ASSERT(bm3d_filter_initialize(short_buffer) == NULL,
              "BM3D init with time_buffer_size < past+future+1 must fail");

  // Patch half-width beyond NLM_HALO_FRAMES (8) would index the frame
  // pointer cache out of bounds and must be rejected.
  Bm3dFilterConfig wide_patch = {
      .spectrum_size = 64,
      .time_buffer_size = 21,
      .patch_size = 18, // half-width 9 > NLM_HALO_FRAMES (8)
      .search_range_freq = 8,
      .search_range_time_past = 16,
      .search_range_time_future = 4,
      .h_parameter = 1.0f,
  };
  TEST_ASSERT(bm3d_filter_initialize(wide_patch) == NULL,
              "BM3D init with patch half-width > NLM_HALO_FRAMES must fail");

  // Exact boundary (half-width == NLM_HALO_FRAMES) stays valid.
  Bm3dFilterConfig edge_patch = {
      .spectrum_size = 64,
      .time_buffer_size = 21,
      .patch_size = 16, // half-width 8 == NLM_HALO_FRAMES
      .search_range_freq = 8,
      .search_range_time_past = 16,
      .search_range_time_future = 4,
      .h_parameter = 1.0f,
  };
  Bm3dFilter* filter = bm3d_filter_initialize(edge_patch);
  TEST_ASSERT(filter != NULL, "patch half-width == NLM_HALO_FRAMES is valid");
  bm3d_filter_free(filter);
}

static void test_bm3d_uniform_and_silence(void) {
  Bm3dFilterConfig config = {
      .spectrum_size = 32,
      .time_buffer_size = 5,
      .patch_size = 4,
      .search_range_freq = 2,
      .search_range_time_past = 2,
      .search_range_time_future = 2,
      .h_parameter = 1.0f,
  };
  Bm3dFilter* filter = bm3d_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "BM3D init should succeed");
  TEST_ASSERT(!bm3d_filter_is_ready(filter), "Not ready before fill");

  float frame[32];
  for (int i = 0; i < 32; i++) {
    frame[i] = 5.0f;
  }
  for (int f = 0; f < 5; f++) {
    bm3d_filter_push_frame(filter, frame);
  }
  TEST_ASSERT(bm3d_filter_is_ready(filter), "Ready after fill");

  float output[32];
  TEST_ASSERT(bm3d_filter_process(filter, output), "Process should succeed");
  for (int i = 0; i < 32; i++) {
    TEST_ASSERT(fabsf(output[i] - 5.0f) < 0.01f,
                "Uniform input must pass through");
  }

  // Silence: valid output, no NaN/Inf.
  float silence[32] = {0.0f};
  for (int f = 0; f < 5; f++) {
    bm3d_filter_push_frame(filter, silence);
  }
  TEST_ASSERT(bm3d_filter_process(filter, output), "Process silence");
  for (int i = 0; i < 32; i++) {
    TEST_ASSERT(isfinite(output[i]), "Silence output must be finite");
    TEST_ASSERT(output[i] >= 0.0f, "Silence output must be >= 0");
  }

  bm3d_filter_reset(filter);
  TEST_ASSERT(!bm3d_filter_is_ready(filter), "Not ready after reset");
  bm3d_filter_free(filter);
}

static void test_bm3d_spike_reduction(void) {
  Bm3dFilterConfig config = {
      .spectrum_size = 64,
      .time_buffer_size = 7,
      .patch_size = 4,
      .paste_block_size = 8,
      .search_range_freq = 8,
      .search_range_time_past = 4,
      .search_range_time_future = 2,
      .h_parameter = 1.0f,
      .num_threads = 1, // single-threaded reference path
  };
  Bm3dFilter* filter = bm3d_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "BM3D init should succeed");

  // Flat background with one isolated spike in the middle frame: the stack
  // mean collapses the outlier residual, so the spike must shrink while the
  // background passes through.
  float flat[64];
  float spiked[64];
  for (int i = 0; i < 64; i++) {
    flat[i] = 2.0f;
    spiked[i] = 2.0f;
  }
  spiked[32] = 20.0f;
  // 4 past frames + spike + 2 future frames: the spiked frame lands exactly
  // on the processed slot (lag = search_range_time_future), so the shrink
  // actually sees it.
  for (int f = 0; f < 4; f++) {
    bm3d_filter_push_frame(filter, flat);
  }
  bm3d_filter_push_frame(filter, spiked);
  for (int f = 0; f < 2; f++) {
    bm3d_filter_push_frame(filter, flat);
  }
  float output[64];
  TEST_ASSERT(bm3d_filter_process(filter, output), "Process should succeed");
  for (int i = 0; i < 64; i++) {
    TEST_ASSERT(isfinite(output[i]), "Output must be finite");
    if (i != 32) {
      TEST_ASSERT(fabsf(output[i] - 2.0f) < 0.5f,
                  "Background must pass through");
    }
  }
  TEST_ASSERT(output[32] < 20.0f, "Isolated spike must be reduced");
  TEST_ASSERT(output[32] >= 2.0f, "Spike must not undershoot background");
  bm3d_filter_free(filter);
}

static void test_bm3d_edge_blocks_scalar_path(void) {
  // Non-8 patch size forces the scalar SSD path; small spectrum forces
  // edge paste blocks through clamped bin lookups (both hit the helpers
  // extracted from bm3d_process_block_range).
  Bm3dFilterConfig config = {
      .spectrum_size = 12,
      .time_buffer_size = 5,
      .patch_size = 3,
      .paste_block_size = 4,
      .search_range_freq = 1,
      .search_range_time_past = 2,
      .search_range_time_future = 1,
      .h_parameter = 1.0f,
      .num_threads = 1,
  };
  Bm3dFilter* filter = bm3d_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "BM3D init should succeed");

  float frame[12];
  for (int i = 0; i < 12; i++) {
    frame[i] = 1.5f + (float)(i % 3) * 0.25f;
  }
  for (int f = 0; f < 5; f++) {
    bm3d_filter_push_frame(filter, frame);
  }
  TEST_ASSERT(bm3d_filter_is_ready(filter), "Ready after fill");

  float output[12];
  TEST_ASSERT(bm3d_filter_process(filter, output),
              "Edge blocks must process without clamping crashes");
  for (int i = 0; i < 12; i++) {
    TEST_ASSERT(isfinite(output[i]), "Output must be finite");
    TEST_ASSERT(output[i] >= 0.0f, "Output must be >= 0");
  }
  bm3d_filter_free(filter);
}

static void test_bm3d_worst_slot_replacement(void) {
  // Many distinct matches exceed BM3D_STACK_MAX so worst-slot replacement
  // runs repeatedly (covers bm3d_worst_slot across a full search).
  Bm3dFilterConfig config = {
      .spectrum_size = 32,
      .time_buffer_size = 9,
      .patch_size = 4,
      .paste_block_size = 8,
      .search_range_freq = 4,
      .search_range_time_past = 4,
      .search_range_time_future = 1,
      .h_parameter = 2.0f, // wide matching so many candidates qualify
      .num_threads = 1,
  };
  Bm3dFilter* filter = bm3d_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "BM3D init should succeed");

  // Distinct random-ish frames: every neighbor is a plausible match.
  float frame[32];
  for (int f = 0; f < 9; f++) {
    for (int i = 0; i < 32; i++) {
      frame[i] = 1.0f + (float)((i * 7 + f * 13) % 5) * 0.3f;
    }
    bm3d_filter_push_frame(filter, frame);
  }
  TEST_ASSERT(bm3d_filter_is_ready(filter), "Ready after fill");

  float output[32];
  TEST_ASSERT(bm3d_filter_process(filter, output),
              "Search with more matches than stack slots must succeed");
  for (int i = 0; i < 32; i++) {
    TEST_ASSERT(isfinite(output[i]), "Output must be finite");
    TEST_ASSERT(output[i] >= 0.0f, "Output must be >= 0");
  }
  bm3d_filter_free(filter);
}

static void test_bm3d_vec8_patch_full_stack(void) {
  // patch == 8 exercises the SIMD fast path (target preload, candidate row
  // loads and vec8 SSD) plus both residual branches of the collaborative
  // shrink when matches land under and over the hard threshold.
  Bm3dFilterConfig config = {
      .spectrum_size = 48,
      .time_buffer_size = 9,
      .patch_size = 8,
      .paste_block_size = 4,
      .search_range_freq = 4,
      .search_range_time_past = 4,
      .search_range_time_future = 1,
      .h_parameter = 1.0f,
      .num_threads = 1,
  };
  Bm3dFilter* filter = bm3d_filter_initialize(config);
  TEST_ASSERT(filter != NULL, "BM3D init should succeed");

  // Mixed content: smooth background plus one sharp peak so some match
  // residuals survive the hard threshold and others collapse to the mean.
  float frame[48];
  for (int f = 0; f < 9; f++) {
    for (int i = 0; i < 48; i++) {
      frame[i] = 2.0f + 0.1f * (float)((i + f) % 7);
    }
    frame[24] = 30.0f; // sharp stationary peak across frames
    bm3d_filter_push_frame(filter, frame);
  }
  TEST_ASSERT(bm3d_filter_is_ready(filter), "Ready after fill");

  float output[48];
  TEST_ASSERT(bm3d_filter_process(filter, output),
              "Vec8 patch path must process");
  for (int i = 0; i < 48; i++) {
    TEST_ASSERT(isfinite(output[i]), "Output must be finite");
    TEST_ASSERT(output[i] >= 0.0f, "Output must be >= 0");
  }
  bm3d_filter_free(filter);
}

int main(void) {
  printf("Running BM3D filter tests...\n\n");
  test_bm3d_lifecycle();
  test_bm3d_geometry_validation();
  test_bm3d_uniform_and_silence();
  test_bm3d_spike_reduction();
  test_bm3d_edge_blocks_scalar_path();
  test_bm3d_worst_slot_replacement();
  test_bm3d_vec8_patch_full_stack();
  printf("\nAll BM3D filter tests passed!\n");
  return 0;
}
