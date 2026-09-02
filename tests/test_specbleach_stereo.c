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

#include "specbleach_stereo.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define BLOCK_SIZE 512
#define CHANNELS 2

static void test_init_and_free(void) {
  printf("Testing stereo init and free...\n");

  TEST_ASSERT(specbleach_stereo_initialize(0, 46.0f, CHANNELS) == NULL,
              "sample_rate == 0 rejected");
  TEST_ASSERT(specbleach_stereo_initialize(44100, 0.0f, CHANNELS) == NULL,
              "frame_size <= 0 rejected");
  TEST_ASSERT(specbleach_stereo_initialize(44100, 46.0f, 0) == NULL,
              "channels == 0 rejected");
  specbleach_stereo_free(NULL);

  specbleach_stereo* stereo =
      specbleach_stereo_initialize(44100, 46.0f, CHANNELS);
  TEST_ASSERT(stereo != NULL, "stereo init");
  TEST_ASSERT(specbleach_stereo_get_channel_count(stereo) == CHANNELS,
              "Channel count roundtrip");
  TEST_ASSERT(specbleach_stereo_get_latency(stereo) > 0, "Latency positive");
  TEST_ASSERT(specbleach_stereo_get_noise_profile_size(stereo) > 0,
              "Profile size positive");
  specbleach_stereo_free(stereo);

  printf("✓ Stereo init and free tests passed\n");
}

static void test_parameter_loading(void) {
  printf("Testing stereo parameter loading...\n");

  specbleach_stereo* stereo =
      specbleach_stereo_initialize(44100, 46.0f, CHANNELS);
  TEST_ASSERT(stereo != NULL, "init for param tests");

  SpecbleachDenoiserParameters params = {0};
  params.learn_noise = SPECBLEACH_LEARN_ALL;
  params.reduction_gain = 0.1f;

  TEST_ASSERT(
      specbleach_stereo_load_parameters(NULL, &params, sizeof(params)) == false,
      "NULL instance rejected");
  TEST_ASSERT(
      specbleach_stereo_load_parameters(stereo, NULL, sizeof(params)) == false,
      "NULL parameters rejected");
  TEST_ASSERT(specbleach_stereo_load_parameters(stereo, &params,
                                                sizeof(params) - 4) == false,
              "Wrong size rejected");
  TEST_ASSERT(
      specbleach_stereo_load_parameters(stereo, &params, sizeof(params)),
      "Valid parameters accepted");

  specbleach_stereo_free(stereo);
  printf("✓ Stereo parameter loading tests passed\n");
}

static void test_process_and_sync(void) {
  printf("Testing stereo processing and profile sync...\n");

  specbleach_stereo* stereo =
      specbleach_stereo_initialize(44100, 46.0f, CHANNELS);
  TEST_ASSERT(stereo != NULL, "init for process tests");

  float input[2][BLOCK_SIZE];
  float output[2][BLOCK_SIZE];
  const float* input_ptrs[CHANNELS] = {input[0], input[1]};
  float* output_ptrs[CHANNELS] = {output[0], output[1]};

  for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
    input[0][s] = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
    input[1][s] = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
  }
  memset(output, 0, sizeof(output));

  TEST_ASSERT(specbleach_stereo_process(NULL, BLOCK_SIZE, input_ptrs,
                                        output_ptrs) == false,
              "NULL instance process rejected");
  TEST_ASSERT(
      specbleach_stereo_process(stereo, 0, input_ptrs, output_ptrs) == false,
      "Zero frames process rejected");

  // Learn phase on noise so profiles get populated per channel
  SpecbleachDenoiserParameters params = {0};
  params.learn_noise = SPECBLEACH_LEARN_ALL;
  params.reduction_gain = 0.1f;
  TEST_ASSERT(
      specbleach_stereo_load_parameters(stereo, &params, sizeof(params)),
      "learn params loaded");

  for (uint32_t block = 0; block < 8; ++block) {
    TEST_ASSERT(
        specbleach_stereo_process(stereo, BLOCK_SIZE, input_ptrs, output_ptrs),
        "process during learn");
  }

  const float intensity = specbleach_stereo_get_transient_intensity(stereo);
  TEST_ASSERT(intensity >= 0.0f && intensity <= 1.0f,
              "transient intensity within range");

  // Reset clears availability again
  TEST_ASSERT(specbleach_stereo_reset_profiles(stereo), "profiles reset");
  TEST_ASSERT(!specbleach_stereo_profile_available_for_channel(stereo, 0, 1),
              "mode 1 unavailable after reset");

  /* Full lifecycle: learn -> stop learning (this finalizes all capture
   * modes) -> every mode available on every channel -> sync succeeds. */
  SpecbleachDenoiserParameters reduce_params = {0};
  reduce_params.learn_noise = SPECBLEACH_LEARN_ALL;
  reduce_params.reduction_gain = 0.1f;
  TEST_ASSERT(specbleach_stereo_load_parameters(stereo, &reduce_params,
                                                sizeof(reduce_params)),
              "re-learn params loaded");
  for (uint32_t block = 0; block < 8; ++block) {
    TEST_ASSERT(
        specbleach_stereo_process(stereo, BLOCK_SIZE, input_ptrs, output_ptrs),
        "process during re-learn");
  }
  reduce_params.learn_noise = SPECBLEACH_LEARN_OFF;
  TEST_ASSERT(specbleach_stereo_load_parameters(stereo, &reduce_params,
                                                sizeof(reduce_params)),
              "reduce params loaded");
  TEST_ASSERT(
      specbleach_stereo_process(stereo, BLOCK_SIZE, input_ptrs, output_ptrs),
      "finalize block processed");

  for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
    TEST_ASSERT(
        specbleach_stereo_profile_available_for_channel(stereo, 0, mode),
        "ch0 mode finalized");
    TEST_ASSERT(
        specbleach_stereo_profile_available_for_channel(stereo, 1, mode),
        "ch1 mode finalized");
  }

  TEST_ASSERT(specbleach_stereo_sync_profiles(stereo),
              "sync succeeds with populated channels");
  for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
    TEST_ASSERT(
        specbleach_stereo_profile_available_for_channel(stereo, 0, mode),
        "ch0 mode available after sync");
    TEST_ASSERT(
        specbleach_stereo_profile_available_for_channel(stereo, 1, mode),
        "ch1 mode available after sync");
  }

  // Per-channel forwarders behave
  const uint32_t profile_size =
      specbleach_stereo_get_noise_profile_size(stereo);
  float* ch0_profile =
      specbleach_stereo_get_noise_profile_for_channel(stereo, 0, 1);
  TEST_ASSERT(ch0_profile != NULL, "ch0 profile pointer valid");
  TEST_ASSERT(
      specbleach_stereo_get_profile_block_count_for_channel(stereo, 0, 1) > 0,
      "ch0 block count positive");

  // Loading ch0's profile into ch1 explicitly must succeed
  TEST_ASSERT(
      specbleach_stereo_load_noise_profile_for_channel(
          stereo, 1, ch0_profile, profile_size,
          specbleach_stereo_get_profile_block_count_for_channel(stereo, 0, 1),
          1),
      "explicit cross-channel profile copy");

  specbleach_stereo_free(stereo);
  printf("✓ Stereo processing and sync tests passed\n");
}

int main(void) {
  test_init_and_free();
  test_parameter_loading();
  test_process_and_sync();

  printf("✅ All specbleach stereo tests passed!\n");
  return 0;
}
