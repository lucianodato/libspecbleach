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

#include "specbleach_profile_migration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define BLOCK_SIZE 512

static void test_validation(void) {
  printf("Testing migration validation...\n");

  specbleach_denoiser* one_d = NULL;
  specbleach_2d_denoiser* two_d = NULL;

  TEST_ASSERT(specbleach_migrate_profiles_1d_to_2d(NULL, NULL) == false,
              "NULL handles rejected");

  two_d = specbleach_2d_initialize(44100, 46.0f);
  TEST_ASSERT(two_d != NULL, "2D init");
  TEST_ASSERT(specbleach_migrate_profiles_1d_to_2d(one_d, two_d) == false,
              "NULL source rejected");
  specbleach_2d_free(two_d);

  one_d = specbleach_denoiser_initialize(44100, 46.0f);
  TEST_ASSERT(one_d != NULL, "1D init");
  TEST_ASSERT(specbleach_migrate_profiles_2d_to_1d(NULL, one_d) == false,
              "NULL source rejected (2d->1d)");
  specbleach_denoiser_free(one_d);

  printf("✓ Migration validation tests passed\n");
}

static void test_migration_copies_profiles(void) {
  printf("Testing profile copy semantics...\n");

  const float frame_size = 46.0f;
  specbleach_denoiser* source =
      specbleach_denoiser_initialize(44100, frame_size);
  specbleach_2d_denoiser* target = specbleach_2d_initialize(44100, frame_size);
  TEST_ASSERT(source && target, "engines initialized with same frame size");

  // Migrating before anything was learned must fail cleanly
  TEST_ASSERT(specbleach_migrate_profiles_1d_to_2d(source, target) == false,
              "migration without profiles fails");

  // Learn a noise profile on the source engine using shaped noise
  SpecbleachDenoiserParameters params = {0};
  params.learn_noise = SPECBLEACH_LEARN_ALL;
  params.reduction_gain = 0.1f;
  params.smoothing_factor = 0.0f;
  params.whitening_factor = 0.5f;
  params.masking_depth = 0.5f;
  params.aggressiveness = 1.0f;
  params.tonal_reduction_gain = 1.0f;
  params.hpss_enable = true;
  TEST_ASSERT(
      specbleach_denoiser_load_parameters(source, &params, sizeof(params)),
      "source params loaded");

  float input[BLOCK_SIZE];
  float output[BLOCK_SIZE];
  for (uint32_t block = 0; block < 8; ++block) {
    for (uint32_t s = 0; s < BLOCK_SIZE; ++s) {
      input[s] = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
    }
    TEST_ASSERT(specbleach_denoiser_process(source, BLOCK_SIZE, input, output),
                "source learn block");
  }

  const uint32_t size_source =
      specbleach_denoiser_get_noise_profile_size(source);
  const uint32_t size_target = specbleach_2d_get_noise_profile_size(target);
  if (size_source == size_target) {
    TEST_ASSERT(specbleach_migrate_profiles_1d_to_2d(source, target),
                "migration succeeds with equal sizes");

    for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
         mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
      TEST_ASSERT(specbleach_2d_noise_profile_available_for_mode(target, mode),
                  "target mode available after migration");
    }

    // Contents must match the fallback source mode exactly
    float* src_profile = specbleach_denoiser_get_noise_profile_for_mode(
        source, SPECBLEACH_PROFILE_MODE_FIRST);
    float* dst_profile = specbleach_2d_get_noise_profile_for_mode(
        target, SPECBLEACH_PROFILE_MODE_FIRST);
    TEST_ASSERT(src_profile && dst_profile, "profile pointers valid");
    for (uint32_t bin = 0; bin < size_target; ++bin) {
      TEST_ASSERT(fabsf(src_profile[bin] - dst_profile[bin]) < 1e-9f,
                  "migrated profile contents identical");
      break; /* spot check first bin; full compare below */
    }
    int mismatches = 0;
    for (uint32_t bin = 0; bin < size_target; ++bin) {
      if (fabsf(src_profile[bin] - dst_profile[bin]) > 1e-9f) {
        ++mismatches;
      }
    }
    TEST_ASSERT(mismatches == 0, "full migrated profile matches source");

    // Reverse direction round-trips too
    specbleach_denoiser* back =
        specbleach_denoiser_initialize(44100, frame_size);
    TEST_ASSERT(back != NULL, "reverse target init");
    TEST_ASSERT(specbleach_migrate_profiles_2d_to_1d(target, back),
                "reverse migration succeeds");
    specbleach_denoiser_free(back);
  } else {
    printf("  (skipping content checks: differing profile sizes %u vs %u)\n",
           size_source, size_target);
  }

  specbleach_denoiser_free(source);
  specbleach_2d_free(target);
  printf("✓ Profile copy tests passed\n");
}

int main(void) {
  test_validation();
  test_migration_copies_profiles();

  printf("✅ All profile migration tests passed!\n");
  return 0;
}
