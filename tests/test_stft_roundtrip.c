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

// STFT identity round-trip: null processing must reconstruct exactly.
// Covers the legacy geometry (512/4x) and the low-latency one (512/8x)
// including the fft-vs-frame term of the OLA normalizer. Kill criterion
// for the 8x low-latency reconfiguration.

#include "../src/shared/stft/stft_processor.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Repo convention: hard failure in every build type (Release defines NDEBUG,
// which silences bare assert()).
#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define SR 48000u
#define FRAME_MS (512.0F * 1000.0F / 48000.0F)

static bool identity_proc(SpectralProcessorHandle handle, float* spectrum) {
  (void)handle;
  (void)spectrum;
  return true;
}

static float max_err_db(const float* in, const float* out, uint32_t n,
                        uint32_t latency, uint32_t skip) {
  float peak = 0.0F;
  for (uint32_t i = latency + skip; i < n; ++i) {
    const float e = fabsf(out[i] - in[i - latency]);
    if (e > peak) {
      peak = e;
    }
  }
  return 20.0F * log10f(peak + 1e-12F);
}

static void run_config(const char* name, uint32_t overlap, ZeroPaddingType pad,
                       uint32_t pad_amount, uint32_t expect_fft,
                       uint32_t expect_hop) {
  StftProcessor* stft = stft_processor_initialize(
      SR, FRAME_MS, overlap, pad, pad_amount, HANN_WINDOW, HANN_WINDOW);
  TEST_ASSERT(stft != NULL, "STFT init must succeed");
  TEST_ASSERT(get_stft_fft_size(stft) == expect_fft, "unexpected FFT size");
  TEST_ASSERT(get_stft_hop_size(stft) == expect_hop, "unexpected hop");
  TEST_ASSERT(get_stft_latency(stft) == 512u, "latency must equal the frame");

  const uint32_t n = 4u * SR;
  const uint32_t flush = SR / 2u;
  float* in = (float*)calloc(n + flush, sizeof(float));
  float* out = (float*)calloc(n + flush, sizeof(float));
  TEST_ASSERT(in != NULL && out != NULL, "buffer alloc must succeed");
  srand(4242);
  for (uint32_t i = 0; i < n; ++i) {
    const float t = (float)i / (float)SR;
    in[i] = 0.4F * sinf(2.0F * (float)M_PI * 1000.0F * t) +
            0.1F * ((float)rand() / (float)RAND_MAX * 2.0F - 1.0F);
  }
  in[SR] += 0.5F; // transient click

  // Odd block sizes to prove chunking independence.
  uint32_t done = 0;
  while (done < n + flush) {
    uint32_t chunk = n + flush - done;
    if (chunk > 300) {
      chunk = 300;
    }
    TEST_ASSERT(stft_processor_run(stft, chunk, &in[done], &out[done],
                                   identity_proc, NULL),
                "STFT run must succeed");
    done += chunk;
  }

  const float err = max_err_db(in, out, n + flush, 512u, 512u);
  printf("%s: round-trip error %.1f dBFS\n", name, err);
  TEST_ASSERT(err < -100.0F, "identity round-trip must reconstruct exactly");
  free(in);
  free(out);
  stft_processor_free(stft);
  stft_processor_free(NULL);
}

int main(void) {
  // Current low-latency geometry.
  run_config("512/4x/512fft", 4u, PAD_TO_VALID_SIZE, 0u, 512u, 128u);
  // Low-latency geometry: 8x overlap, standard padding.
  run_config("512/8x/512fft", 8u, PAD_TO_VALID_SIZE, 0u, 512u, 64u);
  printf("STFT round-trip tests passed.\n");
  return 0;
}
