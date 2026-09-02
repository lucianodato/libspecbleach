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

/**
 * @file denoiser_demo.c
 * @brief Minimal single-channel integration reference for the core API.
 *
 * This is the smallest complete program that denoises audio with
 * libspecbleach. Read this one FIRST; it demonstrates the fundamental
 * lifecycle every integration follows, regardless of how fancy the
 * surrounding application is.
 *
 * ===========================================================================
 * INTEGRATION LIFECYCLE (each step maps to a labeled section in main()):
 * ===========================================================================
 *   1. CREATE     specbleach_denoiser_initialize(sample_rate, frame_size_ms)
 *   2. CONFIGURE  specbleach_denoiser_load_parameters(handle, &params,
 *                                                     sizeof(params))
 *   3. LEARN      process blocks while params.learn_noise is
 *                 SPECBLEACH_LEARN_ALL; the engine accumulates a noise
 *                 profile from whatever you feed it
 *   4. FINALIZE   load parameters again with learn_noise =
 *                 SPECBLEACH_LEARN_OFF. IMPORTANT: profile capture modes are
 *                 only finalized when learning turns OFF, so profiles are
 *                 NOT fully usable until this happens (and one block has
 *                 been processed afterwards).
 *   5. REDUCE     keep calling specbleach_denoiser_process() per block of
 *                 audio for as long as you have audio
 *   6. CLEANUP    specbleach_denoiser_free(handle)
 *
 * ===========================================================================
 * CONTRACTS WORTH KNOWING (all documented in include/specbleach_denoiser.h):
 * ===========================================================================
 * - Block size: the library accepts ANY number of samples per call; it
 *   buffers internally. You do not need to align to frames or hops.
 * - Buffers: input/output are plain float arrays of `number_of_samples`
 *   length per call. Output may safely alias the input buffer.
 * - Parameters: load_parameters() copies everything it needs, including an
 *   optional reduction curve (see reduction_curve_bias ownership docs in
 *   the header). After the call returns you may free or reuse your structs.
 * - Threading: calls on the same handle must not run concurrently with
 *   each other or with _process(). For real-time use, load_parameters()
 *   is allocation-free EXCEPT on the first call after enabling
 *   reduction_curve_enabled.
 * - Latency: query specbleach_denoiser_get_latency() after initialization
 *   and report it to your host/caller if delay compensation matters.
 *
 * Build: cmake -B build -DENABLE_EXAMPLES=ON && cmake --build build
 * Link (downstream apps): pkg-config --cflags --libs specbleach
 *                         or find_package(libspecbleach) ->
 * libspecbleach::libspecbleach
 *
 * For multi-channel processing, click-free engine switching, and noise
 * profile migration between engines, see stereo_denoiser_demo.c which uses
 * the optional extras layer.
 */

#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <sndfile.h>
#include <specbleach_denoiser.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// This is not a deliberate value. The library handles any amount passed
// through a circular buffer
#define BLOCK_SIZE 512
#define NOISE_FRAMES                                                           \
  8 // Amount of frames to capture profile at the beginning of the file (can be
    // anywhere)
#define FRAME_SIZE 50

static bool parse_float_arg(const char* str, float* out, float min_val,
                            float max_val) {
  if (!str || *str == '\0') {
    return false;
  }
  char* endptr = NULL;
  errno = 0;
  float val = strtof(str, &endptr);
  if (errno != 0 || endptr == str || *endptr != '\0') {
    return false;
  }
  if (!isfinite(val) || val < min_val || val > max_val) {
    return false;
  }
  *out = val;
  return true;
}

static bool parse_int_arg(const char* str, int* out, int min_val, int max_val) {
  if (!str || *str == '\0') {
    return false;
  }
  char* endptr = NULL;
  errno = 0;
  long val = strtol(str, &endptr, 10);
  if (errno != 0 || endptr == str || *endptr != '\0') {
    return false;
  }
  if (val < min_val || val > max_val) {
    return false;
  }
  *out = (int)val;
  return true;
}

static bool parse_uint32_arg(const char* str, uint32_t* out, uint32_t min_val,
                             uint32_t max_val) {
  if (!str || *str == '\0') {
    return false;
  }
  char* endptr = NULL;
  errno = 0;
  long val = strtol(str, &endptr, 10);
  if (errno != 0 || endptr == str || *endptr != '\0') {
    return false;
  }
  if (val < 0 || (unsigned long)val < min_val || (unsigned long)val > max_val) {
    return false;
  }
  *out = (uint32_t)val;
  return true;
}

static void print_usage(const char* prog_name) {
  fprintf(stderr, "Usage: %s [options] <noisy input> <denoised output>\n",
          prog_name);
  fprintf(stderr, "Options:\n");
  fprintf(
      stderr,
      "  --reduction <val>          Reduction amount in dB (default: 20.0)\n");
  fprintf(stderr,
          "  --whitening <val>          Whitening factor (default: 50.0)\n");
  fprintf(stderr,
          "  --smoothing <val>          Smoothing factor (default: 0.0)\n");
  fprintf(stderr,
          "  --smoothing-mode <val>     Smoothing mode: 0=temporal, 1=NLM 2D "
          "(default: 0)\n");
  fprintf(
      stderr,
      "  --masking-depth <val>      Masking depth (0.0-1.0, default: 0.5)\n");
  fprintf(stderr,
          "  --steering-response <val>  Steering response (-1.0-1.0, default: "
          "1.0)\n");
  fprintf(stderr,
          "  --adaptive                 Enable adaptive noise estimation\n");
  fprintf(stderr,
          "  --noise-method <val>      Noise estimation method (0-2, default: "
          "0)\n");
  fprintf(stderr,
          "  --frame-size <val>        Frame size in ms (default: 46.0)\n");
  fprintf(stderr,
          "  --learn-frames <val>      Number of learn frames (default: 8)\n");
  fprintf(stderr, "  --help                    Show this help message\n");
}

static void cleanup_resources(SF_INFO* sfinfo, SNDFILE* input_file,
                              SNDFILE* output_file, float* input_buffer,
                              float* output_buffer,
                              specbleach_denoiser* lib_instance) {
  if (input_file) {
    sf_close(input_file);
  }
  if (output_file) {
    sf_close(output_file);
  }
  if (sfinfo) {
    free(sfinfo);
  }
  if (input_buffer) {
    free(input_buffer);
  }
  if (output_buffer) {
    free(output_buffer);
  }
  if (lib_instance) {
    specbleach_denoiser_free(lib_instance);
  }
}

int main(int argc, char** argv) {
  SpecbleachDenoiserParameters parameters = (SpecbleachDenoiserParameters){
      .residual_listen = false,
      .learn_noise = SPECBLEACH_LEARN_ALL, // Learn all modes
      .aggressiveness = 1.0f,              // Use maximum mode for processing
      .reduction_gain = 0.1f,              // 20 dB
      .smoothing_factor = 0.0f,
      .whitening_factor = 0.5f,
      .masking_depth = 0.5f,
      .tonal_reduction_gain = 1.0f,
      .hpss_enable = true};

  static struct option long_options[] = {
      {"reduction", required_argument, 0, 'r'},
      {"whitening", required_argument, 0, 'w'},
      {"smoothing", required_argument, 0, 's'},
      {"smoothing-mode", required_argument, 0, 'S'},
      {"masking-depth", required_argument, 0, 'd'},
      {"steering-response", required_argument, 0, 'l'},
      {"adaptive", no_argument, 0, 'a'},
      {"noise-method", required_argument, 0, 'm'},
      {"frame-size", required_argument, 0, 'f'},
      {"learn-frames", required_argument, 0, 'n'},
      {"help", no_argument, 0, '?'},
      {0, 0, 0, 0}};

  float frame_size_ms = FRAME_SIZE;
  uint32_t learn_frames = NOISE_FRAMES;
  int opt;
  while ((opt = getopt_long(argc, argv, "r:w:s:S:d:l:am:f:n:", long_options,
                            NULL)) != -1) {
    switch (opt) {
      case 'r': {
        float red_db;
        if (!parse_float_arg(optarg, &red_db, 0.0f, INFINITY)) {
          print_usage(argv[0]);
          return 1;
        }
        parameters.reduction_gain = powf(10.0f, -fabsf(red_db) / 20.0f);
        break;
      }
      case 'w': {
        float whitening;
        if (!parse_float_arg(optarg, &whitening, 0.0f, 100.0f)) {
          print_usage(argv[0]);
          return 1;
        }
        parameters.whitening_factor = whitening / 100.0f;
        break;
      }
      case 's': {
        float smoothing;
        if (!parse_float_arg(optarg, &smoothing, 0.0f, 100.0f)) {
          print_usage(argv[0]);
          return 1;
        }
        parameters.smoothing_factor = smoothing / 100.0f;
        break;
      }
      case 'S': {
        int mode;
        if (!parse_int_arg(optarg, &mode, 0, 1)) {
          print_usage(argv[0]);
          return 1;
        }
        parameters.smoothing_mode =
            (mode == 1) ? SB_SMOOTHING_NLM_2D : SB_SMOOTHING_TEMPORAL;
        break;
      }
      case 'd':
        if (!parse_float_arg(optarg, &parameters.masking_depth, 0.0f, 1.0f)) {
          print_usage(argv[0]);
          return 1;
        }
        break;
      case 'l':
        if (!parse_float_arg(optarg, &parameters.aggressiveness, -1.0f, 1.0f)) {
          print_usage(argv[0]);
          return 1;
        }
        break;
      case 'a':
        parameters.adaptive_noise = true;
        break;
      case 'm': {
        int noise_method;
        if (!parse_int_arg(optarg, &noise_method, 0, 2)) {
          print_usage(argv[0]);
          return 1;
        }
        parameters.noise_estimation_method =
            (SpecbleachNoiseEstimationMethod)noise_method;
        break;
      }
      case 'f':
        if (!parse_float_arg(optarg, &frame_size_ms, 0.0001f, INFINITY)) {
          print_usage(argv[0]);
          return 1;
        }
        break;
      case 'n':
        if (!parse_uint32_arg(optarg, &learn_frames, 1, UINT32_MAX)) {
          print_usage(argv[0]);
          return 1;
        }
        break;
      case '?':
      default:
        print_usage(argv[0]);
        return 1;
    }
  }

  if (argc - optind != 2) {
    print_usage(argv[0]);
    return 1;
  }

  const char* input_file_name = argv[optind];
  const char* output_file_name = argv[optind + 1];

  SF_INFO* sfinfo = NULL;
  SNDFILE* input_file = NULL;
  SNDFILE* output_file = NULL;
  float* input_library_buffer = NULL;
  float* output_library_buffer = NULL;
  specbleach_denoiser* lib_instance = NULL;
  int ret = 1;

  do {
    // Allocate memory for SF_INFO
    sfinfo = (SF_INFO*)calloc(1, sizeof(SF_INFO));
    if (!sfinfo) {
      fprintf(stderr, "Error: Failed to allocate memory for SF_INFO\n");
      break;
    }

    // Open input file
    input_file = sf_open(input_file_name, SFM_READ, sfinfo);
    if (!input_file) {
      fprintf(stderr, "Error: Failed to open input file '%s': %s\n",
              input_file_name, sf_strerror(NULL));
      break;
    }

    // Validate audio format
    if (sfinfo->channels != 1) {
      fprintf(stderr,
              "Error: Only mono audio is supported (file has %d "
              "channels)\n",
              sfinfo->channels);
      break;
    }

    // Open output file
    output_file = sf_open(output_file_name, SFM_WRITE, sfinfo);
    if (!output_file) {
      fprintf(stderr, "Error: Failed to open output file '%s': %s\n",
              output_file_name, sf_strerror(NULL));
      break;
    }

    // Allocate buffers
    input_library_buffer = (float*)calloc(BLOCK_SIZE, sizeof(float));
    if (!input_library_buffer) {
      fprintf(stderr, "Error: Failed to allocate input buffer\n");
      break;
    }

    output_library_buffer = (float*)calloc(BLOCK_SIZE, sizeof(float));
    if (!output_library_buffer) {
      fprintf(stderr, "Error: Failed to allocate output buffer\n");
      break;
    }

    // ------------------------------------------------------------------
    // STEP 1: CREATE the engine instance.
    // frame_size_ms controls the STFT window (20-100 ms recommended). It
    // trades frequency resolution against latency. The handle is opaque;
    // pass it back to every specbleach_denoiser_* call from now on.
    // Initialize library instance
    lib_instance = specbleach_denoiser_initialize((uint32_t)sfinfo->samplerate,
                                                  frame_size_ms);
    if (!lib_instance) {
      fprintf(stderr, "Error: Failed to initialize library instance\n");
      break;
    }

    // ------------------------------------------------------------------
    // STEP 2 + 3: CONFIGURE and LEARN.
    //
    // While learn_noise is SPECBLEACH_LEARN_ALL, every processed block
    // updates an internal noise profile instead of applying reduction.
    // Feed it a section containing ONLY the noise you want removed
    // (e.g., the first seconds of a recording before speech starts).
    // NOISE PROFILE LEARN STAGE

    // Load the parameters before doing the denoising or profile learning
    if (!specbleach_denoiser_load_parameters(lib_instance, &parameters,
                                             sizeof(parameters))) {
      fprintf(stderr, "Error: Failed to load parameters\n");
      break;
    }

    // Iterate over some frames (learn_frames) at the beginning of the audio to
    // capture the noise profile
    for (uint32_t i = 0; i < learn_frames; i++) {
      sf_count_t frames_read =
          sf_readf_float(input_file, input_library_buffer, BLOCK_SIZE);
      if (frames_read <= 0) {
        if (frames_read < 0) {
          fprintf(stderr,
                  "Error: Failed to read audio during noise profile learning: "
                  "%s\n",
                  sf_strerror(input_file));
        } else {
          fprintf(stderr,
                  "Warning: End of file reached before capturing noise "
                  "profile\n");
        }
        break;
      }

      // Process the audio to learn the noise profile
      if (!specbleach_denoiser_process(lib_instance, (uint32_t)BLOCK_SIZE,
                                       input_library_buffer,
                                       output_library_buffer)) {
        fprintf(
            stderr,
            "Error: Failed to process audio during noise profile learning\n");
        break;
      }
    }

    // If we broke out of the learn stage due to an error, stop.
    // In adaptive mode, we can proceed even without a pre-learned profile.
    if (!parameters.adaptive_noise &&
        !specbleach_denoiser_noise_profile_available_for_mode(lib_instance,
                                                              1)) {
      fprintf(stderr, "Error: Noise profile was not successfully learned\n");
      break;
    }

    // ------------------------------------------------------------------
    // STEP 4: FINALIZE the profile by turning learning OFF.
    //
    // This reload is not optional bookkeeping: capture modes are only
    // finalized on the learn -> off transition, so reduction quality is
    // undefined until this call happens. Real-time integrations (mixers,
    // plugins) do exactly this when the user releases a "learn" button.
    // NOISE REDUCTION STAGE

    // Turn off noise profile learning to start applying reduction
    parameters.learn_noise = SPECBLEACH_LEARN_OFF;

    // Reload parameters with noise learn off
    if (!specbleach_denoiser_load_parameters(lib_instance, &parameters,
                                             sizeof(parameters))) {
      fprintf(stderr, "Error: Failed to reload parameters\n");
      break;
    }

    // ------------------------------------------------------------------
    // STEP 5: REDUCE. From here on every process() call denoises using
    // the captured profile. In a real-time application this loop is your
    // audio callback: same call shape, any block size. See the threading
    // and latency notes in the file header above.
    // Iterate over the audio to apply denoising
    sf_count_t frames_read;
    while ((frames_read = sf_readf_float(input_file, input_library_buffer,
                                         BLOCK_SIZE)) > 0) {
      // Process the audio
      if (!specbleach_denoiser_process(lib_instance, (uint32_t)BLOCK_SIZE,
                                       input_library_buffer,
                                       output_library_buffer)) {
        fprintf(stderr, "Error: Failed to process audio\n");
        break;
      }

      // Write processed audio
      sf_count_t frames_written =
          sf_writef_float(output_file, output_library_buffer, frames_read);
      if (frames_written != frames_read) {
        fprintf(stderr,
                "Error: Failed to write all frames (wrote %ld of %ld)\n",
                (long)frames_written, (long)frames_read);
        break;
      }
    }

    // Check for read errors
    if (frames_read < 0) {
      fprintf(stderr, "Error: Failed to read audio: %s\n",
              sf_strerror(input_file));
      break;
    }

    // ------------------------------------------------------------------
    // STEP 6: CLEANUP frees all engine memory. Profiles that should survive
    // a restart must be serialized BEFORE this call: copy the arrays from
    // specbleach_denoiser_get_noise_profile_for_mode() for each mode in
    // [SPECBLEACH_PROFILE_MODE_FIRST, SPECBLEACH_PROFILE_MODE_LAST], plus
    // the block counts, and restore them later with
    // specbleach_denoiser_load_noise_profile_for_mode(). See the plugin
    // source (noise-repellent repo) for a complete serialization example.
    // Success
    ret = 0;
  } while (0);

  cleanup_resources(sfinfo, input_file, output_file, input_library_buffer,
                    output_library_buffer, lib_instance);
  return ret;
}
