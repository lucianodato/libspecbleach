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
 * @file stereo_denoiser_demo.c
 * @brief Multi-channel integration reference using the optional extras
 *        layer (specbleach_stereo).
 *
 * This example answers the question: "libspecbleach engines are
 * single-channel — how do I build a stereo/surround denoiser with the
 * conveniences the Noise Repellent plugin is known for?" The answer is the
 * extras layer, a thin orchestration library on top of the core. Nothing in
 * extras does DSP; it only manages engine instances and profile state.
 *
 * ===========================================================================
 * THE EXTRAS MODULES (what problem each solves)
 * ===========================================================================
 * - specbleach_stereo      Owns one core engine per channel. You load
 *                          parameters once and they fan out to every
 *                          channel; process() takes deinterleaved channel
 *                          pointer arrays. Also provides cross-channel
 *                          profile fallback fill and aggregated transient
 *                          reporting.
 *
 * ===========================================================================
 * DATA FLOW IN THIS EXAMPLE (phases are labeled in main())
 * ===========================================================================
 *   PHASE_LEARNING    the group learns the noise profile from the file
 *                     opening (same lifecycle as denoiser_demo.c).
 *   PHASE_REDUCTION   learn turns off (which FINALIZES all capture modes),
 *                     then the group reduces normally.
 *   (--switch-smoothing only) at the file midpoint the smoothing mode flips
 *                     from temporal to 2D NLM by simply reloading
 *                     parameters; the library performs the seamless
 *                     allocation-free crossfade internally and the reported
 *                     latency never changes.
 *
 * Build: cmake -B build -DENABLE_EXAMPLES=ON -DSPECBLEACH_BUILD_EXTRAS=ON \
 *            && cmake --build build
 * Link (downstream apps): pkg-config --cflags --libs specbleach-extras
 *                         or find_package(libspecbleach COMPONENTS extras)
 *                         -> libspecbleach::extras (pulls in the core)
 */

#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <sndfile.h>
#include <specbleach_stereo.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define BLOCK_SIZE 512
#define MAX_CHANNELS 8
#define DEFAULT_LEARN_FRAMES 8
#define DEFAULT_REDUCTION_DB 20.0f
#define DEFAULT_FRAME_SIZE_MS 46.0f

typedef enum DemoPhase {
  PHASE_LEARNING = 0,
  PHASE_REDUCTION = 1,
} DemoPhase;

typedef struct DemoOptions {
  float frame_size_ms;
  float reduction_db;
  uint32_t learn_frames;
  bool switch_smoothing;
} DemoOptions;

static void print_usage(const char* prog_name) {
  fprintf(stderr, "Usage: %s [options] <noisy input> <denoised output>\n",
          prog_name);
  fprintf(stderr, "Options:\n");
  fprintf(stderr,
          "  --reduction <val>     Reduction amount in dB (default: %.1f)\n",
          DEFAULT_REDUCTION_DB);
  fprintf(stderr, "  --frame-size <val>    Frame size in ms (default: %.1f)\n",
          DEFAULT_FRAME_SIZE_MS);
  fprintf(stderr,
          "  --learn-frames <val>  Learn blocks at file start (default: %u)\n",
          DEFAULT_LEARN_FRAMES);
  fprintf(stderr,
          "  --switch-smoothing    Switch temporal -> 2D NLM at the file "
          "midpoint\n");
  fprintf(stderr, "  --help                Show this help message\n");
}

static bool parse_uint_arg(const char* str, uint32_t* out) {
  if (!str || *str == '\0') {
    return false;
  }
  char* endptr = NULL;
  errno = 0;
  const long val = strtol(str, &endptr, 10);
  if (errno != 0 || endptr == str || *endptr != '\0' || val < 1 ||
      (unsigned long)val > UINT32_MAX) {
    return false;
  }
  *out = (uint32_t)val;
  return true;
}

static bool parse_float_arg(const char* str, float* out, const float min_val) {
  if (!str || *str == '\0') {
    return false;
  }
  char* endptr = NULL;
  errno = 0;
  const float val = strtof(str, &endptr);
  if (errno != 0 || endptr == str || *endptr != '\0' || !isfinite(val) ||
      val < min_val) {
    return false;
  }
  *out = val;
  return true;
}

int main(int argc, char** argv) {
  DemoOptions options = {DEFAULT_FRAME_SIZE_MS, DEFAULT_REDUCTION_DB,
                         DEFAULT_LEARN_FRAMES, false};

  static struct option long_options[] = {
      {"reduction", required_argument, 0, 'r'},
      {"frame-size", required_argument, 0, 'f'},
      {"learn-frames", required_argument, 0, 'n'},
      {"switch-smoothing", no_argument, 0, 's'},
      {"help", no_argument, 0, '?'},
      {0, 0, 0, 0}};

  int opt;
  while ((opt = getopt_long(argc, argv, "r:f:n:s", long_options, NULL)) != -1) {
    switch (opt) {
      case 'r':
        if (!parse_float_arg(optarg, &options.reduction_db, 0.0f)) {
          print_usage(argv[0]);
          return 1;
        }
        break;
      case 'f':
        if (!parse_float_arg(optarg, &options.frame_size_ms, 0.0001f)) {
          print_usage(argv[0]);
          return 1;
        }
        break;
      case 'n':
        if (!parse_uint_arg(optarg, &options.learn_frames)) {
          print_usage(argv[0]);
          return 1;
        }
        break;
      case 's':
        options.switch_smoothing = true;
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

  SF_INFO input_info = {0};
  SNDFILE* input_file = sf_open(input_file_name, SFM_READ, &input_info);
  if (!input_file) {
    fprintf(stderr, "Error: cannot open '%s': %s\n", input_file_name,
            sf_strerror(NULL));
    return 1;
  }

  const uint32_t channels = (uint32_t)input_info.channels;
  const uint32_t sample_rate = (uint32_t)input_info.samplerate;
  const uint64_t total_frames = (uint64_t)input_info.frames;

  if (channels > MAX_CHANNELS) {
    fprintf(stderr, "Error: this demo supports up to %u channels\n",
            MAX_CHANNELS);
    sf_close(input_file);
    return 1;
  }

  bool same_file = false;
  if (strcmp(input_file_name, output_file_name) == 0) {
    same_file = true;
  } else {
    struct stat stat_in;
    struct stat stat_out;
    if (stat(input_file_name, &stat_in) == 0 &&
        stat(output_file_name, &stat_out) == 0 &&
        stat_in.st_dev == stat_out.st_dev &&
        stat_in.st_ino == stat_out.st_ino) {
      same_file = true;
    }
  }
  if (same_file) {
    fprintf(stderr, "Error: input and output must be distinct files\n");
    sf_close(input_file);
    return 1;
  }

  SNDFILE* output_file = sf_open(output_file_name, SFM_WRITE, &input_info);
  if (!output_file) {
    fprintf(stderr, "Error: cannot create '%s': %s\n", output_file_name,
            sf_strerror(NULL));
    sf_close(input_file);
    return 1;
  }

  printf("Input: %u channels @ %u Hz, %.1f s\n", channels, sample_rate,
         (double)total_frames / (double)sample_rate);

  specbleach_stereo* group = specbleach_stereo_initialize(
      sample_rate, options.frame_size_ms, channels);
  const uint32_t latency = specbleach_stereo_get_latency(group);

  float* interleaved = calloc((size_t)channels * BLOCK_SIZE, sizeof(float));
  float* channel_data = calloc((size_t)channels * BLOCK_SIZE, sizeof(float));

  bool success = false;
  do {
    if (!group || !interleaved || !channel_data) {
      fprintf(stderr, "Error: initialization failed\n");
      break;
    }

    /* Learn phase parameters */
    SpecbleachDenoiserParameters params = {0};
    params.learn_noise = SPECBLEACH_LEARN_ALL;
    params.reduction_gain = powf(10.0f, -options.reduction_db / 20.0f);
    params.smoothing_mode = SB_SMOOTHING_TEMPORAL;
    params.hpss_enable = true;

    if (!specbleach_stereo_load_parameters(group, &params, sizeof(params))) {
      break;
    }

    uint64_t absolute_frame = 0;
    uint32_t learn_blocks_left = options.learn_frames;
    DemoPhase phase = learn_blocks_left > 0 ? PHASE_LEARNING : PHASE_REDUCTION;
    bool switch_done = false;
    success = true;

    while (success) {
      const sf_count_t read_frames =
          sf_readf_float(input_file, interleaved, BLOCK_SIZE);
      if (sf_error(input_file) != 0) {
        fprintf(stderr, "Error: failed to read input file: %s\n",
                sf_strerror(input_file));
        success = false;
        break;
      }
      if (read_frames == 0) {
        break; /* EOF */
      }

      /* Deinterleave into per-channel processing buffers */
      for (uint32_t ch = 0; ch < channels; ++ch) {
        for (sf_count_t s = 0; s < read_frames; ++s) {
          channel_data[((size_t)ch * BLOCK_SIZE) + (size_t)s] =
              interleaved[((size_t)s * channels) + ch];
        }
      }

      const float* in_ptrs[MAX_CHANNELS] = {0};
      float* out_ptrs[MAX_CHANNELS] = {0};
      for (uint32_t ch = 0; ch < channels; ++ch) {
        in_ptrs[ch] = &channel_data[(size_t)ch * BLOCK_SIZE];
        out_ptrs[ch] =
            &channel_data[(size_t)ch * BLOCK_SIZE]; /* process in place */
      }

      if (phase == PHASE_LEARNING) {
        --learn_blocks_left;
        if (!specbleach_stereo_process(group, (uint32_t)read_frames, in_ptrs,
                                       out_ptrs)) {
          success = false;
          break;
        }
        if (learn_blocks_left == 0) {
          /* Turning learning off finalizes every capture mode */
          params.learn_noise = SPECBLEACH_LEARN_OFF;
          if (!specbleach_stereo_load_parameters(group, &params,
                                                 sizeof(params))) {
            success = false;
            break;
          }
          printf("Noise profile learned after %u blocks\n",
                 options.learn_frames);
          phase = PHASE_REDUCTION;
        }
      } else {
        /* Midpoint: flip the smoothing mode by reloading parameters. The
         * library crossfades internally; profiles and latency are
         * untouched. In a real-time app run the parameter load on a
         * control/message thread, not in the audio callback. */
        if (options.switch_smoothing && !switch_done &&
            absolute_frame >= total_frames / 2U) {
          printf("Switching smoothing mode at t=%.1fs\n",
                 (double)absolute_frame / sample_rate);
          params.smoothing_mode = SB_SMOOTHING_NLM_2D;
          if (!specbleach_stereo_load_parameters(group, &params,
                                                 sizeof(params))) {
            success = false;
            break;
          }
          switch_done = true;
        }

        if (!specbleach_stereo_process(group, (uint32_t)read_frames, in_ptrs,
                                       out_ptrs)) {
          success = false;
          break;
        }
      }

      /* Latency is constant across smoothing mode switches; print it once */
      if (absolute_frame == 0) {
        printf("Latency %u samples (constant for every smoothing mode)\n",
               latency);
      }

      /* Re-interleave processed audio */
      for (uint32_t ch = 0; ch < channels; ++ch) {
        for (sf_count_t s = 0; s < read_frames; ++s) {
          interleaved[((size_t)s * channels) + ch] =
              channel_data[((size_t)ch * BLOCK_SIZE) + (size_t)s];
        }
      }

      if (sf_writef_float(output_file, interleaved, read_frames) !=
          read_frames) {
        success = false;
      }
      absolute_frame += (uint64_t)read_frames;
    }
  } while (0);

  free(interleaved);
  free(channel_data);
  specbleach_stereo_free(group);
  sf_close(input_file);
  sf_close(output_file);

  if (!success) {
    remove(output_file_name);
    return 1;
  }
  printf("Wrote '%s'\n", output_file_name);
  return 0;
}
