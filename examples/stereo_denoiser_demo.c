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
 *        layer (specbleach_stereo + specbleach_transition).
 *
 * This example answers the question: "libspecbleach engines are
 * single-channel — how do I build a stereo/surround denoiser with the
 * conveniences the Noise Repellent plugin is known for?" The answer is the
 * extras layer, a thin orchestration library on top of the core. Nothing in
 * extras does DSP; it only manages engine instances and transitions.
 *
 * ===========================================================================
 * THE THREE EXTRAS MODULES (what problem each solves)
 * ===========================================================================
 * - specbleach_stereo      Owns one core engine per channel. You load
 *                          parameters once and they fan out to every
 *                          channel; process() takes deinterleaved channel
 *                          pointer arrays. Also provides cross-channel
 *                          profile fallback fill and aggregated transient
 *                          reporting.
 * - specbleach_transition  Click-free switching between two ALREADY
 *                          CONFIGURED engine groups. It aligns their
 *                          latencies with an internal delay line, applies an
 *                          equal-power fade over 40 ms, and — when landing
 *                          on the lower-latency engine — ramps its
 *                          compensation delay back out so listeners never
 *                          hear the latency jump. It knows nothing about
 *                          parameters: YOU decide what "switching" means.
 * - specbleach_migrate_profiles_from / specbleach_migrate_profiles_1d_to_2d
 *                          Copies learned noise profiles between groups so
 *                          the new engine keeps the noise the user already
 *                          taught the old one. Without this, every switch
 *                          would sound like the denoiser forgot everything.
 *
 * ===========================================================================
 * DATA FLOW IN THIS EXAMPLE (phases are labeled in main())
 * ===========================================================================
 *   PHASE_LEARNING    primary group learns the noise profile from the file
 *                     opening (same lifecycle as denoiser_demo.c).
 *   PHASE_REDUCTION   learn turns off (which FINALIZES all capture modes),
 *                     then the group reduces normally.
 *   PHASE_TRANSITION  (--switch-engine only) at the file midpoint:
 *                     1. create/prepare the second engine family's group
 *                     2. migrate profiles into it
 *                     3. transition_begin(old_latency, new_latency)
 *                     4. per block: render BOTH groups' wet outputs, hand
 *                        them to specbleach_transition_process(), write the
 *                        blend.
 *   After the fade completes, transition_process() passes the target group
 *   through untouched, so you can keep calling it or drop to single-group
 *   rendering once specbleach_transition_active() returns false.
 *
 * ===========================================================================
 * WHAT A REAL-TIME INTEGRATOR MUST ADD (this demo is offline for brevity)
 * ===========================================================================
 * - Latency reporting: whenever you call transition_begin(), immediately
 *   report specbleach_transition_get_latency() to your host (it equals the
 *   target group's latency from that moment). See the "Latency now ..."
 *   print in main(): that printf stands in for setLatencySamples().
 * - Threading: transition_begin() may allocate its delay line (first use),
 *   and profile migration touches engine state — both belong on a control/
 *   message thread, NOT in the audio callback. process() calls are
 *   real-time safe. For reference, see how the Noise Repellent plugin
 *   queues switch requests to a timer-driven worker with a pause-gate
 *   handshake (Source/PluginProcessor.cpp in the noise-repellent repo).
 * - Parameter changes do NOT require transitions. Transitions are for
 *   swapping ENGINE FAMILIES (spectral <-> 2D NLM) where latency differs.
 *   Plain parameter updates can be loaded directly from the audio thread
 *   (allocation caveat: first load after enabling reduction_curve_enabled).
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
#include <specbleach_transition.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define BLOCK_SIZE 512
#define MAX_CHANNELS 8
#define DEFAULT_LEARN_FRAMES 8
#define DEFAULT_REDUCTION_DB 20.0f
#define DEFAULT_FRAME_SIZE_MS 46.0f

typedef enum DemoPhase {
  PHASE_LEARNING = 0,
  PHASE_REDUCTION = 1,
  PHASE_TRANSITION = 2,
} DemoPhase;

typedef struct DemoOptions {
  float frame_size_ms;
  float reduction_db;
  uint32_t learn_frames;
  bool switch_engine;
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
          "  --switch-engine       Fade from spectral to 2D NLM at the file "
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
  if (errno != 0 || endptr == str || *endptr != '\0' || val < 1) {
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
      {"switch-engine", no_argument, 0, 's'},
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
        options.switch_engine = true;
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

  SNDFILE* output_file = sf_open(output_file_name, SFM_WRITE, &input_info);
  if (!output_file) {
    fprintf(stderr, "Error: cannot create '%s': %s\n", output_file_name,
            sf_strerror(NULL));
    sf_close(input_file);
    return 1;
  }

  const uint32_t channels = (uint32_t)input_info.channels;
  const uint32_t sample_rate = (uint32_t)input_info.samplerate;
  const uint64_t total_frames = (uint64_t)input_info.frames;

  if (channels > MAX_CHANNELS) {
    fprintf(stderr, "Error: this demo supports up to %u channels\n",
            MAX_CHANNELS);
    sf_close(input_file);
    sf_close(output_file);
    return 1;
  }

  printf("Input: %u channels @ %u Hz, %.1f s\n", channels, sample_rate,
         (double)total_frames / (double)sample_rate);

  specbleach_stereo* spectral_group =
      specbleach_stereo_initialize(sample_rate, options.frame_size_ms, channels,
                                   SPECBLEACH_STEREO_ENGINE_SPECTRAL);
  specbleach_stereo* nlm_group =
      options.switch_engine ? specbleach_stereo_initialize(
                                  sample_rate, options.frame_size_ms, channels,
                                  SPECBLEACH_STEREO_ENGINE_NLM_2D)
                            : NULL;
  specbleach_transition* transition =
      specbleach_transition_initialize(sample_rate, BLOCK_SIZE, channels);

  float* interleaved = calloc((size_t)channels * BLOCK_SIZE, sizeof(float));
  float* channel_data = calloc((size_t)channels * BLOCK_SIZE, sizeof(float));
  float* scratch_from = calloc((size_t)channels * BLOCK_SIZE, sizeof(float));

  bool success = false;
  do {
    if (!spectral_group || !transition || !interleaved || !channel_data ||
        !scratch_from || (options.switch_engine && !nlm_group)) {
      fprintf(stderr, "Error: initialization failed\n");
      break;
    }

    /* Learn phase parameters */
    SpecbleachDenoiserParameters params_1d = {0};
    params_1d.learn_noise = SPECBLEACH_LEARN_ALL;
    params_1d.reduction_gain = powf(10.0f, -options.reduction_db / 20.0f);
    params_1d.hpss_enable = true;

    Specbleach2DDenoiserParameters params_2d = {0};
    params_2d.learn_noise = SPECBLEACH_LEARN_ALL;
    params_2d.reduction_gain = powf(10.0f, -options.reduction_db / 20.0f);
    params_2d.hpss_enable = true;

    /* Learn on the primary (spectral) group */
    if (!specbleach_stereo_load_parameters_1d(spectral_group, &params_1d,
                                              sizeof(params_1d))) {
      break;
    }

    uint64_t absolute_frame = 0;
    uint32_t learn_blocks_left = options.learn_frames;
    DemoPhase phase = learn_blocks_left > 0 ? PHASE_LEARNING : PHASE_REDUCTION;
    uint32_t reported_latency = 0;
    success = true;

    while (success) {
      const sf_count_t read_frames =
          sf_readf_float(input_file, interleaved, BLOCK_SIZE);
      if (read_frames <= 0) {
        break; /* EOF */
      }

      /* Deinterleave into per-channel processing buffers */
      for (uint32_t ch = 0; ch < channels; ++ch) {
        for (sf_count_t s = 0; s < read_frames; ++s) {
          channel_data[ch * BLOCK_SIZE + s] = interleaved[s * channels + ch];
        }
      }

      const float* in_ptrs[MAX_CHANNELS] = {0};
      float* out_ptrs[MAX_CHANNELS] = {0};
      float* from_ptrs[MAX_CHANNELS] = {0};
      for (uint32_t ch = 0; ch < channels; ++ch) {
        in_ptrs[ch] = &channel_data[ch * BLOCK_SIZE];
        out_ptrs[ch] = &channel_data[ch * BLOCK_SIZE]; /* process in place */
        from_ptrs[ch] = &scratch_from[ch * BLOCK_SIZE];
      }

      /* -----------------------------------------------------------------
       * PHASE_LEARNING: identical lifecycle to denoiser_demo.c, but one
       * specbleach_stereo call drives every channel at once. Each channel
       * learns its OWN profile from its own input (stereo noise usually
       * differs between ears). */
      if (phase == PHASE_LEARNING) {
        --learn_blocks_left;
        if (!specbleach_stereo_process(spectral_group, (uint32_t)read_frames,
                                       in_ptrs, out_ptrs)) {
          success = false;
          break;
        }
        if (learn_blocks_left == 0) {
          /* Turning learning off finalizes every capture mode */
          params_1d.learn_noise = SPECBLEACH_LEARN_OFF;
          if (!specbleach_stereo_load_parameters_1d(spectral_group, &params_1d,
                                                    sizeof(params_1d))) {
            success = false;
            break;
          }
          printf("Noise profile learned after %u blocks\n",
                 options.learn_frames);
          phase = PHASE_REDUCTION;
        }
        /* Fall through to re-interleave below */
      } else if (phase == PHASE_REDUCTION) {
        /* -----------------------------------------------------------------
         * PHASE_TRANSITION ENTRY (--switch-engine only): the three-step
         * switch recipe. Order matters:
         *   1. load parameters into the target group (learn OFF: it must
         *      reduce, not capture)
         *   2. migrate profiles so the target "remembers" the noise
         *   3. transition_begin() with BOTH groups' latencies; from this
         *      call onward hosts must report the TARGET's latency (see the
         *      transition header for why that is safe even mid-fade).
         * In a real-time app, steps 1-3 run on a control thread with a
         * pause gate around step 2 — see PluginProcessor.cpp. */
        if (options.switch_engine && absolute_frame >= total_frames / 2U) {
          /* Midpoint reached: carry profiles over to the 2D group and
           * start the click-free fade. */
          printf("Switching engines at t=%.1fs (fade %.0f ms)\n",
                 (double)absolute_frame / sample_rate,
                 (double)SPECBLEACH_TRANSITION_FADE_TIME_MS);

          params_2d.learn_noise = SPECBLEACH_LEARN_OFF;
          if (!specbleach_stereo_load_parameters_2d(nlm_group, &params_2d,
                                                    sizeof(params_2d)) ||
              !specbleach_stereo_migrate_profiles_from(nlm_group,
                                                       spectral_group) ||
              !specbleach_transition_begin(
                  transition, specbleach_stereo_get_latency(spectral_group),
                  specbleach_stereo_get_latency(nlm_group))) {
            success = false;
            break;
          }
          phase = PHASE_TRANSITION;
        }

        if (!specbleach_stereo_process(spectral_group, (uint32_t)read_frames,
                                       in_ptrs, out_ptrs)) {
          success = false;
          break;
        }
      } else {
        /* -----------------------------------------------------------------
         * PHASE_TRANSITION: the steady-state pattern while a fade runs.
         * Both groups render their own wet output from the SAME dry input;
         * the transition module aligns latencies and equal-power blends.
         * Note out_ptrs aliases channel_data (in-place) and blended == the
         * target's buffer, which specbleach_transition_process() allows. */
        if (!specbleach_stereo_process(spectral_group, (uint32_t)read_frames,
                                       in_ptrs, from_ptrs) ||
            !specbleach_stereo_process(nlm_group, (uint32_t)read_frames,
                                       in_ptrs, out_ptrs) ||
            !specbleach_transition_process(transition, (uint32_t)read_frames,
                                           (const float**)from_ptrs,
                                           (const float**)out_ptrs, out_ptrs)) {
          success = false;
          break;
        }
      }

      /* -----------------------------------------------------------------
       * LATENCY REPORTING: a real integrator's host-facing duty. While
       * transitioning, the reported latency is the target's (the module
       * aligns the audio internally); afterwards it is simply the active
       * group's. Replace this printf with setLatencySamples()/
       * pw_stream_latency/etc. in your application. */
      const uint32_t current_latency =
          phase == PHASE_TRANSITION
              ? specbleach_transition_get_latency(transition)
              : specbleach_stereo_get_latency(spectral_group);
      if (current_latency != reported_latency) {
        printf("Latency now %u samples\n", current_latency);
        reported_latency = current_latency;
      }

      /* Re-interleave processed audio */
      for (uint32_t ch = 0; ch < channels; ++ch) {
        for (sf_count_t s = 0; s < read_frames; ++s) {
          interleaved[s * channels + ch] = channel_data[ch * BLOCK_SIZE + s];
        }
      }

      sf_writef_float(output_file, interleaved, read_frames);
      absolute_frame += (uint64_t)read_frames;
    }
  } while (0);

  free(interleaved);
  free(channel_data);
  free(scratch_from);
  specbleach_transition_free(transition);
  specbleach_stereo_free(spectral_group);
  specbleach_stereo_free(nlm_group);
  sf_close(input_file);
  sf_close(output_file);

  if (!success) {
    remove(output_file_name);
    return 1;
  }
  printf("Wrote '%s'\n", output_file_name);
  return 0;
}
