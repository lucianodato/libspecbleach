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
 * This is just a simple console application example using the library with the
 * adaptive denoiser processor. It uses libsndfile to read audio, process it
 * with the algorithm and write it to an output file
 */

#include <getopt.h>
#include <sndfile.h>
#include <specbleach_adenoiser.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// This is not a deliberate value. The library handles any amount passed
// through a circular buffer
#define BLOCK_SIZE 512
#define FRAME_SIZE 20

static void print_usage(const char* prog_name) {
  fprintf(stderr, "Usage: %s [options] <noisy input> <denoised output>\n",
          prog_name);
  fprintf(stderr, "Options:\n");
  fprintf(stderr,
          "  --reduction <val>      Reduction amount in dB (default: 20.0)\n");
  fprintf(stderr,
          "  --whitening <val>      Whitening factor (default: 50.0)\n");
  fprintf(stderr, "  --smoothing <val>      Smoothing factor (default: 0.0)\n");
  fprintf(stderr,
          "  --rescale <val>        Noise rescale in dB (default: 6.0)\n");
  fprintf(stderr,
          "  --scaling-type <val>   Noise scaling type (0-2, default: 2)\n");
  fprintf(stderr,
          "  --threshold <val>      Post-filter threshold in dB (default: "
          "-10.0)\n");
  fprintf(stderr,
          "  --method <val>         Noise estimation method (0=Louizou, "
          "1=SPP-MMSE, default: 0)\n");
  fprintf(stderr, "  --help                Show this help message\n");
}

static void cleanup_resources(SF_INFO* sfinfo, SNDFILE* input_file,
                              SNDFILE* output_file, float* input_buffer,
                              float* output_buffer,
                              SpectralBleachHandle lib_instance) {
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
    specbleach_adaptive_free(lib_instance);
  }
}

int main(int argc, char** argv) {
  SpectralBleachParameters parameters =
      (SpectralBleachParameters){.residual_listen = false,
                                 .reduction_amount = 20.F,
                                 .smoothing_factor = 0.F,
                                 .whitening_factor = 50.F,
                                 .noise_scaling_type = 2,
                                 .noise_rescale = 6.F,
                                 .post_filter_threshold = -10.F,
                                 .noise_estimation_method = LOUIZOU_METHOD};

  static struct option long_options[] = {
      {"reduction", required_argument, 0, 'r'},
      {"whitening", required_argument, 0, 'w'},
      {"smoothing", required_argument, 0, 's'},
      {"rescale", required_argument, 0, 'e'},
      {"scaling-type", required_argument, 0, 't'},
      {"threshold", required_argument, 0, 'h'},
      {"method", required_argument, 0, 'm'},
      {"help", no_argument, 0, '?'},
      {0, 0, 0, 0}};

  int opt;
  while ((opt = getopt_long(argc, argv, "r:w:s:e:t:h:m:", long_options,
                            NULL)) != -1) {
    switch (opt) {
      case 'r':
        parameters.reduction_amount = (float)atof(optarg);
        break;
      case 'w':
        parameters.whitening_factor = (float)atof(optarg);
        break;
      case 's':
        parameters.smoothing_factor = (float)atof(optarg);
        break;
      case 'e':
        parameters.noise_rescale = (float)atof(optarg);
        break;
      case 't':
        parameters.noise_scaling_type = atoi(optarg);
        break;
      case 'h':
        parameters.post_filter_threshold = (float)atof(optarg);
        break;
      case 'm':
        parameters.noise_estimation_method = atoi(optarg);
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
  SpectralBleachHandle lib_instance = NULL;
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

    // Initialize library instance
    lib_instance = specbleach_adaptive_initialize((uint32_t)sfinfo->samplerate,
                                                  FRAME_SIZE);
    if (!lib_instance) {
      fprintf(stderr, "Error: Failed to initialize library instance\n");
      break;
    }

    // Load the parameters before doing the denoising
    if (!specbleach_adaptive_load_parameters(lib_instance, parameters)) {
      fprintf(stderr, "Error: Failed to load parameters\n");
      break;
    }

    // Process audio
    sf_count_t frames_read;
    while ((frames_read = sf_readf_float(input_file, input_library_buffer,
                                         BLOCK_SIZE)) > 0) {
      // Process the audio
      if (!specbleach_adaptive_process(lib_instance, (uint32_t)BLOCK_SIZE,
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

    // Success
    ret = 0;
  } while (0);

  cleanup_resources(sfinfo, input_file, output_file, input_library_buffer,
                    output_library_buffer, lib_instance);
  return ret;
}
