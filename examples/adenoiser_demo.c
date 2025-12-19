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

#include <sndfile.h>
#include <specbleach_adenoiser.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// This is not a deliberate value. The library handles any amount passed
// through a circular buffer
#define BLOCK_SIZE 512
#define FRAME_SIZE 20

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
  if (argc != 3) {
    fprintf(stderr, "usage: %s <noisy input> <denoised output>\n", argv[0]);
    return 1;
  }

  const char* input_file_name = argv[1];
  const char* output_file_name = argv[2];

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

    // Configuration of the denoising parameters
    SpectralBleachParameters parameters =
        (SpectralBleachParameters){.residual_listen = false,
                                   .reduction_amount = 20.F,
                                   .smoothing_factor = 0.F,
                                   .whitening_factor = 100.F,
                                   .noise_scaling_type = 2,
                                   .noise_rescale = 2.F,
                                   .post_filter_threshold = -10.F};

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
