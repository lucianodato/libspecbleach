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

#ifndef DSAF_ENGINE_H
#define DSAF_ENGINE_H

#include <stdbool.h>
#include <stdint.h>

#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"

typedef struct SbDsafEngine SbDsafEngine;

typedef struct DsafParameters {
  float reduction_gain;   // Linear gain floor G_min (0..1)
  float smoothing_factor; // Band-gain temporal smoothing (0..1)
  bool adaptive_noise;    // Enable continuous noise tracking
  AdaptiveNoiseEstimationMethod
      noise_estimation_method; // Noise floor tracking method
  float suppression_strength;  // Oversubtraction factor (0..1 mapped)
} DsafParameters;

/**
 * Zero-latency minimum-phase delayless subband adaptive filter engine.
 *
 * The real-time path is a sample-by-sample M-tap FIR convolution. An
 * internal worker thread asynchronously re-synthesizes the FIR
 * coefficients from spectral noise estimates and publishes them through a
 * lock-free double buffer. Algorithmic latency is always zero.
 *
 * Lifecycle: initialize/free are setup-only (allocate, spawn/join the
 * worker). load_parameters and process follow the standard contract
 * (never concurrent with each other on the same instance; process is
 * RT-safe).
 */
SbDsafEngine* dsaf_engine_initialize(uint32_t sample_rate);
void dsaf_engine_free(SbDsafEngine* self);

void dsaf_engine_load_parameters(SbDsafEngine* self,
                                 const DsafParameters* parameters);

bool dsaf_engine_process(SbDsafEngine* self, uint32_t number_of_samples,
                         const float* input, float* output);

uint32_t dsaf_engine_get_latency(SbDsafEngine* self);

#endif
