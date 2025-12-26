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

#ifndef SPECBLEACH_PARAMETERS_H_INCLUDED
#define SPECBLEACH_PARAMETERS_H_INCLUDED

#include <stdbool.h>

typedef struct SpectralBleachParameters {
  /* Sets the processor in listening mode to capture the noise profile.
   * 0 is disabled, 1 will learn all profile types simultaneously.
   * Note: This parameter is ignored by the adaptive denoiser. */
  int learn_noise;

  /* Sets the noise reduction mode to use when learning is disabled.
   * 1 will use the average profile, 2 will use the median profile
   * and 3 will use the max profile.
   * Note: This parameter is ignored by the adaptive denoiser. */
  int noise_reduction_mode;

  /* Enables outputting the residue of the reduction processing. It's either
   * true or false */
  bool residual_listen;

  /* Sets the amount of dBs that the noise will be attenuated. It goes from 0 dB
   * to 40 dB. This controls both the gain attenuation strength and the residual
   * noise mixing level for optimal noise reduction. */
  float reduction_amount;

  /* Percentage of smoothing to apply. Averages the reduction calculation frame
   * per frame so the rate of change is less resulting in less musical noise but
   * if too strong it can blur transient and reduce high frequencies. It goes
   * from 0 to 100 percent */
  float smoothing_factor;

  /* Percentage of whitening that is going to be applied to the residue of the
   * reduction. It modifies the noise floor to be more like white noise. This
   * can help hide musical noise when the noise is colored. It goes from 0 to
   * 100 percent */
  float whitening_factor;

  /* Type of algorithm used to scale noise in order to apply over or under
   * subtraction in different parts of the spectrum while calculating the
   * reduction. 0 is a-posteriori snr scaling using the complete spectrum, 1 is
   * a-posteriori using critical bands, 2 is using masking thresholds and 3 is
   * disabled.
   */
  int noise_scaling_type;

  /* Strength in which the reduction will be applied. It uses the masking
   * thresholds of the signal to determine where in the spectrum the reduction
   * needs to be stronger. This parameter scales how much in each of the
   * frequencies the reduction is going to be applied. It can be a positive dB
   * value in between 0 dB and 12 dB */
  float noise_rescale;

  /* Sets the SNR threshold in dB in which the post-filter will start to blur
   * musical noise. It can be a positive or negative dB value in between -10 dB
   * and 10 dB */
  float post_filter_threshold;

  /* Method used for adaptive noise estimation in the adaptive denoiser.
   * 0: LOUIZOU_METHOD uses minimum statistics (default)
   * 1: SPP_MMSE_METHOD uses Speech Presence Probability with MMSE estimation
   * for lower complexity and unbiased noise tracking. */
  int noise_estimation_method;
} SpectralBleachParameters;

#endif
