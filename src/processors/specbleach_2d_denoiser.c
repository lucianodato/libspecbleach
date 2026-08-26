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

#include "specbleach_2d_denoiser.h"
#include "denoiser2d/spectral_2d_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/stft/stft_processor.h"
#include "shared/utils/general_utils.h"
#include "shared/utils/tonal_detector.h"
#include "specbleach_processor_core.h"
#include <stdlib.h>
#include <string.h>

typedef struct specbleach_2d_denoiser { // NOLINT(readability-identifier-naming)
  uint32_t hop;
  SbProcessorCore* core;
  SpectralProcessorHandle spectral_2d_denoiser;
  float* reduction_curve_copy;
  uint32_t reduction_curve_capacity;
} Sb2DDenoiser;

/**
 * Initializes a 2D spectral denoiser for the specified sample rate and frame size.
 *
 * @param sample_rate Audio sample rate in Hz.
 * @param frame_size Processing frame size.
 * @return A newly initialized denoiser, or `NULL` if the arguments are invalid or initialization fails.
 */
specbleach_2d_denoiser* specbleach_2d_initialize(const uint32_t sample_rate,
                                                 float frame_size) {
  if (sample_rate == 0 || frame_size <= 0.0f) {
    return NULL;
  }

  Sb2DDenoiser* self = (Sb2DDenoiser*)calloc(1U, sizeof(Sb2DDenoiser));
  if (!self) {
    return NULL;
  }

  self->core = sb_processor_core_initialize(
      sample_rate, frame_size, OVERLAP_FACTOR_2D, PADDING_CONFIGURATION_2D,
      ZEROPADDING_AMOUNT_2D, INPUT_WINDOW_TYPE_2D, OUTPUT_WINDOW_TYPE_2D,
      SB_PROCESSOR_CORE_DEFAULT_REAL_SPECTRUM);

  if (!self->core) {
    specbleach_2d_free(self);
    return NULL;
  }

  const uint32_t fft_size = get_stft_fft_size(self->core->stft_processor);
  self->hop = get_stft_hop_size(self->core->stft_processor);

  self->spectral_2d_denoiser = spectral_2d_denoiser_initialize(
      sample_rate, fft_size, OVERLAP_FACTOR_2D, self->core->noise_profile);

  if (!self->spectral_2d_denoiser) {
    specbleach_2d_free(self);
    return NULL;
  }

  return self;
}

/**
 * Releases a 2D denoiser instance and its associated resources.
 * @param instance Denoiser instance to release; may be NULL.
 */
void specbleach_2d_free(specbleach_2d_denoiser* instance) {
  Sb2DDenoiser* self = instance;

  if (!self) {
    return;
  }

  if (self->spectral_2d_denoiser) {
    spectral_2d_denoiser_free(self->spectral_2d_denoiser);
  }
  if (self->core) {
    sb_processor_core_free(self->core);
  }
  if (self->reduction_curve_copy) {
    free(self->reduction_curve_copy);
  }

  free(self);
}

/**
 * Gets the processing latency in samples.
 *
 * @param instance Denoiser instance.
 * @return Total STFT and spectral-denoiser look-ahead latency in samples, or
 *         0 if the instance is invalid.
 */
uint32_t specbleach_2d_get_latency(specbleach_2d_denoiser* instance) {
  Sb2DDenoiser* self = instance;

  if (!self || !self->core || !self->core->stft_processor) {
    return 0;
  }

  // Base STFT latency
  uint32_t stft_latency = get_stft_latency(self->core->stft_processor);

  // Additional NLM look-ahead latency
  uint32_t nlm_latency_frames =
      spectral_2d_denoiser_get_latency_frames(self->spectral_2d_denoiser);
  uint32_t nlm_latency_samples = nlm_latency_frames * self->hop;

  return stft_latency + nlm_latency_samples;
}

/**
 * Processes audio samples through the 2D spectral denoiser.
 *
 * @param instance Denoiser instance to process the audio with.
 * @param number_of_samples Number of samples in the input and output buffers.
 * @param input Input audio samples.
 * @param output Buffer for the processed audio samples.
 * @returns `true` if processing was started, `false` if the instance, sample count, or buffers are invalid.
 */
bool specbleach_2d_process(specbleach_2d_denoiser* instance,
                           const uint32_t number_of_samples, const float* input,
                           float* output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  Sb2DDenoiser* self = instance;

  stft_processor_run(self->core->stft_processor, number_of_samples, input,
                     output, &spectral_2d_denoiser_run,
                     self->spectral_2d_denoiser);

  return true;
}

/**
 * Gets the number of samples in the noise profile.
 * @param instance Denoiser instance.
 * @return Noise profile size, or 0 if the instance is NULL.
 */
uint32_t specbleach_2d_get_noise_profile_size(
    specbleach_2d_denoiser* instance) {
  Sb2DDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_size(self->core) : 0;
}

/**
 * Loads a restored noise profile for the specified processing mode.
 *
 * @param restored_profile Noise-profile data to load.
 * @param profile_size Number of values in the noise profile.
 * @param block_count Number of blocks represented by the profile.
 * @param mode Processing mode associated with the profile.
 * @return `true` if the profile is loaded successfully, `false` otherwise.
 */
bool specbleach_2d_load_noise_profile_for_mode(specbleach_2d_denoiser* instance,
                                               const float* restored_profile,
                                               const uint32_t profile_size,
                                               const uint32_t block_count,
                                               const int mode) {
  Sb2DDenoiser* self = instance;
  return self ? sb_processor_core_load_noise_profile_for_mode(
                    self->core, restored_profile, profile_size, block_count,
                    mode)
              : false;
}

/**
 * Resets the noise profile used by the 2D denoiser.
 *
 * @param instance Denoiser instance whose noise profile is reset.
 * @return `true` if the profile was reset successfully, `false` if the instance is `NULL` or the reset failed.
 */
bool specbleach_2d_reset_noise_profile(specbleach_2d_denoiser* instance) {
  Sb2DDenoiser* self = instance;
  if (!self) {
    return false;
  }
  if (self->spectral_2d_denoiser) {
    spectral_2d_denoiser_reset_noise_profile(self->spectral_2d_denoiser);
  }
  return sb_processor_core_reset_noise_profile(self->core);
}

/**
 * Gets the number of noise-profile blocks available for a profile mode.
 * @param mode Noise-profile mode to query.
 * @returns The number of noise-profile blocks, or 0 if the instance is NULL.
 */
uint32_t specbleach_2d_get_noise_profile_block_count_for_mode(
    specbleach_2d_denoiser* instance, int mode) {
  Sb2DDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_block_count_for_mode(
                    self->core, mode)
              : 0;
}

/**
 * Retrieves the noise profile associated with a profile mode.
 *
 * @param mode Noise-profile mode to retrieve.
 * @return Pointer to the noise profile, or `NULL` if the instance or profile is unavailable.
 */
float* specbleach_2d_get_noise_profile_for_mode(
    specbleach_2d_denoiser* instance, int mode) {
  Sb2DDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_for_mode(self->core, mode)
              : NULL;
}

/**
 * Determines whether a noise profile is available for the specified mode.
 *
 * @param mode Noise-profile mode to query.
 * @return `true` if a noise profile is available for the mode, `false` otherwise.
 */
bool specbleach_2d_noise_profile_available_for_mode(
    specbleach_2d_denoiser* instance, int mode) {
  Sb2DDenoiser* self = instance;
  return self ? sb_processor_core_noise_profile_available_for_mode(self->core,
                                                                   mode)
              : false;
}

/**
 * Loads and applies configuration parameters for the 2D denoiser.
 *
 * @param parameters Parameters to apply, including optional reduction-curve settings.
 * @param parameters_size Size of the parameter structure in bytes.
 * @returns `true` if the parameters are accepted and applied, `false` otherwise.
 */
bool specbleach_2d_load_parameters(
    specbleach_2d_denoiser* instance,
    const Specbleach2DDenoiserParameters* parameters,
    const uint32_t parameters_size) {
  if (!instance || !parameters ||
      parameters_size != sizeof(Specbleach2DDenoiserParameters)) {
    return false;
  }

  Sb2DDenoiser* self = instance;

  const uint32_t profile_size =
      sb_processor_core_get_noise_profile_size(self->core);
  if (parameters->reduction_curve_enabled &&
      parameters->reduction_curve_size != profile_size) {
    return false;
  }

  const float* owned_bias = sb_curve_bias_copy(
      &self->reduction_curve_copy, &self->reduction_curve_capacity,
      profile_size, parameters->reduction_curve_enabled,
      parameters->reduction_curve_bias);
  if (parameters->reduction_curve_enabled && !owned_bias) {
    return false;
  }

  Denoiser2DParameters denoise_parameters =
      sb_denoiser_2d_params_sanitize(parameters);
  denoise_parameters.reduction_curve_bias = owned_bias;

  return load_2d_reduction_parameters(self->spectral_2d_denoiser,
                                      denoise_parameters);
}

/**
 * Gets the current tonal mask.
 *
 * @param instance 2D denoiser instance.
 * @return Pointer to the tonal mask, or `NULL` if the instance is invalid.
 */
const float* specbleach_2d_get_tonal_mask(specbleach_2d_denoiser* instance) {
  Sb2DDenoiser* self = instance;
  return self ? spectral_2d_denoiser_get_tonal_mask(self->spectral_2d_denoiser)
              : NULL;
}

/**
 * Retrieves the detected tonal peak frequencies.
 *
 * @param peak_freqs_hz Buffer receiving peak frequencies in hertz.
 * @param max_peaks Maximum number of peak frequencies to write.
 * @return Number of tonal peaks written, or 0 if the instance is NULL.
 */
uint32_t specbleach_2d_get_tonal_peaks(specbleach_2d_denoiser* instance,
                                       float* peak_freqs_hz,
                                       uint32_t max_peaks) {
  Sb2DDenoiser* self = instance;
  return self ? spectral_2d_denoiser_get_peaks(self->spectral_2d_denoiser,
                                               peak_freqs_hz, max_peaks)
              : 0;
}

/**
 * Detects tonal peaks using the median noise profile when available.
 * @param profile Noise profile to use when no median profile is available.
 * @param profile_size Number of frequency bins in the profile.
 * @param peak_freqs_hz Output array for detected peak frequencies in hertz.
 * @param max_peaks Maximum number of peaks to write.
 * @return Number of detected tonal peaks, or 0 for invalid input.
 */
uint32_t specbleach_2d_get_tonal_peaks_for_profile(
    specbleach_2d_denoiser* instance, const float* profile,
    uint32_t profile_size, float* peak_freqs_hz, uint32_t max_peaks) {
  Sb2DDenoiser* self = instance;
  if (!self || !self->core || !profile || profile_size == 0 ||
      profile_size != sb_processor_core_get_noise_profile_size(self->core)) {
    return 0;
  }
  uint32_t fft_size = get_stft_fft_size(self->core->stft_processor);

  float* median_profile =
      sb_processor_core_get_noise_profile_for_mode(self->core, MEDIAN);
  const float* eval_profile = median_profile ? median_profile : profile;

  return tonal_detector_get_peaks_from_profile(
      eval_profile, profile_size, self->core->sample_rate, fft_size,
      peak_freqs_hz, max_peaks);
}

/**
 * Gets the active noise profile.
 * @param instance Denoiser instance.
 * @return Pointer to the active noise profile, or NULL if instance is NULL.
 */
const float* specbleach_2d_get_active_noise_profile(
    specbleach_2d_denoiser* instance) {
  Sb2DDenoiser* self = instance;
  return self ? spectral_2d_denoiser_get_active_noise_profile(
                    self->spectral_2d_denoiser)
              : NULL;
}

/**
 * Determines whether a transient is currently detected.
 *
 * @param instance Denoiser instance to query.
 * @return `true` if a transient is detected, `false` otherwise.
 */
bool specbleach_2d_is_transient_detected(specbleach_2d_denoiser* instance) {
  Sb2DDenoiser* self = instance;
  return self ? spectral_2d_denoiser_is_transient_detected(
                    self->spectral_2d_denoiser)
              : false;
}

/**
 * Gets the current transient intensity.
 *
 * @param instance 2D denoiser instance.
 * @return The transient intensity, or `0.0f` if the instance is `NULL`.
 */
float specbleach_2d_get_transient_intensity(specbleach_2d_denoiser* instance) {
  Sb2DDenoiser* self = instance;
  return self ? spectral_2d_denoiser_get_transient_intensity(
                    self->spectral_2d_denoiser)
              : 0.0f;
}
