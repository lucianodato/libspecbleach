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

#include "specbleach_denoiser.h"
#include "denoiser/spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/stft/stft_processor.h"
#include "shared/utils/general_utils.h"
#include "shared/utils/tonal_detector.h"
#include "specbleach_processor_core.h"
#include <stdlib.h>
#include <string.h>

typedef struct specbleach_denoiser { // NOLINT(readability-identifier-naming)
  uint32_t hop;
  SbProcessorCore* core;
  SpectralProcessorHandle spectral_denoiser;
  float* reduction_curve_copy;
  uint32_t reduction_curve_capacity;
} SbSpectralDenoiser;

/**
 * Initializes a spectral denoiser for the specified audio format.
 *
 * @param sample_rate Audio sample rate in Hz.
 * @param frame_size STFT frame size.
 * @returns A configured denoiser instance, or `NULL` if initialization fails.
 */
specbleach_denoiser* specbleach_denoiser_initialize(const uint32_t sample_rate,
                                                    float frame_size) {
  SbSpectralDenoiser* self =
      (SbSpectralDenoiser*)calloc(1U, sizeof(SbSpectralDenoiser));
  if (!self) {
    return NULL;
  }

  self->core = sb_processor_core_initialize(
      sample_rate, frame_size, OVERLAP_FACTOR_1D, PADDING_CONFIGURATION_1D,
      ZEROPADDING_AMOUNT_1D, INPUT_WINDOW_TYPE_1D, OUTPUT_WINDOW_TYPE_1D, 0);

  if (!self->core) {
    specbleach_denoiser_free(self);
    return NULL;
  }

  const uint32_t fft_size = get_stft_fft_size(self->core->stft_processor);
  self->hop = get_stft_hop_size(self->core->stft_processor);

  self->spectral_denoiser = spectral_denoiser_initialize(
      self->core->sample_rate, fft_size, OVERLAP_FACTOR_1D,
      self->core->noise_profile);

  if (!self->spectral_denoiser) {
    specbleach_denoiser_free(self);
    return NULL;
  }

  return self;
}

/**
 * Releases a spectral denoiser instance and its owned resources.
 *
 * @param instance Denoiser instance to release; may be NULL.
 */
void specbleach_denoiser_free(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;

  if (!self) {
    return;
  }

  if (self->spectral_denoiser) {
    spectral_denoiser_free(self->spectral_denoiser);
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
 * Gets the total processing latency of the denoiser in samples.
 *
 * @param instance Denoiser instance.
 * @return The combined STFT and spectral-denoiser latency in samples, or 0 if the instance is invalid.
 */
uint32_t specbleach_denoiser_get_latency(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;

  if (!self || !self->core || !self->core->stft_processor) {
    return 0;
  }

  uint32_t stft_latency = get_stft_latency(self->core->stft_processor);
  uint32_t denoiser_latency_frames =
      spectral_denoiser_get_latency_frames(self->spectral_denoiser);
  uint32_t denoiser_latency_samples = denoiser_latency_frames * self->hop;

  return stft_latency + denoiser_latency_samples;
}

/**
 * Processes audio samples through the spectral denoiser.
 *
 * @param instance Denoiser instance to process the audio with.
 * @param number_of_samples Number of samples to process.
 * @param input Input audio samples.
 * @param output Buffer for the processed audio samples.
 * @return `true` if processing is started, `false` if the instance, sample count, or buffers are invalid.
 */
bool specbleach_denoiser_process(specbleach_denoiser* instance,
                                 const uint32_t number_of_samples,
                                 const float* input, float* output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  SbSpectralDenoiser* self = instance;

  stft_processor_run(self->core->stft_processor, number_of_samples, input,
                     output, &spectral_denoiser_run, self->spectral_denoiser);

  return true;
}

/**
 * Gets the number of bins in the configured noise profile.
 *
 * @param instance Denoiser instance.
 * @return Noise-profile size, or 0 if instance is NULL.
 */
uint32_t specbleach_denoiser_get_noise_profile_size(
    specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_size(self->core) : 0;
}

/**
 * Loads a restored noise profile using the specified profile mode.
 *
 * @param restored_profile Noise profile data to load.
 * @param profile_size Number of values in the profile.
 * @param block_count Number of blocks represented by the profile.
 * @param mode Noise profile loading mode.
 * @return `true` if the profile is loaded successfully, `false` otherwise.
 */
bool specbleach_denoiser_load_noise_profile_for_mode(
    specbleach_denoiser* instance, const float* restored_profile,
    const uint32_t profile_size, const uint32_t block_count, const int mode) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_load_noise_profile_for_mode(
                    self->core, restored_profile, profile_size, block_count,
                    mode)
              : false;
}

/**
 * Resets the noise profile maintained by the denoiser.
 *
 * @param instance Denoiser instance whose noise profile should be reset.
 * @return `true` if the profile was reset successfully, `false` if the instance is invalid or the reset failed.
 */
bool specbleach_denoiser_reset_noise_profile(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  if (!self) {
    return false;
  }
  if (self->spectral_denoiser) {
    spectral_denoiser_reset_noise_profile(self->spectral_denoiser);
  }
  return sb_processor_core_reset_noise_profile(self->core);
}

/**
 * Gets the number of noise-profile blocks available for a processing mode.
 * @param instance Denoiser instance.
 * @param mode Processing mode whose noise-profile block count is requested.
 * @returns Number of noise-profile blocks, or 0 if the instance is invalid.
 */
uint32_t specbleach_denoiser_get_noise_profile_block_count_for_mode(
    specbleach_denoiser* instance, int mode) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_block_count_for_mode(
                    self->core, mode)
              : 0;
}

/**
 * Retrieves the noise profile for the specified mode.
 * @param instance Denoiser instance.
 * @param mode Noise-profile mode to retrieve.
 * @return Pointer to the noise profile, or `NULL` if the instance is invalid.
 */
float* specbleach_denoiser_get_noise_profile_for_mode(
    specbleach_denoiser* instance, int mode) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_get_noise_profile_for_mode(self->core, mode)
              : NULL;
}

/**
 * Determines whether a noise profile is available for the specified mode.
 *
 * @param instance Denoiser instance to query.
 * @param mode Noise-profile mode to check.
 * @return `true` if a noise profile is available for the mode, `false` otherwise.
 */
bool specbleach_denoiser_noise_profile_available_for_mode(
    specbleach_denoiser* instance, int mode) {
  SbSpectralDenoiser* self = instance;
  return self ? sb_processor_core_noise_profile_available_for_mode(self->core,
                                                                   mode)
              : false;
}

/**
 * Loads denoiser parameters and applies them to the instance.
 *
 * @param instance Denoiser instance to configure.
 * @param parameters Parameters to load.
 * @param parameters_size Size of the parameters structure.
 * @return `true` if the parameters were loaded successfully, `false` otherwise.
 */
bool specbleach_denoiser_load_parameters(
    specbleach_denoiser* instance,
    const SpecbleachDenoiserParameters* parameters,
    const uint32_t parameters_size) {
  if (!instance || !parameters ||
      parameters_size != sizeof(SpecbleachDenoiserParameters)) {
    return false;
  }

  SbSpectralDenoiser* self = instance;

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

  DenoiserParameters denoise_parameters =
      sb_denoiser_params_sanitize(parameters);
  denoise_parameters.reduction_curve_bias = owned_bias;

  return load_reduction_parameters(self->spectral_denoiser, denoise_parameters);
}

/**
 * Retrieves the current tonal mask.
 *
 * @param instance Denoiser instance.
 * @return Pointer to the tonal mask, or NULL if the instance is invalid.
 */
const float* specbleach_denoiser_get_tonal_mask(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_get_tonal_mask(self->spectral_denoiser)
              : NULL;
}

/**
 * Retrieves the frequencies of currently detected tonal peaks.
 * @param peak_freqs_hz Buffer that receives peak frequencies in hertz.
 * @param max_peaks Maximum number of peak frequencies to write.
 * @returns The number of peak frequencies written.
 */
uint32_t specbleach_denoiser_get_tonal_peaks(specbleach_denoiser* instance,
                                             float* peak_freqs_hz,
                                             uint32_t max_peaks) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_get_peaks(self->spectral_denoiser,
                                            peak_freqs_hz, max_peaks)
              : 0;
}

/**
 * Detects tonal peak frequencies from a noise profile.
 * @param profile Noise profile to analyze.
 * @param profile_size Number of values in the noise profile.
 * @param peak_freqs_hz Output buffer for detected peak frequencies in hertz.
 * @param max_peaks Maximum number of frequencies to write.
 * @returns Number of detected peak frequencies, or zero for invalid input.
 */
uint32_t specbleach_denoiser_get_tonal_peaks_for_profile(
    specbleach_denoiser* instance, const float* profile, uint32_t profile_size,
    float* peak_freqs_hz, uint32_t max_peaks) {
  SbSpectralDenoiser* self = instance;
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
 * Gets the currently active noise profile.
 *
 * @param instance Denoiser instance.
 * @return Pointer to the active noise profile, or NULL if the instance is invalid.
 */
const float* specbleach_denoiser_get_active_noise_profile(
    specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_get_active_noise_profile(
                    self->spectral_denoiser)
              : NULL;
}

/**
 * Determines whether a transient has been detected.
 *
 * @param instance Denoiser instance to query.
 * @return `true` if a transient is detected, `false` otherwise.
 */
bool specbleach_denoiser_is_transient_detected(specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_is_transient_detected(self->spectral_denoiser)
              : false;
}

/**
 * Gets the intensity of the currently detected transient.
 *
 * @param instance Denoiser instance.
 * @returns The transient intensity, or 0.0f if instance is NULL.
 */
float specbleach_denoiser_get_transient_intensity(
    specbleach_denoiser* instance) {
  SbSpectralDenoiser* self = instance;
  return self ? spectral_denoiser_get_transient_intensity(
                    self->spectral_denoiser)
              : 0.0f;
}
