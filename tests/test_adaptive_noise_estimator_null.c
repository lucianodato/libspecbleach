#include <stdio.h>
#include <stdlib.h>

#include "processors/adaptivedenoiser/adaptive_denoiser.h"
#include "processors/denoiser/spectral_denoiser.h"
#include "shared/noise_estimation/adaptive_noise_estimator.h"
#include "shared/noise_estimation/noise_estimator.h"
#include "shared/noise_estimation/noise_profile.h"
#include "shared/post_estimation/noise_floor_manager.h"
#include "shared/post_estimation/postfilter.h"
#include "shared/post_estimation/spectral_whitening.h"
#include "shared/pre_estimation/absolute_hearing_thresholds.h"
#include "shared/pre_estimation/critical_bands.h"
#include "shared/pre_estimation/masking_estimator.h"
#include "shared/pre_estimation/noise_scaling_criterias.h"
#include "shared/pre_estimation/spectral_smoother.h"
#include "shared/pre_estimation/transient_detector.h"
#include "shared/stft/fft_transform.h"
#include "shared/stft/stft_buffer.h"
#include "shared/stft/stft_processor.h"
#include "shared/stft/stft_windows.h"
#include "shared/utils/denoise_mixer.h"
#include "shared/utils/spectral_features.h"
#include "shared/utils/spectral_trailing_buffer.h"
#include "specbleach_adenoiser.h"
#include "specbleach_denoiser.h"

int main(void) {
  printf("Testing NULL pointer safety for all free functions...\n");

  // Core API
  specbleach_free(NULL);
  specbleach_adaptive_free(NULL);

  // Shared Modules
  louizou_estimator_free(NULL);
  spp_mmse_estimator_free(NULL);
  noise_estimation_free(NULL);
  noise_profile_free(NULL);
  spectral_smoothing_free(NULL);
  transient_detector_free(NULL);
  critical_bands_free(NULL);
  masking_estimation_free(NULL);
  absolute_hearing_thresholds_free(NULL);
  noise_scaling_criterias_free(NULL);
  spectral_trailing_buffer_free(NULL);
  spectral_features_free(NULL);
  denoise_mixer_free(NULL);
  spectral_whitening_free(NULL);
  postfilter_free(NULL);
  noise_floor_manager_free(NULL);
  fft_transform_free(NULL);
  stft_processor_free(NULL);
  stft_buffer_free(NULL);
  stft_window_free(NULL);

  // Internal Processors
  spectral_denoiser_free(NULL);
  spectral_adaptive_denoiser_free(NULL);

  printf("✅ NULL pointer safety tests passed!\n");
  return 0;
}
