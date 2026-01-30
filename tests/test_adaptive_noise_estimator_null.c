#include <stdio.h>
#include <stdlib.h>

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

#include "specbleach_denoiser.h"

int main(void) {
  printf("Starting comprehensive module safety and coverage tests...\n");

  uint32_t sample_rate = 44100;
  uint32_t fft_size = 1024;
  uint32_t hop = 256;
  uint32_t overlap_factor = 4;
  uint32_t real_spectrum_size = (fft_size / 2) + 1;

  // 1. Core API
  printf("Testing Core API...\n");
  specbleach_free(NULL);
  SpectralBleachHandle sb = specbleach_initialize(sample_rate, 10.0f);
  if (sb) {
    specbleach_free(sb);
  }

  // 2. Internal Modules
  printf("Testing Internal Modules...\n");

  // Spectral Smoother
  spectral_smoothing_free(NULL);
  SpectralSmoother* ss = spectral_smoothing_initialize(fft_size, FIXED);
  spectral_smoothing_free(ss);

  // Transient Detector
  transient_detector_free(NULL);
  TransientDetector* td = transient_detector_initialize(fft_size);
  transient_detector_free(td);

  // Critical Bands
  critical_bands_free(NULL);
  CriticalBands* cb =
      critical_bands_initialize(sample_rate, fft_size, OPUS_SCALE);
  critical_bands_free(cb);

  // Masking Estimator
  masking_estimation_free(NULL);
  MaskingEstimator* me =
      masking_estimation_initialize(fft_size, sample_rate, POWER_SPECTRUM);
  masking_estimation_free(me);

  // Absolute Hearing Thresholds
  absolute_hearing_thresholds_free(NULL);
  AbsoluteHearingThresholds* aht = absolute_hearing_thresholds_initialize(
      sample_rate, fft_size, POWER_SPECTRUM);
  absolute_hearing_thresholds_free(aht);

  // Adaptive Noise Estimators
  louizou_estimator_free(NULL);
  AdaptiveNoiseEstimator* le =
      louizou_estimator_initialize(real_spectrum_size, sample_rate, fft_size);
  louizou_estimator_free(le);

  spp_mmse_estimator_free(NULL);
  AdaptiveNoiseEstimator* spp =
      spp_mmse_estimator_initialize(real_spectrum_size, sample_rate, fft_size);
  spp_mmse_estimator_free(spp);

  // Noise Profile
  noise_profile_free(NULL);
  NoiseProfile* np = noise_profile_initialize(real_spectrum_size);

  // Noise Estimator (requires noise profile)
  noise_estimation_free(NULL);
  NoiseEstimator* ne = noise_estimation_initialize(fft_size, np);
  noise_estimation_free(ne);

  noise_profile_free(np);

  // Noise Floor Manager
  noise_floor_manager_free(NULL);
  NoiseFloorManager* nfm =
      noise_floor_manager_initialize(fft_size, sample_rate, hop);
  noise_floor_manager_free(nfm);

  // Postfilter
  postfilter_free(NULL);
  PostFilter* pf = postfilter_initialize(fft_size);
  postfilter_free(pf);

  // Spectral Whitening
  spectral_whitening_free(NULL);
  SpectralWhitening* sw = spectral_whitening_initialize(fft_size);
  spectral_whitening_free(sw);

  // STFT Buffer
  stft_buffer_free(NULL);
  StftBuffer* sb_buf = stft_buffer_initialize(fft_size, 0, hop);
  stft_buffer_free(sb_buf);

  // STFT Windows
  stft_window_free(NULL);
  StftWindows* sw_win = stft_window_initialize(fft_size, overlap_factor,
                                               HANN_WINDOW, HANN_WINDOW);
  stft_window_free(sw_win);

  // FFT Transform (internal)
  fft_transform_free(NULL);
  FftTransform* ft = fft_transform_initialize(fft_size, NO_PADDING, 0);
  fft_transform_free(ft);

  // STFT Processor
  stft_processor_free(NULL);
  StftProcessor* stft_p =
      stft_processor_initialize(sample_rate, 10.0f, overlap_factor, NO_PADDING,
                                0, HANN_WINDOW, HANN_WINDOW);
  stft_processor_free(stft_p);

  // Denoise Mixer
  denoise_mixer_free(NULL);
  DenoiseMixer* dm = denoise_mixer_initialize(fft_size, sample_rate, hop);
  denoise_mixer_free(dm);

  // Spectral Trailing Buffer
  spectral_trailing_buffer_free(NULL);
  SpectralTrailingBuffer* stb =
      spectral_trailing_buffer_initialize(real_spectrum_size, 5);
  spectral_trailing_buffer_free(stb);

  // Spectral Features
  spectral_features_free(NULL);
  SpectralFeatures* sfe = spectral_features_initialize(real_spectrum_size);
  spectral_features_free(sfe);

  // Noise Scaling Criteria
  noise_scaling_criterias_free(NULL);
  NoiseScalingCriterias* nsc = noise_scaling_criterias_initialize(
      fft_size, OPUS_SCALE, sample_rate, POWER_SPECTRUM);
  noise_scaling_criterias_free(nsc);

  // Internal Processors
  spectral_denoiser_free(NULL);

  // Dispatcher Coverage (Extra)
  printf("Testing Adaptive Estimator Dispatcher...\n");
  float* noise_profile = (float*)calloc(real_spectrum_size, sizeof(float));

  AdaptiveNoiseEstimator* lou_est =
      louizou_estimator_initialize(real_spectrum_size, sample_rate, fft_size);
  AdaptiveNoiseEstimator* spp_est =
      spp_mmse_estimator_initialize(real_spectrum_size, sample_rate, fft_size);

  adaptive_estimator_set_state(lou_est, noise_profile, LOUIZOU_METHOD);
  adaptive_estimator_set_state(spp_est, noise_profile, SPP_MMSE_METHOD);

  adaptive_estimator_apply_floor(lou_est, noise_profile);
  adaptive_estimator_update_seed(lou_est, noise_profile);

  adaptive_estimator_apply_floor(spp_est, noise_profile);
  adaptive_estimator_update_seed(spp_est, noise_profile);

  louizou_estimator_free(lou_est);
  spp_mmse_estimator_free(spp_est);
  free(noise_profile);

  printf("✅ All safety and basic coverage tests passed!\n");
  return 0;
}
