#include <stdio.h>
#include <stdlib.h>

#include "processors/denoiser/spectral_denoiser.h"

#include "shared/denoiser_logic/core/noise_floor_manager.h"
#include "shared/denoiser_logic/core/noise_profile.h"
#include "shared/denoiser_logic/estimators/adaptive_noise_estimator.h"
#include "shared/denoiser_logic/estimators/noise_estimator.h"
#include "shared/denoiser_logic/processing/spectral_whitening.h"
#include "shared/denoiser_logic/processing/suppression_engine.h"
#include "shared/stft/fft_transform.h"
#include "shared/stft/stft_buffer.h"
#include "shared/stft/stft_processor.h"
#include "shared/stft/stft_windows.h"
#include "shared/utils/absolute_hearing_thresholds.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/masking_estimator.h"
#include "shared/utils/spectral_features.h"
#include "shared/utils/spectral_smoother.h"
#include "shared/utils/transient_detector.h"

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
  MaskingEstimator* me = masking_estimation_initialize(
      fft_size, sample_rate, OPUS_SCALE, POWER_SPECTRUM);
  masking_estimation_free(me);

  // Absolute Hearing Thresholds
  absolute_hearing_thresholds_free(NULL);
  AbsoluteHearingThresholds* aht = absolute_hearing_thresholds_initialize(
      sample_rate, fft_size, POWER_SPECTRUM);
  absolute_hearing_thresholds_free(aht);

  // Adaptive Noise Estimators
  adaptive_estimator_free(NULL);
  AdaptiveNoiseEstimator* mar = adaptive_estimator_initialize(
      real_spectrum_size, sample_rate, fft_size, MARTIN_METHOD);
  adaptive_estimator_free(mar);

  AdaptiveNoiseEstimator* spp = adaptive_estimator_initialize(
      real_spectrum_size, sample_rate, fft_size, SPP_MMSE_METHOD);
  adaptive_estimator_free(spp);

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
  StftWindows* sw_win = stft_window_initialize(
      fft_size, fft_size, fft_size / 4U, HANN_WINDOW, HANN_WINDOW);
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

  // Spectral Features
  spectral_features_free(NULL);
  SpectralFeatures* sfe = spectral_features_initialize(real_spectrum_size);
  spectral_features_free(sfe);

  // Suppression Engine
  suppression_engine_free(NULL);
  SuppressionEngine* se = suppression_engine_initialize(
      real_spectrum_size, sample_rate, OPUS_SCALE, POWER_SPECTRUM);
  suppression_engine_free(se);

  // Internal Processors
  spectral_denoiser_free(NULL);

  // Dispatcher Coverage (Extra)
  printf("Testing Adaptive Estimator Dispatcher...\n");
  float* noise_profile = (float*)calloc(real_spectrum_size, sizeof(float));

  AdaptiveNoiseEstimator* mar_est = adaptive_estimator_initialize(
      real_spectrum_size, sample_rate, fft_size, MARTIN_METHOD);
  AdaptiveNoiseEstimator* spp_est = adaptive_estimator_initialize(
      real_spectrum_size, sample_rate, fft_size, SPP_MMSE_METHOD);

  adaptive_estimator_set_state(mar_est, noise_profile, MARTIN_METHOD);
  adaptive_estimator_set_state(spp_est, noise_profile, SPP_MMSE_METHOD);

  adaptive_estimator_apply_floor(mar_est, noise_profile);
  adaptive_estimator_update_seed(mar_est, noise_profile);

  adaptive_estimator_apply_floor(spp_est, noise_profile);
  adaptive_estimator_update_seed(spp_est, noise_profile);

  adaptive_estimator_free(mar_est);
  adaptive_estimator_free(spp_est);
  free(noise_profile);

  printf("✅ All safety and basic coverage tests passed!\n");
  return 0;
}
