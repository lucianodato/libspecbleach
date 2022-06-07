shared_sources = files(
    'shared/utils/general_utils.c',
    'shared/utils/denoise_mixer.c',
    'shared/utils/spectral_features.c',
    'shared/utils/spectral_utils.c',
    'shared/utils/spectral_trailing_buffer.c',
    'shared/pre_estimation/absolute_hearing_thresholds.c',
    'shared/pre_estimation/masking_estimator.c',
    'shared/pre_estimation/noise_scaling_criterias.c',
    'shared/pre_estimation/transient_detector.c',
    'shared/pre_estimation/spectral_smoother.c',
    'shared/post_estimation/spectral_whitening.c',
    'shared/post_estimation/postfilter.c',
    'shared/noise_estimation/adaptive_noise_estimator.c',
    'shared/noise_estimation/noise_estimator.c',
    'shared/noise_estimation/noise_profile.c',
    'shared/gain_estimation/gain_estimators.c',
    'shared/pre_estimation/critical_bands.c',
    'shared/stft/fft_transform.c',
    'shared/stft/stft_windows.c',
    'shared/stft/stft_buffer.c',
    'shared/stft/stft_processor.c',
)

processors_sources = files(
    'processors/denoiser/spectral_denoiser.c',
    'processors/adaptivedenoiser/adaptive_denoiser.c',
    'processors/specbleach_adenoiser.c',
    'processors/specbleach_denoiser.c',
)

specbleach_sources = [shared_sources, processors_sources]

# Build of the shared object
libspecbleach = library('specbleach',
  sources: specbleach_sources,
  c_args: lib_c_args,
  dependencies: dep,
  include_directories: inc,
  install: true)

# Generate package module configuration
pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries: libspecbleach,
                 version: meson.project_version(),
                 name: 'libspecbleach',
                 description: 'C library for audio noise reduction and other spectral effects')

# Declare libspecbleach dependency
libspecbleach_dep = declare_dependency(include_directories: inc,
  link_with: libspecbleach)