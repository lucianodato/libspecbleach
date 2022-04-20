project('libspecbleach', 'c', version: '0.1.0',default_options: ['c_std=c99'])

#sources to compile
shared_sources = [
    'src/shared/utils/spectral_utils.c',
    'src/shared/utils/spectral_features.c',
    'src/shared/utils/general_utils.c',
    'src/shared/utils/denoise_mixer.c',
    'src/shared/post_estimation/spectral_whitening.c',
    'src/shared/post_estimation/postfilter.c',
    'src/shared/pre_estimation/noise_scaling_criterias.c',
    'src/shared/pre_estimation/critical_bands.c',
    'src/shared/pre_estimation/absolute_hearing_thresholds.c',
    'src/shared/pre_estimation/masking_estimator.c',
    'src/shared/pre_estimation/transient_detector.c',
    'src/shared/pre_estimation/spectral_smoother.c',
    'src/shared/gain_estimation/gain_estimators.c',
    'src/shared/noise_estimation/adaptive_noise_estimator.c',
    'src/shared/noise_estimation/noise_estimator.c',
    'src/shared/noise_estimation/noise_profile.c',
    'src/shared/stft/fft_transform.c',
    'src/shared/stft/stft_windows.c',
    'src/shared/stft/stft_buffer.c',
    'src/shared/stft/stft_processor.c',
]

specbleach_sources = [
    'src/denoiser/spectral_denoiser.c',
    'src/adaptivedenoiser/adaptive_denoiser.c',
    'src/specbleach_adenoiser.c',
    'src/specbleach_denoiser.c',
]

example_sources = [
    'example/adenoiser_demo.c'
]

# public header
inc = include_directories('include')

m_dep = meson.get_compiler('c').find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
dep = [m_dep, fftw_dep]

libspecbleach = library('specbleach',
  shared_sources,
  specbleach_sources,
  include_directories : inc,
  dependencies: dep,
  install: false)

libspecbleach_dep = declare_dependency(include_directories : inc,
  link_with : libspecbleach)

sndfile_dep = meson.get_compiler('c').find_library('sndfile',
                has_headers: ['sndfile.h'],
                required: get_option('enable_example'))

if sndfile_dep.found() and get_option('enable_example')
    executable('adenoiser_demo', example_sources, dependencies : [sndfile_dep,libspecbleach_dep])
endif

