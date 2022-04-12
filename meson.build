project('libspecbleach', 'c', version: '0.1.0',default_options: ['c_std=c99','default_library=static','buildtype=release'])

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
    'src/specbleach.c',
]

# public header
inc = include_directories('include')

m_dep = meson.get_compiler('c').find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
dep = [m_dep, fftw_dep]

#get the host operating system and configure install path and shared object extension
current_os = host_machine.system()
current_arch = host_machine.cpu_family()

#compiler flags
cflags = []

# Add x86_64 optimization where appropriate (not for ARM, not for riscv)
if current_arch != 'aarch64' and current_arch != 'riscv64' and current_os != 'darwin' 
    cflags += ['-ffast-math', '-msse','-msse2','-mfpmath=sse']
endif
# Add osx multiarch flags when appropriate
if current_os == 'darwin'
    cflags += ['-mrecip']
endif
if current_os == 'windows'
    cflags += ['-fno-stack-protector']
endif

libspecbleach = static_library('specbleach',
  shared_sources,
  specbleach_sources,
  include_directories : inc,
  c_args: cflags,
  dependencies: dep,
  install: false)

libspecbleach_dep = declare_dependency(include_directories : inc,
  link_with : libspecbleach)

