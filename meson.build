project('libspecbleach', 'c', version: '0.1.6',default_options: ['c_std=c99'])

# Sources to compile
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

specbleach_headers = [
  'include/specbleach_adenoiser.h',
  'include/specbleach_denoiser.h',
]

example_sources = [
    'example/adenoiser_demo.c'
]

# Dependencies for libspecbleach
m_dep = meson.get_compiler('c').find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
dep = [m_dep, fftw_dep]

# Get the host operating system and cpu architecture
current_os = host_machine.system()
current_arch = build_machine.cpu_family()

# Shared c_args for libraries
lib_c_args = []

# Add default x86 and x86_64 optimizations
if current_arch == 'x86' or current_arch == 'x86_64' and current_os != 'darwin'
    lib_c_args += ['-msse','-msse2','-mfpmath=sse','-ffast-math','-fomit-frame-pointer','-fno-finite-math-only']
endif

# Public header
inc = include_directories('include')
install_headers(specbleach_headers)

# Build of the shared object
libspecbleach = library('specbleach',
  shared_sources,
  specbleach_sources,
  c_args: lib_c_args,
  include_directories : inc,
  dependencies: dep,
  install: true)

# Generate package module configuration
pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : libspecbleach,
                 version: meson.project_version(),
                 name : 'libspecbleach',
                 description : 'C library for audio noise reduction and other spectral effects')

# Declare libspecbleach dependency
libspecbleach_dep = declare_dependency(include_directories : inc,
  link_with : libspecbleach)

# Example compilation
if get_option('enable_example')
  sndfile_dep = meson.get_compiler('c').find_library('sndfile',
                  has_headers: ['sndfile.h'],
                  required: true)

  executable('adenoiser_demo',
    example_sources,
    dependencies : [sndfile_dep,libspecbleach_dep],
    install: true)
endif
