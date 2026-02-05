# Tests configuration
test_data_dir = meson.current_source_dir() / 'test_data'
test_c_args = lib_c_args + ['-DTEST_DATA_DIR="@0@/"'.format(test_data_dir)]
test_link_args = lib_link_args

# Unit tests for utility functions
test_utils = executable('test_utils',
  sources: ['test_utils.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Integration tests for full pipeline
test_integration = executable('test_integration',
  sources: ['test_integration.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Audio regression tests
test_audio_regression = executable('test_audio_regression',
  sources: ['test_audio_regression.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Real audio file tests
test_audio_files = executable('test_audio_files',
  sources: ['test_audio_files.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Audio file regression tests
test_audio_file_regression = executable('test_audio_file_regression',
  sources: ['test_audio_file_regression.c'],
  dependencies: [libspecbleach_dep, m_dep, sndfile_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# FFT tests
test_fft = executable('test_fft',
  sources: ['test_fft.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Spectral Whitening tests
test_spectral_whitening = executable('test_spectral_whitening',
  sources: ['test_spectral_whitening.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# STFT Latency tests
test_stft_latency = executable('test_stft_latency',
  sources: ['test_stft_latency.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Gain Estimation tests
test_gain_estimation = executable('test_gain_estimation',
  sources: ['test_gain_estimation.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Pre-estimation tests
test_pre_estimation = executable('test_pre_estimation',
  sources: ['test_pre_estimation.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Post-filter tests
test_postfilter = executable('test_postfilter',
  sources: ['test_postfilter.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Noise profile tests
test_noise_profile = executable('test_noise_profile',
  sources: ['test_noise_profile.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Specbleach denoiser tests
test_specbleach_denoiser = executable('test_specbleach_denoiser',
  sources: ['test_specbleach_denoiser.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Noise Floor Manager tests
test_noise_floor_manager = executable('test_noise_floor_manager',
  sources: ['test_noise_floor_manager.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Adaptive noise estimator tests
test_adaptive_noise_estimator_null = executable('test_adaptive_noise_estimator_null',
  sources: ['test_adaptive_noise_estimator_null.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: [inc, include_directories('../src')]
)

# NLM filter tests
test_nlm_filter = executable('test_nlm_filter',
  sources: ['test_nlm_filter.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: [inc, include_directories('../src')]
)

# 2D Denoiser wrapper tests
test_specbleach_2d_denoiser = executable('test_specbleach_2d_denoiser',
  sources: ['test_specbleach_2d_denoiser.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)

# Run tests
test('utils', test_utils)
test('fft', test_fft)
test('spectral_whitening', test_spectral_whitening)
test('noise_floor_manager', test_noise_floor_manager)
test('stft_latency', test_stft_latency)
test('integration', test_integration)
test('audio_regression', test_audio_regression)
test('audio_files', test_audio_files)
test('audio_file_regression', test_audio_file_regression)
test('gain_estimation', test_gain_estimation)
test('pre_estimation', test_pre_estimation)
test('postfilter', test_postfilter)
test('noise_profile', test_noise_profile)
test('specbleach_denoiser', test_specbleach_denoiser)
test('specbleach_2d_denoiser', test_specbleach_2d_denoiser)

test('adaptive_noise_estimator_null', test_adaptive_noise_estimator_null)
test('nlm_filter', test_nlm_filter)

test_martin = executable('test_martin',
  sources: ['test_martin_noise_estimator.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)
test('martin', test_martin)

test_brandt = executable('test_brandt',
  sources: ['test_brandt_noise_estimator.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)
test('brandt', test_brandt)

test_spp_mmse = executable('test_spp_mmse',
  sources: ['test_spp_mmse_noise_estimator.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  link_args: test_link_args,
  include_directories: inc
)
test('spp_mmse', test_spp_mmse)

