# Tests configuration
test_c_args = lib_c_args

# Unit tests for utility functions
test_utils = executable('test_utils',
  sources: ['test_utils.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  include_directories: inc
)

# Integration tests for full pipeline
test_integration = executable('test_integration',
  sources: ['test_integration.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  include_directories: inc
)

# Audio regression tests
test_audio_regression = executable('test_audio_regression',
  sources: ['test_audio_regression.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  include_directories: inc
)

# Real audio file tests
test_audio_files = executable('test_audio_files',
  sources: ['test_audio_files.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  include_directories: inc
)

# Audio file regression tests
test_audio_file_regression = executable('test_audio_file_regression',
  sources: ['test_audio_file_regression.c'],
  dependencies: [libspecbleach_dep, m_dep, sndfile_dep],
  c_args: test_c_args,
  include_directories: inc
)

# FFT tests
test_fft = executable('test_fft',
  sources: ['test_fft.c'],
  dependencies: [libspecbleach_dep, m_dep, fftw_dep],
  c_args: test_c_args,
  include_directories: inc
)

# Spectral Whitening tests
test_spectral_whitening = executable('test_spectral_whitening',
  sources: ['test_spectral_whitening.c'],
  dependencies: [libspecbleach_dep, m_dep],
  c_args: test_c_args,
  include_directories: inc
)

# Run tests
test('utils', test_utils)
test('fft', test_fft)
test('spectral_whitening', test_spectral_whitening)
test('integration', test_integration)
test('audio_regression', test_audio_regression)
test('audio_files', test_audio_files)
test('audio_file_regression', test_audio_file_regression)
