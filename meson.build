project('libspecbleach', 'c', version: '0.1.7', default_options: ['c_std=c99'])

# Get the host operating system and cpu architecture
current_os = host_machine.system()
current_arch = build_machine.cpu_family()

# Shared c_args for libraries
lib_c_args = []

# Default x86 and x86_64 optimizations
if current_arch == 'x86' or current_arch == 'x86_64' and current_os != 'darwin'
    lib_c_args += ['-msse','-msse2','-mfpmath=sse','-ffast-math','-fomit-frame-pointer','-fno-finite-math-only']
endif

# Dependencies for libspecbleach
m_dep = meson.get_compiler('c').find_library('m', required: true)
fftw_dep = dependency('fftw3f', required: true)
dep = [m_dep, fftw_dep]

# Public Headers
subdir('include')

# Build library
subdir('src')

# Examples building
if get_option('enable_examples')
  subdir('example')
endif