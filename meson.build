project('libspecbleach', 'c', version: '0.3.0',
        default_options: ['c_std=c17', 'warning_level=2'])

# Note: Meson 0.60.0 or newer is recommended

# Get the host operating system and cpu architecture
current_os = host_machine.system()
current_arch = build_machine.cpu_family()
cc = meson.get_compiler('c')

# Shared c_args and link_args for libraries
lib_c_args = []
lib_link_args = []

# Compiler warning flags
custom_warning_level = get_option('custom_warning_level')
if custom_warning_level >= 1
    lib_c_args += cc.get_supported_arguments([
        '-Wall',
        '-Wextra',
    ])
endif

if custom_warning_level >= 2
    lib_c_args += cc.get_supported_arguments([
        '-Wpedantic',
        '-Wno-unused-parameter',
        '-Wno-unused-variable',
    ])
endif

# Treat warnings as errors if requested
if get_option('treat_warnings_as_errors')
    lib_c_args += ['-Werror']
endif

# Architecture-specific optimizations
# x86 and x86_64 optimizations (excluding macOS)
if (current_arch == 'x86' or current_arch == 'x86_64') and current_os != 'darwin'
    lib_c_args += cc.get_supported_arguments([
        '-msse',
        '-msse2',
        '-mfpmath=sse',
        '-ffast-math',
        '-fomit-frame-pointer',
        '-fno-finite-math-only'
    ])

    if get_option('enable_avx')
        lib_c_args += cc.get_supported_arguments(['-mavx'])
    endif
    if get_option('enable_avx2')
        lib_c_args += cc.get_supported_arguments(['-mavx2', '-mfma'])
    endif
endif

# macOS-specific x86_64 AVX support
if (current_arch == 'x86' or current_arch == 'x86_64') and current_os == 'darwin'
    if get_option('enable_avx')
        lib_c_args += cc.get_supported_arguments(['-mavx'])
    endif
    if get_option('enable_avx2')
        lib_c_args += cc.get_supported_arguments(['-mavx2', '-mfma'])
    endif
endif

# ARM64 optimizations (NEON)
if current_arch == 'aarch64'
    lib_c_args += cc.get_supported_arguments([
        '-mfpu=neon-fp-armv8',
    ])
endif

# Sanitizers (for debug builds)
if get_option('enable_sanitizers')
    if get_option('buildtype') == 'debug'
        sanitizer_args = []
        if get_option('sanitize_address')
            sanitizer_args += ['-fsanitize=address']
        endif
        if get_option('sanitize_undefined')
            sanitizer_args += ['-fsanitize=undefined']
        endif
        # Sanitizer flags must be passed to both compiler AND linker
        lib_c_args += cc.get_supported_arguments(sanitizer_args)
        lib_link_args += cc.get_supported_link_arguments(sanitizer_args)
    else
        warning('Sanitizers are only enabled in debug builds')
    endif
endif

# Dependencies for libspecbleach
m_dep = meson.get_compiler('c').find_library('m', required: true)
fftw_dep = dependency('fftw3f', required: true, static: get_option('static_deps'), include_type: 'system')
openmp_dep = dependency('openmp', required: false)

# Optional dependency for examples and tests
sndfile_dep = dependency('', required: false)
if get_option('enable_examples') or get_option('enable_tests')
  sndfile_dep = dependency('sndfile', required: false, include_type: 'system')
  if not sndfile_dep.found()
    warning('sndfile not found - examples and tests will not be built')
  endif
endif

dep = [m_dep, fftw_dep]
if openmp_dep.found()
  dep += openmp_dep
  lib_c_args += ['-DSB_HAS_OPENMP=1']
else
  lib_c_args += ['-DSB_HAS_OPENMP=0']
endif

# Public Headers
subdir('include')

# Build library
subdir('src')

# Export dependency for parent projects
libspecbleach_dep = declare_dependency(include_directories: [inc, src_inc],
  link_with: libspecbleach)

# Examples building
if get_option('enable_examples') and sndfile_dep.found()
  subdir('examples')
endif

# Tests building
if get_option('enable_tests')
  subdir('tests')
endif

# Code formatting target
run_target('format',
  command: [
    'sh', '-c',
    'find "@SOURCE_ROOT@" -name "*.[ch]" -o -name "*.cpp" -o -name "*.hpp" | grep -v "/build/" | xargs clang-format -i && echo "Code formatting complete!"'
  ]
)
