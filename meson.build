project('libspecbleach', 'c', version: '0.1.7', 
        default_options: ['c_std=c11', 'warning_level=2'])

# Note: Meson 0.60.0 or newer is recommended

# Get the host operating system and cpu architecture
current_os = host_machine.system()
current_arch = build_machine.cpu_family()
cc = meson.get_compiler('c')

# Shared c_args for libraries
lib_c_args = []

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
        lib_c_args += cc.get_supported_arguments(sanitizer_args)
    else
        warning('Sanitizers are only enabled in debug builds')
    endif
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
