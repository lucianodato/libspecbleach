subdir('shared')

subdir('processors')

specbleach_sources = [shared_sources, processors_sources]

src_inc = include_directories('.')

# Add AVX specific compilation
has_avx_support = cc.has_argument('-mavx')
lib_link_whole = []

if has_avx_support and (current_arch == 'x86' or current_arch == 'x86_64')
  avx_sources = files('shared/denoiser_logic/processing/nlm_filter_avx.c')
  
  libspecbleach_avx = static_library('specbleach_avx',
    avx_sources,
    c_args: lib_c_args + ['-mavx'],
    dependencies: dep,
    include_directories: [inc, src_inc]
  )
  
  lib_link_whole += libspecbleach_avx
endif

# Build of the shared object
libspecbleach = library('specbleach',
  sources: specbleach_sources,
  c_args: lib_c_args,
  link_args: lib_link_args,
  link_whole: lib_link_whole,
  dependencies: dep,
  include_directories: [inc, src_inc],
  version: meson.project_version(),
  soversion: '0',
  install: true)

# Generate package module configuration
pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries: libspecbleach,
                 version: meson.project_version(),
                 name: 'libspecbleach',
                 description: 'C library for audio noise reduction and other spectral effects')

# Note: libspecbleach_dep is now declared in the root meson.build