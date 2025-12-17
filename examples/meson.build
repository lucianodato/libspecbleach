sndfile_dep = dependency('sndfile', required: true)

executable('adenoiser_demo',
  sources: 'adenoiser_demo.c',
  dependencies: [sndfile_dep, libspecbleach_dep],
  install: true)

executable('denoiser_demo',
  sources: 'denoiser_demo.c',
  dependencies: [sndfile_dep, libspecbleach_dep],
  install: true)
