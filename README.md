# libspecbleach

C library for audio noise reduction

[![main](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml)
[![development](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml/badge.svg?branch=development)](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml)

## Background

This library is based on the algorithms that were used in noise-repellent. These were extracted into a this standalone library. It was design to be extensible and modular. It uses the concept of a spectral processor which uses a short time Fourier transform (STFT) to process the audio. There are two initial processors in place, one which uses the adaptive part of noise repellent and one that uses the manual capturing profile based denoising. The library could be extended with more spectral processors using any STFT-based algorithm such as de-crackle, de-click and other audio restoration algorithms.

## De-noise algorithms

There several techniques implemented in the library that are being used in the denoisers, such as masking thresholds estimation, onset detectors, etc. All these are being used in conjuntion to improve the very basic spectral substraction algorithm. Most of the papers used are listed in the wiki of the project. Also a block diagram is provided to explain the reduction architecture.

## Build

If you wish to compile yourself and install the library you will need the a C compiling toolchain, Meson build system, ninja compiler, git and fftw3 library.

Installation:

```bash
  git clone https://github.com/lucianodato/noise-repellent.git
  cd noise-repellent
  meson build --buildtype=release --prefix=/usr --libdir=lib (your-os-appropriate-location-fullpath)
  ninja -C build -v
  sudo ninja -C build install
```
