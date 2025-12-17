# libspecbleach

[![build](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml/badge.svg)](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)

C library for audio noise reduction and other spectral effects

## Table of Contents

- [Background](#background)
- [De-noise Algorithms](#de-noise-algorithms)
- [Build](#build)
- [Installation](#installation)
- [Example](#example)
- [Development](#development)
- [Contributing](#contributing)
- [License](#license)

## Background

This library is based on the algorithms that were used in [noise-repellent](https://github.com/lucianodato/noise-repellent). These were extracted into a standalone library to remove the lv2 dependency. It was designed to be extensible and modular. It uses the concept of a spectral processor which itself uses a short time Fourier transform (STFT) to process the audio. There are two initial processors in place, one which uses the adaptive part of noise repellent and one that uses the manual capturing profile based denoising. The library could be extended with more spectral processors using any STFT-based algorithm such as de-crackle, de-click and other audio restoration algorithms.

## De-noise Algorithms

There are several techniques implemented in the library that are being used in the denoisers, such as masking thresholds estimation, onset detectors, etc. All these are being used in conjunction to improve the very basic spectral subtraction algorithm. Most of the papers used are listed in the wiki of the project. Also a block diagram is provided to explain the reduction architecture.

## Build

If you wish to compile yourself and install the library you will need:
- A C compiling toolchain (GCC or Clang)
- [Meson](https://mesonbuild.com/) build system (0.60.0 or newer)
- [Ninja](https://ninja-build.org/) build tool
- [FFTW3](http://www.fftw.org/) library (float version)
- [libsndfile](https://github.com/libsndfile/libsndfile) (optional, for examples)

## Installation

```bash
git clone https://github.com/lucianodato/libspecbleach.git
cd libspecbleach
meson setup build --buildtype=release
meson compile -C build
sudo meson install -C build
```


## Example

Simple console apps examples are provided to demonstrate how to use the library. It needs libsndfile to compile successfully. You can build them with:

```bash
meson setup build --buildtype=release -Denable_examples=true
meson compile -C build
```

### Adaptive noise reduction

```bash
./build/example/adenoiser_demo <input file name> <output file name>
```

### Manual noise reduction

```bash
./build/example/denoiser_demo <input file name> <output file name>
```

It will recognize any libsndfile supported format.

## Development

### Building for Development

For development builds with debugging symbols and warnings:

```bash
meson setup build --buildtype=debug -Dwarning_level=3
meson compile -C build
```

### Code Formatting

The project uses `clang-format` for code formatting. To format the code:

```bash
find . -type f \( -name "*.c" -o -name "*.h" \) ! -path "./build*" | xargs clang-format -i
```

### Running Tests

If tests are enabled:

```bash
meson setup build -Denable_tests=true
meson test -C build
```

## Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute to this project.

## License

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

See [LICENSE](LICENSE) for more details.
