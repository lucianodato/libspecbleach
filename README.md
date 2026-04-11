# libspecbleach

[![build](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml/badge.svg)](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/lucianodato/libspecbleach/branch/main/graph/badge.svg)](https://codecov.io/gh/lucianodato/libspecbleach)
[![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/lucianodato/libspecbleach?utm_source=oss&utm_medium=github&utm_campaign=lucianodato%2Flibspecbleach&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)](https://coderabbit.ai)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)

C library for audio noise reduction and other spectral effects

## Table of Contents

- [Background](#background)
- [De-noise Algorithms](#de-noise-algorithms)
- [Build](#build)
- [Installation](#installation)
- [Usage Examples](#usage-examples)
- [Development](#development)
- [Contributing](#contributing)
- [License](#license)

## Background

This library is a standalone, modular spectral processing engine originally based on [noise-repellent](https://github.com/lucianodato/noise-repellent). It decouples the DSP algorithms from any specific plugin API (like LV2), allowing for easy integration into various applications.

The core architecture uses a unified spectral processor concept. The library currently implements advanced spectral denoising using efficient circular buffering (`SbSpectralCircularBuffer`) and modern STFT processing. It is designed to be extensible, supporting future additions like de-crackle or de-click algorithms.

## De-noise Algorithms

The library implements a sophisticated spectral gating algorithm enhanced by several psychoacoustic and statistical techniques:

### 1. Spectral Gating & Subtraction
The fundamental noise reduction method uses spectral subtraction with proprietary framing and windowing to minimize artifacts.

### 2. 2D Denoising (Time-Frequency Filtering)
A Non-Local Means (NLM) algorithm filters the spectrogram in both time and frequency domains simultaneously. This preserving structural details of the signal while reducing musical noise and "burbling" artifacts often associated with simple spectral subtraction.
*Note: This feature is computationally intensive and requires SIMD optimization (enabled in Release builds).*

### 3. Masking Veto
To preserve transients and prevent over-processing, a psychoacoustic masking model estimates the auditory masking threshold. If the signal components are strong enough to mask the noise naturally, the "veto" system prevents unnecessary noise reduction, preserving the natural character of the audio.

### 4. Tonal Reduction
Specialized handling for tonal noise components allows for more aggressive reduction of static hums and whines without affecting broadband characteristics.

### 5. Transparent Whitening
The whitening feature (noise floor recovery) has been refined to be transparent at 0dB reduction, ensuring that the noise floor texture is natural and consistent with the reduction amount.

### 6. Adaptive Estimation
In addition to manual noise profile capture, the library supports adaptive noise usage for changing noise environments.

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

## Build Options

You can configure the build using `-Doption=value`:

- `enable_examples`: Build example applications (default: `false`). Requires `libsndfile`.
- `enable_tests`: Build unit and integration tests (default: `false`). Requires `libsndfile`.
- `static_deps`: Link internal dependencies (like FFTW3) statically (default: `false`). Useful for creating self-contained libraries.
- `custom_warning_level`: 0-3 (default: `2`). Controls compiler warning verbosity.
- `treat_warnings_as_errors`: Treat compiler warnings as errors (default: `false`).
- `enable_sanitizers`: Enable sanitizers in debug builds (default: `false`).
- `sanitize_address`: Enable AddressSanitizer (default: `false`).
- `sanitize_undefined`: Enable UndefinedBehaviorSanitizer (default: `false`).
- `lv2dir`: Install directory for LV2 bundles (absolute path or relative to prefix) (default: '').

> [!IMPORTANT]
> **Critical Performance Note for Packagers**: The advanced "2D Denoising" (NLM) feature is computationally intensive and relies heavily on SIMD vectorization and function inlining. You **MUST** compile with `--buildtype=release` (or `-O3`) to ensure usability. Debug or unoptimized builds will result in excessive CPU usage and audio dropouts/xruns.

Example for a static build with examples:
```bash
meson setup build -Dstatic_deps=true -Denable_examples=true
meson compile -C build
```

## Usage Examples

Simple console apps examples are provided to demonstrate how to use the library. It needs libsndfile to compile successfully. You can build them with:

```bash
meson setup build --buildtype=release -Denable_examples=true
meson compile -C build
```

### Manual Noise Profile

To process a file using a manually captured noise profile (first N frames):

```bash
./build/example/denoiser_demo --learn-frames 10 input.wav output.wav
```

### Adaptive Mode

To use the adaptive noise estimator:

```bash
./build/example/denoiser_demo --adaptive input.wav output.wav
```

### Full Options

```bash
./build/example/denoiser_demo \
  --adaptive \
  --reduction 20 \
  --whitening 50 \
  --smoothing 0.0 \
  input.wav output.wav
```

It will recognize any libsndfile supported format.

## Development

### Building for Development

For development builds with debugging symbols:

```bash
meson setup build --buildtype=debug
meson compile -C build
```

### Code Formatting

The project uses `clang-format` for code formatting. To format the code:

```bash
meson compile format -C build
```

### Running Tests

If tests are enabled:

```bash
meson setup build -Denable_tests=true
meson test -C build
```

### Coverage

To generate coverage reports locally, you will need `gcovr` or `lcov` installed.

```bash
meson setup build --buildtype=debug -Db_coverage=true
meson compile -C build
meson test -C build
ninja -C build coverage-html
```

The report will be available in `build/meson-logs/coveragereport/index.html`.

## License

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

See [LICENSE](LICENSE) for more details.
