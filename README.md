# libspecbleach

[![build](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml/badge.svg)](https://github.com/lucianodato/libspecbleach/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/lucianodato/libspecbleach/branch/main/graph/badge.svg)](https://codecov.io/gh/lucianodato/libspecbleach)
[![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/lucianodato/libspecbleach?utm_source=oss&utm_medium=github&utm_campaign=lucianodato%2Flibspecbleach&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)](https://coderabbit.ai)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)

C library for audio noise reduction and other spectral effects.

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

This library is a standalone, modular spectral processing engine originally extracted from [noise-repellent](https://github.com/lucianodato/noise-repellent). It decouples DSP algorithms from any specific plugin framework, allowing for clean integration into C/C++ audio software.

The core architecture uses a unified spectral processor concept. The library implements advanced spectral denoising using efficient circular buffering (`SbSpectralCircularBuffer`) and modern STFT processing. It is designed to be extensible, supporting future additions like de-crackle or de-click algorithms.

## De-noise Algorithms

The library implements a sophisticated spectral gating algorithm enhanced by several psychoacoustic and statistical techniques:

### 1. Spectral Gating & Subtraction
The fundamental noise reduction method uses spectral subtraction with proprietary framing and windowing to minimize artifacts.

### 2. 2D Denoising (Time-Frequency Filtering)
A Non-Local Means (NLM) algorithm filters the spectrogram in both time and frequency domains simultaneously. This preserves structural details of the signal while reducing musical noise and "burbling" artifacts often associated with simple spectral subtraction.
*Note: This feature is computationally intensive and benefits from SIMD acceleration (enabled in Release builds).*

### 3. Masking Veto
To preserve transients and prevent over-processing, a psychoacoustic masking model estimates the auditory masking threshold. If the signal components are strong enough to mask the noise naturally, the "veto" system prevents unnecessary noise reduction, preserving the natural character of the audio.

### 4. Tonal Reduction
Specialized handling for tonal noise components allows for more aggressive reduction of static hums and whines without affecting broadband characteristics.

### 5. Transparent Whitening
The whitening feature (noise floor recovery) has been refined to be transparent at 0dB reduction, ensuring that the noise floor texture is natural and consistent with the reduction amount.

### 6. Adaptive Estimation
In addition to manual noise profile capture, the library supports adaptive noise floor estimation using algorithms like SPP-MMSE, Brandt, and Martin Minimum Statistics.

## Build

To compile and install `libspecbleach`, you will need:
- A C compiling toolchain (GCC or Clang supporting C17)
- [CMake](https://cmake.org/) (3.16 or newer)
- [FFTW3](http://www.fftw.org/) library (`libfftw3f`)
- [OpenMP](https://www.openmp.org/) for parallel processing (optional, recommended for NLM 2D denoising)
- [libsndfile](https://github.com/libsndfile/libsndfile) (optional, for test suite and demo tools)

## Installation

```bash
git clone https://github.com/lucianodato/libspecbleach.git
cd libspecbleach
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j4
sudo cmake --install build
```

## Build Options

You can configure the build using `-Doption=VALUE`:

- `BUILD_SHARED_LIBS`: Build shared library instead of static library (default: `OFF`).
- `ENABLE_AVX`: Enable AVX SIMD optimizations on x86_64 architecture (default: `ON`).
- `ENABLE_TESTS`: Build unit, integration, and audio regression test suite (default: `ON`).
- `ENABLE_EXAMPLES`: Build demo executables (default: `ON`, requires `libsndfile`).
- `ENABLE_SANITIZERS`: Enable AddressSanitizer and UndefinedBehaviorSanitizer (default: `OFF`).

> [!IMPORTANT]
> **Critical Performance Note for Packagers**: The advanced "2D Denoising" (NLM) feature is computationally intensive and relies heavily on SIMD vectorization, function inlining, and **multi-core parallelization via OpenMP**. Builds without OpenMP will skip OpenMP-backed NLM parallelization and may run slower.
>
> You **MUST** compile with `-DCMAKE_BUILD_TYPE=Release` (or `-O3`) to ensure usability. Debug or unoptimized builds will result in excessive CPU usage and audio dropouts/xruns.

Example for a static build with tests and examples:
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTS=ON -DENABLE_EXAMPLES=ON
cmake --build build --config Release
```

## Usage Examples

Console demo applications demonstrate library usage. They require `libsndfile` to build:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_EXAMPLES=ON
cmake --build build
```

### Manual Noise Profile

To process an audio file using a manually captured noise profile (first N frames):

```bash
./build/denoiser_demo --learn-frames 10 input.wav output.wav
```

### Adaptive Mode

To use the adaptive noise estimator:

```bash
./build/denoiser_demo --adaptive input.wav output.wav
```

### Full Options

```bash
./build/denoiser_demo \
  --adaptive \
  --reduction 20 \
  --whitening 50 \
  --smoothing 0.0 \
  input.wav output.wav
```

## Development

### Building for Development

For development builds with debug symbols:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_TESTS=ON -DENABLE_EXAMPLES=ON
cmake --build build
```

### Code Formatting

The project uses `clang-format`. Format modified C files before submitting:

```bash
find src include -type f \( -name "*.c" -o -name "*.h" \) | xargs clang-format -i
```

### Running Tests

Run the full CTest suite:

```bash
ctest --test-dir build --output-on-failure
```

Or run individual test executables:

```bash
./build/test_integration
./build/test_audio_file_regression # requires libsndfile
```

## License

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

See [LICENSE](LICENSE) for more details.
