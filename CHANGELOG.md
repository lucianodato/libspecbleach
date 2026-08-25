#Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.2] - Unreleased

### Added
- **Internal Thread Pool**: Added `SbThreadPool` (`src/shared/utils/thread_pool.h`), a fixed-size worker pool with semaphore-based dispatch and static contiguous partitioning, powering multi-threaded NLM 2D smoothing without any external threading runtime. Thread count is configurable per instance via `NlmFilterConfig::num_threads` (default `NLM_NUM_THREADS_DEFAULT`).

### Improved & Refactored
- **OpenMP Removal**: Replaced OpenMP parallelization in the NLM 2D filter with the internal worker pool. Binaries no longer depend on `libomp`/`libgomp`/`libomp140` runtimes, dispatches are deterministic (static partitioning instead of dynamic scheduling), and no threads or locks are created in the audio path.
- **Adaptive Profile Persistence**: Updated standalone adaptive noise estimation to persist learned noise profiles to the noise profile manager so estimated spectral curves remain available when switching to manual mode.
- **Manual Baseline Re-seeding**: Fixed standalone adaptive mode state tracking to reset hybrid initialization state when no manual profile exists, ensuring manual noise profiles captured while adaptive mode is active immediately morph and seed the baseline floor.

### Fixed


## [0.3.1] - 2026-08-07

### Added
- **`pkg-config` Support**: Added `specbleach.pc.in` template to generate and install `specbleach.pc` into system pkg-config directories for downstream dependency lookup.
- **System Installation Targets**: Added standard CMake `install()` targets for header files (`${
  CMAKE_INSTALL_INCLUDEDIR}`), shared/static libraries (`${
  CMAKE_INSTALL_LIBDIR}`), and CMake export targets.
- **Configurable FFTW3 Dependency**: Added `USE_SYSTEM_FFTW` option (default `ON`) allowing downstream builders to link system-installed `libfftw3f` or automatically fetch and compile static FFTW3 via `FetchContent`.
- **Autoresearch Framework**: Integrated automated optimization & research loop scripts (`autoresearch/`).

### Improved & Refactored
- **NLM 2D Denoising**: Resolved index bounds and array access inconsistencies in the Non-Local Means time-frequency smoothing algorithm.
- **Masking Veto**: Refined psychoacoustic masking veto model for improved transient preservation and noise estimation thresholding.
- **Tonal Masking & Detector**: Enhanced tonal component calculation and peak detection accuracy.
- **Parameter Mappings**: Simplified parameter scaling by replacing `remap_percentage_log_like_unity` with linear scaling for smoothing factors.
- **Spectral Whitening & OpenMP**: Optimized spectral whitening inner loops and updated OpenMP target detection.

### Changed
- **Default Library Type**: Changed default `BUILD_SHARED_LIBS` setting to `ON` to conform with Linux distribution packaging guidelines.
- **Standardized Install Paths**: Integrated `GNUInstallDirs` to ensure installation paths adhere to target OS and distribution filesystem standards.
- **Cleaned Dependency Resolution**: Removed OS-specific hardcoded search paths (`/opt/homebrew`, `/usr/local`) from OpenMP and `libsndfile` detection routines to prevent interference with isolated sandbox builds.

## [0.3.0] - 2026-07-31

### Changed
- **Build System Migration**: Migrated primary build system from Meson to CMake 3.16+.
- Updated CTest test suite configuration and options (`ENABLE_TESTS`, `ENABLE_EXAMPLES`, `ENABLE_AVX`, `ENABLE_SANITIZERS`).
- Improved CMake target export (`libspecbleach::libspecbleach`) for easy downstream inclusion via `add_subdirectory` or `find_package`.

## [0.2.1] - 2026-01-28
### Fixed
- **Critical**: Fixed NLM filter ring buffer indexing bug. The filter was using past search range instead of future look-ahead for target frame alignment, causing significant output latency.
### Added
- Added `whitening_factor` support to `Spectral2DDenoiser`, replacing usage of `noise_rescale`.

## [0.2.0] - 2026-01-28

### Added
- SOVERSION support for proper shared library versioning (SOVERSION 0)
- Versioning and ABI stability documentation
- This CHANGELOG to track project changes
- Regression test suite for audio file processing
- Code coverage reporting with codecov integration
- Demo applications with configurable command-line parameters
- Comprehensive unit tests for all core modules
- `libsndfile` as optional dependency

### Changed
- Made demo applications accept command-line arguments for all processing parameters
- Improved STFT processor input latency calculation
- Enhanced memory initialization and error handling across all modules
- Library now generates versioned shared objects (libspecbleach.so.0.2.0, libspecbleach.so.0)

### Fixed
- STFT input latency bug causing incorrect delay calculations
- Memory initialization issues causing CI failures under MALLOC_PERTURB_
- Various memory leaks and uninitialized memory access patterns
- `M_PI` definition issues in test suite
- Code formatting consistency across the codebase

## [0.1.0] - Initial Release

### Added
- Core spectral processing framework with STFT
- Adaptive noise reduction processor
- Profile-based noise reduction processor
- Spectral gating, whitening, and other audio effects
- FFT-based frequency analysis utilities
- Meson build system with comprehensive options
- Basic demo applications
- LGPL 2.1 licensing
- GitHub Actions CI/CD pipeline

---

## Versioning Guidelines

This project follows **Semantic Versioning**:

- **MAJOR version** (X.y.z): Incompatible API/ABI changes
  - When MAJOR changes, SOVERSION will be incremented
  - Breaking changes to public API or ABI
  
- **MINOR version** (x.Y.z): New features, backward-compatible
  - New functionality added without breaking existing code
  - New functions or capabilities
  
- **PATCH version** (x.y.Z): Backward-compatible bug fixes
  - Bug fixes that don't change API/ABI
  - Documentation updates
  - Internal improvements

### ABI Stability Promise

- Libraries with the same **SOVERSION** are binary compatible
- Current **SOVERSION: 0** (pre-1.0 release, API may still evolve)
- SOVERSION will be incremented to 1 when API reaches stability (1.0.0 release)
- After 1.0.0, SOVERSION will match MAJOR version number
