# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
