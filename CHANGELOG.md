# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Modernized build system with C11 standard support
- Compiler warning flags (-Wall, -Wextra, -Wpedantic) with configurable levels
- ARM64 (aarch64) optimization support
- Build options for sanitizers (AddressSanitizer, UndefinedBehaviorSanitizer)
- Multi-platform CI/CD testing (Ubuntu, macOS x86_64, macOS ARM64)
- Code formatting with clang-format
- EditorConfig for consistent editor settings
- Separate lint workflow for code quality checks
- CHANGELOG.md following Keep a Changelog format
- CONTRIBUTING.md with development guidelines
- Improved README with table of contents and badges
- Updated BUILD_INSTRUCTIONS.md with modern Meson syntax

### Changed
- Upgraded C standard from C99 to C11
- Updated GitHub Actions workflows to use Ubuntu 22.04
- Updated GitHub Actions checkout action from v2 to v4
- Modernized Meson build syntax (using `meson setup` instead of `meson build`)
- Fixed README clone URL (was pointing to noise-repellent)
- Fixed missing `#include <stdio.h>` in example files

### Fixed
- Missing include directives in example demos
- Incorrect repository URL in README installation instructions
- Outdated Meson build commands in documentation

## [0.1.7] - 2022-XX-XX

### Added
- Initial release of libspecbleach
- Spectral denoiser processor
- Adaptive denoiser processor
- Example applications demonstrating library usage

[Unreleased]: https://github.com/lucianodato/libspecbleach/compare/v0.1.7...HEAD
[0.1.7]: https://github.com/lucianodato/libspecbleach/releases/tag/v0.1.7

