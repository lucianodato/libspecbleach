# Test Data and Reference Files

This directory contains test data and infrastructure for validating the libspecbleach audio processing library.

## Test Data Files

- `Speech.wav` - Original test audio file containing speech with background noise
- `Speech_denoised.wav` - Reference output from spectral denoising algorithm
- `Speech_adaptive_denoised.wav` - Reference output from adaptive denoising algorithm

## Reference File Generation

The reference denoised files are generated from the original `Speech.wav` using the library's denoising algorithms. These files serve as regression tests to ensure the denoising produces consistent results.

### Regenerating Reference Files

If you need to update the reference files (e.g., after algorithm changes), use the provided script:

```bash
# Generate reference files automatically
./tests/generate_reference_files.sh
```

Or manually:

```bash
# Build the examples
meson setup build -Denable_examples=true
meson compile -C build

# Generate reference files
./build/examples/denoiser_demo tests/test_data/Speech.wav tests/test_data/Speech_denoised.wav
./build/examples/adenoiser_demo tests/test_data/Speech.wav tests/test_data/Speech_adaptive_denoised.wav
```

### Why Reference Files Aren't Committed

- **Large binary files**: Audio files are large and bloat the repository
- **Regenerable**: Can be recreated from the original audio using the algorithms
- **Version control**: Focus on committing the source audio, not processed outputs
- **CI/CD**: Reference files are generated during testing if needed

## Running Tests

```bash
# Run all tests
meson test -C build

# Run specific test suites
meson test -C build --suite utils
meson test -C build --suite integration
meson test -C build --suite audio_regression
meson test -C build --suite audio_files
```

## Test Categories

- **Unit Tests** (`test_utils`): Validate individual mathematical functions
- **Integration Tests** (`test_integration`): Test full denoising pipeline
- **Audio Regression** (`test_audio_regression`): Synthetic audio processing validation
- **File Validation** (`test_audio_files`): Ensure test data accessibility
