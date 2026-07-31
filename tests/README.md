# Test Data and Reference Files

This directory contains test data and infrastructure for validating the libspecbleach audio processing library.

## Test Data Files

- `Speech.wav` - Original test audio file containing speech with background noise
- `Speech_denoised.wav` - Reference output from spectral denoising algorithm
- `Speech_adaptive_denoised.wav` - Reference output from adaptive denoising algorithm

## Reference File Generation

The reference denoised files are generated from the original `Speech.wav` using the library's denoising algorithms. These files serve as regression tests to ensure the denoising produces consistent results.

### Regenerating Reference Files

If you need to update the reference files (e.g., after algorithm changes), use the provided script. 

> [!NOTE]
> This script is a manual utility and is **not** part of the automated CTest suite. It is used to establish the "ground truth" for the regression tests.

```bash
# Generate reference files automatically with canonical parameters
./tests/generate_reference_files.sh
```

Or manually using the demo applications:

```bash
# Build the examples
cmake -B build -DENABLE_EXAMPLES=ON
cmake --build build

# Generate reference files with canonical parameters
./build/denoiser_demo --reduction 20.0 --whitening 50.0 --smoothing 0.0 --masking-depth 0.5 --masking-elasticity 0.1 tests/test_data/Speech.wav tests/test_data/Speech_denoised.wav
./build/denoiser_demo --reduction 20.0 --whitening 50.0 --smoothing 0.0 --masking-depth 0.5 --masking-elasticity 0.1 --learn-avg 3 tests/test_data/Speech.wav tests/test_data/Speech_adaptive_denoised.wav
```

### Why Reference Files Aren't Committed

- **Large binary files**: Audio files are large and bloat the repository.
- **Regenerable**: They can be recreated from the original audio using the algorithms.
- **Version control**: Focus on committing the source audio, not processed outputs.
- **CI/CD**: Reference files are generated during testing if needed.

## Running Tests

```bash
# Configure and build tests
cmake -B build -DENABLE_TESTS=ON
cmake --build build

# Run all tests
cd build && ctest --output-on-failure

# Run specific test suites
ctest -R test_utils --output-on-failure
ctest -R test_integration --output-on-failure
ctest -R test_audio_regression --output-on-failure
ctest -R test_audio_files --output-on-failure
ctest -R test_audio_file_regression --output-on-failure
```

## Test Categories

- **Unit Tests** (`test_utils`, `test_fft`, etc.): Validate individual mathematical functions
- **Integration Tests** (`test_integration`): Test full denoising pipeline
- **Audio Regression** (`test_audio_regression`): Synthetic audio processing validation
- **File Validation** (`test_audio_files`): Ensure test data accessibility
- **Audio File Regression** (`test_audio_file_regression`): Compare library output against reference `.wav` files bit-for-bit.
