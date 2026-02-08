#!/bin/bash

# Script to generate reference denoised audio files for testing
# Usage: ./tests/generate_reference_files.sh

set -e

echo "Generating reference denoised audio files..."

# Check if we're in the right directory
if [ ! -f "meson.build" ]; then
    echo "Error: Run this script from the project root directory"
    exit 1
fi

# Check if test data exists
if [ ! -f "tests/test_data/Speech.wav" ]; then
    echo "Error: tests/test_data/Speech.wav not found"
    echo "Make sure the test audio file is present"
    exit 1
fi

# Build examples to ensure latest changes are included
echo "Building examples..."
if [ ! -d "build" ]; then
    meson setup build -Denable_examples=true
else
    meson configure build -Denable_examples=true
fi
meson compile -C build

# Generate reference files with explicit parameters to ensure consistency
echo "Generating spectral denoiser reference..."
./build/examples/denoiser_demo \
    --reduction 20.0 \
    --whitening 50.0 \
    --smoothing 0.0 \
    --masking-depth 0.5 \
    --masking-elasticity 0.1 \
    tests/test_data/Speech.wav tests/test_data/Speech_denoised.wav

echo "Generating adaptive denoiser reference..."
./build/examples/denoiser_demo \
    --reduction 20.0 \
    --whitening 50.0 \
    --smoothing 0.0 \
    --masking-depth 0.5 \
    --masking-elasticity 0.1 \
    --adaptive \
    --frame-size 20.0 \
    --learn-frames 0 \
    tests/test_data/Speech.wav tests/test_data/Speech_adaptive_denoised.wav

echo "Reference files generated successfully!"
echo "Files created:"
echo "  - tests/test_data/Speech_denoised.wav"
echo "  - tests/test_data/Speech_adaptive_denoised.wav"

