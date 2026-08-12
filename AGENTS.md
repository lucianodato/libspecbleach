# libspecbleach Agent Context

💡 *Foundational mandates and architectural rules specifically for `libspecbleach`.*

## Foundational Mandates

1. **DSP Integrity:** Maintain complete determinism in the processing path.
2. **SIMD Awareness & FTZ/DAZ:** 
   - The 2D Denoising (NLM) filter is sensitive to optimizations. Verify SIMD auto-vectorization or explicit vector instructions on changes.
   - For real-time safety, always enable FTZ (Flush-To-Zero) and DAZ (Denormals-Are-Zero) using `sb_simd_enable_ftz_daz()` and restore state with `sb_simd_restore_state(...)` during heavy SIMD processing loops.
3. **`sb_sel8` Mask Contract:** The mask argument for `sb_sel8` must be normalized (any non-zero lane treated as true) to guarantee consistent behavior across SSE, AVX, ARM NEON, and Scalar backends.
4. **Reference Generation:** Regenerate reference audio files (`./tests/generate_reference_files.sh`) only when algorithm intentional changes occur.

## Project Structure

- **Core Logic:** `src/shared/` — Core math operations.
- **Processors:** `src/processors/` — Orchestrates shared modules into the public API.
- **SIMD Utilities:** `src/shared/utils/simd_utils.h` — Explicit SIMD abstractions.
- **Memory Management:** Uses a handle-based pattern initialized via `specbleach_initialize`.

## Build & Test Commands

```bash
# Build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j4

# Run Audio Regression Tests
cd build && ctest -R test_audio_file_regression --output-on-failure
```
## Architectural Notes

- **Dynamic AVX Dispatch**: The NLM filter uses runtime dynamic dispatch. `nlm_filter.c` initializes `process_fn` check with `__builtin_cpu_supports("avx")` and maps to `nlm_filter_process_avx` (compiled under the `libspecbleach_avx` target with `-mavx`) or `nlm_filter_process_generic`.
- **Boundary Safety**: NLM patches (e.g., 8x8 or 4x4) must check frequency and time bounds, falling back to safe clamping paths on frequency boundaries.
- **Circular Buffering**: The library relies on `SbSpectralCircularBuffer` for efficient time-frequency access.
- **STFT**: Framing and windowing are standardized across the library (see `configurations.h`).
- **Memory**: The library uses a handle-based pattern. Memory is allocated during `specbleach_initialize` and must be strictly managed.

## C Coding Standards

When editing any source file (`.c` or `.h`) in `libspecbleach` or `noise-repellent`, follow these standards:

### 1. Portable OpenMP Loop Variable Declarations
- **Loop Counters Outside For Statements**: Always declare OpenMP loop counter variables **before** the `#pragma omp parallel for` directive (e.g., `uint32_t i = 0; #pragma omp parallel for ... for (i = 0; ...)`). Never declare loop variables directly inside the `for` statement (`for (uint32_t i = 0; ...)`). MSVC's OpenMP parser enforces C89/C90 variable scoping rules and will fail with `error C3015: initialization in OpenMP 'for' statement has improper form`.

### 2. Clean Unused Includes
- Review the `#include` block at the top of the file.
- Verify if any symbol (function, struct, define, typedef) from that header is actually used.
- **If unused, REMOVE the `#include` statement.**

### 2. Enforce LGPL Header
- Every source file must start with the following standard license header:
```c
/*
libspecbleach - A spectral processing library

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
```

### 3. Testing Conventions
- **Unit Tests**: Maintain a 1:1 mapping between source modules and unit test files (e.g., `my_module.c` -> `test_my_module.c`).
- **Integration Tests**: Ensure the filename contains the keyword `integration` (e.g., `test_integration_audio_flow.c`).

### 4. Code Formatting & Quality (Mandatory)
- **Always Lint Files Before Finishing**: Run `rtk clang-format -i <modified_files>` on every modified or newly created `.c` and `.h` file before completing tasks so CI lint checks pass.
- Run `clang-tidy` to check for static analysis warnings or style violations.

### 5. CMake Build Standards
- **Reusing Build Folders**: Do **NOT** delete the build directory (`rm -rf build`) unless it is corrupted.
- **Reconfiguration**: Use `cmake -B build [options]` to reconfigure the build directory while preserving existing object files.

---

## Real-Time DSP Programming Guidelines

These guidelines are mandatory for writing code that executes in the real-time audio thread:

### 1. Prohibited Operations in the Audio Callback
The audio thread runs with strict system deadlines. Unbounded execution times will cause audible dropouts/glitches.
*   **No Memory Allocation**: Never call `malloc`, `free`, `calloc`, `realloc`, or `new`/`delete`. All dynamic buffers must be pre-allocated during initialization.
*   **No Locks or Mutexes**: Never use blocking synchronization primitives like `pthread_mutex_lock` or `std::mutex`. Use lock-free Single-Producer-Single-Consumer (SPSC) ring buffers or atomic variables instead.
*   **No File I/O or System Calls**: Do not perform file reads/writes, network sockets, or console printing (like `printf` or standard logging). Log events to a lock-free queue and consume them on a low-priority thread.

### 2. Denormal Number Prevention
Extremely small floating-point values near zero can trigger costly FPU exceptions and massive CPU spikes.
*   **Flush-to-Zero (FTZ) & Denormals-are-Zero (DAZ)**: Enable hardware-level FTZ/DAZ during heavy math loops.
*   **DC Bias**: Alternatively, if a DC bias is proposed to prevent values from decaying to denormal thresholds, it must be defined as a named tuning constant in `src/shared/configurations.h` (do not use hardcoded magic values) and requires algorithm-specific regression validation to verify its effect on filter response and check for unwanted DC offset.

### 3. Inter-Thread Communication
*   **Atomics**: Use atomic types (C11 `stdatomic.h` or C++ `std::atomic`) for thread-safe controls (bypass, levels, parameters).
*   **SPSC Ring Buffers**: Transfer blocks of data, spectral states, or events using thread-safe, lock-free ring buffers. Never resize ring buffers inside the audio callback.

### 4. Performance Optimizations
*   **Memory Alignment**: Align buffers to SIMD boundaries (16 or 32 bytes) to ensure compiler auto-vectorization is optimal.
*   **Reciprocal Multiplication**: Use multiplications by reciprocals ($1.0 / x$) instead of division operations inside hot processing loops.

