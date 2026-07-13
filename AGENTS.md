# libspecbleach Agent Context

This file contains foundational mandates and architectural context for Gemini CLI when working on `libspecbleach`.

## Foundational Mandates

1. **DSP Integrity**: Never introduce non-deterministic logic or blocking calls (locks, I/O, malloc) into the processing path (functions like `specbleach_process` or anything called within it).
2. **Scientific Rigor**: All tuning constants MUST be defined in `src/shared/configurations.h`. Do not use magic numbers in implementation files.
3. **Regression Testing**: After any change to DSP logic, you MUST run the audio regression suite:
   - `meson test -C build --suite audio_file_regression`
4. **SIMD Awareness & FTZ/DAZ**: The "2D Denoising" (NLM) feature is extremely sensitive to optimization. Always verify that changes don't break SIMD auto-vectorization or explicit vector instructions.
   - **FTZ/DAZ**: For real-time safety, always enable FTZ (Flush-To-Zero) and DAZ (Denormals-Are-Zero) using `sb_simd_enable_ftz_daz()` and restore with `sb_simd_restore_state(...)` during heavy SIMD processing loops.
   - **sb_sel8 Mask Contract**: The mask argument for `sb_sel8` must be normalized (any non-zero lane treated as true) to guarantee consistent behavior across SSE, AVX, ARM NEON, and Scalar backends.
5. **No Reverts**: Do not revert changes unless they break the build or fail tests significantly.

## Project Structure & Workflow

- **Core Logic**: Located in `src/shared/`. This is where the math happens.
- **Processors**: `src/processors/` orchestrates shared modules into the public API.
- **SIMD Utilities**: Explicit SIMD abstractions are centralized in [simd_utils.h](file:///Users/luciandato/repos/libspecbleach/src/shared/utils/simd_utils.h).
- **Reference Generation**: If an algorithm change is intentional and the reference audio needs updating, use:
  - `./tests/generate_reference_files.sh`
- **Build**: Use `meson setup build --buildtype=release` for performance testing. Use `debug` for development.

## Architectural Notes

- **Dynamic AVX Dispatch**: The NLM filter uses runtime dynamic dispatch. `nlm_filter.c` initializes `process_fn` check with `__builtin_cpu_supports("avx")` and maps to `nlm_filter_process_avx` (compiled under the `libspecbleach_avx` target with `-mavx`) or `nlm_filter_process_generic`.
- **Boundary Safety**: NLM patches (e.g., 8x8 or 4x4) must check frequency and time bounds, falling back to safe clamping paths on frequency boundaries.
- **Circular Buffering**: The library relies on `SbSpectralCircularBuffer` for efficient time-frequency access.
- **STFT**: Framing and windowing are standardized across the library (see `configurations.h`).
- **Memory**: The library uses a handle-based pattern. Memory is allocated during `specbleach_initialize` and must be strictly managed.

## C Coding Standards

When editing any source file (`.c` or `.h`) in `libspecbleach` or `noise-repellent`, follow these standards:

### 1. Clean Unused Includes
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

### 4. Code Formatting & Quality
- Ensure `clang-format` is applied.
- Run `clang-tidy` to check for static analysis warnings or style violations.

### 5. Meson Build Standards
- **Reusing Build Folders**: Do **NOT** delete the build directory (`rm -rf build`) unless it is corrupted.
- **Reconfiguration**: Use `meson setup build --reconfigure [options]` to reconfigure the build directory while preserving existing object files.

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
*   **DC Bias**: Alternatively, add a tiny bias (e.g., `1.0e-18f`) to internal filter feedback states to prevent values from decaying to denormal thresholds.

### 3. Inter-Thread Communication
*   **Atomics**: Use atomic types (C11 `stdatomic.h` or C++ `std::atomic`) for thread-safe controls (bypass, levels, parameters).
*   **SPSC Ring Buffers**: Transfer blocks of data, spectral states, or events using thread-safe, lock-free ring buffers. Never resize ring buffers inside the audio callback.

### 4. Performance Optimizations
*   **Memory Alignment**: Align buffers to SIMD boundaries (16 or 32 bytes) to ensure compiler auto-vectorization is optimal.
*   **Reciprocal Multiplication**: Use multiplications by reciprocals ($1.0 / x$) instead of division operations inside hot processing loops.
