---
name: implement_dsp_processor
description: Instructions for implementing a new audio processor in libspecbleach.
---

# Implement DSP Processor Skill

This skill guides you through adding a new audio processing algorithm to `libspecbleach`.

## 1. Core Implementation
Location: `src/processors/<name>/`

1.  **Header (`.h`)**:
    *   Define the configuration struct (e.g., `<Name>Config`).
    *   Define the state struct (opaque if possible).
    *   Declare `initialize`, `process`, `reset`, and `free`.
2.  **Source (`.c`)**:
    *   Implement the lifecycle and processing logic.
    *   **Optimization**: Use `configurations.h` constants. Guard SIMD code with `#ifdef __ARM_NEON` (etc).
    *   **Memory**: Ensure `free` handles NULL.

## 2. Public Wrapper
Location: `src/processors/specbleach_<name>.c`

1.  Create a standard wrapper matching `specbleach_denoiser` style.
2.  Use `SpectralBleachHandle` (void*).
3.  Implement `specbleach_<name>_initialize`, `_process`, `_free`.
4.  Implement `specbleach_<name>_get_latency` (critical for DAW compensation).

## 3. Build System
Location: `src/meson.build`

*   Add the new source files to the `sources` list.

## 4. Testing
Location: `tests/`

1.  **Unit Test**: Create `tests/test_<name>.c` covering the core logic.
2.  **Integration Test**: Add a case to `tests/test_integration.c` using the wrapper.
3.  **Meson**: Add the new test executable to `tests/meson.build`.
4.  **Performance**: Always verify performance with `--buildtype=release`. Debug builds disable critical SIMD/inlining optimizations.
