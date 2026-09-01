# Contributing to libspecbleach

We welcome contributions to this project. We aim to make contributing as straightforward and transparent as possible, including:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Improving documentation

## Development Workflow

We use **CMake** for building.

1. **Clone the repo:**
   ```bash
   git clone https://github.com/lucianodato/libspecbleach.git
   cd libspecbleach
   ```

2. **Setup build:**
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_TESTS=ON -DENABLE_EXAMPLES=ON -DENABLE_SANITIZERS=ON
   ```

3. **Compile and Test:**
   ```bash
   cmake --build build -j4
   cd build && ctest --output-on-failure
   ```

4. **Format Code:**
   We use `clang-format`. Please format your code before submitting:
   ```bash
   find src include -type f \( -name "*.c" -o -name "*.h" \) | xargs clang-format -i
   ```

## Pull Requests

1. Fork the repo and create your branch from `main`.
2. If you've added code that should be tested, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes (`ctest`).
5. Make sure your code lints and formats cleanly.
6. Submit the pull request.

## Coding Guidelines

- Follow the existing code style (enforced by `clang-format`)
- Use C17 standard
- Write clear, documented code
- **Testing Conventions:**
  - **Unit Tests:** Maintain a 1:1 mapping between source modules and unit test files (e.g., `my_module.c` -> `test_my_module.c`).
  - **Integration Tests:** Ensure filenames contain the word `integration` (e.g., `test_integration.c`).
- Add tests for new functionality
- Update documentation for API changes
- **Clean up unused includes:** Remove any `#include` directives that are not actually used in the source file.
- **License Header:** Ensure every source file begins with the standard LGPL license comment.

### Safety and Robustness

1. **Robust Freeing:** Every `_free` function MUST include a NULL check: `if (!self) return;`.
2. **Atomic Initialization:** Initialization functions should handle partial allocation failures gracefully. If a sub-allocation fails, all previously allocated resources within that function must be freed before returning `NULL`.
3. **FFT Buffer Sizing:** Always use the correct buffer size for FFT operations. Remember that for real-to-complex FFTs of size `N`, the real spectrum size is `N/2 + 1`.
4. **Verification:** Always run your changes through static analysis (`clang-tidy`, `cppcheck`) and with sanitizers enabled (`-DENABLE_SANITIZERS=ON`).

## License

By contributing, you agree that your contributions will be licensed under the LGPL-2.1 License.

## Adding a New Processor

libspecbleach is designed to grow STFT-based processors (denoisers, dereverberators, enhancers, ...). The architecture keeps that cheap — here is the verified recipe and the contracts that keep it true.

### What you get for free (do not reimplement)

- **STFT engine**: `stft_processor_run()` works like `qsort()` — it takes any
  `bool spectral_processing(SpectralProcessorHandle, float* spectrum)`
  callback and handles framing, windowing, overlap-add, and latency around
  it (`src/shared/stft/stft_processor.h`). Your processor is one function.
- **Noise profile storage & modes**: `SbProcessorCore` bundles the STFT
  engine with multi-mode noise profile state, including serialization
  accessors used by plugins (`specbleach_processor_core.{c,h}`).
- **Orchestration layer** (`extras/`, opt-in): stereo groups and profile
  migration are processor-agnostic at the group level. A new engine family
  gets them by adding entries, not logic.

### Recipe

1. **Spectral callback + internal params**: implement
   `bool my_processor_run(SpectralProcessorHandle self, float* spectrum)`
   plus a private parameters struct in `src/processors/my_processor/`.
   Sanitize/clamp user values at this boundary; never trust public inputs in
   the DSP path (see `sb_denoiser_params_sanitize`).
2. **Public header**: `include/specbleach_my_processor.h` following
   `specbleach_denoiser.h` as the template:
   - Own opaque handle type: `typedef struct specbleach_my_processor
     specbleach_my_processor;` — never share `void*` handles across
     families; distinct types are what make cross-engine misuse a compile
     error.
   - Own `SpecbleachMyProcessorParameters` struct; new structs never touch
     existing ones, so old binaries keep loading.
   - Every declaration annotated `SPECBLEACH_API`; document threading,
     ownership, and RT-safety per function.
   - `load_parameters(handle, const Params*, sizeof)` — always pointer+size.
   - Shared enums live in `specbleach_common.h` only if genuinely common;
     numeric values are contract.
 3. **Wrapper**: `src/processors/specbleach_my_processor.c` modeled on
    `specbleach_denoiser.c`: private state struct tagged with the public
    type name, own the STFT + noise profile lifetime, delegate to your
    processing callback, own any copied-in data (see `curve_bias_copy`)
    and free it in `_free`.
 4. **Extras wiring** (so integrators get stereo groups immediately):
    - Extend `extras/src/specbleach_stereo.c` to forward the new engine's
      per-channel calls when needed.
5. **Tests**: 1:1 unit test per module (`tests/test_my_processor.c`),
   registered in CMakeLists; add regression audio references only for
   intentional algorithm changes.
6. **Gates that must stay green**: full ctest (audio regression is
   bit-identical unless you changed DSP intentionally), clang-format,
   clang-tidy (workflow lints `src include extras`), cppcheck.

### When to refactor extras

Extras dispatches per engine family with explicit switches (~17 branches).
That is deliberate at two families. If you are adding the **third**, budget
time to collapse those switches into an ops table (function-pointer struct
per family) so family count grows without touching dispatch sites. Before
then, do not pre-abstract.
