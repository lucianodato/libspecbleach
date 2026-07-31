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
