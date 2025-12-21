# Contributing to libspecbleach

We welcome contributions to this project. We aim to make contributing as straightforward and transparent as possible, including:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Improving documentation

## Development Workflow

We use **Meson** and **Ninja** for building.

1. **Clone the repo:**
   ```bash
   git clone https://github.com/lucianodato/libspecbleach.git
   cd libspecbleach
   ```

2. **Setup build:**
   ```bash
   meson setup build --buildtype=debug -Denable_sanitizers=true
   ```

3. **Compile and Test:**
   ```bash
   meson compile -C build
   ```

4. **Format Code:**
   We use `clang-format`. Please format your code before submitting:
   ```bash
   meson compile format -C build
   ```

## Development Build Options

For development purposes, additional features may be enabled:

```bash
# Enable tests and examples for development
meson setup build --buildtype=debug -Denable_tests=true -Denable_examples=true -Denable_sanitizers=true

# Run tests
meson test -C build

# Build examples
meson compile -C build
```

## Pull Requests

1. Fork the repo and create your branch from `main`.
2. If you've added code that should be tested, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code lints.
6. Submit the pull request.

## Coding Guidelines

- Follow the existing code style (enforced by clang-format)
- Use C17 standard
- Write clear, documented code
- Add tests for new functionality
- Update documentation for API changes

## License

By contributing, you agree that your contributions will be licensed under its LGPL-2.1 License.
