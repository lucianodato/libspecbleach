# Contributing to libspecbleach

Thank you for your interest in contributing to libspecbleach! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

- Be respectful and considerate of others
- Welcome newcomers and help them learn
- Focus on constructive feedback

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone https://github.com/your-username/libspecbleach.git`
3. Create a branch for your changes: `git checkout -b feature/your-feature-name`
4. Make your changes
5. Test your changes
6. Submit a pull request

## Development Setup

### Prerequisites

- C compiler (GCC or Clang)
- Meson build system (0.60.0 or newer)
- Ninja build tool
- FFTW3 library (float version)
- libsndfile (optional, for examples)

### Building for Development

```bash
# Configure debug build with maximum warnings
meson setup build --buildtype=debug -Dwarning_level=3

# Compile
meson compile -C build

# Run tests (if available)
meson test -C build
```

## Code Style

### Formatting

The project uses `clang-format` for code formatting. Before submitting a pull request, ensure your code is properly formatted:

```bash
# Format all C source and header files
find . -type f \( -name "*.c" -o -name "*.h" \) ! -path "./build*" | xargs clang-format -i
```

The project uses the LLVM style with C-specific adjustments. See [.clang-format](.clang-format) for the complete configuration.

### General Guidelines

- Use 2 spaces for indentation (not tabs)
- Maximum line length: 80 characters
- Use meaningful variable and function names
- Add comments for complex algorithms or non-obvious code
- Follow the existing code style in the file you're editing

### C Standard

The project uses C11. Please ensure your code compiles with C11 standard enabled.

## Commit Messages

Write clear, descriptive commit messages:

- Use the imperative mood ("Add feature" not "Added feature")
- Keep the first line under 72 characters
- Provide additional context in the body if needed
- Reference issues/PRs when applicable

Example:
```
Add ARM64 optimization support

Enable NEON intrinsics for ARM64 architectures to improve
performance on Apple Silicon and other ARM-based systems.

Fixes #123
```

## Testing

- Test your changes on multiple platforms if possible
- Ensure the library builds without warnings (at warning_level=2)
- Test with both debug and release builds
- If adding new features, consider adding tests

## Pull Request Process

1. **Update Documentation**: If you're adding features or changing behavior, update the README and CHANGELOG.md
2. **Run Formatting**: Ensure code is formatted with clang-format
3. **Check Build**: Verify the project builds successfully
4. **Write Tests**: Add tests for new functionality if applicable
5. **Create PR**: Submit a pull request with a clear description

### PR Checklist

- [ ] Code follows the project's style guidelines
- [ ] Code is properly formatted with clang-format
- [ ] Build succeeds without warnings (at default warning level)
- [ ] Documentation is updated (README, CHANGELOG, etc.)
- [ ] Tests pass (if applicable)
- [ ] Commit messages are clear and descriptive

## Reporting Issues

When reporting bugs or requesting features:

- Use a clear, descriptive title
- Provide steps to reproduce (for bugs)
- Include information about your environment (OS, compiler, etc.)
- Include relevant code snippets or error messages
- Check if the issue already exists before creating a new one

## Questions?

If you have questions about contributing, feel free to:
- Open an issue for discussion
- Check existing issues and pull requests
- Review the codebase to understand patterns

Thank you for contributing to libspecbleach!

