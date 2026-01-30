---
name: development_guidelines
description: "Code formatting, style, and safety rules for libspecbleach"
---

# libspecbleach Development Guidelines

## Code Formatting and Workflow

**CRITICAL**: Always run `ninja -C build format` after modifying any C/C++ files (.c, .h, .cpp, .hpp).

- CI/CD pipeline includes clang-format checks.
- Formatting must be done before committing or creating pull requests.

```bash
# Format all C/C++ files
ninja -C build format

# Format individual files
clang-format -i path/to/file.c
```

## Build System Rules

**Use Meson for everything** - no Makefiles or separate scripts.

Common commands:
```bash
meson setup build -Denable_examples=true
meson compile -C build
meson test -C build
```

## Code Style Guidelines

### C Code Standards
- LLVM clang-format style (2-space indentation, 80-char limit).
- **Mathematical Constants**: Get `M_PI` from `src/shared/configurations.h`.

### Safety and Robustness
- **Initialization & Cleanup**: Always implement matching `_initialize` and `_free`.
- **NULL Checks**: Every `_free` function **MUST** start with `if (!self) return;`.
- **Atomic Initialization**: Check every allocation in `_initialize`. If any fail, call `_free` and return `NULL`.
- **Allocation**: Use `calloc` for zero-initialized structures.

### FFT Buffer Management
- Uses FFTW Half-Complex (R2HC) format.
- Ensure buffers are sized to `fft_size`, not just `real_spectrum_size`, when mirroring is needed.

### Public API Design
- Keep internal configurations, macros, and enums out of `include/specbleach/`.
- Use primitive types (`int`, `float`, `bool`) in public parameter structures.

## Commit Workflow

1. Make changes.
2. Run `ninja -C build format`.
3. Test compilation: `meson compile -C build`.
4. Run tests: `meson test -C build`.
5. Commit with clear message.

## Quality Assurance

### Testing Requirements
- Unit tests must pass: `meson test -C build`.
- Audio quality verified with examples.
- No new compiler warnings.

### Code Review Checklist
- [ ] Code is properly formatted (`ninja -C build format`)
- [ ] Tests pass
- [ ] Documentation updated for API changes
- [ ] No new compiler warnings
- [ ] Follows established code style

## Performance Considerations
- SIMD acceleration through FFTW.
- Memory-aligned allocations for optimal performance.
- Optimized for low-latency real-time processing.
