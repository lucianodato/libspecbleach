---
name: meson_build
description: Best practices and commands for building, testing, and installing libspecbleach using Meson.
---

# Meson Build Skill for libspecbleach

This skill provides instructions on how to use the Meson build system for the `libspecbleach` project.

## Core Rules

1.  **Use a Dedicated Build Directory**: Always use a directory named `build` in the project root.
2.  **Avoid Redundant Build Folders**: Never create folders like `build_new`, `build2`, etc.
3.  **Reconfigure if Needed**: If the `build` directory exists but you need to change options, use:
    ```bash
    meson setup build --reconfigure [options]
    ```

## Common Workflows

### Standard Build
```bash
meson setup build
meson compile -C build
```

### Build with Tests (Recommended for development)
```bash
meson setup build -Denable_tests=true
meson compile -C build
meson test -C build
```

### Build for Code Coverage
```bash
meson setup build -Denable_tests=true -Db_coverage=true
meson compile -C build
meson test -C build
# Then use lcov or similar to generate reports from build/
```

### Debug Build with Sanitizers
```bash
meson setup build -Dbuildtype=debug -Denable_sanitizers=true -Dsanitize_address=true -Dsanitize_undefined=true
meson compile -C build
```

## Project Options

Refer to `meson_options.txt` for the full list. Key options:

| Option | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `enable_tests` | boolean | `false` | Enable building tests |
| `enable_examples` | boolean | `false` | Enable building example applications |
| `custom_warning_level` | integer | `2` | Compiler warning level (0-3) |
| `static_deps` | boolean | `false` | Link dependencies (fftw3, etc.) statically |
| `enable_sanitizers` | boolean | `false` | Enable sanitizers in debug builds |

## Tips
- Use `meson configure build` to see all current configuration options and their values.
- If you encounter build issues, try `ninja -C build clean` before recompiling, or as a last resort, `rm -rf build` and start over.
