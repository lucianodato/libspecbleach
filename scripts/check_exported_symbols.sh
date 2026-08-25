#!/usr/bin/env bash
#
# Verifies that a libspecbleach shared library exports ONLY public API
# symbols. Prints one [PASS] line when clean; prints offending symbols and
# exits non-zero otherwise.
#
# Usage: check_exported_symbols.sh <path-to-shared-library>

set -u

LIB="${1:?Usage: check_exported_symbols.sh <shared-library>}"

if [ ! -f "$LIB" ]; then
  echo "[FAIL] $LIB does not exist" >&2
  exit 1
fi

OS="$(uname -s)"
case "$OS" in
  Darwin)
    RAW_SYMBOLS="$(nm -gUj "$LIB" 2>/dev/null)"
    ;;
  Linux)
    RAW_SYMBOLS="$(nm -D --defined-only "$LIB" 2>/dev/null | awk '{print $NF}')"
    ;;
  *)
    echo "[SKIP] symbol export check unsupported on $OS" >&2
    exit 0
    ;;
esac

if [ -z "$RAW_SYMBOLS" ]; then
  echo "[FAIL] could not read dynamic symbols from $LIB" >&2
  exit 1
fi

# Public API surface: every SPECBLEACH_API-annotated function.
ALLOWED_RE='^_?(specbleach_denoiser_|specbleach_2d_|specbleach_stereo_|specbleach_transition_|specbleach_migrate_)'

# Standard linker-provided symbols that always appear in the export table.
STANDARD_RE='^__(mh_execute_header|chkstk|chkstk_ms|isa_available|x86_features)|^(__bss_start|_edata|_end|_init|_fini)$'

LEAKS="$(printf '%s\n' "$RAW_SYMBOLS" | grep -vE "$ALLOWED_RE" | grep -vE "$STANDARD_RE" || true)"

if [ -n "$LEAKS" ]; then
  echo "[FAIL] $LIB leaks non-public symbols:" >&2
  printf '%s\n' "$LEAKS" >&2
  exit 1
fi

echo "[PASS] $(basename "$LIB"): exported symbol surface clean"
