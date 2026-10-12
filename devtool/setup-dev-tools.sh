#!/usr/bin/env bash
# Install the tools devtool/check.sh uses for tests and formatting:
# GoogleTest (unit tests in src/unit_test) and clang-format (incl. git-clang-format).
#
# Usage: devtool/setup-dev-tools.sh
#   Linux (apt): uses sudo.  macOS: uses Homebrew.
set -euo pipefail

info() { echo "==> $*"; }

case "$(uname -s)" in
    Linux)
        command -v apt-get >/dev/null || { echo "error: apt-get not found; install googletest and clang-format manually" >&2; exit 1; }
        if [[ $EUID -ne 0 ]] && ! sudo -n true 2>/dev/null && [[ ! -t 0 ]]; then
            echo "error: sudo needs a password but there is no terminal; run this script in your own terminal" >&2
            exit 1
        fi
        info "installing libgtest-dev clang-format (sudo)"
        sudo apt-get update
        sudo apt-get install -y libgtest-dev clang-format
        ;;
    Darwin)
        command -v brew >/dev/null || { echo "error: Homebrew not found" >&2; exit 1; }
        info "installing googletest clang-format (brew)"
        brew install googletest clang-format
        ;;
    *)
        echo "error: unsupported OS $(uname -s)" >&2; exit 1 ;;
esac

info "checking"
ok=1
if command -v clang-format >/dev/null; then
    clang-format --version
else
    echo "NG: clang-format not found"; ok=0
fi
if command -v git-clang-format >/dev/null; then
    echo "git-clang-format: $(command -v git-clang-format)"
else
    echo "NG: git-clang-format not found (devtool/check.sh needs it)"; ok=0
fi
gtest_header=""
for d in /usr/include /usr/local/include "$(brew --prefix 2>/dev/null || echo /nonexistent)/include"; do
    [[ -f "$d/gtest/gtest.h" ]] && gtest_header="$d/gtest/gtest.h" && break
done
if [[ -n "$gtest_header" ]]; then
    echo "gtest: $gtest_header"
else
    echo "NG: gtest/gtest.h not found"; ok=0
fi

((ok)) || exit 1
cat <<'EOF'

done. devtool/check.sh now builds src/unit_test and runs ctest
(the build directory is reconfigured automatically on the next run).
EOF
