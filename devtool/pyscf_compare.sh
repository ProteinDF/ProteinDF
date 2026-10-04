#!/usr/bin/env bash
# Wrapper script to execute pyscf_compare.py with the configured PySCF Python environment.
#
# Usage: devtool/pyscf_compare.sh [options passed to pyscf_compare.py]
#
# Configuration:
#   Python environment and test paths can be set via environment variables or in
#   $(git rev-parse --git-common-dir)/regress.conf:
#     PYSCF_PYTHON        Path to Python binary with PySCF installed
#     PROTEINDF_TEST_DIR  Path to ProteinDF_test repository

set -uo pipefail

top=$(git rev-parse --show-toplevel) || exit 1
common_dir=$(git rev-parse --git-common-dir) || exit 1

env_pyscf_python="${PYSCF_PYTHON:-}"
env_test_dir="${PROTEINDF_TEST_DIR:-}"

conf_file="$common_dir/regress.conf"
if [[ -f "$conf_file" ]]; then
    # shellcheck disable=SC1090
    source "$conf_file"
fi

python_bin="${env_pyscf_python:-${PYSCF_PYTHON:-}}"
export PROTEINDF_TEST_DIR="${env_test_dir:-${PROTEINDF_TEST_DIR:-}}"

if [[ -z "$python_bin" ]]; then
    if command -v python3 >/dev/null 2>&1 && python3 -c "import pyscf" >/dev/null 2>&1; then
        python_bin="python3"
    else
        echo "ERROR: PYSCF_PYTHON is not set, and python3 in PATH does not have pyscf installed." >&2
        echo "Please set PYSCF_PYTHON in environment or $conf_file." >&2
        exit 1
    fi
fi

if [[ ! -x "$python_bin" && "$python_bin" != "python3" ]]; then
    echo "ERROR: PYSCF_PYTHON ($python_bin) is not executable." >&2
    exit 1
fi

exec "$python_bin" "$top/devtool/pyscf_compare.py" "$@"
