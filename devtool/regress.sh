#!/usr/bin/env bash
# Regression test runner using ProteinDF_test.
# Builds and installs ProteinDF, prepares Python environment, and runs tests.
#
# Usage: devtool/regress.sh [options]
#   --suite <suite>       Test suite to run (default: serial_dev)
#   --entries <list>      Comma-separated list of test entries to run (default: all in suite)
#   --build-dir <dir>     CMake build directory (default: build-regress)
#   --install-dir <dir>   Install prefix (default: <build-dir>/install)
#   --no-build            Skip build and install (use existing binaries)
#   --format <fmt>        Reference format: auto|h5 (default: auto)
#                           auto: same as h5 (reference.h5)
#                           h5:   require reference.h5 (SKIP if missing)
#   --energy-threshold <val>
#                         Total energy comparison threshold passed to
#                         pdf-test-h5.py (default: derived from
#                         fl_Userinput's threshold-energy clamped to max(val, 1e-8))
#   --keep-work-dir <dir> Copy each entry's work directory (fl_Userinput,
#                           pdfresults.h5, logs, ...) into <dir>/<entry>/
#                           after the run, instead of discarding it with the
#                           temporary workspace. Useful for inspecting a run
#                           or harvesting pdfresults.h5 as a new reference.h5.
#   --full                Print every entry (default: PASS rows are hidden from the
#                           progress lines and the summary table)
#   -j, --jobs <N>        Parallel build jobs (default: nproc)
#   -h, --help            Show this help message
#
# Configuration:
#   Repository paths must be configured via environment variables or in
#   $(git rev-parse --git-common-dir)/regress.conf (KEY="value" format).
#   Environment variables take precedence over regress.conf.
#
#   PROTEINDF_TEST_DIR    Path to ProteinDF_test (required)
#   PROTEINDF_PYTOOLS_DIR Path to ProteinDF_pytools (required)
#   PROTEINDF_BRIDGE_DIR  Path to ProteinDF_bridge (required)
#   PDF_CMAKE_ARGS        Extra options passed to cmake configure (optional)
#   REGRESS_VENV_DIR      Path to shared python venv (default: <git-common-dir>/regress-venv)
#   REGRESS_LOGS_DIR      Path to log directory (default: <git-common-dir>/regress-logs)
#

set -uo pipefail

show_help() {
    sed -n '2,/^set -uo pipefail/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'
    exit 0
}

suite="serial_dev"
entries_arg=""
build_dir="build-regress"
install_dir=""
do_build=1
jobs=$(nproc 2>/dev/null || echo 4)
ref_format="auto"
energy_threshold_opt=""
keep_work_dir=""
full=0

while (($#)); do
    case "$1" in
        --suite)
            suite="$2"; shift 2 ;;
        --entries)
            entries_arg="$2"; shift 2 ;;
        --build-dir)
            build_dir="$2"; shift 2 ;;
        --install-dir)
            install_dir="$2"; shift 2 ;;
        --no-build)
            do_build=0; shift ;;
        --format)
            ref_format="$2"; shift 2 ;;
        --energy-threshold)
            energy_threshold_opt="$2"; shift 2 ;;
        --keep-work-dir)
            keep_work_dir="$2"; shift 2 ;;
        --full)
            full=1; shift ;;
        -j|--jobs)
            jobs="$2"; shift 2 ;;
        -h|--help)
            show_help ;;
        *)
            echo "Unknown option: $1" >&2
            show_help ;;
    esac
done

case "$ref_format" in
    auto|h5) ;;
    *)
        echo "ERROR: --format must be one of: auto, h5 (SQLite db is no longer supported)" >&2
        exit 1 ;;
esac

top=$(git rev-parse --show-toplevel) || exit 1
cd "$top"
common_dir=$(git rev-parse --git-common-dir) || exit 1

if [[ -z "$install_dir" ]]; then
    # Make install_dir absolute under build_dir
    if [[ "$build_dir" = /* ]]; then
        install_dir="$build_dir/install"
    else
        install_dir="$top/$build_dir/install"
    fi
fi

# Save environment variables to give them precedence over regress.conf
env_test_dir="${PROTEINDF_TEST_DIR:-}"
env_pytools_dir="${PROTEINDF_PYTOOLS_DIR:-}"
env_bridge_dir="${PROTEINDF_BRIDGE_DIR:-}"
env_venv_dir="${REGRESS_VENV_DIR:-}"
env_logs_dir="${REGRESS_LOGS_DIR:-}"

conf_file="$common_dir/regress.conf"
if [[ -f "$conf_file" ]]; then
    # shellcheck disable=SC1090
    source "$conf_file"
fi

test_repo="${env_test_dir:-${PROTEINDF_TEST_DIR:-}}"
pytools_repo="${env_pytools_dir:-${PROTEINDF_PYTOOLS_DIR:-}}"
bridge_repo="${env_bridge_dir:-${PROTEINDF_BRIDGE_DIR:-}}"
venv_dir="${env_venv_dir:-${REGRESS_VENV_DIR:-$common_dir/regress-venv}}"
logs_dir="${env_logs_dir:-${REGRESS_LOGS_DIR:-$common_dir/regress-logs}}"

missing_vars=()
[[ -z "$test_repo" ]] && missing_vars+=("PROTEINDF_TEST_DIR")
[[ -z "$pytools_repo" ]] && missing_vars+=("PROTEINDF_PYTOOLS_DIR")
[[ -z "$bridge_repo" ]] && missing_vars+=("PROTEINDF_BRIDGE_DIR")

if ((${#missing_vars[@]} > 0)); then
    echo "ERROR: Missing required configuration variables: ${missing_vars[*]}" >&2
    echo >&2
    echo "Please set them via environment variables or in $conf_file." >&2
    echo "Example $conf_file:" >&2
    echo "  PROTEINDF_TEST_DIR=\"/path/to/ProteinDF_test\"" >&2
    echo "  PROTEINDF_PYTOOLS_DIR=\"/path/to/ProteinDF_pytools\"" >&2
    echo "  PROTEINDF_BRIDGE_DIR=\"/path/to/ProteinDF_bridge\"" >&2
    exit 1
fi

mkdir -p "$logs_dir"

if [[ -n "$keep_work_dir" ]]; then
    mkdir -p "$keep_work_dir"
    keep_work_dir=$(cd "$keep_work_dir" && pwd)
fi

echo "=== ProteinDF Regression Test ==="
echo "Worktree:    $top"
echo "Branch:      $(git branch --show-current) ($(git rev-parse --short HEAD))"
echo "Suite:       $suite"
echo "Build dir:   $build_dir"
echo "Install dir: $install_dir"
echo "Venv:        $venv_dir"
echo "Logs dir:    $logs_dir"
echo

# -----------------------------------------------------------------------------
# 1. Python environment (uv venv)
# -----------------------------------------------------------------------------
echo "--- 1. Setting up Python environment ---"
t_py_start=$(date +%s)

venv_ready=0
if [[ -x "$venv_dir/bin/python" && -x "$venv_dir/bin/pdf-archive-h5.py" && -x "$venv_dir/bin/pdf-test-h5.py" ]]; then
    if "$venv_dir/bin/python" -c "import proteindf_tools as p; import proteindf_bridge as b; import numpy as np; import h5py" >/dev/null 2>&1; then
        # a venv made by an older version of this script carries sitecustomize.py patches
        sp_dir=$("$venv_dir/bin/python" -c "import site; print(site.getsitepackages()[0])")
        [[ -e "$sp_dir/sitecustomize.py" ]] || venv_ready=1
    fi
fi

if ((venv_ready)); then
    echo "Using existing Python environment at $venv_dir"
else
    echo "Creating virtual environment at $venv_dir (Python 3.12)..."
    if ! command -v uv >/dev/null 2>&1; then
        echo "ERROR: 'uv' command not found." >&2
        exit 1
    fi

    if [[ ! -d "$bridge_repo" || ! -d "$pytools_repo" ]]; then
        echo "ERROR: Source repositories not found:" >&2
        echo "  ProteinDF_bridge: $bridge_repo" >&2
        echo "  ProteinDF_pytools: $pytools_repo" >&2
        exit 1
    fi

    uv venv --clear "$venv_dir" --python 3.12 > "$logs_dir/venv-create.log" 2>&1 || {
        echo "Failed to create uv venv. See $logs_dir/venv-create.log" >&2
        exit 1
    }

    echo "Installing ProteinDF_bridge and dependencies..."
    uv pip install --python "$venv_dir/bin/python" "$bridge_repo" > "$logs_dir/venv-bridge-install.log" 2>&1 || {
        echo "Failed to install ProteinDF_bridge. See $logs_dir/venv-bridge-install.log" >&2
        exit 1
    }

    echo "Installing ProteinDF_pytools..."
    uv pip install --python "$venv_dir/bin/python" "$pytools_repo" > "$logs_dir/venv-pytools-install.log" 2>&1 || {
        echo "Failed to install ProteinDF_pytools. See $logs_dir/venv-pytools-install.log" >&2
        exit 1
    }

    echo "Environment setup complete."
fi

t_py_end=$(date +%s)
echo "Python environment checked/prepared in $((t_py_end - t_py_start))s."
echo

# Ensure venv tools (such as yml2mpac) are available during build
export PATH="$venv_dir/bin:$PATH"

# -----------------------------------------------------------------------------
# 2. Build and install ProteinDF
# -----------------------------------------------------------------------------
if ((do_build)); then
    echo "--- 2. Building and installing ProteinDF ---"
    t_build_start=$(date +%s)
    mkdir -p "$build_dir"

    gen=()
    if [[ ! -f "$build_dir/CMakeCache.txt" ]]; then
        ninja --version >/dev/null 2>&1 && gen=(-G Ninja)
    fi

    # Ensure -O2 optimization flag is present
    cmake_cxx_flags="-O2 -DNDEBUG"

    echo "Configuring CMake..."
    # shellcheck disable=SC2086
    if ! cmake -S . -B "$build_dir" "${gen[@]}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS="$cmake_cxx_flags" \
        -DCMAKE_INSTALL_PREFIX="$install_dir" \
        ${PDF_CMAKE_ARGS:-} > "$build_dir/regress-configure.log" 2>&1; then
        echo "CMake configure failed. See $build_dir/regress-configure.log" >&2
        tail -n 30 "$build_dir/regress-configure.log" >&2
        exit 1
    fi

    echo "Building (jobs: $jobs)..."
    if ! cmake --build "$build_dir" -j "$jobs" > "$build_dir/regress-build.log" 2>&1; then
        echo "Build failed. See $build_dir/regress-build.log" >&2
        grep -E "error|Error [0-9]" "$build_dir/regress-build.log" | sort -u | head -30 >&2
        exit 1
    fi

    echo "Installing to $install_dir..."
    if ! cmake --install "$build_dir" > "$build_dir/regress-install.log" 2>&1; then
        echo "Install failed. See $build_dir/regress-install.log" >&2
        tail -n 30 "$build_dir/regress-install.log" >&2
        exit 1
    fi

    t_build_end=$(date +%s)
    echo "Build & install completed in $((t_build_end - t_build_start))s."
    echo
else
    echo "--- 2. Skipping build (--no-build specified) ---"
    echo
fi

# Verify installation
if [[ ! -x "$install_dir/bin/PDF.x" ]]; then
    echo "ERROR: $install_dir/bin/PDF.x not found or not executable." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# 3. Resolve test entries
# -----------------------------------------------------------------------------
suite_dir="$test_repo/$suite"
if [[ ! -d "$suite_dir" ]]; then
    echo "ERROR: Suite directory not found: $suite_dir" >&2
    exit 1
fi

declare -a entries=()
if [[ -n "$entries_arg" ]]; then
    IFS=',' read -r -a raw_entries <<< "$entries_arg"
    for e in "${raw_entries[@]}"; do
        # trim whitespace
        e_trimmed=$(echo "$e" | xargs)
        if [[ -n "$e_trimmed" ]]; then
            entries+=("$e_trimmed")
        fi
    done
else
    # Find all directories containing fl_Userinput
    while IFS= read -r dir; do
        entry_name=$(basename "$dir")
        entries+=("$entry_name")
    done < <(find "$suite_dir" -mindepth 1 -maxdepth 1 -type d -exec test -f '{}/fl_Userinput' ';' -print | sort)
fi

if ((${#entries[@]} == 0)); then
    echo "ERROR: No test entries found in $suite_dir" >&2
    exit 1
fi

echo "--- 3. Running Regression Tests (${#entries[@]} entries in '$suite') ---"

# Set up test execution environment
export PDF_HOME="$install_dir"
export PATH="$venv_dir/bin:$install_dir/bin:$PATH"

# Temporary workspace for running tests
tmp_workspace=$(mktemp -d -t pdf-regress-XXXXXX)
cleanup() {
    rm -rf "$tmp_workspace"
}
trap cleanup EXIT INT TERM

declare -a result_status=()
declare -a result_time=()
declare -a result_fmt=()
declare -a result_iter=()
declare -a result_calc_te=()
declare -a result_std_te=()
declare -a result_diff_te=()
declare -a result_thr=()
declare -a result_notes=()

total_tests=${#entries[@]}
overall_pass=1
t_test_start=$(date +%s)

for idx in "${!entries[@]}"; do
    entry="${entries[$idx]}"
    src_entry_dir="$suite_dir/$entry"
    num=$((idx + 1))

    if [[ ! -d "$src_entry_dir" ]]; then
        printf "[%2d/%2d] %-22s %-4s (0s) [entry directory does not exist]\n" "$num" "$total_tests" "$entry" "FAIL"
        result_status+=("NOT_FOUND")
        result_time+=("0s")
        result_fmt+=("-")
        result_iter+=("-")
        result_calc_te+=("N/A")
        result_std_te+=("N/A")
        result_diff_te+=("N/A")
        result_thr+=("-")
        result_notes+=("entry directory does not exist")
        overall_pass=0
        continue
    fi

    entry_work_dir="$tmp_workspace/$entry"
    rm -rf "$entry_work_dir"
    cp -r "$src_entry_dir" "$entry_work_dir"

    log_file="$logs_dir/${entry}.log"
    : > "$log_file"

    t0=$(date +%s%N 2>/dev/null || date +%s)

    # Step A: Setup
    (
        cd "$entry_work_dir"
        pdf-clean.sh all >> "$log_file" 2>&1
        pdf-setup.sh >> "$log_file" 2>&1
        if [[ -f ./pre_pdf.sh ]]; then
            ./pre_pdf.sh >> "$log_file" 2>&1
        fi
    )

    # Step B: Run ProteinDF (serial)
    pdf_exit=0
    (
        cd "$entry_work_dir"
        PDF.x >> "$log_file" 2>&1
    ) || pdf_exit=$?

    # Reference: reference.h5 (SKIP when missing)
    ref_h5="$entry_work_dir/reference.h5"
    ref_format_used=""
    if [[ -f "$ref_h5" ]]; then
        ref_format_used="h5"
    fi

    # Step C: Archive results (HDF5)
    archive_exit=0
    calc_file="$entry_work_dir/pdfresults.h5"
    if [[ $pdf_exit -eq 0 ]]; then
        (
            cd "$entry_work_dir"
            pdf-archive-h5.py >> "$log_file" 2>&1
        ) || archive_exit=$?
    fi

    # Step D: Check results and compare with the reference file
    test_exit=0
    ref_file=""
    [[ -n "$ref_format_used" ]] && ref_file="$ref_h5"

    calc_te="N/A"
    std_te="N/A"
    diff_te="N/A"
    iter_str="-"
    thr_str="-"
    entry_thr=""
    note=""

    # Extract Total Energy and iterations via proteindf_tools
    # (PdfParam_H5), same quantity pdf-test-h5.py compares.
    if [[ -f "$calc_file" || -f "$ref_file" ]]; then
        c_file_arg=""
        s_file_arg=""
        [[ -f "$calc_file" ]] && c_file_arg="$calc_file"
        [[ -f "$ref_file" ]] && s_file_arg="$ref_file"

        py_res=$("$venv_dir/bin/python" - "$c_file_arg" "$s_file_arg" << 'PYEOF'
import sys
import proteindf_tools as p

calc_path = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else None
std_path = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] else None

c_te, s_te, d_te, i_str = "N/A", "N/A", "N/A", "-"
c_iter, s_iter = None, None


def load(path):
    obj = p.PdfParam_H5()
    obj.open(path)
    itr = obj.iterations
    return itr, obj.get_total_energy(itr)


if calc_path:
    try:
        c_iter, val_c = load(calc_path)
        if val_c is not None:
            c_te = f"{val_c:.10f}"
    except Exception:
        pass

if std_path:
    try:
        s_iter, val_s = load(std_path)
        if val_s is not None:
            s_te = f"{val_s:.10f}"
    except Exception:
        pass

if c_iter is not None and s_iter is not None:
    if c_iter == s_iter:
        i_str = str(c_iter)
    else:
        i_str = f"{c_iter}/{s_iter}"
elif c_iter is not None:
    i_str = str(c_iter)
elif s_iter is not None:
    i_str = f"-/{s_iter}"

if c_te != "N/A" and s_te != "N/A":
    try:
        diff = float(c_te) - float(s_te)
        d_te = f"{diff:+.2e}"
    except Exception:
        pass

print(f"{c_te}\t{s_te}\t{d_te}\t{i_str}")
PYEOF
)
        IFS=$'\t' read -r calc_te std_te diff_te iter_str <<< "$py_res"
    fi

    if [[ $pdf_exit -ne 0 ]]; then
        test_status="FAIL"
        note="PDF.x exited with $pdf_exit"
        overall_pass=0
    elif [[ $archive_exit -ne 0 || ! -f "$calc_file" ]]; then
        test_status="FAIL"
        note="pdf-archive-h5.py failed"
        overall_pass=0
    elif [[ -z "$ref_format_used" ]]; then
        test_status="SKIP"
        note="no reference"
    else
        entry_thr=""
        if [[ -n "$energy_threshold_opt" ]]; then
            entry_thr="$energy_threshold_opt"
            thr_str="$energy_threshold_opt"
        elif [[ "$ref_format_used" = "h5" ]]; then
            userinput_file="$entry_work_dir/fl_Userinput"
            if [[ -f "$userinput_file" ]]; then
                parsed_thr=$("$venv_dir/bin/python" - "$userinput_file" << 'PYEOF'
import sys
import re

path = sys.argv[1]
pattern = re.compile(r"^\s*convergence/threshold-energy\s*=\s*([^\s#;]+)", re.IGNORECASE)
val = None
try:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = pattern.match(line)
            if m:
                val = m.group(1)
                break
    if val is not None:
        clamped = max(float(val), 1e-8)
        m_str, exp_str = f"{clamped:.8e}".split("e")
        m_str = m_str.rstrip("0").rstrip(".")
        exp = int(exp_str)
        print(f"{m_str}e{exp}")
except Exception:
    pass
PYEOF
)
                if [[ -n "$parsed_thr" ]]; then
                    entry_thr="$parsed_thr"
                    thr_str="$parsed_thr"
                else
                    thr_str="default"
                fi
            else
                thr_str="default"
            fi
        else
            thr_str="default"
        fi

        (
            cd "$entry_work_dir"
            thr_args=()
            [[ -n "$entry_thr" ]] && thr_args=(--energy-threshold "$entry_thr")
            echo "energy threshold: $thr_str" >> "$log_file"
            pdf-test-h5.py -v "${thr_args[@]}" "$calc_file" "$ref_file" >> "$log_file" 2>&1
        ) || test_exit=$?

        if [[ $test_exit -eq 0 ]]; then
            test_status="PASS"
        else
            test_status="FAIL"
            # Extract mismatch messages from log_file
            mismatch_line=$(grep -E "test: .* != " "$log_file" | head -n 1)
            if [[ -n "$mismatch_line" ]]; then
                note=$(echo "$mismatch_line" | sed 's/.*test: //')
            else
                note="discrepancy detected"
            fi
            overall_pass=0
        fi
    fi

    if [[ -n "$keep_work_dir" ]]; then
        rm -rf "${keep_work_dir:?}/${entry:?}"
        mkdir -p "$keep_work_dir/$entry"
        cp -a "$entry_work_dir/." "$keep_work_dir/$entry/"
    fi

    t1=$(date +%s%N 2>/dev/null || date +%s)
    # calculate elapsed seconds
    if [[ ${#t0} -gt 10 ]]; then
        elapsed_sec=$(awk "BEGIN {printf \"%.2fs\", ($t1 - $t0) / 1000000000}")
    else
        elapsed_sec="$((t1 - t0))s"
    fi

    fmt_str="${ref_format_used:--}"

    if ((full)) || [[ "$test_status" != "PASS" ]]; then
        printf "[%2d/%2d] %-22s %-4s (%s)%s\n" "$num" "$total_tests" "$entry" "$test_status" "$elapsed_sec" "${note:+ [$note]}"
    fi

    result_status+=("$test_status")
    result_time+=("$elapsed_sec")
    result_fmt+=("$fmt_str")
    result_iter+=("$iter_str")
    result_calc_te+=("$calc_te")
    result_std_te+=("$std_te")
    result_diff_te+=("$diff_te")
    result_thr+=("$thr_str")
    result_notes+=("$note")
done

t_test_end=$(date +%s)
total_elapsed=$((t_test_end - t_test_start))

# -----------------------------------------------------------------------------
# 4. Summary Table
# -----------------------------------------------------------------------------
echo
echo "=== Regression Test Summary ==="
echo "Suite: $suite | Total: $total_tests | Time: ${total_elapsed}s"
echo

printf "%-4s %-22s %-8s %-4s %-6s %-18s %-18s %-11s %-8s %s\n" "No." "Entry" "Status" "Fmt" "Iter" "Calc TE (a.u.)" "Std TE (a.u.)" "Diff" "Thr" "Note"
printf '%0.s-' {1..120}
echo

for idx in "${!entries[@]}"; do
    entry="${entries[$idx]}"
    num=$((idx + 1))
    status="${result_status[$idx]}"
    time_str="${result_time[$idx]}"
    fmt="${result_fmt[$idx]}"
    istr="${result_iter[$idx]}"
    cte="${result_calc_te[$idx]}"
    ste="${result_std_te[$idx]}"
    dte="${result_diff_te[$idx]}"
    thr="${result_thr[$idx]}"
    nt="${result_notes[$idx]}"

    # compact view: PASS rows are folded into the count below (--full shows them)
    if ((!full)) && [[ "$status" == "PASS" ]]; then
        continue
    fi
    printf "%2d.  %-22s %-8s %-4s %-6s %-18s %-18s %-11s %-8s %s\n" \
        "$num" "$entry" "$status" "$fmt" "$istr" "$cte" "$ste" "$dte" "$thr" "$nt"
done

printf '%0.s-' {1..120}
echo
if ((!full)); then
    echo "(PASS rows are hidden; use --full to list all $total_tests entries)"
fi

# Print detailed diffs for failed tests
failed_count=0
skipped_count=0
passed_count=0
for idx in "${!entries[@]}"; do
    if [[ "${result_status[$idx]}" = "PASS" ]]; then
        ((passed_count++))
    elif [[ "${result_status[$idx]}" = "SKIP" ]]; then
        ((skipped_count++))
    else
        ((failed_count++))
    fi
done

if ((failed_count > 0)); then
    echo
    echo "=== Failed Tests Details ($failed_count test(s)) ==="
    for idx in "${!entries[@]}"; do
        if [[ "${result_status[$idx]}" != "PASS" && "${result_status[$idx]}" != "SKIP" ]]; then
            entry="${entries[$idx]}"
            echo "--- [$entry] log excerpt ($logs_dir/${entry}.log) ---"
            grep -E "(ERROR|CRITICAL|failed|Assertion|not consistent|!=)" "$logs_dir/${entry}.log" | head -n 20 || true
            echo
        fi
    done
fi

echo "All test logs saved to: $logs_dir/"

if ((failed_count == 0)); then
    if ((skipped_count > 0)); then
        echo "Result: ALL TESTS PASSED ($passed_count passed, $skipped_count skipped, ${total_elapsed}s)"
    else
        echo "Result: ALL TESTS PASSED ($passed_count passed, ${total_elapsed}s)"
    fi
    exit 0
else
    echo "Result: $failed_count TEST(S) FAILED (${total_elapsed}s)"
    exit 1
fi
