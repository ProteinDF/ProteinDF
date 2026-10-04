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
#   -j, --jobs <N>        Parallel build jobs (default: nproc)
#   -h, --help            Show this help message
#
# Environment variables:
#   PDF_CMAKE_ARGS        Extra options passed to cmake configure
#   PROTEINDF_TEST_DIR    Path to ProteinDF_test (default: ~/work/dev/pdf-dev/ProteinDF_test)
#   PROTEINDF_PYTOOLS_DIR Path to ProteinDF_pytools (default: ~/work/dev/pdf-dev/ProteinDF_pytools)
#   PROTEINDF_BRIDGE_DIR  Path to ProteinDF_bridge (default: ~/work/dev/pdf-dev/ProteinDF_bridge)
#   REGRESS_VENV_DIR      Path to shared python venv (default: $(git rev-parse --git-common-dir)/regress-venv)
#   REGRESS_LOGS_DIR      Path to log directory (default: $(git rev-parse --git-common-dir)/regress-logs)
#
# Python Environment Compatibility Patches (sitecustomize.py in venv):
#   - PdfArchive alias:
#       proteindf_tools exposes PdfArchive_Sqlite3, but scripts reference
#       proteindf_tools.PdfArchive. Added alias to prevent AttributeError.
#   - Force dict conversion:
#       PdfParamObject.set_by_raw_data expects force items to be indexed by integer.
#       Some SQLite DBs store force records as dicts, causing KeyError: 0.
#       Patched to convert dict {0:x, 1:y, 2:z} to list [x, y, z].
#   - numpy<2:
#       proteindf_bridge calls numpy.ndarray.tostring() removed in NumPy 2.0.
#       Pinned to numpy<2 during installation.

set -uo pipefail

show_help() {
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
}

suite="serial_dev"
entries_arg=""
build_dir="build-regress"
install_dir=""
do_build=1
jobs=$(nproc 2>/dev/null || echo 4)

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
        -j|--jobs)
            jobs="$2"; shift 2 ;;
        -h|--help)
            show_help ;;
        *)
            echo "Unknown option: $1" >&2
            show_help ;;
    esac
done

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

test_repo="${PROTEINDF_TEST_DIR:-$HOME/work/dev/pdf-dev/ProteinDF_test}"
pytools_repo="${PROTEINDF_PYTOOLS_DIR:-$HOME/work/dev/pdf-dev/ProteinDF_pytools}"
bridge_repo="${PROTEINDF_BRIDGE_DIR:-$HOME/work/dev/pdf-dev/ProteinDF_bridge}"
venv_dir="${REGRESS_VENV_DIR:-$common_dir/regress-venv}"
logs_dir="${REGRESS_LOGS_DIR:-$common_dir/regress-logs}"

mkdir -p "$logs_dir"

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
# 1. Build and install ProteinDF
# -----------------------------------------------------------------------------
if ((do_build)); then
    echo "--- 1. Building and installing ProteinDF ---"
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
    echo "--- 1. Skipping build (--no-build specified) ---"
    echo
fi

# Verify installation
if [[ ! -x "$install_dir/bin/PDF.x" ]]; then
    echo "ERROR: $install_dir/bin/PDF.x not found or not executable." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# 2. Python environment (uv venv)
# -----------------------------------------------------------------------------
echo "--- 2. Setting up Python environment ---"
t_py_start=$(date +%s)

venv_ready=0
if [[ -x "$venv_dir/bin/python" && -x "$venv_dir/bin/pdf-archive.py" && -x "$venv_dir/bin/pdf-test.py" ]]; then
    if "$venv_dir/bin/python" -c "import proteindf_tools as p; import proteindf_bridge as b; import numpy as np; assert hasattr(p, 'PdfArchive'); assert hasattr(np.ndarray, 'tostring')" >/dev/null 2>&1; then
        venv_ready=1
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

    uv venv "$venv_dir" --python 3.12 > "$logs_dir/venv-create.log" 2>&1 || {
        echo "Failed to create uv venv. See $logs_dir/venv-create.log" >&2
        exit 1
    }

    echo "Installing ProteinDF_bridge and dependencies..."
    # numpy<2 is required for numpy.ndarray.tostring() compatibility used by proteindf_bridge
    uv pip install --python "$venv_dir/bin/python" "numpy<2" "$bridge_repo" > "$logs_dir/venv-bridge-install.log" 2>&1 || {
        echo "Failed to install ProteinDF_bridge. See $logs_dir/venv-bridge-install.log" >&2
        exit 1
    }

    echo "Installing ProteinDF_pytools..."
    uv pip install --python "$venv_dir/bin/python" "$pytools_repo" > "$logs_dir/venv-pytools-install.log" 2>&1 || {
        echo "Failed to install ProteinDF_pytools. See $logs_dir/venv-pytools-install.log" >&2
        exit 1
    }

    # Workaround: sitecustomize.py for Python 3.12 / pytools compatibility without touching repos
    sp_dir=$("$venv_dir/bin/python" -c "import site; print(site.getsitepackages()[0])")
    cat << 'EOF' > "$sp_dir/sitecustomize.py"
import proteindf_tools
import proteindf_tools.pdfarchive_sqlite3
import proteindf_tools.pdfparam_object

# Patch PdfArchive alias
if not hasattr(proteindf_tools, 'PdfArchive'):
    proteindf_tools.PdfArchive = proteindf_tools.PdfArchive_Sqlite3
if not hasattr(proteindf_tools.pdfarchive_sqlite3, 'PdfArchive'):
    proteindf_tools.pdfarchive_sqlite3.PdfArchive = proteindf_tools.pdfarchive_sqlite3.PdfArchive_Sqlite3

# Patch force dict KeyError: 0 in set_by_raw_data
_orig_set_by_raw_data = proteindf_tools.pdfparam_object.PdfParamObject.set_by_raw_data
def _patched_set_by_raw_data(self, odict):
    if "force" in odict:
        force_dat = odict["force"]
        if isinstance(force_dat, list):
            new_force = []
            for item in force_dat:
                if isinstance(item, dict):
                    new_force.append([item.get(0, item.get('0', 0.0)),
                                      item.get(1, item.get('1', 0.0)),
                                      item.get(2, item.get('2', 0.0))])
                else:
                    new_force.append(item)
            odict["force"] = new_force
    return _orig_set_by_raw_data(self, odict)
proteindf_tools.pdfparam_object.PdfParamObject.set_by_raw_data = _patched_set_by_raw_data
EOF

    echo "Environment setup complete."
fi

t_py_end=$(date +%s)
echo "Python environment checked/prepared in $((t_py_end - t_py_start))s."
echo

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
declare -a result_iter=()
declare -a result_calc_te=()
declare -a result_std_te=()
declare -a result_diff_te=()
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
        result_iter+=("-")
        result_calc_te+=("N/A")
        result_std_te+=("N/A")
        result_diff_te+=("N/A")
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

    # Step C: Archive results to SQLite db
    archive_exit=0
    if [[ $pdf_exit -eq 0 ]]; then
        (
            cd "$entry_work_dir"
            pdf-archive.py >> "$log_file" 2>&1
        ) || archive_exit=$?
    fi

    # Step D: Check results and compare with pdfresults_std.db
    test_exit=0
    std_db="$entry_work_dir/pdfresults_std.db"
    calc_db="$entry_work_dir/pdfresults.db"

    calc_te="N/A"
    std_te="N/A"
    diff_te="N/A"
    iter_str="-"
    note=""

    # Extract Total Energy and iterations via proteindf_tools.PdfArchive
    if [[ -f "$calc_db" || -f "$std_db" ]]; then
        c_db_arg=""
        s_db_arg=""
        [[ -f "$calc_db" ]] && c_db_arg="$calc_db"
        [[ -f "$std_db" ]] && s_db_arg="$std_db"

        py_res=$("$venv_dir/bin/python" - "$c_db_arg" "$s_db_arg" << 'PYEOF'
import sys
import proteindf_tools as p

calc_path = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else None
std_path = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] else None

c_te, s_te, d_te, i_str = "N/A", "N/A", "N/A", "-"
c_iter, s_iter = None, None

if calc_path:
    try:
        arc_c = p.PdfArchive(calc_path)
        c_iter = arc_c.iterations
        val_c = arc_c.get_total_energy(c_iter)
        if val_c is not None:
            c_te = f"{val_c:.10f}"
    except Exception:
        pass

if std_path:
    try:
        arc_s = p.PdfArchive(std_path)
        s_iter = arc_s.iterations
        val_s = arc_s.get_total_energy(s_iter)
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
    elif [[ $archive_exit -ne 0 || ! -f "$calc_db" ]]; then
        test_status="FAIL"
        note="pdf-archive.py failed"
        overall_pass=0
    elif [[ ! -f "$std_db" ]]; then
        test_status="SKIP"
        note="no reference"
    else
        (
            cd "$entry_work_dir"
            pdf-test.py "$calc_db" "$std_db" >> "$log_file" 2>&1
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

    t1=$(date +%s%N 2>/dev/null || date +%s)
    # calculate elapsed seconds
    if [[ ${#t0} -gt 10 ]]; then
        elapsed_sec=$(awk "BEGIN {printf \"%.2fs\", ($t1 - $t0) / 1000000000}")
    else
        elapsed_sec="$((t1 - t0))s"
    fi

    printf "[%2d/%2d] %-22s %-4s (%s)%s\n" "$num" "$total_tests" "$entry" "$test_status" "$elapsed_sec" "${note:+ [$note]}"

    result_status+=("$test_status")
    result_time+=("$elapsed_sec")
    result_iter+=("$iter_str")
    result_calc_te+=("$calc_te")
    result_std_te+=("$std_te")
    result_diff_te+=("$diff_te")
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

printf "%-4s %-22s %-8s %-6s %-18s %-18s %-11s %s\n" "No." "Entry" "Status" "Iter" "Calc TE (a.u.)" "Std TE (a.u.)" "Diff" "Note"
printf '%0.s-' {1..105}
echo

for idx in "${!entries[@]}"; do
    entry="${entries[$idx]}"
    num=$((idx + 1))
    status="${result_status[$idx]}"
    time_str="${result_time[$idx]}"
    istr="${result_iter[$idx]}"
    cte="${result_calc_te[$idx]}"
    ste="${result_std_te[$idx]}"
    dte="${result_diff_te[$idx]}"
    nt="${result_notes[$idx]}"

    printf "%2d.  %-22s %-8s %-6s %-18s %-18s %-11s %s\n" \
        "$num" "$entry" "$status" "$istr" "$cte" "$ste" "$dte" "$nt"
done

printf '%0.s-' {1..105}
echo

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
        echo "Result: ALL TESTS PASSED ($skipped_count skipped, ${total_elapsed}s)"
    else
        echo "Result: ALL TESTS PASSED (${total_elapsed}s)"
    fi
    exit 0
else
    echo "Result: $failed_count TEST(S) FAILED (${total_elapsed}s)"
    exit 1
fi
