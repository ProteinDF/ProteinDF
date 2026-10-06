#!/usr/bin/env bash
# Definition-of-done checks for a working branch. Run it in the branch's worktree.
# Both the implementing agent and the reviewer run this and paste its output verbatim.
#
# Usage: devtool/check.sh [--base <ref>] [--build-dir <dir>] [--no-build] [--verbose]
#   --base       ref to diff against (default: develop)
#   --build-dir  CMake build directory (default: build-check)
#   --no-build   only run the diff checks
#   --verbose    print the whole ctest output (default: the last lines of the ctest
#                output and the failed tests; the full output is in <build-dir>/check-ctest.log)
# Extra CMake options can be passed via PDF_CMAKE_ARGS (e.g. "-DUSE_HDF5=on").
set -uo pipefail

base=develop
build_dir=build-check
do_build=1
verbose=0
while (($#)); do
    case "$1" in
        --base) base=$2; shift 2 ;;
        --build-dir) build_dir=$2; shift 2 ;;
        --no-build) do_build=0; shift ;;
        --verbose) verbose=1; shift ;;
        *) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
    esac
done

top=$(git rev-parse --show-toplevel) || exit 1
cd "$top"
merge_base=$(git merge-base "$base" HEAD) || exit 1

declare -a results=()
failed=0
record() { # <name> <PASS|FAIL|SKIP> [note]
    results+=("$(printf '%-14s %-4s %s' "$1" "$2" "${3:-}")")
    [[ "$2" == FAIL ]] && failed=1
    return 0
}

echo "branch: $(git branch --show-current)  HEAD: $(git rev-parse --short HEAD)  base: $base ($(git rev-parse --short "$merge_base"))"
echo "commits:"
git log --oneline "$merge_base..HEAD" | sed 's/^/  /'
mapfile -t changed < <(git diff --name-only --diff-filter=d "$merge_base" HEAD)
mapfile -t changed_src < <(printf '%s\n' "${changed[@]}" | grep -E '\.(h|hh|hpp|c|cc|cpp|cxx)$')
echo "changed files: ${#changed[@]} (C/C++: ${#changed_src[@]})"
[[ -n "$(git status --porcelain --untracked-files=no)" ]] \
    && echo "WARNING: uncommitted changes exist; they are built but not part of the diff checks"

# 1. whitespace errors in the branch diff
echo; echo "--- whitespace (git diff --check)"
if git diff --check "$merge_base" HEAD; then
    record whitespace PASS
else
    record whitespace FAIL
fi

# 2. clang-format on changed lines only (the existing code is not fully formatted)
echo; echo "--- clang-format (changed lines)"
if ((${#changed_src[@]} == 0)); then
    record clang-format SKIP "no C/C++ changes"
elif ! command -v git-clang-format >/dev/null || ! command -v clang-format >/dev/null; then
    record clang-format SKIP "clang-format not installed"
else
    out=$(git clang-format --diff "$merge_base" -- "${changed_src[@]}" 2>&1)
    if [[ -z "$out" || "$out" == *"no modified files"* || "$out" == *"did not modify"* ]]; then
        record clang-format PASS
    else
        echo "$out"
        record clang-format FAIL "run: git clang-format $merge_base"
    fi
fi

if ((do_build)); then
    # 3. build
    echo; echo "--- build ($build_dir)"
    log="$build_dir/check-build.log"
    gen=()
    if [[ ! -f "$build_dir/CMakeCache.txt" ]]; then
        mkdir -p "$build_dir"
        ninja --version >/dev/null 2>&1 && gen=(-G Ninja)
    fi
    # configure every run so that newly installed packages (e.g. GTest) are picked up
    # shellcheck disable=SC2086
    if cmake -S . -B "$build_dir" "${gen[@]}" ${PDF_CMAKE_ARGS:-} >"$build_dir/check-configure.log" 2>&1; then
        configured=1
    else
        configured=0
        tail -30 "$build_dir/check-configure.log"
        record configure FAIL "see $build_dir/check-configure.log"
    fi
    if ((configured)); then
        if cmake --build "$build_dir" -j "$(nproc)" >"$log" 2>&1; then
            record build PASS
        else
            grep -E "error|Error [0-9]" "$log" | sort -u | head -40
            record build FAIL "see $log"
        fi
        # warnings in the files this branch touched (informational)
        if ((${#changed_src[@]})); then
            n=$(grep -F -f <(printf '%s\n' "${changed_src[@]}") "$log" | grep -c "warning:")
            grep -F -f <(printf '%s\n' "${changed_src[@]}") "$log" | grep "warning:" | sort -u | head -20
            record warnings "$([[ $n == 0 ]] && echo PASS || echo FAIL)" \
                "$n warning(s) in changed files (only files recompiled this run)"
        fi

        # 4. tests
        echo; echo "--- tests (ctest)"
        ntests=$(ctest --test-dir "$build_dir" -N 2>/dev/null | sed -n 's/^Total Tests: //p')
        if [[ -z "$ntests" || "$ntests" == 0 ]]; then
            record tests SKIP "no tests registered (GTest not found?)"
        else
            ctest_log="$build_dir/check-ctest.log"
            ctest --test-dir "$build_dir" --output-on-failure -j "$(nproc)" >"$ctest_log" 2>&1
            ctest_status=$?
            if ((verbose)); then
                cat "$ctest_log"
            elif ((ctest_status == 0)); then
                tail -n 6 "$ctest_log"
            else
                # failed tests with their output, then the ctest summary
                grep -E '^\[  FAILED  \]|Failed|\*\*\*' "$ctest_log" | sort -u | head -40
                tail -n 15 "$ctest_log"
            fi
            if ((ctest_status == 0)); then
                record tests PASS "$ntests test(s)"
            else
                record tests FAIL "see $ctest_log"
            fi
        fi
    fi
fi

echo; echo "=== summary"
printf '%s\n' "${results[@]}"
exit $failed
