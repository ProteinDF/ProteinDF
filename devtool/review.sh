#!/usr/bin/env bash
# Reviewer-side verification of a working branch (run by Claude, from any worktree).
# Rebuilds the branch from scratch, runs check.sh with PDF_HOME unset, optionally
# runs MPI tests on 4 processes, and saves the full output for the review record.
#
# Usage: devtool/review.sh [--keep-build] [--mpi <gtest filter>] [--np <n>] <branch>
#   --keep-build  keep <worktree>/build-review afterwards (default: removed)
#   --mpi         also run pdf-xtest.MPI with this --gtest_filter on --np processes
#                 and count the per-process results
#   --np          number of MPI processes (default: 4)
#
# Output is saved under <git-common-dir>/review-logs/.
set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
keep=0
mpi_filter=""
np=4
while (($#)); do
    case "$1" in
        --keep-build) keep=1; shift ;;
        --mpi) mpi_filter=$2; shift 2 ;;
        --np) np=$2; shift 2 ;;
        -*) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
        *) break ;;
    esac
done
(($# == 1)) || { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
branch=$1

wt=$("$here/flow.sh" path "$branch") || exit 1
logdir="$(git rev-parse --path-format=absolute --git-common-dir)/review-logs"
mkdir -p "$logdir"
log="$logdir/$(date +%Y%m%d-%H%M%S)_${branch//\//-}.log"
build="$wt/build-review"

run() {
    echo "branch:   $branch"
    echo "worktree: $wt"
    echo
    echo "--- commits (develop..$branch)"
    git -C "$wt" log --oneline "develop..$branch"
    echo
    echo "--- diffstat"
    git -C "$wt" diff --stat "develop...$branch"
    echo
    echo "--- worktree status (uncommitted / untracked)"
    git -C "$wt" status --short
    echo

    rm -rf "$build"
    echo "--- check.sh (fresh build, PDF_HOME unset)"
    # use the reviewer's (develop) check.sh so an edited copy on the branch cannot weaken it
    (cd "$wt" && env -u PDF_HOME "$here/check.sh" --build-dir build-review)
    check_status=$?

    mpi_status=0
    if [[ -n "$mpi_filter" ]]; then
        echo
        echo "--- MPI: mpirun -np $np pdf-xtest.MPI --gtest_filter='$mpi_filter'"
        local exe tmp out
        exe=$(find "$build" -name pdf-xtest.MPI -type f | head -1)
        if [[ -z "$exe" ]]; then
            echo "pdf-xtest.MPI not built"
            mpi_status=1
        else
            # run in a scratch dir: the binary writes output.log to its cwd
            tmp=$(mktemp -d)
            out=$(cd "$tmp" && timeout 600 mpirun -np "$np" "$exe" --gtest_filter="$mpi_filter" 2>&1)
            mpi_status=$?
            rm -rf "$tmp"
            out=$(sed 's/\x1b\[[0-9;]*m//g' <<<"$out")
            echo "mpirun exit: $mpi_status"
            echo "processes reporting PASSED: $(grep -c '^\[  PASSED  \]' <<<"$out") / $np"
            echo "FAILED lines: $(grep -c '^\[  FAILED  \]' <<<"$out")"
            grep -E '^\[  (PASSED|FAILED)  \]' <<<"$out" | sort | uniq -c
        fi
    fi

    if ((keep)); then
        echo; echo "kept $build"
    else
        rm -rf "$build"
    fi
    echo
    echo "=== review: check.sh exit $check_status${mpi_filter:+, mpi exit $mpi_status}"
    ((check_status == 0 && mpi_status == 0))
}

run 2>&1 | tee "$log"
status=${PIPESTATUS[0]}
echo "log: $log"
exit "$status"
