#!/usr/bin/env bash
# Summarize the task files and the working branches.
# Task files are private: <develop worktree>/private/tasks/TASK_*.md, where private/ is a
# symlink to the private notes repository (override the directory with PDF_TASK_DIR).
#
# Usage: devtool/tasks.sh [--all]
#   default: tasks that are not merged yet, then the working branches
#   --all:   include merged tasks
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
task_dir=${PDF_TASK_DIR:-$("$here/flow.sh" path develop)/private/tasks}
[[ -d "$task_dir" ]] || { echo "error: no task directory $task_dir" >&2; exit 1; }
all=0
[[ "${1:-}" == --all ]] && all=1

printf '| %s | %s | %s | %s |\n' TASK Branch 状態 レビュー回数
printf '|---|---|---|---|\n'
for f in "$task_dir"/TASK_*.md; do
    [[ -e "$f" ]] || continue
    state=$(grep -m1 -oP '^- \*\*状態\*\*:\s*\K.*' "$f" || echo "?")
    ((all)) || [[ "$state" != マージ済み* ]] || continue
    branch=$(grep -m1 -oP '\*\*Branch\*\*:\s*`\K[^`]+' "$f" || echo "?")
    reviews=$(grep -c '^## レビュー結果' "$f" || true)
    printf '| %s | `%s` | %s | %s |\n' "$(basename "$f" .md)" "$branch" "$state" "$reviews"
done
echo
"$here/flow.sh" list
