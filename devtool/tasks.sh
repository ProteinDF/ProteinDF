#!/usr/bin/env bash
# Summarize the task files (doc/tasks/TASK_*.md) and the working branches.
#
# Usage: devtool/tasks.sh [--all]
#   default: tasks that are not merged yet, then the working branches
#   --all:   include merged tasks
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
top=$(git rev-parse --show-toplevel)
all=0
[[ "${1:-}" == --all ]] && all=1

printf '| %s | %s | %s | %s |\n' TASK Branch 状態 レビュー回数
printf '|---|---|---|---|\n'
for f in "$top"/doc/tasks/TASK_*.md; do
    [[ -e "$f" ]] || continue
    state=$(grep -m1 -oP '^- \*\*状態\*\*:\s*\K.*' "$f" || echo "?")
    ((all)) || [[ "$state" != マージ済み* ]] || continue
    branch=$(grep -m1 -oP '\*\*Branch\*\*:\s*`\K[^`]+' "$f" || echo "?")
    reviews=$(grep -c '^## レビュー結果' "$f" || true)
    printf '| %s | `%s` | %s | %s |\n' "$(basename "$f" .md)" "$branch" "$state" "$reviews"
done
echo
"$here/flow.sh" list
