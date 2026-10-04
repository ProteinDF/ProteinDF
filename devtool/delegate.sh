#!/usr/bin/env bash
# Hand a TASK file to an implementing agent (agy or GitHub Copilot CLI).
# The agent works only in the branch's own worktree (created if missing).
#
# Usage: devtool/delegate.sh [-a agy|copilot] [-i] [-c] [-W] [-b <branch>] [-m <model>] <TASK.md> [extra instructions]
#   -a  agent (default: $PDF_AGENT or agy)
#   -i  interactive session (permissions are asked as usual)
#       default is non-interactive: all tool permissions are auto-approved
#   -c  continue the agent's most recent session in the worktree
#       (e.g. when it stopped to ask something); extra instructions are sent as the reply
#   -W  if the agent stops on its usage limit, wait until the limit resets and
#       continue the session once (-c). Without -W the script exits with 75.
#   -b  branch (default: the "Branch:" line in the TASK file)
#   -m  model passed to the agent CLI
#
# Output is saved under <git-common-dir>/agent-logs/.
# Exit status: the agent's status, or 75 when it stopped on its usage limit.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
agent=${PDF_AGENT:-agy}
interactive=0
continue_session=0
wait_quota=0
branch=""
model=""
while getopts "a:icWb:m:h" opt; do
    case "$opt" in
        a) agent=$OPTARG ;;
        i) interactive=1 ;;
        c) continue_session=1 ;;
        W) wait_quota=1 ;;
        b) branch=$OPTARG ;;
        m) model=$OPTARG ;;
        *) sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
(($# >= 1)) || { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

task=$(realpath "$1"); shift
extra="$*"
[[ -f "$task" ]] || { echo "error: no such file: $task" >&2; exit 1; }
command -v "$agent" >/dev/null || { echo "error: '$agent' not found in PATH" >&2; exit 1; }

if [[ -z "$branch" ]]; then
    # e.g. "- **Branch**: `fix/gcc15-build`"
    branch=$(grep -m1 -oP '\*\*Branch\*\*:\s*`\K[^`]+' "$task" || true)
    [[ -n "$branch" ]] || { echo "error: no **Branch** line in $task; use -b" >&2; exit 1; }
fi

if git show-ref --verify --quiet "refs/heads/$branch"; then
    wt=$("$here/flow.sh" path "$branch")
else
    wt=$("$here/flow.sh" start "$branch" | tail -1)
fi

logdir="$(git rev-parse --path-format=absolute --git-common-dir)/agent-logs"
mkdir -p "$logdir"
log="$logdir/$(date +%Y%m%d-%H%M%S)_${branch//\//-}_${agent}.log"

read -r -d '' prompt <<EOF || true
あなたはProteinDFの実装担当です。
- 作業ディレクトリ: ${wt}
- 作業ブランチ: ${branch}(このブランチ以外にコミットしないこと)
- タスク指示書: ${task}

まず ${wt}/AGENTS.md を読み、次にタスク指示書を読んで、その指示に従って実装してください。
タスク指示書の末尾に「レビュー結果」がある場合は、最新のレビュー結果の修正依頼に対応してください。
タスク指示書と、作業ディレクトリ以外のファイルは編集しないこと。マージ・push・ブランチ削除はしないこと。
タスク指示書の範囲内のファイル編集・コマンド実行・コミットは承認済みなので、確認を求めずに進めること。判断が必要な点があれば、完了報告の「判断が必要な点」に書くこと。
最後に、AGENTS.md の「完了報告」の形式で報告を出力してください。
${extra}
EOF
if ((continue_session)); then
    prompt=${extra:-続けてください。}
fi

echo "agent:    $agent ($([[ $interactive == 1 ]] && echo interactive || echo non-interactive, auto-approve))"
echo "branch:   $branch"
echo "worktree: $wt"
echo "task:     $task"
echo "log:      $log"
echo

cmd=("$agent" --add-dir "$(dirname "$task")")
((continue_session)) && cmd+=(--continue)
[[ -n "$model" ]] && cmd+=(--model "$model")
case "$agent" in
    agy) auto_approve=--dangerously-skip-permissions ;;
    copilot) auto_approve=--allow-all-tools ;;
    *) echo "error: unknown agent '$agent' (agy|copilot)" >&2; exit 1 ;;
esac

set +e
if ((interactive)); then
    # keep the terminal attached; the session itself is not logged
    echo "(interactive session; not logged)" >"$log"
    (cd "$wt" && exec "${cmd[@]}" -i "$prompt")
    status=$?
else
    (cd "$wt" && exec "${cmd[@]}" -p "$prompt" "$auto_approve") 2>&1 | tee "$log"
    status=${PIPESTATUS[0]}
fi
set -e

echo
echo "=== $branch after the agent run (exit $status)"
git -C "$wt" log --oneline "develop..$branch"
git -C "$wt" status --short
echo

# usage limit (agy: "RESOURCE_EXHAUSTED ... Resets in 2h46m23s")
if ((!interactive)) && grep -qE 'RESOURCE_EXHAUSTED|quota reached' "$log"; then
    reset=$(grep -m1 -oP 'Resets in \K[0-9hms]+' "$log" || true)
    secs=0
    [[ "$reset" =~ ([0-9]+)h ]] && secs=$((secs + BASH_REMATCH[1] * 3600))
    [[ "$reset" =~ ([0-9]+)m ]] && secs=$((secs + BASH_REMATCH[1] * 60))
    [[ "$reset" =~ ([0-9]+)s ]] && secs=$((secs + BASH_REMATCH[1]))
    echo "usage limit reached (resets in ${reset:-unknown}, at $(date -d "+${secs} sec" +%H:%M))"
    resume_msg="利用上限のエラーで中断したので再開してください。worktreeのコミット・未コミットの変更を引き継いで、タスク指示書の残りの作業を行い、AGENTS.md の「完了報告」の形式で報告してください。"
    if ((wait_quota)) && [[ -n "$reset" ]]; then
        echo "waiting $((secs + 60))s, then continuing the session"
        sleep $((secs + 60))
        resume=("$0" -a "$agent" -c -b "$branch")
        [[ -n "$model" ]] && resume+=(-m "$model")
        exec "${resume[@]}" "$task" "$resume_msg"
    fi
    echo "resume later: devtool/delegate.sh -c $task  (or switch: -a copilot)"
    exit 75
fi

echo "next: ask Claude to review, e.g. \"$branch をレビューして\" (log: $log)"
exit "$status"
