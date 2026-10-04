#!/usr/bin/env bash
# git-flow style branch operations, with one git worktree per working branch.
#
# Usage:
#   devtool/flow.sh start  feature/<name>|fix/<name>|chore/<name>   # from develop
#   devtool/flow.sh start  release/<YYYY.M.PATCH>                   # from develop, bumps version
#   devtool/flow.sh start  hotfix/<YYYY.M.PATCH>                    # from main, bumps version
#   devtool/flow.sh finish <branch> [--keep]                        # merge --no-ff (+ tag)
#   devtool/flow.sh path   <branch>                                 # print worktree path
#   devtool/flow.sh list
#
# Worktrees are created next to the worktree that has `develop` checked out
# (named after the branch, e.g. feature-foo). Override with PDF_WORKTREE_ROOT.
# Nothing is pushed; push manually when the user says so.
set -euo pipefail

DEVELOP=develop
MAIN=$(git config --get gitflow.branch.master || echo main)

die() { echo "error: $*" >&2; exit 1; }
info() { echo "==> $*"; }

# Path of the worktree that has <branch> checked out (empty if none).
wt_of_branch() {
    git worktree list --porcelain | awk -v ref="refs/heads/$1" '
        /^worktree /{p=substr($0, 10)} $0=="branch "ref{print p; exit}'
}

worktree_root() {
    if [[ -n "${PDF_WORKTREE_ROOT:-}" ]]; then
        echo "$PDF_WORKTREE_ROOT"; return
    fi
    local dev_wt
    dev_wt=$(wt_of_branch "$DEVELOP")
    [[ -n "$dev_wt" ]] || die "'$DEVELOP' is not checked out in any worktree; set PDF_WORKTREE_ROOT"
    dirname "$dev_wt"
}

require_clean() {
    [[ -z "$(git -C "$1" status --porcelain --untracked-files=no)" ]] \
        || die "worktree $1 has uncommitted changes"
}

branch_type() { echo "${1%%/*}"; }

check_version() {
    [[ "$1" =~ ^[0-9]{4}\.[1-9][0-9]?\.[0-9]+$ ]] || die "version must be YYYY.M.PATCH (got '$1')"
}

# Rewrite PROJECT_VERSION_* in CMakeLists.txt and commit (as done in 1781e97).
bump_version() {
    local wt=$1 ver=$2 major minor rev
    IFS=. read -r major minor rev <<<"$ver"
    sed -i -E \
        -e "s/^(set\(PROJECT_VERSION_MAJOR )\"[^\"]*\"/\1\"$major\"/" \
        -e "s/^(set\(PROJECT_VERSION_MINOR )\"[^\"]*\"/\1\"$minor\"/" \
        -e "s/^(set\(PROJECT_VERSION_REVISION )\"[^\"]*\"/\1\"$rev\"/" \
        "$wt/CMakeLists.txt"
    if git -C "$wt" diff --quiet -- CMakeLists.txt; then
        info "version already $ver"
    else
        git -C "$wt" commit -q -m "style: update release version" -- CMakeLists.txt
        info "bumped version to $ver"
    fi
}

cmd_start() {
    local branch=${1:-} type name base wt
    [[ "$branch" == */* ]] || die "usage: flow.sh start <type>/<name>"
    type=$(branch_type "$branch"); name=${branch#*/}
    case "$type" in
        feature|fix|chore|release) base=$DEVELOP ;;
        hotfix) base=$MAIN ;;
        *) die "unknown branch type '$type' (feature|fix|chore|release|hotfix)" ;;
    esac
    [[ "$type" == release || "$type" == hotfix ]] && check_version "$name"
    git check-ref-format --branch "$branch" >/dev/null || die "invalid branch name '$branch'"
    git show-ref --verify --quiet "refs/heads/$branch" && die "branch '$branch' already exists"

    wt="$(worktree_root)/${branch//\//-}"
    [[ -e "$wt" ]] && die "$wt already exists"

    info "creating $branch from $base at $wt"
    git worktree add -q -b "$branch" "$wt" "$base"
    [[ "$type" == release || "$type" == hotfix ]] && bump_version "$wt" "$name"
    echo "$wt"
}

# Merge <branch> into <target> with --no-ff inside the worktree that has
# <target> checked out (a temporary worktree is used if there is none).
merge_into() {
    local target=$1 src=$2 wt tmp=""
    wt=$(wt_of_branch "$target")
    if [[ -z "$wt" ]]; then
        tmp=$(mktemp -d); wt=$tmp/wt
        git worktree add -q "$wt" "$target"
    fi
    require_clean "$wt"
    info "merging $src into $target (in $wt)"
    if ! git -C "$wt" merge --no-ff --no-edit "$src"; then
        die "merge conflict in $wt; resolve it there and commit, then re-run finish"
    fi
    if [[ -n "$tmp" ]]; then
        git worktree remove "$wt"; rm -rf "$tmp"
    fi
}

cmd_finish() {
    local branch=${1:-} keep=0 type name wt
    [[ -n "$branch" ]] || die "usage: flow.sh finish <branch> [--keep]"
    [[ "${2:-}" == --keep ]] && keep=1
    git show-ref --verify --quiet "refs/heads/$branch" || die "no branch '$branch'"
    type=$(branch_type "$branch"); name=${branch#*/}
    wt=$(wt_of_branch "$branch")
    if [[ -n "$wt" ]]; then
        require_clean "$wt"
        # git worktree remove refuses untracked files; check before merging
        if ((!keep)) && [[ -n "$(git -C "$wt" status --porcelain)" ]]; then
            git -C "$wt" status --short
            die "worktree $wt has untracked files; remove them (or use --keep) and re-run finish"
        fi
    fi

    case "$type" in
        feature|fix|chore)
            merge_into "$DEVELOP" "$branch"
            ;;
        release|hotfix)
            check_version "$name"
            git show-ref --verify --quiet "refs/tags/$name" && die "tag '$name' already exists"
            merge_into "$MAIN" "$branch"
            git tag -a "$name" -m "$name" "$MAIN"
            info "tagged $name"
            merge_into "$DEVELOP" "$name"
            ;;
        *) die "unknown branch type '$type'" ;;
    esac

    if (( keep )); then
        info "kept branch $branch${wt:+ and worktree $wt}"
    else
        [[ -n "$wt" ]] && git worktree remove "$wt" && info "removed worktree $wt"
        git branch -d "$branch"
    fi
    if [[ "$type" == release || "$type" == hotfix ]]; then
        info "done. push when approved: git push github $DEVELOP $MAIN $name"
    else
        info "done. push when approved: git push github $DEVELOP"
    fi
}

cmd_path() {
    local wt
    wt=$(wt_of_branch "${1:?usage: flow.sh path <branch>}")
    [[ -n "$wt" ]] || die "'$1' is not checked out in any worktree"
    echo "$wt"
}

cmd_list() {
    local ref wt
    git worktree list
    echo
    for ref in $(git for-each-ref --format='%(refname:short)' \
            refs/heads/feature refs/heads/fix refs/heads/chore refs/heads/release refs/heads/hotfix); do
        wt=$(wt_of_branch "$ref")
        printf '%-40s %3s commits ahead of %s  %s\n' "$ref" \
            "$(git rev-list --count "$DEVELOP..$ref")" "$DEVELOP" "${wt:-(no worktree)}"
    done
}

case "${1:-}" in
    start)  shift; cmd_start "$@" ;;
    finish) shift; cmd_finish "$@" ;;
    path)   shift; cmd_path "$@" ;;
    list)   shift; cmd_list ;;
    *) sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
