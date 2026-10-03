# 実装担当エージェント(agy / GitHub Copilot)向けルール

このリポジトリでは、**実装はagyまたはGitHub Copilot、レビューはClaude、マージの承認はユーザー**が担当する。
全体の流れは`CONTRIBUTING.md`を参照。ここには実装担当が守るべきことだけを書く。

## 作業場所とブランチ(MUST)

- 作業はタスク指示書(`doc/tasks/TASK_*.md`)で指定されたブランチの**専用worktree**で行う。
  worktreeは`devtool/flow.sh start <branch>`(または`devtool/delegate.sh`)が作る
  (例: `~/orca/workspaces/ProteinDF/fix-gcc15-build`)。
- **`develop`のworktree(Claudeの作業場所)や他のworktreeのファイルを編集しない。** タスク指示書も読むだけにする(レビュー結果はClaudeが書き足す)。
- そのブランチ以外にコミットしない。`develop`・`main`へのマージ、push、ブランチやworktreeの削除、`git stash`はしない。
- `git status`と`git branch --show-current`で、作業前に場所とブランチを確認する。

## 実装

- タスク指示書の「対象」「完了の定義」に従う。指示書に書かれていない設計判断が必要になったら、勝手に決めずに完了報告の「判断が必要な点」に書いて止まってよい。
- スコープ外の変更(無関係なリファクタリング、整形のみの差分、ファイルの大量移動)をしない。
- 既存コードのスタイルに合わせる(`.clang-format`、`.editorconfig`)。変更した行だけを整形し、既存ファイル全体を整形し直さない。
- エラーを握りつぶさない(黙ってデフォルト値に置き換えない)。コンパイラの警告・エラーを`-fpermissive`や`-Wno-*`で抑制して済ませない(指示書で許可された場合を除く)。
- コミットは意味のある単位で分け、メッセージは既存の履歴に合わせる(`feat: ...`、`fix: ...`、`test: ...`、`docs: ...`、`style: ...`、英語)。

## 完了の確認

worktreeで次を実行し、成功することを確認する。

```bash
devtool/check.sh
```

`check.sh`は、ブランチの差分の空白エラー、変更行の`clang-format`、ビルド(`build-check/`)、`ctest`を実行する。
SKIPになった項目(GTestやclang-formatが未インストールなど)は、そのまま報告すればよい。

## 完了報告

最後に次の形式で報告する。**数字やコマンドの出力は要約・転記せず、実際の出力をそのまま貼る**
(過去に、存在しないテストファイルや実際と違うテスト件数が報告されたことがある)。

```markdown
## 完了報告: <branch>

### コミット
<git log --oneline develop..HEAD の出力>

### 実施内容
<完了の定義の各項目に対して何をしたか。項目番号を対応させる>

### 確認結果
<devtool/check.sh の summary 以降の出力をそのまま>
<指示書が求めた追加の確認(実行例・計測値など)とその出力>

### 判断が必要な点・未対応の点
<なければ「なし」>
```
