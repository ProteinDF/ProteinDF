# Contributing

## ブランチモデル(git-flow)

- **`main`**: リリース済みの状態だけを反映する。直接コミットしない。`release/*`・`hotfix/*`のマージだけを受け取り、マージごとにバージョンタグ(例: `2025.3.0`、`v`なし)を打つ。
- **`develop`**: 日常開発の統合ブランチ。
- **`feature/*`・`fix/*`・`chore/*`**: `develop`から切り、レビュー承認後に`develop`へ`--no-ff`でマージする。
- **`release/YYYY.M.PATCH`**: `develop`から切ってリリースを準備し、`main`(タグ付き)と`develop`へマージする。
- **`hotfix/YYYY.M.PATCH`**: `main`から切って緊急修正し、`main`(タグ付き)と`develop`へマージする。

バージョンは`YYYY.M.PATCH`(例: `2025.3.0`)で、`CMakeLists.txt`の`PROJECT_VERSION_MAJOR/MINOR/REVISION`に書く。

git-flowのCLIは使わず、`devtool/flow.sh`で操作する。**作業ブランチごとに専用のgit worktreeを作る**(`develop`のworktreeの隣。例: `~/orca/workspaces/ProteinDF/feature-foo`)。

```bash
devtool/flow.sh start feature/foo        # developから切ってworktreeを作る
devtool/flow.sh start release/2026.10.0  # developから切り、バージョンを書き換えてコミット
devtool/flow.sh start hotfix/2025.3.1    # mainから切り、バージョンを書き換えてコミット
devtool/flow.sh finish feature/foo       # developへ--no-ffマージし、worktreeとブランチを削除
devtool/flow.sh finish release/2026.10.0 # mainへマージしてタグを打ち、タグをdevelopへマージ
devtool/flow.sh list                     # 作業ブランチとworktreeの一覧
```

`flow.sh`はpushしない。pushは別途行う(`git push github develop`など)。

## 役割分担(AIエージェントとの開発)

| 担当 | 役割 | ルール |
|---|---|---|
| ユーザー | 方針・設計判断、マージの承認、push | |
| Claude | タスク指示書の作成、レビュー、承認後のマージ | `CLAUDE.md` |
| agy / GitHub Copilot | 実装(作業ブランチの専用worktreeで) | `AGENTS.md` |

Claudeと実装担当はworktreeを分けるので、同じ作業ディレクトリでブランチを切り替え合って変更が混ざる事故は起きない。

### 流れ

1. **指示書**: ユーザーとClaudeが相談し、Claudeが`doc/tasks/TASK_<name>.md`を書いて`develop`にコミットする。指示書には作業ブランチ名(`**Branch**:`の行)、対象、完了の定義を書く。
2. **実装**: 実装担当に指示書を渡す。
   ```bash
   devtool/delegate.sh doc/tasks/TASK_<name>.md            # agy(非対話、権限は自動承認)
   devtool/delegate.sh -a copilot doc/tasks/TASK_<name>.md # GitHub Copilot CLI
   devtool/delegate.sh -i doc/tasks/TASK_<name>.md         # 対話モード(権限は都度確認)
   ```
   ブランチとworktreeがなければ作る。非対話モードの出力は`.git/agent-logs/`に保存される。
3. **レビュー**: Claudeがworktreeで`/code-review`と`devtool/check.sh`を実行し、結果を指示書の末尾に追記する。要修正なら手順2に戻る(同じコマンドで、実装担当は最新のレビュー結果に対応する)。
4. **マージ**: 収束したらユーザーが承認し、Claudeが`devtool/flow.sh finish <branch>`を実行する。

## 完了の確認

テストと整形チェックには GoogleTest と clang-format が必要。初回に次を実行する(Linuxはapt+sudo、macOSはHomebrew)。

```bash
devtool/setup-dev-tools.sh
```

```bash
devtool/check.sh   # 差分の空白エラー、変更行のclang-format、ビルド、ctest
```
