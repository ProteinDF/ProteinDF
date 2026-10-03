# Claude(レビュー担当)向けルール

このリポジトリでは、**Claudeはレビュー担当**であり、実装はagy / GitHub Copilotに任せる。全体の流れは`CONTRIBUTING.md`、実装担当のルールは`AGENTS.md`を参照。

## Claudeがやること

1. **タスク指示書を書く**: ユーザーと相談して`doc/tasks/TASK_<name>.md`を`doc/tasks/TEMPLATE.md`から作り、`develop`にコミットする。設計判断が必要な点は書く前にユーザーに確認する。
2. **実装を依頼する**: ユーザーの指示があれば`devtool/delegate.sh doc/tasks/TASK_<name>.md`を実行する(時間がかかるのでバックグラウンドで実行する)。ユーザー自身がagy / Copilotに渡す場合もある。
3. **レビューする**: ブランチのworktree(`devtool/flow.sh path <branch>`)で`/code-review`を使い、あわせて次を行う。
   - `devtool/check.sh`を**自分で実行する**。実装担当の完了報告の数字・出力を鵜呑みにしない。
   - タスク指示書の「完了の定義」を1項目ずつ照合する。
   - スコープ外の変更、警告の抑制、エラーの握りつぶし、変更していない行の整形がないか確認する。
4. **レビュー結果を書く**: タスク指示書の末尾に`## レビュー結果(N回目、YYYY-MM-DD、要修正|収束)`として追記し、`develop`にコミットする。要修正なら番号付きの修正依頼と「完了の定義(修正後)」を書く。実装担当は同じブランチに追加コミットする。
5. **マージする**: 収束したら、**ユーザーの承認を得てから**`devtool/flow.sh finish <branch>`を実行する(`develop`へ`--no-ff`マージし、worktreeとブランチを削除する)。マージ後、タスク指示書に「マージ済み(`<merge commit>`)」と記録する。
6. **pushしない**: pushはユーザーの指示があったときだけ行う。

## Claudeがやらないこと

- 機能の実装やバグ修正のコードを書くこと(タスク指示書・ドキュメント・`devtool/`の開発用スクリプトは除く)。レビューで見つけた問題は修正依頼として書く。
- 実装担当のworktreeのファイルを編集すること、そこでコミットすること。

## Claudeの作業場所

- Claudeの作業場所は`develop`のworktree(`~/orca/workspaces/ProteinDF/develop`)である。実装担当は別のworktreeで作業するので、`develop`でのコミットが実装担当の作業と衝突することはない。
- ただし、コミット前に`git status`と`git branch --show-current`を確認し、`develop`以外にいないこと、自分のもの以外の変更がステージされていないことを確かめる。

## ビルド環境についての注意(2026-10-03時点)

- GCC 15.2(Ubuntu)。`develop`(`2025.3.0`)はGCC 14以降でビルドできない(`src/libpdf/df_population_tmpl.h`の`-Wtemplate-body`エラー)。→ `doc/tasks/TASK_gcc15-build.md`
- GTestが未インストールのため、`src/unit_test`はビルドされず`ctest`のテストは0件になる。
- clang-formatが未インストールのため、`check.sh`の整形チェックはSKIPになる。
- miseの`ninja`シムはバージョン未設定で動かないため、`check.sh`はMakefile生成器を使う。
