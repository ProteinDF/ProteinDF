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

- GCC 15.2(Ubuntu 26.04)。GCC 14以降でビルドできない問題は`fix/gcc15-build`で修正済み。
- GTest・clang-formatは`devtool/setup-dev-tools.sh`で入る。未インストールの間は`check.sh`のtests・clang-formatがSKIPになる。
- GTest 1.13以降(Ubuntu 26.04の`libgtest-dev`は1.17)はC++17以上が必要。aptのGTestはCMakeの設定で`cxx_std_17`を持つので、本体がC++14でも`pdf-xtest`だけはC++17でビルドされる(ソースからビルドしたGTestをFindGTestのモジュール方式で使うと`#error C++ versions less than C++17 are not supported`になる)。既定をC++17に上げる変更は`TASK_cxx17.md`。
- `xtest`は316件中2件が失敗する(`TlDenseSymmetricMatrix_Lapack.multiplication_MV`・`multiplication_VMV`、例外`basic_string: construction from null is not valid`)。`xtest.mpi`は成功。xtestの実行に約3分かかる。
- miseの`ninja`シムはバージョン未設定で動かないため、`check.sh`はMakefile生成器を使う。
