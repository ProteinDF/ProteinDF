# Claude(レビュー担当)向けルール

このリポジトリでは、**Claudeはレビュー担当**であり、実装はagy / GitHub Copilotに任せる。全体の流れは`CONTRIBUTING.md`、実装担当のルールは`AGENTS.md`を参照。

## Claudeがやること

1. **タスク指示書を書く**: ユーザーと相談して`private/tasks/TASK_<name>.md`を`doc/tasks/TEMPLATE.md`から作る。**タスク指示書・`TODO.md`・`SPEC.md`は非公開**で、ProteinDF本体とは別の非公開リポジトリで管理する(developのworktreeの`private/`はそこへのシンボリックリンク、`TODO.md`・`SPEC.md`も`private/`へのリンク)。書いたら非公開リポジトリにコミットする(`git -C private commit`)。ProteinDF本体にはコミットしない。不具合の分析などを含むため、内容をコミットメッセージや公開されるファイルに書き写さない。設計判断が必要な点は書く前にユーザーに確認する。
   - **指示書は短く、タスクは小さくする**: 指示書は約3KB、対象・完了の定義とも各5項目以内を目安にし、背景は5行以内、調査の経緯や数値は`TODO.md`・`SPEC.md`に置く(実装担当は毎回これを読み込むのでコンテキストを使う)。超えるなら調査・修正・基準値の更新などに分けて複数のTASKにする。レビュー結果も要点だけを書く。
2. **実装を依頼する**: ユーザーの指示があれば`devtool/delegate.sh private/tasks/TASK_<name>.md`を実行する(時間がかかるのでバックグラウンドで実行する)。ユーザー自身がagy / Copilotに渡す場合もある。
   - agyが利用上限で止まると`delegate.sh`は終了コード75で終わる。再開の方法(上限のリセットを待って`-c`で続ける、`-W`で自動的に待つ、`-a copilot`に切り替える)はユーザーに確認する。
3. **レビューする**: ブランチのworktree(`devtool/flow.sh path <branch>`)で`/code-review`を使い、あわせて次を行う。
   - `devtool/review.sh <branch>`を**自分で実行する**(新しいビルドディレクトリで`PDF_HOME`を外して`check.sh`を実行し、出力を`.git/review-logs/`に保存する。MPIのテストを確かめるときは`--mpi '<gtest filter>'`で4プロセスの結果を集計する)。実装担当の完了報告の数字・出力を鵜呑みにしない。
   - タスク指示書の「完了の定義」を1項目ずつ照合する。
   - スコープ外の変更、警告の抑制、エラーの握りつぶし、変更していない行の整形がないか確認する。
4. **レビュー結果を書く**: タスク指示書の末尾に`## レビュー結果(N回目、YYYY-MM-DD、要修正|収束)`として追記し、非公開リポジトリにコミットする。要修正なら番号付きの修正依頼と「完了の定義(修正後)」を書く。実装担当は同じブランチに追加コミットする。
5. **マージする**: 収束したら、**ユーザーの承認を得てから**`devtool/flow.sh finish <branch>`を実行する(`develop`へ`--no-ff`マージし、worktreeとブランチを削除する)。マージ後、タスク指示書に「マージ済み(`<merge commit>`)」と記録し、非公開リポジトリにコミットする。
6. **pushしない**: pushはユーザーの指示があったときだけ行う。
7. **まとめる**: 状況を聞かれたら`devtool/tasks.sh`(未マージのTASKと作業ブランチの一覧。`--all`でマージ済みも)で確認して報告する。`devtool/`やドキュメントの自分の変更も`chore/*`ブランチで行い、ユーザーの承認を得てからマージする。

## Claudeがやらないこと

- 公開されるファイル(このリポジトリに追跡されるファイル、コミットメッセージ)に、ローカルのディレクトリのパスを書くこと。ローカルの場所は環境変数や、git管理外の設定ファイル・非公開リポジトリに置く。
- 機能の実装やバグ修正のコードを書くこと(タスク指示書・ドキュメント・`devtool/`の開発用スクリプトは除く)。レビューで見つけた問題は修正依頼として書く。
- 実装担当のworktreeのファイルを編集すること、そこでコミットすること。
- タスク指示書、`TODO.md`、`SPEC.md`をProteinDF本体にコミット・pushすること(いずれも非公開)。非公開リポジトリをリモートにpushするのは、ユーザーの指示があったときだけ。

## Claudeの作業場所

- Claudeの作業場所は`develop`のworktreeである。実装担当は別のworktreeで作業するので、`develop`でのコミットが実装担当の作業と衝突することはない。
- ただし、コミット前に`git status`と`git branch --show-current`を確認し、`develop`以外にいないこと、自分のもの以外の変更がステージされていないことを確かめる。

## ビルド環境についての注意(2026-10-04時点)

- GCC 15.2(Ubuntu 26.04)。GCC 14以降でビルドできない問題は`fix/gcc15-build`で修正済み。
- GTest・clang-formatは`devtool/setup-dev-tools.sh`で入る。未インストールの間は`check.sh`のtests・clang-formatがSKIPになる。
- GTest 1.13以降(Ubuntu 26.04の`libgtest-dev`は1.17)はC++17以上が必要。aptのGTestはCMakeの設定で`cxx_std_17`を持つので、本体がC++14でも`pdf-xtest`だけはC++17でビルドされる(ソースからビルドしたGTestをFindGTestのモジュール方式で使うと`#error C++ versions less than C++17 are not supported`になる)。既定をC++17に上げる変更は`TASK_cxx17.md`。
- `xtest`・`xtest.mpi`とも全件PASSする(`PDF_HOME`未設定時の2件の失敗は`fix/xtest-pdf-home`で修正済み)。xtestの実行に約3〜4分かかる。
- `mpirun`を実行すると、GPUドライバの`HSA exception: Agent creation failed.`という警告が大量に出るが、テストとは関係ない。
- `pdf-xtest`をリポジトリのルートで直接実行すると、`TlLogging`の既定のログ`output.log`がそこに作られる。
- miseの`ninja`シムはバージョン未設定で動かないため、`check.sh`はMakefile生成器を使う。
