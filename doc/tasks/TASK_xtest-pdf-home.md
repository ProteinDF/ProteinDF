# TASK: xtestの2件の失敗(PDF_HOME未設定時の例外)の修正

- **Branch**: `fix/xtest-pdf-home`
- **作成**: 2026-10-03(Claude)
- **状態**: 未着手

> C++17でビルドした`xtest`は316件中2件が失敗する。
>
> ```text
> [  FAILED  ] TlDenseSymmetricMatrix_Lapack.multiplication_MV
> [  FAILED  ] TlDenseSymmetricMatrix_Lapack.multiplication_VMV
> C++ exception with description "basic_string: construction from null is not valid" thrown in the test body.
> ```
>
> **原因(Claude確認済み)**: 2件のテスト(`src/unit_test/tl_dense_symmetric_matrix_lapack_test.cc`の255行目〜)は、`TlSystem::getEnv("PDF_HOME")`でデータファイル(`data/unit_test/{M,V,MV,VMV}.mat`、リポジトリに含まれている)の場所を決めている。`TlSystem::getEnv`(`src/libpdftl/TlSystem.cpp:58`)は`std::string ans(std::getenv(key.c_str()));`で、環境変数が未設定だと`std::getenv`が返すNULLから`std::string`を作る。これは未定義動作で、GCC 15のlibstdc++では例外になる。`PDF_HOME`にリポジトリのルートを設定して実行すると、2件ともPASSすることを確認した。
>
> 本体のコード(`PdfKeyword.cpp`・`Fl_Db_Basis.cpp`・`DfInitialGuessHarris.cpp`)は`std::getenv`を直接使っていて、`TlSystem::getEnv`を呼んでいるのはこのテストだけである。

## 役割分担・ブランチ運用(MUST)

- 実装はagy / GitHub Copilot、レビューはClaudeが担当する。`AGENTS.md`を必ず読むこと。
- 上記のブランチの専用worktreeで作業する。`develop`へは自分でマージしない。
- この指示書は編集しない(レビュー結果はClaudeが追記する)。
- **`chore/cxx17`(`TASK_cxx17.md`)がdevelopにマージされてから着手する。** C++14のままではGTestがビルドできないため。

## 対象

1. **`TlSystem::getEnv`**: 環境変数が未設定のときに未定義動作にならないようにする。未設定の場合の戻り値(空文字列にする、など)は実装担当が決め、docコメントに書く。
2. **ctestでの`PDF_HOME`**: `ctest`(および`devtool/check.sh`)を実行したときに、`PDF_HOME`を手で設定しなくても2件がデータファイルを見つけられるようにする。例: `src/unit_test/CMakeLists.txt`で、`PDF_HOME`が未設定ならテストの`ENVIRONMENT`プロパティにソースディレクトリを設定する。ユーザーが`PDF_HOME`を設定している場合の扱い(優先するかどうか)も決めて、完了報告に書く。
3. **黙ってPASSしないようにする**: いまの2件のテストは、データファイルが読めなかった場合に、行数0のまま比較のループを1回も回らずにPASSしてしまう可能性がある。読み込みに失敗した場合や、行列のサイズが期待と違う場合にテストが失敗するようにする(`load`の戻り値の確認、サイズの`ASSERT_EQ`など)。`load`が失敗をどう伝えるか(戻り値・例外)を調べてから書くこと。

## 完了の定義

1. `PDF_HOME`を設定せずに(`env -u PDF_HOME`)`devtool/check.sh`を実行し、build・testsがPASSになる(`xtest`の316件すべてと`xtest.mpi`)。gtestの出力の`[  PASSED  ]`の行をそのまま完了報告に貼る。
2. `PDF_HOME`を存在しないディレクトリ(例: `/nonexistent`)にして`pdf-xtest --gtest_filter='TlDenseSymmetricMatrix_Lapack.multiplication_*'`を直接実行すると、2件が(例外ではなく、読み込み失敗として)FAILすることを確認し、その出力を貼る。対象2で「ユーザー設定を優先しない」と決めた場合は、代わりにデータファイルを一時的に読めなくして確認する。
3. `TlSystem::getEnv`の未設定時の挙動を確認するテストを`src/unit_test`に追加する(既存のテストファイルの置き方に合わせる)。
4. `devtool/check.sh`が通る(clang-formatを含む)。
5. `AGENTS.md`の形式で完了報告を出す。
