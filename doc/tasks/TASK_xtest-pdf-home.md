# TASK: xtestの2件の失敗(PDF_HOME未設定時の例外)の修正

- **Branch**: `fix/xtest-pdf-home`
- **作成**: 2026-10-03(Claude)
- **状態**: レビュー中(収束、マージ承認待ち)

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

## レビュー結果(1回目、2026-10-03、要修正)

`fix/xtest-pdf-home`(`aef2ac5`・`5f7e677`・`588b718`・`0032e54`)をレビューした。`TlSystem::getEnv`の修正(`aef2ac5`)、`TlSystem::getEnv`のテスト(`588b718`)、ctestへの`PDF_HOME`の設定とテストの読み込み・サイズの確認(`0032e54`)は妥当である。ユーザーが設定した`PDF_HOME`を優先する設計も、`data/unit_test`がインストール先にも置かれる(`data/unit_test/CMakeLists.txt`)ので問題ない。

### 修正依頼

1. **【実バグ・スコープ外】`5f7e677`(`TlDenseSymmetricMatrixObject::load`・`TlDenseGeneralMatrixObject::load`の`throw;`を`return false;`に変更)を取り消す。** 変更前は、行列ファイルが開けない・形式が不正な場合にプログラムが異常終了していた(catch節の外の`throw;`による`std::terminate`)。変更後は`false`が返るが、本体(`src/libpdf`・`src/pdf`・`src/tools`)で`load`を呼んでいる箇所のうち、戻り値を確認しているのは約13か所で、約269か所は戻り値を無視している(Claudeがgrepで数えた。行列以外の`load`も含む概数)。たとえば`DfTotalEnergy.h:389`の`rho.load(...)`や`DfObject::getPInMatrix`の`P.load(path)`は、ファイルが読めなくても空の行列のまま計算を続けることになる。これは`AGENTS.md`の「エラーを握りつぶさない」に反し、「エラーで止まる」から「黙って誤った結果を出す」への悪化である。また、ライブラリの挙動変更はこのタスクの対象外である。
   - `git revert 5f7e677`で取り消す(履歴を書き換えない)。
   - テスト側は、`load`の前に`ASSERT_TRUE(TlFile::isExistFile(path))`でファイルの存在を確認する。これで`PDF_HOME`が誤っている場合は、異常終了せずに読み込み失敗としてFAILになる。`ASSERT_TRUE(M.load(...))`とサイズの確認はそのまま残してよい。
   - 形式が不正なファイルの場合は(変更前と同じく)異常終了するが、それでよい。`throw;`を適切な例外に置き換える改善は、必要なら別のタスクにする(本体の多数の呼び出し元に影響するため)。

### 完了の定義(修正後)

1. 上記に対応し、同じブランチに追加コミットする。
2. `env -u PDF_HOME devtool/check.sh`でbuild・testsがPASSする(`[  PASSED  ]`の行を貼る)。
3. `PDF_HOME=/nonexistent`で`pdf-xtest --gtest_filter='TlDenseSymmetricMatrix_Lapack.multiplication_*'`を直接実行し、2件が異常終了せずにFAILすることを確認し、出力を貼る。
4. `git diff develop...HEAD -- src/libpdftl/tl_dense_symmetric_matrix_object.cc src/libpdftl/tl_dense_general_matrix_object.cc`の出力が空であること(本体の`load`が変わっていないこと)を貼る。

## レビュー結果(2回目、2026-10-03、収束)

修正コミット`aeb9019`(`5f7e677`のrevert)・`fe5b294`(読み込み前のファイル存在確認)を確認した。修正依頼1に対応済み。Claudeが次を実際に確認した。

- `git diff develop...HEAD -- src/libpdftl/tl_dense_symmetric_matrix_object.cc src/libpdftl/tl_dense_general_matrix_object.cc`が空(本体の`load`は変わっていない)。
- `env -u PDF_HOME devtool/check.sh --build-dir build-review`(新しいビルドディレクトリ): whitespace・clang-format・build・warnings・testsすべてPASS(ctestの`xtest`・`xtest.mpi`とも100%)。
- `PDF_HOME=/nonexistent`で`pdf-xtest --gtest_filter='TlDenseSymmetricMatrix_Lapack.multiplication_*'`を実行すると、2件が異常終了せずに`TlFile::isExistFile`のアサーションでFAILする。

**残っている点(別タスクの候補)**: `TlDenseSymmetricMatrixObject::load`・`TlDenseGeneralMatrixObject::load`は、ファイルが開けない・形式が不正な場合にcatch節の外の`throw;`で`std::terminate`する。エラーで止まるという点では安全だが、適切な例外に置き換えるかどうかは、本体の多数の呼び出し元に関わるので別途判断する。
