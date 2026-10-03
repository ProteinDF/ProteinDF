# TASK: ScaLAPACK版のcatch節の外の`throw;`(3か所)を例外にする

- **Branch**: `fix/scalapack-throw`
- **作成**: 2026-10-03(Claude)
- **状態**: マージ済み(`f7d610f`)

> `TASK_matrix-load-exception.md`と同じ問題(catch節の外の`throw;`は再送出する例外がないので`std::terminate`で強制終了する)が、ScaLAPACK版の次の3か所にある。これを理由を持った例外を投げるように直す(ユーザー依頼、2026-10-03)。
>
> | ファイル | 場所 | 条件 |
> |---|---|---|
> | `src/libpdftl/tl_dense_general_matrix_impl_scalapack.cc` | `TlDenseGeneralMatrix_ImplScalapack::inverse()` | `pdgetri_`の`info != 0` |
> | 同上 | 同上 | `pdgetrf_`の`info != 0`(特異行列など) |
> | `src/libpdftl/tl_dense_symmetric_matrix_impl_scalapack.cc` | `TlDenseSymmetricMatrix_ImplScalapack(const TlDenseGeneralMatrix_ImplScalapack&)` | 行数と列数が違う |
>
> **Claudeが確認した前提**:
> - 本体に例外を握りつぶす`catch`はない(`TASK_matrix-load-exception.md`参照)。例外にしても、捕まえなければ今までどおりプログラムは止まり、理由が表示されるようになる。
> - `pdgetrf_`・`pdgetri_`の`info`と行列の次元は全プロセスで同じ値になるので、全プロセスが同じように例外を投げる。プロセスごとに挙動が分かれる心配はない(違うと考える根拠があれば完了報告に書くこと)。
> - `inverse()`は`IPIV`・`WORK`・`IWORK`を`new[]`で確保していて、いまの`throw;`の経路では解放されない。例外を捕まえられるようになるとメモリリークが問題になる。
> - MPI版のテスト(`xtest.mpi`、`src/unit_test/CMakeLists.txt`の`PDF_XTEST_MPI_SOURCES`)に`TlDenseGeneralMatrix_Scalapack.inverse`がある(`tl_dense_general_matrix_scalapack_test.cc:428`)。

## 役割分担・ブランチ運用(MUST)

- 実装はagy / GitHub Copilot、レビューはClaudeが担当する。`AGENTS.md`を必ず読むこと。
- 上記のブランチの専用worktreeで作業する。`develop`へは自分でマージしない。
- この指示書は編集しない(レビュー結果はClaudeが追記する)。
- **`fix/matrix-load-exception`(`TASK_matrix-load-exception.md`)がdevelopにマージされてから着手する。** 例外の型(`std::runtime_error`か、その派生クラスか)は、そちらで決めたものに合わせる。

## 対象

1. 上記3か所の`throw;`を、`TASK_matrix-load-exception.md`で決めた例外の型を投げるように変える。メッセージには、関数名と理由(`pdgetrf_`/`pdgetri_`と`info`の値、または行数・列数)を含める。既存の出力(`std::cout`・`log_.critical`)は残してよい。
2. `inverse()`の`IPIV`・`WORK`・`IWORK`・`WORK_SIZE`・`IWORK_SIZE`を`std::vector`に置き換え、例外の経路でもリークしないようにする。計算の手順(`pdgetrf_` → ワークサイズの問い合わせ → `pdgetri_`)は変えない。
3. 上記以外の`throw;`(catch節の中の正しい再送出)は変更しない。

## 完了の定義

1. `src/unit_test`のMPI版テストに次を追加する(既存の`tl_dense_general_matrix_scalapack_test.cc`・`tl_dense_symmetric_matrix_scalapack_test.cc`の書き方に合わせる)。
   - 特異行列(例: すべて0の行がある行列)の`inverse()`が例外になる(`EXPECT_THROW`)。
   - 行数と列数が違う一般行列から対称行列を作ると例外になる。
2. 既存の`TlDenseGeneralMatrix_Scalapack.inverse`が従来どおりPASSする。
3. `xtest.mpi`を、既定のプロセス数(4)で実行してPASSする。テストの途中で一部のプロセスだけが止まる(ハングする)ことがないことを確認する。
4. `devtool/check.sh`が通る(clang-formatを含む)。
5. `AGENTS.md`の形式で完了報告を出す。

## レビュー結果(1回目、2026-10-03、収束)

`fix/scalapack-throw`(`df2e8be`・`90484d3`)をレビューした。対象の3か所だけが`std::runtime_error`(関数名と`info`の値、または行数・列数を含む)に置き換わり、`inverse()`の作業配列はすべて`std::vector`になった。計算の手順(`pdgetrf_` → ワークサイズの問い合わせ → `pdgetri_`)は変わっていない。catch節の中の`throw;`は変更されていない。Claudeが次を実際に確認した。

- `env -u PDF_HOME devtool/check.sh --build-dir build-review`(新しいビルドディレクトリ): whitespace・clang-format・build・warnings・testsすべてPASS(ctestの`xtest`・`xtest.mpi`とも100%)。
- `mpirun -np 4 pdf-xtest.MPI --gtest_filter='*throws*:TlDenseGeneralMatrix_Scalapack.inverse'`: `mpirun`の終了コード0。4プロセスすべてで3件ともOK(`[  PASSED  ] 3 tests`が4回、FAILEDは0件)。一部のプロセスだけが止まることはなかった。
- `main_MPI.cpp`は各プロセスが`RUN_ALL_TESTS()`の結果をそれぞれ返すので、0番以外のプロセスの失敗も`mpirun`の終了コードに現れる(ctestで検出される)。
- `pdgetrf_`の`INFO`はScaLAPACKの仕様で全プロセス共通の値(global output)なので、全プロセスが同じように例外を投げる。

**補足**: このマシンで`mpirun`を実行すると、GPUドライバの`HSA exception: Agent creation failed.`という警告が大量に出るが、テストとは関係ない。

ユーザー承認のうえ、2026-10-03にdevelopへマージした(`f7d610f`)。
