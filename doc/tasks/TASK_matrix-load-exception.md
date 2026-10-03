# TASK: 行列ファイル読み込み(load)の失敗を例外にする

- **Branch**: `fix/matrix-load-exception`
- **作成**: 2026-10-03(Claude)
- **状態**: レビュー中(収束、マージ承認待ち)

> `TlDenseGeneralMatrixObject::load`(`src/libpdftl/tl_dense_general_matrix_object.cc`、3か所)と`TlDenseSymmetricMatrixObject::load`(`src/libpdftl/tl_dense_symmetric_matrix_object.cc`、2か所)は、ファイルが開けない・形式が不正・未対応の形式の場合に、ログを出したあと**catch節の外で`throw;`**している。catch節の外の`throw;`は再送出する例外がないので、`std::terminate`が呼ばれてプログラムが強制終了する(例外は投げられない)。これを、理由を持った例外を投げるように直す(ユーザー依頼、2026-10-03)。
>
> **Claudeが確認した前提**:
> - 本体(`src/libpdf`・`src/libpdftl`・`src/pdf`・`src/tools`)に`catch`は4か所しかなく、すべて`bad_alloc`または`catch (...)`での再送出である。例外を握りつぶす箇所はない。したがって、例外にしても、捕まえなければ今までどおりプログラムは止まる(`terminate called after throwing ... what(): <理由>`のように理由が出るようになる)。
> - 本体で`load`の戻り値を確認しているのは約13か所、無視しているのは約269か所(行列以外の`load`を含む概数)。`fix/xtest-pdf-home`のレビューで、`return false`にする案は「黙って空の行列で計算を続ける」ことになるため却下した(`TASK_xtest-pdf-home.md`のレビュー結果1回目)。**`return false`にはしないこと。**
> - 今回の対象外だが、同じく「catch節の外の`throw;`」が`tl_dense_general_matrix_impl_scalapack.cc`(逆行列の`pdgetri_`・`pdgetrf_`失敗時、2か所)と`tl_dense_symmetric_matrix_impl_scalapack.cc`(コンストラクタの次元不一致、1か所)にある。catch節の中の`throw;`(`tl_dense_general_matrix_impl_lapack.cc`の2か所、`tl_dense_vector_impl_lapack.cc`の1か所)は正しい再送出なので触らない。

## 役割分担・ブランチ運用(MUST)

- 実装はagy / GitHub Copilot、レビューはClaudeが担当する。`AGENTS.md`を必ず読むこと。
- 上記のブランチの専用worktreeで作業する。`develop`へは自分でマージしない。
- この指示書は編集しない(レビュー結果はClaudeが追記する)。

## 対象

1. 上記2つの`load`の5か所の`throw;`を、`std::runtime_error`(またはその派生クラス)を投げるように変える。メッセージには、ファイルパスと理由(開けない・形式が不正・未対応の形式)を含める。既存のログ出力(`log_.critical`)は残してよい。
   - 派生クラスを新しく作るかどうかは実装担当が決めてよい。作る場合は`src/libpdftl`に置き、理由を完了報告に書く。
2. `load`のdocコメント(ヘッダー)に、失敗時に例外を投げることと、戻り値`bool`が何を意味するか(例外以外で`false`を返す経路が残るかどうか)を書く。**例外を投げる経路と`false`を返す経路の両方が残る場合は、どの条件でどちらになるかを完了報告に表で書くこと。**
3. `load`を呼んで戻り値を確認している本体の箇所(約13か所)について、例外になったことで挙動が変わる箇所(例: 「ファイルがなければ`false`で別の処理」を前提にしている箇所)がないか調べ、結果を一覧にして完了報告に書く。挙動が変わる箇所があった場合は、修正せずに「判断が必要な点」に書いて止まってよい。
4. 上記以外の`throw;`(上の前提に書いたscalapackの3か所)は変更しない。

## 完了の定義

1. `src/unit_test`に、`load`の失敗が例外になることのテストを追加する(既存のテストファイルの置き方に合わせる)。
   - 存在しないファイル → 例外(`EXPECT_THROW`)。general・symmetricの両方。
   - 形式が不正なファイル(例: テスト内で一時ファイルに数バイトのゴミを書く)→ 例外。general・symmetricの両方。
   - 正しいファイルは従来どおり読める(既存のテストで確認されていればそれでよい。どのテストかを完了報告に書く)。
2. 既存の`TlDenseSymmetricMatrix_Lapack.multiplication_MV`・`multiplication_VMV`(読み込み前に`TlFile::isExistFile`で存在を確認している)は変更しなくてよい。
3. 存在しないファイルを`load`する小さなプログラム(またはテストの外での実行)で、プログラムが`terminate`するときに理由のメッセージが表示されることを確認し、その出力を完了報告に貼る(テスト用のコードはコミットしなくてよい)。
4. `devtool/check.sh`が通る(clang-formatを含む)。
5. `AGENTS.md`の形式で完了報告を出す。

## レビュー結果(1回目、2026-10-03、収束)

`fix/matrix-load-exception`(`592fca5`・`5ad74a1`)をレビューした。対象の5か所の`throw;`だけが`std::runtime_error`(パスと理由を含む)に置き換わっており、`return false`にはなっていない。scalapackの3か所は変更されていない。派生クラスを作らない判断も妥当。Claudeが次を実際に確認した。

- `env -u PDF_HOME devtool/check.sh --build-dir build-review`(新しいビルドディレクトリ): whitespace・clang-format・build・warnings・testsすべてPASS(ctestの`xtest`・`xtest.mpi`とも100%)。
- 追加テスト(`throwsOnNonExistentFile`・`throwsOnCorruptedFile`)はtyped testのテンプレートに追加され、Lapack・Eigen・Eigen_FP32の一般行列・対称行列すべてで実行される。形式不正のテストは一時ファイルを作り、終わったら消している。
- 呼び出し元の調査結果のうち`TlMatrixCache.h:235`(事前に`isExistFile`で確認)を読んで確かめた。

**残っている点(別タスクの候補)**:
1. `TlDenseSymmetricMatrixObject::load`の`default:`(RLHD以外の形式)は、例外を投げず、**`resize(row)`で0埋めされた行列のまま`false`を返す**(もともとある経路で、今回の対象外)。戻り値を確認しない呼び出し元では、0の行列で計算が黙って続く。同じ関数の`row != col`の場合も、ログを出すだけで処理を続ける。
2. 存在しないファイルを読むと、メッセージは「cannot open matrix file」ではなく「illegal matrix format」になる(ヘッダーの読み取りで先に失敗するため)。誤解を招くが、止まること自体は正しい。
3. `pdf-xtest`をリポジトリのルートで直接実行すると、`TlLogging`の既定のログファイル`output.log`がそこに作られる(このブランチのworktreeにも残っていた)。
