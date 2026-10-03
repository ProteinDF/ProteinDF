# TASK: 対称行列のloadで、未対応の形式・正方でない行列を例外にする

- **Branch**: `fix/symmetric-load-format`
- **作成**: 2026-10-03(Claude)
- **状態**: レビュー中(収束、マージ承認待ち)

> `TlDenseSymmetricMatrixObject::load`(`src/libpdftl/tl_dense_symmetric_matrix_object.cc`)には、失敗しても止まらない経路が2つ残っている(`TASK_matrix-load-exception.md`のレビュー結果1回目で発見、ユーザー依頼 2026-10-03)。
>
> 1. **未対応の形式**: ヘッダーの`matrixType`が`RLHD`以外(`RSFD`・`CSFD`・`RUHD`など、`tl_matrix_object.h`の`MatrixType`)の場合、`switch`の`default:`でログを出すだけで、**`resize(row)`で0埋めされた行列のまま`false`を返す**。本体の`load`の呼び出しのほとんどは戻り値を確認しないので、0の行列で計算が黙って続く。
> 2. **正方でない行列**: `row != col`の場合も、ログ(`illegal format`)を出すだけで`resize(row)`して先に進む。
>
> **Claudeが確認した前提**:
> - ファイルが開けない場合・ヘッダーが不正な場合は、`TASK_matrix-load-exception.md`ですでに`std::runtime_error`を投げるようになっている。この2つの経路も同じ扱いにそろえる。
> - ヘッダーの検査(`TlMatrixUtils::getHeaderInfo`)はファイルサイズが行列のサイズと一致するかまで見ているので、ファイルが途中で切れている場合はすでに例外になる。
> - `TASK_matrix-load-exception.md`の調査では、戻り値を確認して「`false`なら別の処理をする」呼び出し元はなかった。
> - **他の形式からの変換はしない。** 例えば正方の`RSFD`ファイルを対称行列として読む(下三角を取る)ことは、今回はしない(実際に必要になった場合に別途判断する)。

## 役割分担・ブランチ運用(MUST)

- 実装はagy / GitHub Copilot、レビューはClaudeが担当する。`AGENTS.md`を必ず読むこと。
- 上記のブランチの専用worktreeで作業する。`develop`へは自分でマージしない。
- この指示書は編集しない(レビュー結果はClaudeが追記する)。
- `pdf-xtest`をリポジトリのルートで直接実行すると`output.log`が作られる。終わったら削除すること。

## 対象

1. `row != col`の場合は、`resize`する前に`std::runtime_error`を投げる。メッセージにはパスと行数・列数を含める。
2. `default:`(`RLHD`以外の形式)では`std::runtime_error`を投げる。メッセージにはパスと`matrixType`の値を含める。
3. これで`false`を返す経路がなくなるなら、ヘッダーのdocコメント(`tl_dense_symmetric_matrix_object.h`)を`TlDenseGeneralMatrixObject::load`と同じ書き方(「成功時はtrue。失敗時は例外を投げるため、falseを返す経路はない」)に直す。戻り値の型`bool`は変えない。
4. **回帰の調査**: この変更で、今まで(0の行列のまま)動いていた本体のコードが例外で止まるようになる可能性がある。次を調べて、結果を一覧にして完了報告に書く。
   - `src/libpdf`・`src/pdf`・`src/tools`で、対称行列の型(`TlDenseSymmetricMatrix_*`)に`load`しているファイルのパスのうち、一般行列の型(`TlDenseGeneralMatrix_*`)で`save`しているパスと同じもの(同じパス関数、例: `DfObject::get*Path`)があるか。
   - あれば、その箇所と、一般行列として保存したものを対称行列として読む意図があるかどうか(コードとコメントから読み取れる範囲で)。
   - 該当する箇所があった場合は、修正せずに「判断が必要な点」に書いて止まってよい。

## 完了の定義

1. `src/unit_test`の`dense_symmetric_matrix_test_template.h`(`TASK_matrix-load-exception.md`で`throwsOnNonExistentFile`・`throwsOnCorruptedFile`を追加したファイル)に、次のテストを追加する。
   - 一般行列(正方、`RSFD`)として`save`したファイルを対称行列で`load`すると例外になる。
   - 一般行列(正方でない)として`save`したファイルを対称行列で`load`すると例外になる。
   - いずれも、例外のあとで対称行列が0埋めの行列に変わっていない(`load`前の状態のまま)かどうかを確認し、どちらになるかを完了報告に書く(どちらでもよいが、挙動を明らかにする)。
2. 既存の対称行列の`save` → `load`のテスト(`doesSaveAndLoad`など)が従来どおりPASSする。
3. `devtool/check.sh`が通る(clang-formatを含む)。
4. `AGENTS.md`の形式で完了報告を出す。

## レビュー結果(1回目、2026-10-04、収束)

`fix/symmetric-load-format`(`1d492da`・`6187c3b`)をレビューした(agyは利用上限で一度中断し、2026-10-04 02:15に`delegate.sh -c`で再開した)。`row != col`は`resize`の前に、`RLHD`以外の形式は`default:`で、それぞれパスと値を含む`std::runtime_error`を投げるようになり、`false`を返す経路はなくなった。docコメントも一般行列と同じ書き方になった。Claudeが次を実際に確認した。

- `env -u PDF_HOME devtool/check.sh --build-dir build-review`(新しいビルドディレクトリ): whitespace・clang-format・build・warnings・testsすべてPASS(ctestの`xtest`・`xtest.mpi`とも100%)。
- 追加テスト(`throwsOnGeneralSquareMatrixFile`・`throwsOnGeneralNonSquareMatrixFile`)はLapack・Eigen・Eigen_FP32の対称行列すべてで実行される。
- 回帰の調査(一般行列として保存したパスを対称行列で読む箇所は0件)のうち、`src/tools/main_mat_show.cpp`がヘッダーの`matrixType`で対称行列・一般行列を振り分けてから`load`していることを読んで確かめた。

**残っている軽微な点(対応不要)**: 正方の一般行列のファイルを読んだ場合は、`resize(row)`のあとで例外を投げるので、行列のサイズは`load`前から変わる(既存の要素は保持され、広がった部分は0)。正方でない場合は`resize`の前に投げるので変わらない。例外で止まるので黙って計算が続くことはないが、失敗時に状態を変えないようにしたい場合は、形式の確認を`resize`の前に移せばよい。
