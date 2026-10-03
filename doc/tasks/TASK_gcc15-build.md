# TASK: GCC 14以降でビルドできない問題の修正

- **Branch**: `fix/gcc15-build`
- **作成**: 2026-10-03(Claude)
- **状態**: マージ済み(`0e6357a`)

> `develop`(`2025.3.0`)は、GCC 15.2(Ubuntu 15.2.0-16ubuntu1)で`libpdf`のビルドに失敗し、`pdf`本体までビルドできない。原因は`src/libpdf/df_population_tmpl.h`のクラステンプレート`DfPopulation_tmpl`が、自分(と基底クラス`DfObject`)にないメンバ関数を呼んでいること。GCC 14から、インスタンス化されないテンプレート本体のこの種のエラーが`-Wtemplate-body`として既定でエラーになった。同名の関数は非テンプレート版の`src/libpdf/DfPopulation.h`にあり、テンプレート化の途中で移し忘れたものと思われる。

## 役割分担・ブランチ運用(MUST)

- 実装はagy / GitHub Copilot、レビューはClaudeが担当する。`AGENTS.md`を必ず読むこと。
- 上記のブランチの専用worktreeで作業する。`develop`へは自分でマージしない。
- この指示書は編集しない(レビュー結果はClaudeが追記する)。

## 再現

```bash
devtool/check.sh
```

Claudeが確認したエラー(`cmake --build -- -k`で全件を出したもの。重複を除くとこの3件):

```text
src/libpdf/df_population_tmpl.h:111:15: error: 'class DfPopulation_tmpl<SymmetricMatrix, Vector>' has no member named 'getGrossAtomPop'; did you mean 'getGrossOrbPop'? [-Wtemplate-body]
src/libpdf/df_population_tmpl.h:157:11: error: 'class DfPopulation_tmpl<SymmetricMatrix, Vector>' has no member named 'calcPop' [-Wtemplate-body]
src/libpdf/df_population_tmpl.h:159:39: error: 'class DfPopulation_tmpl<SymmetricMatrix, Vector>' has no member named 'getSumOfNucleiCharges' [-Wtemplate-body]
```

このヘッダーを含む`df_population_{lapack,eigen,scalapack}.cc`、`DfInitialGuess.cpp`、`DfInitialGuessHarris.cpp`、`DfInitialGuessHarris_Parallel.cpp`のコンパイルが失敗する。

## 対象

1. `DfPopulation_tmpl`に、足りないメンバ関数(`getGrossAtomPop`、`calcPop`、`getSumOfNucleiCharges`)を、`DfPopulation.h`の実装をもとに追加する。`DfPopulation`の派生クラス(`DfPopulation_Parallel`など)が`calcPop`を`virtual`で上書きしている点にも注意し、テンプレート版で同じ扱いが必要かどうかを調べて完了報告に書く。
2. `-Wno-template-body`や`-fpermissive`でエラーを抑制して済ませないこと。
3. 上記3件を直したあとに別のコンパイルエラー・リンクエラーが出た場合は、同じ方針で直してよい。ただし、`df_population_tmpl.h`以外の修正が5ファイルを超える、または計算結果が変わりうる修正になる場合は、そこで止めて完了報告の「判断が必要な点」に書くこと。

## 完了の定義

1. `devtool/check.sh`の`build`がPASSになる(`pdf`本体を含む全ターゲットがビルドできる)。
2. 追加したメンバ関数が、`DfPopulation.h`の対応する関数と同じ計算をしていることを、完了報告で関数ごとに対応づけて説明する(どの行をどう移したか)。
3. ビルドの警告のうち、変更したファイルで新たに出たものがない(`check.sh`の`warnings`)。
4. `devtool/check.sh`が通る(GTest・clang-formatが未インストールのためのSKIPは可)。
5. `AGENTS.md`の形式で完了報告を出す。

## レビュー結果(1回目、2026-10-03、収束)

`fix/gcc15-build`(`3d78d56`)をレビューした。変更は`src/libpdf/df_population_tmpl.h`の1ファイル(+82行)のみで、警告の抑制はしていない。Claudeが次を実際に確認した。

- `devtool/check.sh --build-dir build-review`(新しいビルドディレクトリで最初からビルド): whitespace PASS、build PASS(`PDF.x`・`PPDF.x`まで生成)、変更ファイルの警告0件。clang-format・testsは未インストールのためSKIP。
- 追加された`setNucleiCharges`・`calcPop`・`getGrossAtomPop`・`getSumOfNucleiCharges`は、`DfPopulation.cpp`・`DfPopulation.h`の対応する関数と計算が同じ(違いは、行列型をクラスのテンプレート引数から取ることと、OpenMPのcritical名だけ)。
- コンストラクタで`setNucleiCharges()`を呼ぶ変更(未初期化だった`nucleiCharges_`の初期化)について: 既存の実行経路で`DfPopulation_Lapack`・`DfPopulation_Eigen`を使っているのは初期推定(`DfInitialGuess.h`・`DfInitialGuessHarris.h`)の`getSumOfElectrons()`だけで、`nucleiCharges_`を使わない。増えるのは`Fl_Geometry`の構築1回だけで、既存の計算結果は変わらない。

**残っている軽微な点(対応不要)**:
1. `calcPop`は`virtual`だが、テンプレート版の派生クラスで上書きしているものはない(`DfPopulation`との一貫性のためで、害はない)。
2. 今回コンパイルされるようになった`getReport`・`getAtomPopulation`は、現時点でテンプレート版からは呼ばれていない。テストがないため、これらの数値は未検証である。テンプレート版を実際に使うようになったときに、`DfPopulation`の出力と比較すること。

ユーザー承認のうえ、2026-10-03にdevelopへマージした(`0e6357a`)。
