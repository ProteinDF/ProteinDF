# 変更履歴

## 2026.10.0

### 計算結果が変わる変更

- **UKS の GGA(BLYP・B3LYP など)**: 交差スピン項で、逆スピンの密度勾配を使うように修正した。UKS の BLYP・B3LYP の全エネルギーが変わる(閉殻の計算には影響しない)。
- **B3LYP の局所相関(VWN)**: Gaussian・libxc と同じ定義(VWN の RPA パラメータ、von Barth–Hedin の一次補間)にした。開殻の B3LYP の全エネルギーが変わる(O2 で約 3×10⁻⁴ Hartree、エネルギーが上がる方向)。閉殻の計算は変わらない。旧定義は残していない。
- **Becke88 の grid-free 法**: β スピンの微分に誤ってα スピンの量を使っていたのを修正した(UKS・ROKS の grid-free 法のみ)。

### 入力の変更

- `xc-potential` を省略すると、warning を出して `b3lyp` を使う(従来は `unknown XC functional` で異常終了した)。
- `xc_functional` が空のときに未定義動作になる問題を修正した。

### ビルドの変更

- 既定の C++ 標準を C++17 にした。GCC 14 以降でビルドできるようにした。
- `data/pdfkwd.mpac` をビルド時に `data/pdfkwd.yaml` から生成するようにした。**ビルドに `ProteinDF_bridge` の `yml2mpac` が必要になった**(PATH に置くか、CMake の `-DYML2MPAC_EXECUTABLE=` で指定する)。

### そのほか

- 行列ファイルの読み込みに失敗したときや、ScaLAPACK の行列操作でエラーが起きたときは、異常終了せずに `std::runtime_error` を投げるようにした。
- LAPACK 版の対称行列が、必要な約 2 倍のメモリを確保していた問題を修正した(N=10000 で約 780 MB → 約 400 MB)。
- 開発用ツールを追加した(`devtool/`: 差分の検査、回帰テスト、PySCF との比較など。使い方は `CONTRIBUTING.md`)。
