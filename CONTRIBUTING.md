# Contributing

## ブランチモデル(git-flow)

- **`main`**: リリース済みの状態だけを反映する。直接コミットしない。`release/*`・`hotfix/*`のマージだけを受け取り、マージごとにバージョンタグ(例: `2025.3.0`、`v`なし)を打つ。
- **`develop`**: 日常開発の統合ブランチ。
- **`feature/*`・`fix/*`・`chore/*`**: `develop`から切り、レビュー承認後に`develop`へ`--no-ff`でマージする。
- **`release/YYYY.M.PATCH`**: `develop`から切ってリリースを準備し、`main`(タグ付き)と`develop`へマージする。
- **`hotfix/YYYY.M.PATCH`**: `main`から切って緊急修正し、`main`(タグ付き)と`develop`へマージする。

バージョンは`YYYY.M.PATCH`(例: `2025.3.0`)で、`CMakeLists.txt`の`PROJECT_VERSION_MAJOR/MINOR/REVISION`に書く。

作業ブランチは`develop`から切り(`hotfix/*`は`main`から)、`--no-ff`でマージする。作業ごとに専用のgit worktreeを使うと、ブランチの切り替え合いによる変更の混在を避けられる。

```bash
git worktree add -b feature/foo ../feature-foo develop   # developから切ってworktreeを作る
git switch develop && git merge --no-ff feature/foo      # レビュー後にdevelopへマージ
git worktree remove ../feature-foo && git branch -d feature/foo
```

`main`へのマージとタグ付け(`release/*`・`hotfix/*`)は、メンテナーが行う。

## 完了の確認

テストと整形チェックには GoogleTest と clang-format が必要。初回に次を実行する(Linuxはapt+sudo、macOSはHomebrew)。

```bash
devtool/setup-dev-tools.sh
```

ビルドには`data/pdfkwd.yaml`から`pdfkwd.mpac`を生成するため`ProteinDF_bridge`の`yml2mpac`が必要です(`regress-venv`等にインストールしてPATHを通すか`-DYML2MPAC_EXECUTABLE`で指定)。

```bash
devtool/check.sh   # 差分の空白エラー、変更行のclang-format、ビルド、ctest
devtool/regress.sh # 計算結果の回帰テスト(ProteinDF_testのserial_devスイート)
```

`regress.sh`の実行には隣接リポジトリ(`ProteinDF_test`・`ProteinDF_pytools`・`ProteinDF_bridge`)が必要。`$(git rev-parse --git-common-dir)/regress.conf`(または環境変数)で場所を指定する。
特定のテストのみを実行する場合は `--entries` を指定する(例: `devtool/regress.sh --entries O2_UB3LYP,N2_UB3LYP`)。
