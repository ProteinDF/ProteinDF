# TASK: 既定のC++規格をC++14からC++17に上げる

- **Branch**: `chore/cxx17`
- **作成**: 2026-10-03(Claude)
- **状態**: 未着手

> `CMakeLists.txt`は`CMAKE_CXX_STANDARD`の既定を14にしている。GoogleTest 1.13以降はC++17以上を要求するため(Ubuntu 26.04の`libgtest-dev`は1.17)、このままでは`src/unit_test`がビルドできない(`gtest-port.h:273: #error C++ versions less than C++17 are not supported.`)。既定をC++17に上げる(ユーザー承認済み、2026-10-03)。
>
> Claudeが事前に確認したこと: `-DCMAKE_CXX_STANDARD=17`を指定すると、GCC 15.2で本体(`PDF.x`・`PPDF.x`)とテスト(`pdf-xtest`・`pdf-xtest.MPI`)がエラーなくビルドできる。`xtest`は316件中2件が失敗する(`TlDenseSymmetricMatrix_Lapack.multiplication_MV`・`multiplication_VMV`、例外`basic_string: construction from null is not valid`)。この2件は今回のスコープ外で、別のTASKで扱う。

## 役割分担・ブランチ運用(MUST)

- 実装はagy / GitHub Copilot、レビューはClaudeが担当する。`AGENTS.md`を必ず読むこと。
- 上記のブランチの専用worktreeで作業する。`develop`へは自分でマージしない。
- この指示書は編集しない(レビュー結果はClaudeが追記する)。

## 前提

GTestとclang-formatがインストールされていること(`devtool/setup-dev-tools.sh`、ユーザーが実行する)。`check.sh`のtestsがSKIPになる場合は、インストールされていない。その場合は完了の定義2・3の代わりに、その旨を完了報告に書く。

## 対象

1. `CMakeLists.txt`の`CMAKE_CXX_STANDARD`の既定を`17`にする。あわせて`CMAKE_CXX_STANDARD_REQUIRED ON`と`CMAKE_CXX_EXTENSIONS OFF`を設定する。ただし`-DCMAKE_CXX_STANDARD=...`で外から指定された値は、これまでどおり優先されるようにする。
2. README・ドキュメント(`doc/source`)にC++の要件が書かれていれば、C++17に更新する(書かれていなければ追加しなくてよい)。
3. C++17に上げたことで出るようになった警告・エラーがあれば直す。ただし`xtest`の上記2件の失敗は直さない。
4. `.vscode/c_cpp_properties.json`の`cppStandard`(現在`c++11`)は、`c++17`にしてよい(任意)。

## 完了の定義

1. 新しいビルドディレクトリで`devtool/check.sh --build-dir build-cxx17`を実行し、buildがPASSになる。configureのログ(`build-cxx17/check-configure.log`)にGTestが見つかったことが出ている。
2. `ctest`で`xtest.mpi`がPASSし、`xtest`の失敗が上記の2件だけである(gtestの出力の`[  FAILED  ]`の行をそのまま完了報告に貼る)。
3. ビルドログの警告の総数を、C++14(`develop`、`-DCMAKE_CXX_STANDARD=14`)とC++17(このブランチ)で比べ、両方の数と、C++17で増えた警告の種類を完了報告に書く(どちらも新しいビルドディレクトリで最初からビルドすること)。
4. `-DCMAKE_CXX_STANDARD=14`を指定すると、C++14でビルドされること(外からの指定が優先されること)を確認する(`compile_commands.json`やビルドログの`-std=`を確認する)。
5. `AGENTS.md`の形式で完了報告を出す。
