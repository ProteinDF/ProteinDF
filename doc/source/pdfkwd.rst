* pdf_param_path

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "pdf_param_path"
  "説明", "パラメータファイルのパス"
  "デフォルト値", "pdfparam.mpac"
  "形式", ".*"


* step_control

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "step_control"
  "説明", "計算ステップの指定"
  "デフォルト値", "create"
  "形式", "((create|integral|guess|scf|force)\s*)+"


* comment

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "comment"
  "説明", "コメント"
  "デフォルト値", ""
  "形式", ".*"


* linear_algebra_package

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "linear_algebra_package"
  "説明", "線形ライブラリの指定"
  "デフォルト値", "lapack"
  "形式", "(lapack|scalapack)"


* scalapack_block_size

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scalapack_block_size"
  "説明", "ScaLAPACKブロックサイズの指定"
  "デフォルト値", "64"
  "形式", "(\d+)"


* save_distributed_matrix_to_local_disk

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "save_distributed_matrix_to_local_disk"
  "説明", "大域分散行列を各ノードで分散して保存"
  "デフォルト値", "no"
  "形式", "(yes|no)"


* local_disk_path

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "local_disk_path"
  "説明", "ローカルディスクへの保存に使用するパス"
  "デフォルト値", "/tmp"
  "形式", ".*"


* memory_size

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "memory_size"
  "説明", "使用可能なメモリサイズ"
  "デフォルト値", "1GB"
  "形式", "\d+\s*(MB|GB)"


* use_mapfile

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "use_mapfile"
  "説明", "仮想ディスクを使用する"
  "デフォルト値", "no"
  "形式", "(yes|no)"


* mapfile_size

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "mapfile_size"
  "説明", "仮想ディスクのサイズ"
  "デフォルト値", "auto"
  "形式", "(\d+\s*(MB|GB)|auto)"


* mapfile_basename

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "mapfile_basename"
  "説明", "仮想ディスクファイルの基本ファイル名"
  "デフォルト値", "pdf"
  "形式", ".*"


* work_on_disk

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "work_on_disk"
  "説明", "いつでも仮想ディスクを使用する"
  "デフォルト値", "no"
  "形式", "(yes|no)"


* parallel_processing_type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "parallel_processing_type"
  "説明", "並列アルゴリズムの選択"
  "デフォルト値", "divide_and_conquer"
  "形式", "(divide_and_conquer|DC|master_slave|MS)"


* cleanup

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "cleanup"
  "説明", "中間ファイルの削除"
  "デフォルト値", "no"
  "形式", "(yes|no)"


* show_keyword

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "show_keyword"
  "説明", "キーワードを表示する"
  "デフォルト値", "no"
  "形式", "(yes|no)"


* show_input

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "show_input"
  "説明", "入力パラメータを表示する"
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* show_coordinates

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "show_coordinates"
  "説明", "分子座標を表示する"
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* show_orbital_basis

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "show_orbital_basis"
  "説明", "基底関数情報の表示形式の指定"
  "デフォルト値", "gamess"
  "形式", "(gamess|amoss|none)"


* guess

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "guess"
  "説明", "guessの指定"
  "デフォルト値", "rho"
  "形式", "(rho|file_rho|lcao|density|core|huckel|harris)"


* cut-value

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "cut-value"
  "説明", "積分カットオフ値"
  "デフォルト値", "1.0E-10"
  "形式", "(double > 0)"


* block-size

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "block-size"
  "説明", "積分ブロックサイズ"
  "デフォルト値", "1024000"
  "形式", "(integer > 0)"


* charge-extrapolate-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "charge-extrapolate-number"
  "説明", "ダミー電荷を段階的に外挿する回数"
  "デフォルト値", "1"
  "形式", "(integer >= 1)"


* orbital-overlap-correspondence

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "orbital-overlap-correspondence"
  "説明", "軌道の対応付けを行う"
  "デフォルト値", "off"
  "形式", "on/off"


* orbital-overlap-correspondence-first

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "orbital-overlap-correspondence-first"
  "説明", "軌道の対応付けを初回から行う"
  "デフォルト値", "off"
  "形式", "on/off"


* orbital-overlap-correspondence-method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "orbital-overlap-correspondence-method"
  "説明", "軌道対応付けの方法"
  "デフォルト値", "mo-overlap"
  "形式", "mo-overlap/mo-projection"


* orbital-overlap-correspondence-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "orbital-overlap-correspondence-number"
  "説明", ""
  "デフォルト値", "3"
  "形式", "integer"


* summary

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "summary"
  "説明", "計算サマリーを表示する"
  "デフォルト値", "none"
  "形式", "(none|convergence|every-scf)"


* analyze_population

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "analyze_population"
  "説明", "形式電荷の計算を行う"
  "デフォルト値", "none"
  "形式", "(none|convergence|every-scf)"


* max-iteration

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "max-iteration"
  "説明", "SCF繰り返し計算の最大値"
  "デフォルト値", "100"
  "形式", "(integer)"


* method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method"
  "説明", "計算方法"
  "デフォルト値", "rks"
  "形式", "(rks|uks|roks)"


* method/rks/occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/rks/occlevel"
  "説明", "rks計算における占有軌道"
  "デフォルト値", ""
  "形式", "(array of integer >= 0)"


* method/rks/electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/rks/electrons"
  "説明", "rks計算のおける電子数"
  "デフォルト値", ""
  "形式", "(integer >= 2)"


* method/uks/alpha_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/uks/alpha_occlevel"
  "説明", "uks計算におけるalpha電子の占有軌道"
  "デフォルト値", ""
  "形式", "(array of integer >= 0)"


* method/uks/alpha_electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/uks/alpha_electrons"
  "説明", "uks計算におけるalpha電子の数"
  "デフォルト値", ""
  "形式", "(integer >= 1)"


* method/uks/beta_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/uks/beta_occlevel"
  "説明", "uks計算におけるbeta電子の占有軌道"
  "デフォルト値", ""
  "形式", "(array of integer >= 0)"


* method/uks/beta_electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/uks/beta_electrons"
  "説明", "uks計算におけるbeta電子の数"
  "デフォルト値", ""
  "形式", "(integer >= 1)"


* method/roks/closed_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/roks/closed_occlevel"
  "説明", "roks計算における閉殻軌道の指定"
  "デフォルト値", ""
  "形式", "(array of integer >= 1)"


* method/uks/close_electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/uks/close_electrons"
  "説明", "roks計算における閉殻軌道の電子数"
  "デフォルト値", ""
  "形式", "(array of integer >= 0)"


* method/roks/open_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/roks/open_occlevel"
  "説明", "roks計算における開殻軌道の指定"
  "デフォルト値", ""
  "形式", "(array of integer >= 1)"


* method/uks/beta_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "method/uks/beta_occlevel"
  "説明", "uks計算におけるbeta電子の占有軌道"
  "デフォルト値", ""
  "形式", "(array of integer >= 0)"


* save_diff_density_matrix

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "save_diff_density_matrix"
  "説明", ""
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* use_matrix_cache

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "use_matrix_cache"
  "説明", ""
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* force_loading_from_disk

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "force_loading_from_disk"
  "説明", ""
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* show_cache_report

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "show_cache_report"
  "説明", "キャッシュ情報の表示"
  "デフォルト値", "no"
  "形式", "(yes|no)"


* disk-utilization

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "disk-utilization"
  "説明", ""
  "デフォルト値", "no"
  "形式", "yes/no"


* update_method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "update_method"
  "説明", ""
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* orbital-independence-threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "orbital-independence-threshold"
  "説明", ""
  "デフォルト値", "0.007"
  "形式", "(real >= 0)"


* convergence/type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "convergence/type"
  "説明", "収束の対象となるパラメータ"
  "デフォルト値", "density"
  "形式", "(fock|density|dcoef)"


* convergence/threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "convergence/threshold"
  "説明", "収束の閾値"
  "デフォルト値", "1e-3"
  "形式", "(real > 0)"


* convergence/threshold-energy

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "convergence/threshold-energy"
  "説明", "エネルギーにおける収束条件"
  "デフォルト値", "1e-4"
  "形式", "(real > 0)"


* scf-acceleration

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration"
  "説明", "収束加速法の選択"
  "デフォルト値", "damping"
  "形式", "(damping|anderson|diis)"


* scf-acceleration/damping/damping-factor

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration/damping/damping-factor"
  "説明", "damping法におけるダンピング係数"
  "デフォルト値", "0.85"
  "形式", "(real)"


* scf-acceleration/daming/start-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration/daming/start-number"
  "説明", "damping法を開始するイテレーション"
  "デフォルト値", "0"
  "形式", "(int)"


* scf-acceleration/damping/damping-type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration/damping/damping-type"
  "説明", "damping法を適用する物理量"
  "デフォルト値", "density"
  "形式", "(fock|density_matrix|density|dcoef)"


* scf-acceleration/anderson/damping-factor

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration/anderson/damping-factor"
  "説明", "anderson法におけるダンピング係数"
  "デフォルト値", "0.50"
  "形式", "(real)"


* scf-acceleration/anderson/start-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration/anderson/start-number"
  "説明", "anderson法を適用開始するイテレーション"
  "デフォルト値", "3"
  "形式", "(int)"


* level-shift

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "level-shift"
  "説明", "レベルシフト法を適用する"
  "デフォルト値", "no"
  "形式", "(yes|no)"


* level-shift/start-iteration

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "level-shift/start-iteration"
  "説明", "レベルシフト法を開始するイテレーション"
  "デフォルト値", "1"
  "形式", "(int)"


* level-shift/ls-closed-mo

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "level-shift/ls-closed-mo"
  "説明", "閉殻軌道に使用するレベルシフト値"
  "デフォルト値", "0.00"
  "形式", "(real)"


* level-shift/ls-open-mo

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "level-shift/ls-open-mo"
  "説明", "開殻軌道に使用するレベルシフト値"
  "デフォルト値", "0.00"
  "形式", "(real)"


* level-shift/ls-virtual-mo

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "level-shift/ls-virtual-mo"
  "説明", "空軌道に使用するレベルシフト値"
  "デフォルト値", "0.00"
  "形式", "(real)"


* scf-acceleration/diis/number-of-diis

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration/diis/number-of-diis"
  "説明", ""
  "デフォルト値", "3"
  "形式", "(integer)"


* scf-acceleration/diis/start-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration/diis/start-number"
  "説明", "DIIS法を開始するイテレーション"
  "デフォルト値", "3"
  "形式", "(integer >= 0)"


* scf-acceleration/diis/start-extrapolation

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-acceleration/diis/start-extrapolation"
  "説明", ""
  "デフォルト値", "6"
  "形式", "(integer >= 0)"


* xc-potential

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "xc-potential"
  "説明", "交換相関ポテンシャルの指定"
  "デフォルト値", "svwn~"
  "形式", "(svwn~|svwn|blyp|b3lyp)"


* grid_free

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "grid_free"
  "説明", "grid-free法を使用する"
  "デフォルト値", "false"
  "形式", "(yes|no)"


* xc-potential/grid-type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "xc-potential/grid-type"
  "説明", "使用するグリッドの選択"
  "デフォルト値", "sg-1"
  "形式", "fine/medium-fine/medium/coarse/sg-1"


* xc-update

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "xc-update"
  "説明", ""
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* xc-density-threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "xc-density-threshold"
  "説明", ""
  "デフォルト値", "1.0E-16"
  "形式", "real"


* TEI-integral-driven

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "TEI-integral-driven"
  "説明", ""
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* J_engine

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "J_engine"
  "説明", "クーロン項の計算方法の選択"
  "デフォルト値", "RI_J"
  "形式", "(conventional|RI_J|CD)"


* K_engine

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "K_engine"
  "説明", "Fock交換項の計算方法の選択"
  "デフォルト値", "conventional"
  "形式", "(conventional|RI_K|CD)"


* CDAM_tau

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "CDAM_tau"
  "説明", "CDAM法におけるτ値"
  "デフォルト値", "1.0E-10"
  "形式", "(real)"


* CD_epsilon

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "CD_epsilon"
  "説明", "CD法におけるε値"
  "デフォルト値", "1.0E-4"
  "形式", "(real)"


* scf-memory-saving

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "scf-memory-saving"
  "説明", ""
  "デフォルト値", "no"
  "形式", "yes/no"


* geometry/cartesian/input

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "geometry/cartesian/input"
  "説明", "カーテシアン座標による#分子座標の指定"
  "デフォルト値", "nil"
  "形式", "nil/stored"


* geometry/cartesian/unit

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "geometry/cartesian/unit"
  "説明", ""
  "デフォルト値", "a.u."
  "形式", "au/a.u./angstrom"


* basis-set/orbital

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "basis-set/orbital"
  "説明", ""
  "デフォルト値", "nil"
  "形式", "nil/stored"


* basis-set/density-auxiliary

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "basis-set/density-auxiliary"
  "説明", ""
  "デフォルト値", "nil"
  "形式", "nil/stored"


* basis-set/exchange-auxiliary

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "basis-set/exchange-auxiliary"
  "説明", ""
  "デフォルト値", "nil"
  "形式", "nil/stored"


* xc_density_threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "xc_density_threshold"
  "説明", ""
  "デフォルト値", "1.0E-16"
  "形式", "real"


* geometry

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "geometry"
  "説明", ""
  "デフォルト値", "cartesian"
  "形式", "cartesian/file"


* coordinates

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "coordinates"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* basis_sets

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "basis_sets"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* basis_sets_j

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "basis_sets_j"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* basis_sets_k

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "basis_sets_k"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* independent-orbital-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "independent-orbital-number"
  "説明", ""
  "デフォルト値", "0"
  "形式", "(integer > 0)"


* xc-potential/gxalpha/alpha-value

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "xc-potential/gxalpha/alpha-value"
  "説明", ""
  "デフォルト値", "0.7"
  "形式", "real"


* myu-nyu-zero

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "myu-nyu-zero"
  "説明", ""
  "デフォルト値", "no"
  "形式", "(yes|no)"


* guess/nsp-ppq

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "guess/nsp-ppq"
  "説明", ""
  "デフォルト値", "nil"
  "形式", "nil/stored"


* guess/sp-ppq

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "guess/sp-ppq"
  "説明", ""
  "デフォルト値", "nil"
  "形式", "nil/stored"


* guess/trans-angle-threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "guess/trans-angle-threshold"
  "説明", ""
  "デフォルト値", "1.0 1.5 20 30"
  "形式", "nil/stored"


* guess/make-myu-nyu

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "guess/make-myu-nyu"
  "説明", ""
  "デフォルト値", "meth0"
  "形式", "meth0/meth1/meth2/meth3/meth4"


* guess/vct-normalize

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "guess/vct-normalize"
  "説明", ""
  "デフォルト値", "ON OFF OFF"
  "形式", ""


* guess/part-normalize

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "guess/part-normalize"
  "説明", ""
  "デフォルト値", "nil"
  "形式", "nil/stored"


* guess/user-vector

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "guess/user-vector"
  "説明", ""
  "デフォルト値", "nil"
  "形式", "(nil|stored)"


* num_of_atoms

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "num_of_atoms"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* num_of_dummy_atoms

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "num_of_dummy_atoms"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* num_of_AOs

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "num_of_AOs"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* num_of_MOs

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "num_of_MOs"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* num_of_auxCDs

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "num_of_auxCDs"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* num_of_auxXCs

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "num_of_auxXCs"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* TE

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "TE"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* debug/file_warning

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/file_warning"
  "説明", ""
  "デフォルト値", "yes"
  "形式", "(yes|no)"


* debug/save_J

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/save_J"
  "説明", ""
  "デフォルト値", "no"
  "形式", "(yes|no)"


* debug/save_K

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/save_K"
  "説明", ""
  "デフォルト値", "no"
  "形式", "(yes|no)"


* debug/save_Fxc_pure

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/save_Fxc_pure"
  "説明", ""
  "デフォルト値", "no"
  "形式", "(yes|no)"


* debug/save_forces

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/save_forces"
  "説明", ""
  "デフォルト値", "no"
  "形式", "(yes|no)"


* cutoff_distribution

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "cutoff_distribution"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* length_scale_parameter

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "length_scale_parameter"
  "説明", ""
  "デフォルト値", "1"
  "形式", ""


* control

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "control"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* debug/eri/exact_J

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/eri/exact_J"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* cutoff_density

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "cutoff_density"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* cutoff_epsilon3

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "cutoff_epsilon3"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* debug/eri/exact_K

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/eri/exact_K"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* new_engine

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "new_engine"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* debug/eri/output_K

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/eri/output_K"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* debug/eri/output_J

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "debug/eri/output_J"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* num_of_iterations

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "num_of_iterations"
  "説明", ""
  "デフォルト値", ""
  "形式", ""


* stat

.. csv-table::
  :widths: 20,80
  :stub-columns: 1
  
  "キーワード", "stat"
  "説明", ""
  "デフォルト値", ""
  "形式", ""



