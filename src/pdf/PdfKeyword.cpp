#include <iostream>
#include <cassert>
#include "PdfKeyword.h"
#include "TlUtils.h"

PdfKeyword::PdfKeyword()
{
    this->initialize();
}


PdfKeyword::~PdfKeyword()
{
}


TlSerializeData PdfKeyword::getDefault() const
{
    TlSerializeData data;

    KeywordListType::const_iterator pEnd = this->kwdList_.end();
    for (KeywordListType::const_iterator p = this->kwdList_.begin(); p != pEnd; ++p) {
        const std::string keyword = p->keyword;
        const std::string value = p->defaultValue;
        data["model"][keyword] = value;
    }

    return data;
}


void PdfKeyword::initialize()
{
    this->kwdList_.clear();
    KeywordInfo item;

    // MAIN ==============================================================
    item.keyword          = "pdf_param_path";
    item.explanation      = "the file path for ProteinDF parameters.";
    item.explanationJ     = "";
    item.defaultValue     = "pdfparam.mpac";
    item.syntax           = "*";
    item.type             = KWD_DEBUG;
    this->kwdList_.push_back(item);

    item.keyword          = "step_control";
    item.explanation      = "Job scheme of the calculation.";
    item.explanationJ     = "計算手順を指示する。createは初期化, integralは初期積分, guessは初期値作成, scfはSCF繰り返し計算の処理を表す。";
    item.defaultValue     = "create";
    item.syntax           = "(create|integral|guess|scf)+";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "comment";
    item.explanation      = "Comment of the calcculation. ";
    item.explanationJ     = "コメントを指定する。";
    item.defaultValue     = "";
    item.syntax           = "nil/(none)";
    item.type             = KWD_DEFAULT;
    
    this->kwdList_.push_back(item);

    item.keyword          = "restart";
    item.explanation      = "The calculation is restarted with the previous result.";
    item.explanationJ     = "再計算(リスタート)を行う(yes)。";
    item.defaultValue     = "no";
    item.syntax           = "no/yes";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "linear_algebra_package";
    item.explanation      = "linear algebra package";
    item.explanationJ     = "行列演算ライブラリを指定する。";
    item.defaultValue     = "lapack";
    item.syntax           = "(lapack|scalapack)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "scalapack_block_size";
    item.explanation      = "set block size for ScaLAPACK";
    item.explanationJ     = "ScaLAPACKのブロックサイズを指定する。";
    item.defaultValue     = "64";
    item.syntax           = "integer";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "save_distributed_matrix_to_local_disk";
    item.explanation      = "";
    item.explanationJ     = "分散行列をディスクに分散保存する";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "local_disk_path";
    item.explanation      = "";
    item.explanationJ     = "分散行列を分散保持するパスを指定する";
    item.defaultValue     = "/tmp";
    item.syntax           = "";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);
    
    item.keyword          = "memory_size";
    item.explanation      = "";
    item.explanationJ     = "プロセスあたりの使用可能なメモリサイズを指定する。";
    item.defaultValue     = "1GB";
    item.syntax           = "";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "use_mapfile";
    item.explanation      = "";
    item.explanationJ     = "メモリが不足時はディスクを補助記憶領域として利用する";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "mapfile_size";
    item.explanation      = "";
    item.explanationJ     = "補助記憶領域のサイズを指定する";
    item.defaultValue     = "auto";
    item.syntax           = "";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "mapfile_basename";
    item.explanation      = "";
    item.explanationJ     = "補助記憶領域のパスを指定する";
    item.defaultValue     = "/tmp";
    item.syntax           = "";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "work_on_disk";
    item.explanation      = "";
    item.explanationJ     = "計算時にメモリに余裕があっても必ずディスク領域を用いる(yes)。";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "parallel_processing_type";
    item.explanation      = "computing method in parallel calculation";
    item.explanationJ     = "並列計算の計算方法を指定する。";
    //item.defaultValue     = "divide_and_conquer";
    item.defaultValue     = "MS";
    item.syntax           = "(divide_and_conquer|DC|master_slave|MS)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "cleanup";
    item.explanation      = "clean up intermediate files";
    item.explanationJ     = "";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);
    
    // INPUT =============================================================
    item.keyword          = "show_keyword";
    item.explanation      = "print keyword list";
    item.explanationJ     = "ログファイルにキーワードのリストを出力する(yes)。";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "show_input";
    item.explanation      = "print input perameters";
    item.explanationJ     = "ログファイルに入力データを出力する(yes)。";
    item.defaultValue     = "yes";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "show_coordinates";
    item.explanation      = "print coordinates";
    item.explanationJ     = "ログファイルに座標データを出力する(yes)。";
    item.defaultValue     = "yes";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "show_orbital_basis";
    item.explanation      = "orbital basis output type";
    item.explanationJ     = "基底関数の情報を指定のフォーマットで出力する。";
    item.defaultValue     = "gamess";
    item.syntax           = "(gamess|amoss|none)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    // GUESS =============================================================
    item.keyword          = "scf-start-guess";
    item.explanation      = "star scf from guess lcao and occupation files, or rho and occupation (fl_Userinput) ";
    item.explanationJ     = "SCF 計算における初期値を指定する。";
    item.defaultValue     = "rho";
    item.syntax           = "(rho|file_rho|lcao|density|core|huckel|harris)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    // INTEGRAL ==========================================================
    item.keyword          = "cut-value";
    item.explanation      = "cut off threshold for integrals";
    item.explanationJ     = "積分計算のカットオフ値を指定する。";
    item.defaultValue     = "1.0E-10";
    item.syntax           = "(double > 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "block-size";
    item.explanation      = "block size for integrals";
    item.explanationJ     = "積分計算における計算粒度を指定する。";
    item.defaultValue     = "1024000";
    item.syntax           = "(integer > 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    // SCF ===============================================================
    item.keyword          = "charge-extrapolate-number";
    item.explanation      = "";
    item.explanationJ     = "点電荷を段階的に挿入するときに使用する。";
    item.defaultValue     = "0";
    item.syntax           = "(integer >= 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "orbital-overlap-correspondence";
    item.explanation      = "";
    item.explanationJ     = "軌道対応のチェックを行う。";
    item.defaultValue     = "off";
    item.syntax           = "on/off";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "orbital-overlap-correspondence-first";
    item.explanation      = "skip 1st orbital correspondence to 0th orbital in DfDmatrix";
    item.explanationJ     = "";
    item.defaultValue     = "off";
    item.syntax           = "on/off";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    // 射影演算子法を使うために必要なキーワード orbital-overlap-correspondence-method
    // orbital-overlap-correspondence-method = mo-overlap    : 軌道の重なりの対応
    // orbital-overlap-correspondence-method = mo-projection : 射影演算子法
    // ただし、最初の数iterationを射影演算子法、
    // DfDmatrix.cxxのMO_OVERLAP_ITERで指定した数以降のiterationでは軌道の重なりの対応
    // を使用する場合は、mo-overlapを指定。

    item.keyword          = "orbital-overlap-correspondence-method";
    item.explanation      = "orbital correspondence method";
    item.explanationJ     = "mo-overlap法もしくはmo-projection法を利用するかどうかを指定する。";
    item.defaultValue     = "mo-overlap";
    item.syntax           = "mo-overlap/mo-projection";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "orbital-overlap-correspondence-number";
    item.explanation      = "maxmum iteration number by orbital projection method";
    item.explanationJ     = "";
    item.defaultValue     = "3";
    item.syntax           = "integer";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "summary";
    item.explanation      = "print summary data";
    item.explanationJ     = "計算結果の要約を表示する。";
    item.defaultValue     = "none";
    item.syntax           = "(none|convergence|every-scf)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "analyze_population";
    item.explanation      = "analize population";
    item.explanationJ     = "電荷解析を行う。";
    item.defaultValue     = "none";
    item.syntax           = "(none|convergence|every-scf)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "max-iteration";
    item.explanation      = "";
    item.explanationJ     = "最大繰り返し回数を指定する。";
    item.defaultValue     = "100";
    item.syntax           = "(integer)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "method";
    item.explanation      = "";
    item.explanationJ     = "計算方法を指定する。";
    item.defaultValue     = "nsp";
    item.syntax           = "sp/nsp/roks";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "method/nsp/occlevel";
    item.explanation      = "";
    item.explanationJ     = "占有軌道のレベルを指定する。";
    item.defaultValue     = "";
    item.syntax           = "(array of integer >= 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "RKS/electrons";
    item.explanation      = "";
    item.explanationJ     = "電子数を指定する。";
    item.defaultValue     = "";
    item.syntax           = "(integer >= 2)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "method/sp/alpha-spin-occlevel";
    item.explanation      = "";
    item.explanationJ     = "占有軌道のレベルを指定する。";
    item.defaultValue     = "";
    item.syntax           = "(array of integer >= 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "UKS/alphaElectrons";
    item.explanation      = "";
    item.explanationJ     = "電子数を指定する。";
    item.defaultValue     = "";
    item.syntax           = "(integer >= 1)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "method/sp/beta-spin-occlevel";
    item.explanation      = "";
    item.explanationJ     = "占有軌道のレベルを指定する。";
    item.defaultValue     = "";
    item.syntax           = "(array of integer >= 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "UKS/betaElectrons";
    item.explanation      = "";
    item.explanationJ     = "電子数を指定する。";
    item.defaultValue     = "";
    item.syntax           = "(integer >= 1)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "save_diff_density_matrix";
    item.explanation      = "";
    item.explanationJ     = "差電子密度行列をファイルに保存する。";
    item.defaultValue     = "yes";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);
    
    item.keyword          = "use_matrix_cache";
    item.explanation      = "";
    item.explanationJ     = "行列をキャッシュする。";
    item.defaultValue     = "yes";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "show_cache_report";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "new_engine";
    item.explanation      = "using new engine(2011)";
    item.explanationJ     = "";
    item.defaultValue     = "yes";
    item.syntax           = "(yes|no)";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    // for ROKS
    item.keyword          = "method/roks/electron-number";
    item.explanation      = "electron numbers on ROKS calculation";
    item.explanationJ     = "占有軌道のレベルを指定する。";
    item.defaultValue     = "";
    item.syntax           = "(integer >= 1)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "method/roks/electron-number-alpha";
    item.explanation      = "alpha electron numbers on ROKS calculation";
    item.explanationJ     = "電子数を指定する。";
    item.defaultValue     = "";
    item.syntax           = "(integer >= 1)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "method/roks/electron-number-beta";
    item.explanation      = "beta electron numbers on ROKS calculation";
    item.explanationJ     = "電子数を指定する。";
    item.defaultValue     = "";
    item.syntax           = "(integer >= 1)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "method/roks/closed-shell";
    item.explanation      = "closed shell orbital on ROKS calculation";
    item.explanationJ     = "閉殻軌道のレベルを指定する。";
    item.defaultValue     = "";
    item.syntax           = "(array of integer >= 1)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "method/roks/open-shell";
    item.explanation      = "open shell orbital on ROKS calculation";
    item.explanationJ     = "開殻軌道のレベルを指定する。";
    item.defaultValue     = "";
    item.syntax           = "(array of integer >= 1)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "disk-utilization";
    item.explanation      = "";
    item.explanationJ     = "DISK法を用いる。";
    item.defaultValue     = "no";
    item.syntax           = "yes/no";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "orbital-independence-threshold";
    item.explanation      = "";
    item.explanationJ     = "基底関数の線形独立性を判定するための閾値を指定する。";
    item.defaultValue     = "0.007";
    item.syntax           = "(real >= 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

//     item.keyword          = "convergence";
//     item.explanation      = "";
//     item.explanationJ     = "";
//     item.defaultValue     = "this";
//     item.syntax           = "this";
//     item.type             = KWD_DEFAULT;
//     this->kwdList_.push_back(item);

    item.keyword          = "convergence/type";
    item.explanation      = "";
    item.explanationJ     = "収束判定に用いる物理量の種類（全エネルギー以外）を指定する。";
    item.defaultValue     = "density";
    item.syntax           = "(fock|density|dcoef)";
    item.type             = KWD_DEFAULT;;
    this->kwdList_.push_back(item);

    item.keyword          = "convergence/threshold";
    item.explanation      = "";
    item.explanationJ     = "指定された物理量における収束判定のための閾値を指定する。";
    item.defaultValue     = "1e-3";
    item.syntax           = "(real > 0)";
    item.type             = KWD_DEFAULT;;
    this->kwdList_.push_back(item);

    item.keyword          = "convergence/threshold-energy";
    item.explanation      = "";
    item.explanationJ     = "収束判定における全エネルギー用のしきい値を指定する。";
    item.defaultValue     = "1e-4";
    item.syntax           = "(real > 0)";
    item.type             = KWD_DEFAULT;;
    this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration";
    item.explanation      = "";
    item.explanationJ     = "SCF計算における収束加速法を指定する。";
    item.defaultValue     = "anderson";
    item.syntax           = "(damping|anderson|diis)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration/damping/damping-factor";
    item.explanation      = "";
    item.explanationJ     = "damping 法における damping 割合を指定する。";
    item.defaultValue     = "0.85";
    item.syntax           = "(real)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration/daming/start-number";
    item.explanation      = "start iteration number of damping";
    item.explanationJ     = "damping法の開始SCF回数を指定する。";
    item.defaultValue     = "0";
    item.syntax           = "(int)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration/damping/damping-type";
    item.explanation      = "";
    item.explanationJ     = "damping法を適応する物理量を指定する。";
    item.defaultValue     = "density";
    item.syntax           = "(fock|density_matrix|density|dcoef)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration/anderson/damping-factor";
    item.explanation      = "";
    item.explanationJ     = "anderson法でのdamping係数を指定する。";
    item.defaultValue     = "0.50";
    item.syntax           = "(real)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration/anderson/start-number";
    item.explanation      = "start iteration number of Anderson's converge";
    item.explanationJ     = "anderson法の開始SCF回数を指定する。";
    item.defaultValue     = "3";
    item.syntax           = "(int)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "level-shift";
    item.explanation      = "add level shift to F' matrix";
    item.explanationJ     = "level-shift 法を使用する。";
    item.defaultValue     = "off";
    item.syntax           = "off/no";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "level-shift/start-iteration";
    item.explanation      = "level shift is carried out if start-iteration >= SCF iteration";
    item.explanationJ     = "レベルシフトの開始SCF回数を指定する。";
    item.defaultValue     = "1";
    item.syntax           = "(int)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "level-shift/ls-closed-mo";
    item.explanation      = "level shift value for closed mo";
    item.explanationJ     = "閉殻軌道に加えるレベルシフト値を指定する。a.u. 単位で指定する。";
    item.defaultValue     = "0.00";
    item.syntax           = "(real)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "level-shift/ls-open-mo";
    item.explanation      = "level shift value for open mo";
    item.explanationJ     = "開殻軌道に加えるレベルシフト値を指定する。";
    item.defaultValue     = "0.00";
    item.syntax           = "(real)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "level-shift/ls-virtual-mo";
    item.explanation      = "level shift value for virtural mo";
    item.explanationJ     = "空軌道に加えるレベルシフト値を指定する。";
    item.defaultValue     = "0.00";
    item.syntax           = "(real)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

//     item.keyword          = "level-shift/delta-group-closed";
//     item.explanation      = "level shift value between closed group";
//     item.explanationJ     = "";
//     item.defaultValue     = "0.00";
//     item.syntax           = "(real)";
//     item.type             = KWD_DEFAULT;
//     this->kwdList_.push_back(item);

//     item.keyword          = "level-shift/delta-group-open";
//     item.explanation      = "level shift value between open group";
//     item.explanationJ     = "";
//     item.defaultValue     = "0.00";
//     item.syntax           = "(real)";
//     item.type             = KWD_DEFAULT;
//     this->kwdList_.push_back(item);

//     item.keyword          = "level-shift/delta-group-virtual";
//     item.explanation      = "level shift value between virtual group";
//     item.explanationJ     = "";
//     item.defaultValue     = "0.00";
//     item.syntax           = "(real)";
//     item.type             = KWD_DEFAULT;
//     this->kwdList_.push_back(item);

//     item.keyword          = "scf-acceleration/diis/type";
//     item.explanation      = "";
//     item.explanationJ     = "";
//     item.defaultValue     = "pulay";
//     item.syntax           = "pulay/sellers";
//     item.type             = KWD_DEFAULT;
//     this->kwdList_.push_back(item);

//     item.keyword          = "scf-acceleration/diis/property";
//     item.explanation      = "";
//     item.explanationJ     = "";
//     item.defaultValue     = "fock";
//     item.syntax           = "fock/density";
//     item.type             = KWD_DEFAULT;
//     this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration/diis/number-of-diis";
    item.explanation      = "";
    item.explanationJ     = "diis 外挿に使用する過去の行列またはベクトル数を指定する。";
    item.defaultValue     = "3";
    item.syntax           = "(integer)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration/diis/start-number";
    item.explanation      = "";
    item.explanationJ     = "diis 外挿の準備を始める SCF 回数を指定する。";
    item.defaultValue     = "3";
    item.syntax           = "(integer >= 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "scf-acceleration/diis/start-extrapolation";
    item.explanation      = "";
    item.explanationJ     = "diis 外挿を始める SCF 回数を指定する。";
    item.defaultValue     = "6";
    item.syntax           = "(integer >= 0)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "xc-potential";
    item.explanation      = "XC-potential type";
    item.explanationJ     = "交換相関ポテンシャルを指定する。";
    item.defaultValue     = "svwn~";
    item.syntax           = "(svwn~|svwn|blyp|b3lyp)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "xc-potential/grid-type";
    item.explanation      = "Grid type selection for grid-used methods";
    item.explanationJ     = "グリッドを指定する。";
    item.defaultValue     = "sg-1";
    item.syntax           = "fine/medium-fine/medium/coarse/sg-1";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "xc-update";
    item.explanation      = "whether incDFT method is executed or not";
    item.explanationJ     = "incDFT法を用いる。";
    item.defaultValue     = "yes";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "xc-density-threshold";
    item.explanation      = "cutoff threshold of density value on each grid";
    item.explanationJ     = "各グリッド上の電子密度の閾値を指定する。";
    item.defaultValue     = "1.0E-16";
    item.syntax           = "real";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);
    
    item.keyword          = "TEI-integral-driven";
    item.explanation      = "method of two-electron integral";
    item.explanationJ     = "2電子積分の計算方法を指定する。";
    item.defaultValue     = "yes";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "RI_J";
    item.explanation      = "using RI method on building J matrix";
    item.explanationJ     = "クーロン項をRIにて計算する。";
    item.defaultValue     = "yes";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "RI_K";
    item.explanation      = "using RI method on building K matrix";
    item.explanationJ     = "Fockの交換項をRIにて計算する。";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);
    
    item.keyword          = "scf-memory-saving";
    item.explanation      = "yes : use less memories and increase the times of calculation of three-index integrals";
    item.explanationJ     = "";
    item.defaultValue     = "no";
    item.syntax           = "yes/no";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "geometry/cartesian/input";
    item.explanation      = "Input of nuclear geometry by cartesian.";
    item.explanationJ     = "核座標をカーテシアン座標で入力する。";
    item.defaultValue     = "nil";
    item.syntax           = "nil/stored";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "geometry/cartesian/unit";
    item.explanation      = "Unit of inputted nuclear geometry. Default is a.u.";
    item.explanationJ     = "カーテシアン座標の単位を指定する。";
    item.defaultValue     = "a.u.";
    item.syntax           = "au/a.u./angstrom";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "basis-set/orbital";
    item.explanation      = "Input of basiss set for orbitals.";
    item.explanationJ     = "基底関数を指定する。";
    item.defaultValue     = "nil";
    item.syntax           = "nil/stored";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "basis-set/density-auxiliary";
    item.explanation      = "Input of basis set for electron density.";
    item.explanationJ     = "クーロン項用の補助基底関数を指定する。";
    item.defaultValue     = "nil";
    item.syntax           = "nil/stored";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "basis-set/exchange-auxiliary";
    item.explanation      = "Input of basi set for exchange potential.";
    item.explanationJ     = "交換相関項用の補助基底関数を指定する。";
    item.defaultValue     = "nil";
    item.syntax           = "nil/stored";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

    item.keyword          = "xc_density_threshold";
    item.explanation      = "";
    item.explanationJ     = "数値積分の各グリッド上の密度の閾値を指定する";
    item.defaultValue     = "1.0E-16";
    item.syntax           = "real";
    item.type             = KWD_DEFAULT;
    this->kwdList_.push_back(item);

//     item.keyword          = "incDFT";
//     item.explanation      = "";
//     item.explanationJ     = "incDFT法を用いる。";
//     item.defaultValue     = "yes";
//     item.syntax           = "(yes|no)";
//     item.type             = KWD_DEFAULT;
//     this->kwdList_.push_back(item);

    // internal ================================================================
    item.keyword          = "geometry";
    item.explanation      = "set geometry.";
    item.explanationJ     = "";
    item.defaultValue     = "cartesian";
    item.syntax           = "cartesian/file";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "coordinates";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "";
    item.syntax           = "";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "basis_set";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "";
    item.syntax           = "";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);
    
    item.keyword          = "basis_set_auxD";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "";
    item.syntax           = "";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "basis_set_auxXC";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "";
    item.syntax           = "";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "independent-orbital-number";
    item.explanation      = "number of linearly independent orbitals";
    item.explanationJ     = "";
    item.defaultValue     = "0";
    item.syntax           = "(integer > 0)";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "xc-potential/gxalpha/alpha-value";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "0.7";
    item.syntax           = "real";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "myu-nyu-zero";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "guess/nsp-ppq";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "nil";
    item.syntax           = "nil/stored";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "guess/sp-ppq";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "nil";
    item.syntax           = "nil/stored";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "guess/trans-angle-threshold";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "1.0 1.5 20 30";
    item.syntax           = "nil/stored";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "guess/make-myu-nyu";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "meth0";
    item.syntax           = "meth0/meth1/meth2/meth3/meth4";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "guess/vct-normalize";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "ON OFF OFF";
    item.syntax           = "";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "guess/part-normalize";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "nil";
    item.syntax           = "nil/stored";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "guess/user-vector";
    item.explanation      = "";
    item.explanationJ     = "";
    item.defaultValue     = "nil";
    item.syntax           = "(nil|stored)";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "debug/file_warning";
    item.explanation      = "";
    item.explanationJ     = "ファイルが存在しない場合に警告をログに表示する";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    item.keyword          = "debug/save_J";
    item.explanation      = "";
    item.explanationJ     = "クーロン成分(J)行列を保存する";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);
    
    item.keyword          = "debug/save_forces";
    item.explanation      = "save the force of each term.";
    item.explanationJ     = "";
    item.defaultValue     = "no";
    item.syntax           = "(yes|no)";
    item.type             = KWD_HIDDEN;
    this->kwdList_.push_back(item);

    // alias
    this->setupAliasList();
}


void PdfKeyword::setupAliasList()
{
    // old_keyword -> new_keyword
    this->kwdAlias_["step-control"] = "step_control";
    this->kwdAlias_["linear-algebra-package"] = "linear_algebra_package";
    this->kwdAlias_["method/nsp/electron-number"]  = "RKS/electrons";
    this->kwdAlias_["method/sp/alpha-elec-number"] = "UKS/alphaElectrons";
    this->kwdAlias_["method/sp/beta-elec-number"]  = "UKS/betaElectrons";

    this->kwdAlias_["parallel-processing-type"]  = "parallel_processing_type";
}


void PdfKeyword::convertAlias(TlSerializeData* pData)
{
    assert(pData != NULL);

    TlSerializeData::MapIterator p = pData->beginMap();
    TlSerializeData::MapIterator pEnd = pData->endMap();
    while (p != pEnd) {
        const std::string keyword = p->first.getStr();
        const TlSerializeData value = p->second;

        AliasContainerType::const_iterator it = this->kwdAlias_.find(keyword);
        if (it != this->kwdAlias_.end()) {
            const std::string newKeyword = it->second;

            if (pData->hasKey(newKeyword) != true) {
                (*pData)[newKeyword] = value;
                pData->erase(keyword);
            } else {
                std::cerr << TlUtils::format("[WARN] alias keyword \"%s\" could not overwrite to \"%s\".",
                                             keyword.c_str(), newKeyword.c_str())
                          << std::endl;
            }
            p = pData->beginMap();
        } else {
            ++p;
        }
    }
}


std::string PdfKeyword::getCSV(bool showHiddenItem) const
{
    std::string output = "keyword, explanation, default, syntax\n";

    const int numOfItems = this->kwdList_.size();
    for (int i = 0; i < numOfItems; ++i) {
        const KeywordInfo& item = this->kwdList_[i];
        if (((item.type & KWD_HIDDEN) != 0) &&
            (showHiddenItem == false)) {
            continue;
        }
        output += TlUtils::format("\"%s\", \"%s\", \"%s\", \"%s\"\n",
                                  item.keyword.c_str(),
                                  item.explanation.c_str(),
                                  item.defaultValue.c_str(),
                                  item.syntax.c_str());
    }

    return output;
}


std::string PdfKeyword::getCSV_J(bool showHiddenItem) const
{
    std::string output = "keyword, explanation, default, syntax\n";

    const int numOfItems = this->kwdList_.size();
    for (int i = 0; i < numOfItems; ++i) {
        const KeywordInfo& item = this->kwdList_[i];
        if (((item.type & KWD_HIDDEN) != 0) &&
            (showHiddenItem == false)) {
            continue;
        }
        output += TlUtils::format("\"%s\", \"%s\", \"%s\", \"%s\"\n",
                                  item.keyword.c_str(),
                                  item.explanationJ.c_str(),
                                  item.defaultValue.c_str(),
                                  item.syntax.c_str());
    }

    return output;
}


TlSerializeData PdfKeyword::getSerializeData() const
{
    TlSerializeData data;

    KeywordListType::const_iterator itEnd = this->kwdList_.end();
    for (KeywordListType::const_iterator it = this->kwdList_.begin(); it != itEnd; ++it) {
        TlSerializeData property;
        property["description"] = it->explanation;
        //property["description_jp"] = it->explanationJ;
        property["default"] = it->defaultValue;
        property["syntax"] = it->syntax;
        property["hidden"] = (it->type != KWD_DEFAULT);

        TlSerializeData item;
        item["keyword"] = it->keyword;
        item["property"] = property;
        
        data.pushBack(item);
    }
    
    return data;
}


void PdfKeyword::checkInputParam(const TlSerializeData& param) const
{
    if (param.getType() != TlSerializeData::MAP) {
        std::cerr << TlUtils::format("[WARN] input param type mismatch. file: %s, line: %d",
                                     __FILE__, __LINE__)
                  << std::endl;
    }

    TlSerializeData::MapConstIterator mapItEnd = param.endMap();
    for (TlSerializeData::MapConstIterator mapIt = param.beginMap(); mapIt != mapItEnd; ++mapIt) {
        const std::string keyword = mapIt->first.getStr();
        const bool hasKeyword = this->hasKeyword(keyword);
        if (hasKeyword == false) {
            std::cerr << TlUtils::format("[WARN] unknown keyword: %s", keyword.c_str())
                      << std::endl;
        }
    }
}


bool PdfKeyword::hasKeyword(const std::string& keyword) const
{
    if (this->kwdDb_.empty() == true) {
        this->makeDB();
    }

    KeywordDbType::const_iterator it = this->kwdDb_.find(keyword);
    return (it != this->kwdDb_.end());
}


void PdfKeyword::makeDB() const
{
    this->kwdDb_.clear();

    const std::size_t numOfKwds = this->kwdList_.size();
    for (std::size_t i = 0; i < numOfKwds; ++i) {
        const KeywordInfo& info = this->kwdList_[i];
        const std::string& keyword = info.keyword;

        this->kwdDb_[keyword] = info;
    }

    if (this->kwdDb_.size() != numOfKwds) {
        std::cerr << TlUtils::format("[WARN] size mismatch: file: %s, line: %s.", __FILE__, __LINE__)
                  << std::endl;
    }
}

