#ifndef DFINITIALGUESS_RHO_H
#define DFINITIALGUESS_RHO_H

#include <vector>
#include <string>

#include "DfObject.h"
#include "TlVector.h"
#include "TlSerializeData.h"

#define VER1_0    1   // Normal Auxiliary type;
#define VER1_1    0   // Add one S_orbital type;

/** 初期密度および初期交換相関ポテンシャルの合成を行うクラス
 */
class DfInitialguess : public DfObject {
public:
    DfInitialguess(TlSerializeData* pPdfParam);
    ~DfInitialguess();

    int dfGusMain();
    void Print();

private:
    void densityFitOnly();

private:
    /** 計算に必要な数値を算出したり、fl_Globalinput から値を入力する（計算の前準備ルーチン）
     */
    int prepare();

    /** 各原子が何番目の展開番号から始まっているかを算出する
     */
    int prepare2();

    /// このクラス使用中に必要なメモリを確保する
    int getMemory();

    /// 計算に必要なデータを読み込む
    int readData();

    /** 原子固有のデータをセットする。 prepare2()でセットしていない他のEachADのデータメンバを登録する。
     */
    int setStruct();

    /** 本クラスでのキーワード、オプションなどを読み込む
     */
    int setKeyword();

    /** 初期値合成を効率的に行うためのテーブルを作成する
     */
    int makeStepTable();

    /** 現在はダミーメンバ関数
     *
     *  将来的にDB に登録されている補助関数と違う補助関数を用いて計算を行う場合に
     *  基底の変換を行って初期値を合成しなければならない。その時のための関数。
     */
    int setGDBdatafile();

    /** 現在はダミーメンバ関数
     *
     *  基底変換の実行を行うメンバ関数として用意。
     */
    int CorrespOrbital();

    /** fl_GuessDB/fl_Gdb_atomに入っている原子のデータを構造体に持ち込む
     */
    int setAtomCoef();

    /** 計算を行うメイン関数
     */
    int calcMain();

    /** 密度行列の入力があった場合に値を読み込み、バイナリの形でファイルに出力する
     */
    void readPpq();

    /** DfScf::Densityfitonly()を呼び出す
     */
    void callDensityfitting();

    /** 密度行列の入力があった場合に交換相関ポテンシャルの合成を行う
     */
    void makeXcpotential();

    /** 初期値を原子単位から合成する
     */
    void EXE_ATOM(int);

    /** 現在はダミーメンバ関数
     *
     *  結合方向を考慮した原子単位の合成を行うために用意。
     */
    void Consider_Connect() {};

    /** 初期値を原子単位以外の分子DB から合成する
     */
    void EXE_MOLECULE(int);

    /** ユーザが定義した部分的ベクトルファイルを展開係数ベクトルにはめ込む操作を行う
     */
    void CoverVector();

    // Optimal Angle Serch ;
    /** 最適な回転角度を検索するためのメンバ関数
     *
     *  修正パウエル法に基づく。
     */
    void powell(double*,
                double*,double*,double*,
                double*,double*,double*,int);

    /** 回転角度θx、θy、θz の値を±180 度の間の角になるように調整する
     */
    void AngleTrans(double*,int);

    /** 一変数検索を行う
     *
     *  成功と失敗のルーチン法に基づく
     */
    double rsf(double*,double,double*,
               double*,double*,double*,
               double*,double*,double*,int);

    /** DB 分子の結合原子の座標を
     *  ｙ軸周りにθy、ｘ軸周りにθx、ｚ軸周りにθz だけ回転させた時に、
     *  計算分子の座標との距離の誤差の和を返り値として戻す。
     */
    double AfinFunction(double*,
                        double*,double*,double*,
                        double*,double*,double*,int);

    // Function for getting molecular density.
    /** 展開係数を変換するための回転マトリクスを決定する
     */
    int makeAfinmatrix(int);

    /** makeAfinmatrix で決定した回転マトリクスを利用して展開係数の変換操作を行う
     */
    int transOrbital(int);

    int calcCoefDumy1();
    int calcCoefDumy2();

    /** 利用するDB 分子と計算分子の原子の並びが同じであるかどうか調べる
     */
    void CheckDBFormat(int);

    /** NSP 計算の場合の展開係数規格化メインルーチン
     */
    int VectorNormalyze();

    /** NSP 計算の場合に展開係数ベクトルrhoの規格化を行う
     */
    void VctRouNorm();

    /** NSP 計算の場合に展開係数ベクトルmyuの規格化を行う
     */
    void VctMyuNorm();

    /** NSP 計算の場合に展開係数ベクトルnyuの規格化を行う
     */
    void VctNyuNorm();

    /** SP 計算の場合の展開係数規格化メインルーチン
     */
    int VectorNormalyze2();

    /** SP 計算の場合に展開係数ベクトルrhoの規格化を行う
     */
    void VctRouNorm2();

    /** SP 計算の場合に展開係数ベクトルmyuの規格化を行う
     */
    void VctMyuNorm2();

    /** SP 計算の場合に展開係数ベクトルnyuの規格化を行う
     */
    void VctNyuNorm2();

    /** 合成した初期値を０と１の番号のファイルにバイナリの形で出力する
     */
    int GusOutput();

    /** テンポラルな出力関数
     */
    int GusError(const std::string& str);

    /** rhoの値を設定する
     */
    void setCoefRoudata();

    /** ベクトルの規格化を行う
     */
    void VctNormfromIGuess2();

    //---------------------------------------------------
    // Other functions
    /** ２つのベクトルが作る角度を計算
     */
    double calc_angle(double,double,double,double,double,double);

    /** 等価原子の入れ替えを行う
     */
    void swap(int,int,double*,double*,double*);
    //---------------------------------------------------
    //   Matrix Libraly Member functions

    /** 行列の積を行う
     */
    void multiply_afinmat(double*,double*,double*);
    void multiply_afinmat(std::vector<double>& ,double*,double*);

    /** 行列のコピーを行う
     */
    void dainyu_afinmat(double*,double*);
    void dainyu_afinmat(std::vector<double>&, double*);
    void dainyu_afinmat(double*, std::vector<double>&);

    /** ベクトルのコピーを行う
     */
    void dainyu_vector(double*,double*,int);

    /** 単位行列を作る
     */
    void initialyze_AfinMatrix(double*);
    void initialyze_AfinMatrix(std::vector<double>&);

    /** ゼロ行列を作る
     */
    void clear_AfinMatrix(double*);

    /** 行列を表示する
     */
    void display_afinmatrix(double*);

private:
    int tempflag;
    int tempflagw1;

private:
    //TlSerializeData* pPdfParam_;

private:
    // ---------- Common vaiable in DfInitialguess ----------
    enum { MaxGusMemory  = 10000000};
    enum { MaxAtomName   = 10 };
    enum { MaxAuxSetName = 80 }; // MaxNumber of Auxiliary SetName.
    enum { MaxLabelChar  = 80 }; // MaxNumber of Label.
    enum { MaxSPDterm    =  4 }; // Vetcotr Number for routerm (and myu,nyu).
    enum { aboutAuxTerm  = 80 }; // Average Number of AuxTerm for one Atom.
    enum { MaxFname      = 256 }; // MaxNumber of filename in DfInitialguess.
    enum { LineBuf       = 80 };
    enum { MaxCharcterLength = 300 };//MaxNumber of keyword_rigth_value(char*);
    enum { MaxEqualPairCount = 100 };  // by class Fl_GuessDb
    enum { MaxEqualPair      = 1000 };  // by class Fl_GuessDb
    enum { MaxCntNum         = 10 };   // by class Fl_GuessDb
    enum { MaxStepNum = 1000 }; // Number of MaxDbabel Count.
    // This Number is for StepNumber.
    enum { MaxPartialNum = 20 }; // Max Number of Partial appointment for all.
    enum { Max1DbAtomNum = 100 }; // MaxAtomNum for 1DB's Molecule.

private:
    //--for scf data---
    enum SCFtype { NSP , SP }; // NSP==0 ,SP==1 ;
    SCFtype scftype;

    struct Auxtype {
        std::string AtomKind;   /// 原子名
        std::string AuxSetName; /// 補助関数セット名
        std::vector<int> rDbShellNum; /// rhoのSnum、Pnum、Dnum、SPDnumの数字
        std::vector<int> mDbShellNum; /// myuのSnum、Pnum、Dnum、SPDnumの数字を保存
        std::vector<int> nDbShellNum; /// nyuのSnum、Pnum、Dnum、SPDnumの数字を保存
        std::string Label1;     /// ラベル1
        std::string Label2;     /// ラベル2
        std::vector<double> alpha;    /// rhoの補助関数の広がりを表す指数
        std::vector<double> gamma;    /// myuの補助関数の広がりを表す指数
        std::vector<double> gamma2;   /// nyuの補助関数の広がりを表す指数
        std::vector<int>    routerm;  /// rhoの展開項数 s, p, d型関数のそれぞれの展開項数
        std::vector<int>    myuterm;  /// myuの展開項数 s, p, d型関数のそれぞれの展開項数
        std::vector<int>    nyuterm;  /// nyuの展開項数 s, p, d型関数のそれぞれの展開項数
        //   ↓
        // example) for routerm.
        // [0]=total=[1]*1 + [2]*3 + [3]*5 ;
        //            s       p       d
        std::vector<double> rou; /// rhoの展開係数
        std::vector<double> myu; /// myuの展開係数
        std::vector<double> nyu; /// nyuの展開係数
    };
    //typedef struct Auxtype* GusSet; // Declaration for Auxtype_Object.

private:
    /// 初期値合成のためにどのDBファイルを使用するかを指定する型
    enum Unit_Method {
        fl_Atom,
        fl_Namino,
        fl_Pepamino,
        fl_Camino,
        fl_Molec,
        fl_User
    };

private:
    std::string GDBfile0;  // for Atom.
    std::string GDBfile1;  // for Namino.
    std::string GDBfile2;  // for Pepamino.
    std::string GDBfile3;  // for Camino.
    std::string GDBfile4;  // for public Molecular.
    std::string GDBfile5;  // for User Molecular.

    std::string GDBfile1_2;// for Namino.
    std::string GDBfile2_2;// for Pepamino.
    std::string GDBfile3_2;// for Camino.
    std::string GDBfile4_2;// for public Molecular.
    std::string GDBfile5_2;// for User Molecular.

    /// 初期値合成のためのデータベースのファイル名を指定する
    void setGuessDBfile();

private:
    // --------- Struct for each all_Atom data. ----------
    struct eachAtomData {
        int   Rou_StartNum;    // StartNumber of each Atom's Auxiliary(rou).
        int   Myu_StartNum;    // StartNumber of each Atom's Auxiliary(myu).
        int   Nyu_StartNum;    // StartNumber of each Atom's Auxiliary(nyu).
        // up3 value define in prepare2().
        int   AuxTypeNumber;   // Number of Auxtype(struct).
        // The Auxtype_Data's Number corresopnd to
        // Atom,Label2,Auxiliary's Sets.
        // This value defile in setStruct().
        Unit_Method  Unitfile ;    // making Unit for each Atom.
    };

    // Declaration for eachAtomData_Object.
//   typedef  struct eachAtomData*   EachAtomData;
//   EachAtomData   EachAD;
    std::vector<eachAtomData> EachAD;

    // -----------------------------------------------------------

    // ================ Keyword and option ==================
public:
    enum  KeyOnOff  { ON , OFF };
//   typedef enum KeyOnOff* Key_ONorOFF;

public:
    enum  HowMyuNyu { Meth0 , Meth1 , Meth2 , Meth3 , Meth4 };
    // Meth0 : Make from AtomDB.
    // Meth1 : Apporoximation Myu = (Rou)^(-3) in Expresion;
    // Meth2 : Apporoximation(Meth0) for All Space.(fitting);
    // Meth3 : make from converge file for grand state calc.
    // Meth4 : Content is empty;
    // Meth2,Meth4 is dumy function now.
public:
    enum  WhichRMN  { Rou , Myu , Nyu }; // Rou==0 ; Myu==1 ; Nyu==2 ;

private:
    std::vector<KeyOnOff>  rouNormOnOff; // On/Off : Partial Normalize for each atom.
    std::vector<KeyOnOff>  myuNormOnOff; // On/Off : Partial Normalize for each atom.
    std::vector<KeyOnOff>  nyuNormOnOff; // On/Off : Partial Normalize for each atom.

private:
    struct  KeyContent {
        int      Num;       // Appointment Number of Partial data.
        WhichRMN  RMN;       // Rou or Myu or Nyu;
        int      From1;     // Appointment of From_AtomNumber;
        int      To1;       // Appointment of To___AtomNumber;
        int      From2;     // Appointment of From_AtomNumber;
        int      To2;       // Appointment of To___AtomNumber;
//     char*     Filename;  // Filename for UserDefine_File;
        std::string     Filename;  // Filename for UserDefine_File;
        double    ElecNum;   // Electron Number for UserDefineNum;
    };

    // Declaration for struct KeyContent_Object.
//   typedef  struct  KeyContent*   Keyword;

    int       UserVctRouNum;      // Number of Partial UserDefine Vector;
    int       UserVctMyuNum;      // Number of Partial UserDefine Vector;
    int       UserVctNyuNum;      // Number of Partial UserDefine Vector;
    int       PartialRouNormNum;  // Number of Partial_Nomalize_Vector;
    int       PartialMyuNormNum;  // Number of Partial_Nomalize_Vector;
    int       PartialNyuNormNum;  // Number of Partial_Nomalize_Vector;

    std::string PpqFname     ;      // Appointment : Matrix_Input_FileName;
    std::string AlphaPpqFname;      // Appointment : Matrix_Input_FileName;
    std::string BetaPpqFname ;      // Appointment : Matrix_Input_FileName;

    HowMyuNyu  MethodMyuNyu;       // Method Appointment : making Myu and Nyu;
    std::vector<KeyContent> UserVctRou;         // Appointment : Partial UserDefine Vector;
    std::vector<KeyContent> UserVctMyu;         // Appointment : Partial UserDefine Vector;
    std::vector<KeyContent> UserVctNyu;         // Appointment : Partial UserDefine Vector;
    KeyOnOff   VctRouNormalize;    // Normalize for All Vector(Rou);
    KeyOnOff   VctMyuNormalize;    // Normalize for All Vector(Myu);
    KeyOnOff   VctNyuNormalize;    // Normalize for All Vector(Nyu);
    std::vector<KeyContent> PartialRouNormalize; // Normalize for Partial Vector(Rou);
    std::vector<KeyContent> PartialMyuNormalize; // Normalize for Partial Vector(Rou);
    std::vector<KeyContent> PartialNyuNormalize; // Normalize for Partial Vector(Rou);

    // ============== struct for make StepTable =================
private:
    struct Step {
        Unit_Method Unitfile;
//     char*       DbLabel;
        std::string       DbLabel;
        int         AtomNum;
        int        FirstAtomNum;
        KeyOnOff    UserVctRouExist;
        KeyOnOff    UserVctMyuExist;
        KeyOnOff    UserVctNyuExist;
    };
    //typedef  struct  Step*    StepTable;

    //StepTable  GuessStep;
    std::vector<Step> GuessStep;

    int STEPnumber;

    DfInitialguess::KeyOnOff   serchUserVctExist(WhichRMN,int,int);

    //====== Option : Threshold Angle Gosa for Trans Density=============
private:
    double  FirstAveGosa;
    double  FirstEachGosa;
    double  LastAveGosa;
    double  LastEachGosa;

    // -------------- for data member ---------------------
private:
    int     AtomNum;     // Total atom number.
    int     AtomKindNum; // Number of Atomtype by AtomKind and Lbel2.
    int     AtomKindNumInclDummy; // Number of Atomtype including Dummy Atom.
    int     DumyAtomNum; // Number of Dummy Atom.
    int    ElectronNum; // Total electron Num.
    int    AlphaSpinNum;// Alpha Spin Num.
    int    BetaSpinNum; // Beta  Spin Num.
    double  nAlpha;      // for Xcpotfitting. ex) alpha=0.7

    TlVector CoefRou;
    TlVector CoefMyu;
    TlVector CoefNyu;

    int    MaxTermRou;  // Sum of TermRou.
    int    MaxTermMyu;  // Sum of TermMyu.
    int    MaxTermNyu;  // Sum of TermNyu.

    //GusSet  Atomtype;    // declaration for structure Auxtype* .
    std::vector<Auxtype>  Atomtype;    // declaration for structure Auxtype* .

    TlVector w1Rou, w2Rou;  // Work vector for Rou.
    TlVector w1Myu, w2Myu;  // Work vector for Myu.
    TlVector w1Nyu, w2Nyu;  // Work vector for Nyu.

    static int w1RouCount;
    static int w1MyuCount;
    static int w1NyuCount;

    // (only for SP(UHF) )
    TlVector CoefRouAlpha; // Vector of Initial density for RouAlpha.
    TlVector CoefRouBeta ; // Vector of Initial density for RouBeta.
    TlVector CoefMyuAlpha; // Vector of Initial density for MyuAlpha.
    TlVector CoefMyuBeta;  // Vector of Initial density for MyuBeta.
    TlVector CoefNyuAlpha; // Vector of Initial density for NyuAlpha.
    TlVector CoefNyuBeta;  // Vector of Initial density for NyuBeta.

public:
    struct AtomData {
        std::string auxname;
        std::string atomname;
        double x;
        double y;
        double z;

        std::vector<double> afinmat;
        std::vector<double> CoefRou;
        std::vector<double> CoefMyu;
        std::vector<double> CoefNyu;

        int Routerm;
        int Myuterm;
        int Nyuterm;

        std::vector<int> rouDbShellNum;
        std::vector<int> myuDbShellNum;
        std::vector<int> nyuDbShellNum;
        int ConectNum;
        std::vector<int> Conection;
    };

    struct AminoData {
        int Res1AtomNum;
        std::vector<AtomData> Atomdata;
        std::vector<int> Order;
        int equalPairNum;
        std::vector<int> equalPaircount;
        std::vector<int> equalPair;
        int peptidoTotAtmNum;
        std::vector<int> peptidoAtomnumber;
        int TotalRouterm;
        int TotalMyuterm;
        int TotalNyuterm;
    };

    AminoData  InpAmino;
    AminoData  DbAmino;
};

#endif // DFINITIALGUESS_RHO_H
