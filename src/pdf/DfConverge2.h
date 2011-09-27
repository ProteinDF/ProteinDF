#ifndef DFCONVERGE2_H
#define DFCONVERGE2_H

#include "DfObject.h"
#include "TlMatrix.h"
#include "TlVector.h"

/// DIIS 法の処理を行うクラス
class DfConverge2 : public DfObject {
public:
    DfConverge2(TlSerializeData* pPdfParam, int num_iter);
    virtual ~DfConverge2();

    bool DfConv2Main();

private:
    /** エラーベクトルの計算を行う
     */
    void calculate_ErrorVector(const RUN_TYPE runType, int iteration, TlMatrix& A);

    /** 前回のB 行列を読み込む
     */
    void read_PreviousBmatrix(const std::string& type, int iteration, TlMatrix& X);

    /** B 行列の更新を行う
     */
    void update_Bmatrix(const std::string& type, int iteration, TlMatrix& B, TlMatrix& E);

    /** Fock行列にDIIS法の処理を行う
     */
    void interpolate_KhonSham(const std::string& type, int iteration, TlVector& c);

    /** DIIS 処理を行うための係数を求める
     */
    void solve_DIIS(TlMatrix& B, TlVector& c);

    /** 係数を求めるための連立1 次方程式を解く
     */
    void gauss(TlMatrix& B, TlVector& X, TlVector& L);


private:
//     int number_iteration;

//     int number_mo_basis;
//     int number_ao_basis;

//     std::string scftype;

    // DIIS
    int diis_dimension;           /// diis 外挿に使用する行列またはベクトル数を指定するキーワード
    int diis_property;            /// ユーザ設定のscf-acceleration/diis/property キーワード
    int diis_start_number;        /// diis 外挿の準備を始める SCF 回数を指定するキーワード
    int diis_start_extrapolation; /// diis 外挿を始める SCF 回数を指定するキーワード
    std::string diis_type;        /// ユーザ設定のscf-acceleration/diis/type キーワード
    int diis_actual_dimension;    /// diis 処理に入ってからのiteration 数
};

#endif // DFCONVERGE2_H


