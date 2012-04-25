#ifndef DFCONVCHECK_H
#define DFCONVCHECK_H

#include <string>
#include "DfObject.h"

// used for DfConvcheck class
struct THRESHOLD {
    double den;
    double ene;
};


/// 収束判定を行うクラス
class DfConvcheck : public DfObject {
public:
    DfConvcheck(TlSerializeData* pPdfParam, int num_iter);
    virtual ~DfConvcheck();

    void DfConvcheckMain();
    int judge() {
        return converged_flag;
    }

protected:
    template<class SymmetricMatrixType>
    void main(int iteration);

    virtual void showResults();

    /** 全エネルギーのずれを計算する
     */
    double dev_total_energy(int iteration);

    /** 密度行列展開係数の最大のずれを計算する
     */
    double dev_cd_coefficient(RUN_TYPE runType, int iteration);

    /** 交換相関ポテンシャル係数の最大のずれを計算する
     */
    double dev_xc_coefficient(RUN_TYPE runType, int iteration);

    /** X-alpha法の時の展開係数の最大のずれを計算する
     */
    double dev_xa_coefficient(RUN_TYPE runType, int iteration);

protected:
    /** Fock行列要素の最大のずれを計算する
     */
    template<class SymmetricMatrixType>
    double dev_kohn_sham_matrix(const RUN_TYPE, int iteration);
    
    /** 密度行列要素の最大のずれを計算する
     */
    template<class SymmetricMatrixType>
    double dev_density_matrix(RUN_TYPE runType, int iteration);
    
    template<class SymmetricMatrixType>
    double dev_standard_dev_cd(RUN_TYPE runType, int iteration);

    
protected:
    std::string convergence_type;

    double dev_sd; /// 密度行列要素全体の前回と今回の値のずれの平均
    double dev_dm; /// 密度行列要素の前回と今回の値のずれの最大値
    double dev_te; /// 全エネルギーの前回と今回の値のずれ
    double dev_ks; /// フォック行列要素の前回と今回の値のずれの最大値
    double dev_cd; /// 密度行列展開係数の前回と今回の値のずれの最大値
    double dev_xc; /// 交換相関ポテンシャル展開係数の前回と今回の値のずれの最大値
    double dev_xa; /// X-alpha法の時の展開係数の前回と今回の値のずれの最大値
    double threshold_cri;     /// 値のずれに対して許容するしきい値
    double threshold_cri_ene; /// エネルギーの値のずれに対して許容するしきい値

    double Flimit;   // フォック行列要素の値のずれに対して許容するしきい値
    double Dlimit;   // 密度行列展開係数の値のずれに対して許容するしきい値
    double DMlimit;  // 密度行列要素の値のずれに対して許容するしきい値
    double DElimit;  // 全エネルギーの値のずれに対して許容するしきい値

    int converged_flag; /// 収束判定に関するフラッグ。収束時には0が入る

    double dev_sd_a;
    double dev_sd_b;
    double dev_dm_a;
    double dev_dm_b;
    double dev_dm_c;
    double dev_dm_o;
    double dev_ks_a;
    double dev_ks_b;
    double dev_cd_a;
    double dev_cd_b;
    double dev_xc_a;
    double dev_xc_b;
    double dev_xa_a;
    double dev_xa_b;
};


template<class SymmetricMatrixType>
void DfConvcheck::main(const int iteration)
{
    this->dev_te = dev_total_energy(iteration);

    switch (this->m_nMethodType) {
    case METHOD_RKS:
        {
            if (this->J_engine_ == J_ENGINE_RI_J) {
                this->dev_sd = this->dev_standard_dev_cd<SymmetricMatrixType>(RUN_RKS, iteration);
            }
            this->dev_dm = this->dev_density_matrix<SymmetricMatrixType>(RUN_RKS, iteration);
            this->dev_ks = this->dev_kohn_sham_matrix<SymmetricMatrixType>(RUN_RKS, iteration);
            if (this->m_bIsXCFitting == true) {
                this->dev_cd = dev_cd_coefficient(RUN_RKS, iteration);
                //this->dev_xc = dev_xc_coefficient  (type, iteration );
                //this->dev_xa = dev_xa_coefficient  (type, iteration );
            }
        }
        break;

    case METHOD_UKS:
        {
            if (this->J_engine_ == J_ENGINE_RI_J) {
                this->dev_sd_a = dev_standard_dev_cd<SymmetricMatrixType>(RUN_UKS_ALPHA, iteration);
                this->dev_sd_b = dev_standard_dev_cd<SymmetricMatrixType>(RUN_UKS_BETA,  iteration);
            }

            this->dev_dm_a = dev_density_matrix<SymmetricMatrixType>(RUN_UKS_ALPHA, iteration);
            this->dev_ks_a = dev_kohn_sham_matrix<SymmetricMatrixType>(RUN_UKS_ALPHA, iteration);
            this->dev_cd_a = dev_cd_coefficient(RUN_UKS_ALPHA, iteration);
            this->dev_xc_a = dev_xc_coefficient(RUN_UKS_ALPHA, iteration);
            this->dev_xa_a = dev_xa_coefficient(RUN_UKS_ALPHA, iteration);

            this->dev_dm_b = dev_density_matrix<SymmetricMatrixType>(RUN_UKS_BETA,  iteration);
            this->dev_ks_b = dev_kohn_sham_matrix<SymmetricMatrixType>(RUN_UKS_BETA,  iteration);
            this->dev_cd_b = dev_cd_coefficient(RUN_UKS_BETA,  iteration);
            this->dev_xc_b = dev_xc_coefficient(RUN_UKS_BETA,  iteration);
            this->dev_xa_b = dev_xa_coefficient(RUN_UKS_BETA,  iteration);
            
            this->dev_sd = (dev_sd_a > dev_sd_b) ? dev_sd_a : dev_sd_b;
            this->dev_dm = (dev_dm_a > dev_dm_b) ? dev_dm_a : dev_dm_b;
            this->dev_ks = (dev_ks_a > dev_ks_b) ? dev_ks_a : dev_ks_b;
            this->dev_cd = (dev_cd_a > dev_cd_b) ? dev_cd_a : dev_cd_b;
            this->dev_xc = (dev_xc_a > dev_xc_b) ? dev_xc_a : dev_xc_b;
            this->dev_xa = (dev_xa_a > dev_xa_b) ? dev_xa_a : dev_xa_b;
        }
        break;

    case METHOD_ROKS:
        {
            this->dev_sd_a = dev_standard_dev_cd<SymmetricMatrixType>(RUN_UKS_ALPHA,  iteration);
            this->dev_sd_b = dev_standard_dev_cd<SymmetricMatrixType>(RUN_UKS_BETA,   iteration);
            this->dev_dm_c = dev_density_matrix<SymmetricMatrixType>(RUN_ROKS_CLOSE, iteration);
            this->dev_dm_o = dev_density_matrix<SymmetricMatrixType>(RUN_ROKS_OPEN,  iteration);
            this->dev_ks = this->dev_kohn_sham_matrix<SymmetricMatrixType>(RUN_ROKS,           iteration);
            this->dev_cd_a = dev_cd_coefficient(RUN_UKS_ALPHA,  iteration);
            this->dev_xc_a = dev_xc_coefficient(RUN_UKS_ALPHA,  iteration);
            this->dev_xa_a = dev_xa_coefficient(RUN_UKS_ALPHA,  iteration);
            this->dev_cd_b = dev_cd_coefficient(RUN_UKS_BETA,   iteration);
            this->dev_xc_b = dev_xc_coefficient(RUN_UKS_BETA,   iteration);
            this->dev_xa_b = dev_xa_coefficient(RUN_UKS_BETA,   iteration);
            
            this->dev_sd = (dev_sd_a > dev_sd_b) ? dev_sd_a : dev_sd_b;
            this->dev_dm = (dev_dm_c > dev_dm_o) ? dev_dm_c : dev_dm_o;
            this->dev_cd = (dev_cd_a > dev_cd_b) ? dev_cd_a : dev_cd_b;
            this->dev_xc = (dev_xc_a > dev_xc_b) ? dev_xc_a : dev_xc_b;
            this->dev_xa = (dev_xa_a > dev_xa_b) ? dev_xa_a : dev_xa_b;
        }
        break;

    default:
        abort();
        break;
    }
}


template<class SymmetricMatrixType>
double DfConvcheck::dev_density_matrix(const RUN_TYPE runType, const int iteration)
{
    // density matrix convergence
    SymmetricMatrixType prevP = DfObject::getPpqMatrix<SymmetricMatrixType>(runType, iteration -1);
    SymmetricMatrixType P = DfObject::getPpqMatrix<SymmetricMatrixType>(runType, iteration);

    // get maximum deviation
    P -= prevP;
    index_type i_max, j_max;
    const double deviation_value = P.getMaxAbsoluteElement(&i_max, &j_max);

    return deviation_value;
}


// kohn-sham matrix convergence
template<class SymmetricMatrixType>
double DfConvcheck::dev_kohn_sham_matrix(const RUN_TYPE runType, const int iteration)
{
    assert(iteration >= 2);
    SymmetricMatrixType prevF = DfObject::getFpqMatrix<SymmetricMatrixType>(runType, iteration -1);
    SymmetricMatrixType F = DfObject::getFpqMatrix<SymmetricMatrixType>(runType, iteration);

    F -= prevF;
    const double deviation_value = F.getMaxAbsoluteElement();

    return deviation_value;
}


template<class SymmetricMatrixType>
double DfConvcheck::dev_standard_dev_cd(const RUN_TYPE runType, const int iteration)
{
    assert(iteration >= 2);

    const index_type numOfAuxDens = this->m_nNumOfAux;
    double standard_deviation = 0.0;

    // cd coefficient convergence
    TlVector prevRho = DfObject::getRho<TlVector>(runType, iteration -1);
    TlVector rho = DfObject::getRho<TlVector>(runType, iteration);
    assert(prevRho.getSize() == this->m_nNumOfAux);
    assert(rho.getSize() == this->m_nNumOfAux);

    // calc DeltaRou
    TlVector deltaRho = rho - prevRho;
        
    SymmetricMatrixType Sab2 = DfObject::getSab2Matrix<SymmetricMatrixType>();
    assert(Sab2.getNumOfRows() == numOfAuxDens);
    assert(Sab2.getNumOfCols() == numOfAuxDens);
    TlVector v(numOfAuxDens);
    for (index_type p = 0; p < numOfAuxDens; ++p) {
        TlVector Sab2_v = Sab2.getRowVector(p);
        Sab2_v.dot(deltaRho);
        v += Sab2_v;
    }
        
    // calc DeltaRho * (Sab2 * DeltaRho^dagger)
    v.dot(deltaRho);
    const double IntegSquDeltaRho = v.sum();
    
    // calc standard_deviation
    standard_deviation = std::sqrt(std::fabs(IntegSquDeltaRho)) / (double)this->m_nNumOfAux;
    
    this->log_.info(TlUtils::format("IntegSquDeltaRho         is %14.4le", IntegSquDeltaRho));
    this->log_.info(TlUtils::format("standard deviation of cd is %14.4le", standard_deviation));

    return standard_deviation;
}

#endif // DFCONVCHECK_H
