// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef DFERIENGINE_H
#define DFERIENGINE_H

#include <cassert>
#include <vector>
#include <map>
#include <set>
#include <bitset>
#include <cmath>

#include "DfEngineObject.h"
#include "TlPosition.h"
#include "TlMath.h"
#include "TlAngularMomentumVector.h"
#include "TlAngularMomentumVectorSet.h"
#include "TlUtils.h"
#include "TlMath.h"
#include "TlOrbitalInfoObject.h"
// #include "TlTime.h"
#include "TlLogging.h"

/// 仕様
///
/// 変数名・式は主にJ. Phys. Chem. 94, 5564-5572 (1990) に基づく。
///
/// 変数説明:
/// (unnormalized) primitive Cartesian Gaussian Function
///   psi(r) = (x-Ax)^ax * (y-Ay)^ay * (z-Az)^az * exp(-alpha*(r-A)^2)
///   [a = (a_x, a_y, a_z); angular momentum vector] (a_x + a_y + a_z; angular momentum)
///   [A; position vector]
///   [alpha; exponent]
/// 
/// 制限値:
/// - angular momentum L
///   L = a + b + c + d + a_bar + b_bar + c_bar + d_bar
///   [a, b, c, d, a_bar, b_bar, c_bar, d_bar; angular momentum]
///   d型(=2)と2回微分(=2)をサポート =>  L_max = 2 * 4 + 2 = 10 あれば十分。
///   f型(=3)と1回微分(=1)をサポート =>  L_max = 3 * 4 + 1 = 13 あれば十分。
///   コード中は L_MAX で定義される。
///

/// query上限値
// #define A_BAR_MAX (2) // 2階微分
// #define B_BAR_MAX (2) // 2階微分
// #define A_MAX (2) // s=0, p=1, d=2
// #define B_MAX (2) // s=0, p=1, d=2

/// 本プログラムで扱うことのできる angular momentum の最大値
#define ERI_L_MAX ((ERI_A_BAR_MAX + ERI_B_BAR_MAX + ERI_A_MAX + ERI_B_MAX)*2)

/// 縮約数(4中心)の積の最大値 
#define ERI_KPKQ_MAX (100)

/// 最大Angular Momentumのvector数
#define ERI_NUM_OF_AMVS (((ERI_L_MAX +1) +1)*((ERI_L_MAX +1) +2) / 2)

/// [0](m)から[L](m)までの総要素: (L+1)*(L+2)*(L+3)/6)
/// 基数が0であるため+1が必要-> (L+1)
#define ERI_NUM_OF_R_KINDS (((ERI_L_MAX +1) +1)*((ERI_L_MAX +1) +2)*((ERI_L_MAX +1) +3)/6)

/// [0](0)~[L](m)
#define ERI_NUM_OF_RM_KINDS (ERI_NUM_OF_R_KINDS * ((ERI_L_MAX +1) +1))

/// a^の最大数 (1回微分までサポート)
#define ERI_A_BAR_MAX (1)
/// b^の最大数 (1回微分までサポート)
#define ERI_B_BAR_MAX (1)
/// aの最大数 (aのf(=3)とbのf,そしてその1回微分(+1)までサポート)
#define ERI_A_MAX (7)
/// bの最大数 (bのf(=3)をサポート)
#define ERI_B_MAX (3)
/// pの最大数
#define ERI_P_MAX (7)
/// a'の最大数
#define ERI_A_PRIME_MAX (1)
/// b'の最大数
#define ERI_B_PRIME_MAX (7)
/// p'の最大数
#define ERI_P_PRIME_MAX (7)

/// 制限値を調べるときは CHECK_MAX_COUNTを立てること。
// #define CHECK_MAX_COUNT
#define ERI_NR_DASH_SIZE (31360 +1)
#define ERI_NUM_OF_ERI_STATES ((ERI_A_BAR_MAX +1) * (ERI_B_BAR_MAX +1)  \
                               * (ERI_A_MAX +1) * (ERI_B_MAX +1) * (ERI_P_MAX +1) \
                               * (ERI_A_PRIME_MAX +1) * (ERI_B_PRIME_MAX +1) * (ERI_P_PRIME_MAX +1))
//#define ERI_NUM_OF_AMVS (108 +1)
#define ERI_MAX_BATCH (5880 +1)

class DfEriEngine : public DfEngineObject {
public:
    typedef int index_type;
    
    struct AngularMomentum2 {
        AngularMomentum2(int ibar =0, int jbar =0, int i =0, int j =0)
            : a_bar(ibar), b_bar(jbar), a(i), b(j) {
            assert(ibar <= ERI_A_BAR_MAX);
            assert(jbar <= ERI_B_BAR_MAX);
            assert(i <= ERI_A_MAX);
            assert(j <= ERI_B_MAX);
        }

        int sum() const {
            return a_bar + b_bar + a + b;
        };

        int index() const {
            return ((a_bar * (ERI_B_BAR_MAX +1) +b_bar) * (ERI_A_MAX +1) + a) * (ERI_B_MAX +1) + b;
        }

        static int maxIndex() {
            return (ERI_A_BAR_MAX +1) * (ERI_B_BAR_MAX +1) * (ERI_A_MAX +1) * (ERI_B_MAX +1);
        }
        
    public:
        // 8bitで-31 ~ +32 まで OK
        int a_bar; // grad i
        int b_bar; // grad j
        int a;
        int b;
    };

    // for pGTO ----------------------------------------------------------------
    class PGTO_Pair {
    public:
        // default constructer for container template
        PGTO_Pair()
            : alpha2_(0.0), beta2_(0.0), zeta2_(0.0),
              sigma_P_(0.0), P_(0.0, 0.0, 0.0), U_P_(0.0) {
        };

        PGTO_Pair(const int a, const double coefA, const double expA,
                  const int b, const double coefB, const double expB,
                  const TlPosition& A,
                  const TlPosition& B)
            : alpha2_(2.0 * expA), beta2_(2.0 * expB) {
            const static double U_COEF = std::sqrt(8.0 * M_PI * M_PI * M_PI);
            this->zeta2_ = this->alpha2_ + this->beta2_;
            this->sigma_P_ = 1.0 / this->zeta2_;
            this->P_ = (this->alpha2_ * A + this->beta2_ * B) * this->sigma_P_;
            this->U_P_ = U_COEF * std::pow(this->sigma_P_, (a + b + 1.5))
                * coefA * coefB * std::exp(-2.0 * expA * expB * this->sigma_P_ * A.squareDistanceFrom(B));
        }

    public:
        double alpha2() const {
            return this->alpha2_;
        }

        double beta2() const {
            return this->beta2_;
        }
        
        double zeta2() const {
            return this->zeta2_;
        }

        double sigma_P() const {
            return this->sigma_P_;
        }

        TlPosition P() const {
            return this->P_;
        }

        double U_P() const {
            return this->U_P_;
        }
            
    private:
        // double alpha_;
        // double beta_;
        double alpha2_; // = alpha * 2
        double beta2_;  // = beta * 2
        double zeta2_; // = alpha *2 + beta *2
        double sigma_P_; // eq.(20)
        TlPosition P_;   // eq.(21)
        double U_P_;     // eq.(22)
    };
    typedef std::vector<PGTO_Pair> PGTO_Pairs;

    // PGTO_Pairsを比較するためのfunctor
    struct PGTO_sort_functor_cmp {
        bool operator()(const PGTO_Pair& a, const PGTO_Pair& b) const {
            return (std::fabs(a.U_P()) > std::fabs(b.U_P()));
        }
    };

    // for cGTO ----------------------------------------------------------------
    struct CGTO_Pair {
    public:
        CGTO_Pair(const TlPosition& ab, const PGTO_Pairs& ps) : AB(ab), PS(ps) {
        }
        
    public:
        TlPosition AB;
        PGTO_Pairs PS;
    };

    // Elementary 4-Center Quantities ------------------------------------------
private:
    struct E4 {
    public:
        E4() {
        }
        
        E4(const TlPosition& P, const TlPosition& Q,
             double sigma_P, double sigma_Q,
             double U_P, double U_Q) {
            this->R = Q - P;
            this->R2 = this->R.squareDistanceFrom();
            this->_2theta2 = 1.0 / (sigma_P + sigma_Q);
            this->_2T = this->_2theta2 * this->R2;
            this->U = U_P * U_Q;
        }

    public:
        TlPosition R;
        double R2;
        double _2theta2;
        double _2T;
        double U;
    };
    

    // for Contract ------------------------------------------------------------
private:
    struct ContractScale {
    public:
        explicit ContractScale(const int ap =0, const int bp =0, const int pp =0)
            : a_prime(ap), b_prime(bp), p_prime(pp) {
        }

        int index() const {
            const int v = ((a_prime * ERI_B_PRIME_MAX) + b_prime) * ERI_P_PRIME_MAX + p_prime;
            return v;
        }

    public:
        // int a_prime : 6; // -31 ~ +32 まで OK
        // int b_prime : 6;
        // int p_prime : 6;
        int a_prime;
        int b_prime;
        int p_prime;
    };

    struct ContractScale_cmp {
        bool operator()(const ContractScale& rhs1, const ContractScale& rhs2) const {
            return (rhs1.index() < rhs2.index());
        }
    };
    
    typedef std::set<ContractScale, ContractScale_cmp> ContractScalesSet;
    typedef std::vector<ContractScale> ContractScalesVector;    

    
    struct ContractState {
    public:
        explicit ContractState(const TlAngularMomentumVector& R = TlAngularMomentumVector(),
                               const int ap =0, const int bp =0, const int pp =0,
                               const int cp =0, const int dp =0, const int qp =0)
            : r_index(R.index()), r(R.angularMomentum()),
              a_prime(ap), b_prime(bp), p_prime(pp),
              c_prime(cp), d_prime(dp), q_prime(qp),
              isCached_index_(false), cached_index_(0) {
            assert(this->r < (1 << 4));
            assert(this->r_index < (1 << 6));
            
            assert(ap < (1 << 3));
            assert(bp < (1 << 3));
            assert(pp < (1 << 3));
            assert(cp < (1 << 3));
            assert(dp < (1 << 3));
            assert(qp < (1 << 3));
        }

        ContractState(const ContractState& rhs)
            : r_index(rhs.r_index), r(rhs.r),
              a_prime(rhs.a_prime), b_prime(rhs.b_prime), p_prime(rhs.p_prime),
              c_prime(rhs.c_prime), d_prime(rhs.d_prime), q_prime(rhs.q_prime),
              isCached_index_(rhs.isCached_index_), cached_index_(rhs.cached_index_) {
        }

    public:
        void setABP(const int ap, const int bp, const int pp) {
            this->a_prime = ap;
            this->b_prime = bp;
            this->p_prime = pp;
            this->isCached_index_ = false;
        }

        void setCDQ(const int cp, const int dp, const int qp) {
            this->c_prime = cp;
            this->d_prime = dp;
            this->q_prime = qp;
            this->isCached_index_ = false;
        }

        void setR(const TlAngularMomentumVector& R) {
            this->r_index = R.index();
            this->r = R.angularMomentum();
            this->isCached_index_ = false;
        }

        int getCprime() const {
            return this->c_prime;
        }
        int getDprime() const {
            return this->d_prime;
        }
        int getQprime() const {
            return this->q_prime;
        }
        
        std::string debugOut() const {
            return TlUtils::format("[%u(%u); (%u %u %u|%u %u %u)]",
                                   r, r_index,
                                   a_prime, b_prime, p_prime,
                                   c_prime, d_prime, q_prime);
        }

        // TODO: hotspot
        int index() const {
            if (this->isCached_index_ != true) {
                this->cached_index_ =
                    (((((( this->a_prime)*(ERI_B_PRIME_MAX +1)
                         + this->b_prime)*(ERI_P_PRIME_MAX +1)
                         + this->p_prime)*(ERI_A_PRIME_MAX +1)
                         + this->c_prime)*(ERI_B_PRIME_MAX +1)
                         + this->d_prime)*(ERI_P_PRIME_MAX +1)
                         + this->q_prime)* ERI_NUM_OF_AMVS;

                const int r = this->r;
                this->cached_index_ += (r)*(r+1)*(r+2)/6 + this->r_index;
                this->isCached_index_ = true;
            }

            return this->cached_index_;
        }
        
    private:
        // unsigned int r_index : 16; // -32767 ~ +32768 までOK
        // unsigned int r       :  6; // -31 ~ +32 まで OK
        // unsigned int a_prime :  6; // -7 ~ +8 まで OK
        // unsigned int b_prime :  6;
        // unsigned int p_prime :  6;
        // unsigned int c_prime :  6;
        // unsigned int d_prime :  6;
        // unsigned int q_prime :  6;
        int r_index;
        int r;
        int a_prime;
        int b_prime;
        int p_prime;
        int c_prime;
        int d_prime;
        int q_prime;

    private:
        // mutable bool isCached_index_abp_;
        // mutable std::size_t cached_index_abp_;
        // mutable bool isCached_index_cdq_;
        // mutable std::size_t cached_index_cdq_;
        // mutable bool isCached_index_r_;
        // mutable std::size_t cached_index_r_;
        mutable bool isCached_index_;
        mutable int cached_index_;
    };

    struct nR_dash {
    public:
        nR_dash(const ContractState& inCS = ContractState(),
                const std::vector<double>& inValues = std::vector<double>())
            : cs(inCS), values(inValues) {
        }

    public:
        ContractState cs;
        std::vector<double> values;
    };
    typedef std::vector<nR_dash> nR_dash_Type;

    // recursive relation state ------------------------------------------------
private:
    struct ERI_State {
    public:
        explicit ERI_State(int a_bar_in =0, int b_bar_in =0,
                           int a_in =0, int b_in =0, int p_in =0,
                           int a_prime_in =0, int b_prime_in =0, int p_prime_in =0)
            : a_bar(a_bar_in), b_bar(b_bar_in), a(a_in), b(b_in), p(p_in),
            a_prime(a_prime_in), b_prime(b_prime_in), p_prime(p_prime_in) {
        }

        // for debug
        std::string debugOut() const {
            return TlUtils::format("(%d %d, %d %d %d, %d %d %d|",
                                   a_bar, b_bar, a, b, p,
                                   a_prime, b_prime, p_prime);
        }

        int index() const {
            const int index =
                ((((((( this->a_bar)  *(ERI_B_BAR_MAX +1)
                      + this->b_bar)  *(ERI_A_MAX +1)
                      + this->a)      *(ERI_B_MAX +1)
                      + this->b)      *(ERI_P_MAX +1)
                      + this->p)      *(ERI_A_PRIME_MAX +1)
                      + this->a_prime)*(ERI_B_PRIME_MAX +1)
                      + this->b_prime)*(ERI_P_PRIME_MAX +1)
                      + this->p_prime;
            assert(index < ERI_NUM_OF_ERI_STATES);
            return index;
        }
        
    public:
        int a_bar;
        int b_bar;
        int a;
        int b;
        int p;
        int a_prime;
        int b_prime;
        int p_prime;
    };

    // struct ERI_State_cmp {
    //     bool operator()(const ERI_State& rhs1, const ERI_State& rhs2) const {
    //         return (rhs1.index() < rhs2.index());
    //     }
    // };

    typedef std::vector<std::vector<std::vector<double> > > EriDataType;

        
    // construct & destruct ----------------------------------------------------
public:
    DfEriEngine();
    virtual ~DfEriEngine();

    // double getElapseCalcTime() const;

    static CGTO_Pair getCGTO_pair(const TlOrbitalInfoObject& orbInfo,
                                  const index_type shellIndexP,
                                  const index_type shellIndexQ = -1,
                                  const double threshold = 1.0E-20);
    static CGTO_Pair getCGTO_pair(const TlOrbitalInfoObject& orbInfo1,
                                  const TlOrbitalInfoObject& orbInfo2,
                                  const index_type shellIndex1,
                                  const index_type shellIndex2 = -1,
                                  const double threshold = 1.0E-20);
    
    virtual void calc(const int diff1, const TlOrbitalInfoObject& orbInfo1, const index_type shell1,
                      const int diff2, const TlOrbitalInfoObject& orbInfo2, const index_type shell2,
                      const int diff3, const TlOrbitalInfoObject& orbInfo3, const index_type shell3,
                      const int diff4, const TlOrbitalInfoObject& orbInfo4, const index_type shell4);

    virtual double value(const index_type index) const;
    
    void calcGrad(const AngularMomentum2& qAB, const AngularMomentum2& qCD,
                  const CGTO_Pair& IJ, const CGTO_Pair& KL);

private:
    void initialize();
    
    void calc0(const AngularMomentum2& qAB, const AngularMomentum2& qCD,
               const CGTO_Pair& IJ, const CGTO_Pair& KL);

    void calc(const AngularMomentum2& qAB, const AngularMomentum2& qCD);

    void calcGrad(const AngularMomentum2& qAB, const AngularMomentum2& qCD);

    void calcGrad_sub(const AngularMomentum2& qAB, const AngularMomentum2& qCD);
    
    void copyResultsToOutputBuffer(const AngularMomentum2& qAB,
                                   const AngularMomentum2& qCD,
                                   double* pOutput);
    
    void calcE4CQ();

    // [0]^(m)テーブルを作成する ===============================================
    void calc0m();

    // [r]^(0)テーブルを作成する ===============================================
    void calcR0();
    void calcRM(const TlAngularMomentumVector& r, const int m);
    int indexRM(const TlAngularMomentumVector& amv, const int m) const;
    int initiativeRM(const TlAngularMomentumVector& amv) const;

    // contract ================================================================
    ContractScalesVector choice(const AngularMomentum2& AB);
    ContractScalesVector choice(const AngularMomentum2& AB1, const AngularMomentum2& AB2);
    
    void choice(const int a_bar, const int b_bar,
                const int a, const int b, const int p,
                const int a_prime, const int b_prime, const int p_prime,
                ContractScalesSet* pContractList);
    // int index_contract(const int a_prime, const int b_prime, const int p_prime) const;

    ContractScalesVector transContractScales_SetToVector(const ContractScalesSet& contractScales);

    void contract(const AngularMomentum2& qAB,
                  const AngularMomentum2& qCD,
                  const ContractScalesVector& bra_contractScales,
                  const ContractScalesVector& ket_contractScales);
    void contract_bra(const AngularMomentum2& qAB,
                      const TlAngularMomentumVector& r,
                      const int a_prime, const int b_prime, const int p_prime,
                      const int nR_dash_index);
    // void contract_ket(const AngularMomentum2& qCD,
    //                   const ContractState& cs, const std::vector<double>& KQ_values);
    void get_contract_ket_coef_numerators(const int c_prime,
                                          const int d_prime,
                                          const int zeta_exp,
                                          std::vector<double>* pCoefNumerators);
    // void contract_ket(const ContractState& cs,
    //                   const std::vector<double>& coef_numerators,
    //                   const std::vector<double>& KQ_values);
    
    // calc PQ
    void calcPQ(const AngularMomentum2& qAB,
                const AngularMomentum2& qCD,
                const ContractScalesVector& bra_contractScales,
                const ContractScalesVector& ket_contractScales);

    /// calcPQ内で利用するcs_indexを返す
    /// 
    /// 目的は高速化。cs_indexの格納順序はcalcPQ内の多重ループに従う。
    std::vector<int> 
    get_csindex_for_calcPQ(const int angularMomentumP,
                           const AngularMomentum2& qAB,
                           const int angularMomentumQ,
                           const AngularMomentum2& qCD,
                           const ContractScalesVector& bra_contractScales,
                           const ContractScalesVector& ket_contractScales);

    void transpose(const AngularMomentum2& qAB,
                   const AngularMomentum2& qCD,
                   const ContractScalesVector& ket_contractScales);
    
    // ERI ---------------------------------------------------------------------
    void calcERI(const int a_bar, const int b_bar,
                 const int a, const int b, const int p,
                 const int a_prime, const int b_prime, const int p_prime,
                 EriDataType* pERI);

    int index(const TlAngularMomentumVector& a_bar, const TlAngularMomentumVector& b_bar,
              const TlAngularMomentumVector& a, const TlAngularMomentumVector& b,
              const TlAngularMomentumVector& p) const;

    /// eq.44
    void ERI_EQ44(const ERI_State eriState, EriDataType* pERI);

    /// optimized EQ44
    ///
    /// ERI_EQ44() in case of "a^=b^=0"
    void ERI_EQ44_00xxx(const ERI_State eriState, EriDataType* pERI);

    /// eq.45
    void ERI_EQ45(const ERI_State eriState, EriDataType* pERI);

    /// optimized EQ45
    ///
    /// ERI_EQ45() in case of "a^=b^=b=p=0"
    void ERI_EQ45_00x00(const ERI_State eriState, EriDataType* pERI);

    /// optimized EQ45
    ///
    /// ERI_EQ45() in case of "a^=b^=b=0"
    void ERI_EQ45_00x0x(const ERI_State eriState, EriDataType* pERI);

    void ERI_EQ46(const ERI_State eriState, EriDataType* pERI); // 非常用
    void ERI_EQ43(const ERI_State eriState, EriDataType* pERI);
    void ERI_EQ47(const ERI_State eriState, EriDataType* pERI);

    // transform 6D to 5D
    void transform6Dto5D(const AngularMomentum2& qAB,
                         const AngularMomentum2& qCD,
                         double* pOutput);
    void transform6Dto5D_i(const int I_, const int J_, const int K_, const int L_,
                           const int J, const int K, const int L,
                           const double* pInput, double* pOutput);
    void transform6Dto5D_j(const int I_, const int J_, const int K_, const int L_,
                           const int I, const int K, const int L,
                           const double* pInput, double* pOutput);
    void transform6Dto5D_k(const int I_, const int J_, const int K_, const int L_,
                           const int I, const int J, const int L,
                           const double* pInput, double* pOutput);
    void transform6Dto5D_l(const int I_, const int J_, const int K_, const int L_,
                           const int I, const int J, const int K,
                           const double* pInput, double* pOutput);

    void transform10Fto7F_i(const int I_, const int J_, const int K_, const int L_,
                           const int J, const int K, const int L,
                           const double* pInput, double* pOutput);
    void transform10Fto7F_j(const int I_, const int J_, const int K_, const int L_,
                           const int I, const int K, const int L,
                           const double* pInput, double* pOutput);
    void transform10Fto7F_k(const int I_, const int J_, const int K_, const int L_,
                           const int I, const int J, const int L,
                           const double* pInput, double* pOutput);
    void transform10Fto7F_l(const int I_, const int J_, const int K_, const int L_,
                           const int I, const int J, const int K,
                           const double* pInput, double* pOutput);

    void compD(const AngularMomentum2& qAB,
               const AngularMomentum2& qCD);
    
public:
    double* WORK;
    double* WORK_A;
    double* WORK_B;
    double* WORK_C;
    double* WORK_D;
    
private:
    static const int OUTPUT_BUFFER_SIZE;
    
    // E1 = (1, 0, 0), (0, 1, 0), (0, 0, 1);
    static const TlAngularMomentumVector E1_[3];

    // E2 = (2, 0, 0), (0, 2, 0), (0, 0, 2);
    static const TlAngularMomentumVector E2_[3];

    static const double INV_SQRT3;
    
    TlPosition AB_;
    TlPosition CD_;
    
    int sumOfAngularMomentums_; // sum of angular momentums
    PGTO_Pairs bra_;
    PGTO_Pairs ket_;
    
    /// Gm(T)を格納するテーブル
    double* pGmT_;
    
    /// 長さを持つE4CQ構造体のリスト
    E4* pE4_;

    /// [0]^(m) のテーブル
    /// m の最大値は L_MAX
    /// インデックスは index = m * KPKQ_MAX + KPKQ;
    /// すなわち、KPKQの方向が隣同士。
    double** p0M_;
    
    // for [r](0) --------------------------------------------------------------
    /// [r]^(m) のテーブル
    double** pRM_;
    
    /// [r]^(m) をすでに計算したかどうかを保持するテーブル
    std::bitset<ERI_NUM_OF_RM_KINDS> calcdRM_;

    // for contract ------------------------------------------------------------
    std::map<int, ContractScalesVector> choice_tbl_;
    
    double* pContractBraCoef_;
    
    //n_R_dash_Type n_R_dash_;
    nR_dash_Type nR_dash_;

    /// abp(r)cdq を格納
    /// indexはContractState::index()で求める。
    double* p_abpRcdq_; 

    // for calcPQ
    std::map<int, std::vector<int> > csindeces_forCalcPQ_;
    
    // for recursive relations
    EriDataType ERI_bra_;
    EriDataType ERI_ket_;
    
    // std::map<ERI_State, bool, ERI_State_cmp> isCalcdERI_;
    std::vector<int> isCalcdERI_;
    int ERI_batch_;

    /// 6D->5D変換用バッファ
    double* pTransformBuf_;

    //
    // TlTime time_calc_all_;
    TlLogging& log_;

#ifdef CHECK_MAX_COUNT
    int maxSizeOf_sumOfAngularMomentums_;
    int max_a_bar_;
    int max_b_bar_;
    int max_a_;
    int max_b_;
    int max_p_;
    int max_a_prime_;
    int max_b_prime_;
    int max_p_prime_;
    int maxSizeOf_nR_dash_;
    int maxNumOfAMVs_;
    int maxERI_batch_;
#endif //
};


inline double DfEriEngine::value(const index_type index) const
{
    return this->WORK[index];
}


#endif // DFERIENGINE
