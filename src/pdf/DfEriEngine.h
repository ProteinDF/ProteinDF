#ifndef DFERIENGINE_H
#define DFERIENGINE_H

#include <cassert>
#include <vector>
#include <map>
#include <set>
#include <bitset>
#include <cmath>

#include "TlPosition.h"
#include "TlMath.h"
#include "TlAngularMomentumVector.h"
#include "TlAngularMomentumVectorSet.h"
#include "TlUtils.h"
#include "TlMath.h"
#include "TlOrbitalInfoObject.h"

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

/// 本プログラムで扱うことのできる angular momentum の最大値
#define ERI_L_MAX (9 +1)

/// 縮約数(4中心)の積の最大値 
#define ERI_KPKQ_MAX (100)

/// 最大Angular Momentumのvector数
#define ERI_NUM_OF_AMVS ((ERI_L_MAX+1)*(ERI_L_MAX+2) / 2)

/// [0](m)から[L](m)までの総要素: (L+1)*(L+2)*(L+3)/6)
/// 基数が0であるため+1が必要-> (L+1)
#define ERI_NUM_OF_R_KINDS ((ERI_L_MAX +1)*(ERI_L_MAX +2)*(ERI_L_MAX +2)/6)

/// [0](0)~[L](m)
#define ERI_NUM_OF_RM_KINDS (ERI_NUM_OF_R_KINDS * (ERI_L_MAX +1))

/// a^の最大数 (d(=2)の微分までサポート)
#define ERI_A_BAR_MAX (3 +1)
/// b^の最大数 (d(=2)の微分までサポート)
#define ERI_B_BAR_MAX (3 +1)
/// aの最大数 (aのd(=2)とbのd,そしてその2回微分(+2)までサポート)
#define ERI_A_MAX (7 +1)
/// bの最大数 (bのd(=2)をサポート)
#define ERI_B_MAX (3 +1)
/// pの最大数
#define ERI_P_MAX (7 +1)
/// a'の最大数
#define ERI_A_PRIME_MAX (3 +1)
/// b'の最大数
#define ERI_B_PRIME_MAX (7 +1)
/// p'の最大数
#define ERI_P_PRIME_MAX (7 +1)

/// 制限値を調べるときは CHECK_MAX_COUNTを立てること。
//#define CHECK_MAX_COUNT
#define ERI_NR_DASH_SIZE (6820 +1)
#define ERI_NUM_OF_ERI_STATES (ERI_A_BAR_MAX * ERI_B_BAR_MAX \
                               * ERI_A_MAX * ERI_B_MAX * ERI_P_MAX \
                               * ERI_A_PRIME_MAX * ERI_B_PRIME_MAX * ERI_P_PRIME_MAX)
//#define ERI_NUM_OF_AMVS (108 +1)
#define ERI_MAX_BATCH (1456 +1)

class DfEriEngine {
public:
    typedef int index_type;
    
    struct Query {
        Query(int ibar =0, int jbar =0, int i =0, int j =0)
            : a_bar(ibar), b_bar(jbar), a(i), b(j) {
            assert(ibar < ERI_A_BAR_MAX);
            assert(jbar < ERI_B_BAR_MAX);
            assert(i < ERI_A_MAX);
            assert(j < ERI_B_MAX);
        }

        int sum() const {
            return a_bar + b_bar + a + b;
        };
        
    public:
        // int a_bar : 8; // grad i
        // int b_bar : 8; // grad j
        // int a : 8;
        // int b : 8;
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
            : alpha(0.0), beta(0.0), sigma_P(0.0), P(0.0, 0.0, 0.0), U_P(0.0) {
        };

        PGTO_Pair(const int a, const double coefA, const double expA,
                  const int b, const double coefB, const double expB,
                  const TlPosition& A,
                  const TlPosition& B)
            : alpha(expA), beta(expB) {
            const static double U_COEF = std::sqrt(8.0 * M_PI * M_PI * M_PI);
            //const double ab = expA + expB;
            this->sigma_P = 1.0 / (2.0 * (this->alpha + this->beta));
            this->P = (this->alpha * A + this->beta * B) * 2.0 * this->sigma_P;
            this->U_P = U_COEF * std::pow(this->sigma_P, (a + b + 1.5))
                * coefA * coefB * std::exp(-2.0 * this->alpha * this->beta * this->sigma_P * A.squareDistanceFrom(B));
        }

    public:
        double alpha;
        double beta;
        double sigma_P; // eq.(20)
        TlPosition P;   // eq.(21)
        double U_P;     // eq.(22)
    };
    typedef std::vector<PGTO_Pair> PGTO_Pairs;

    // PGTO_Pairsを比較するためのfunctor
    struct PGTO_sort_functor_cmp {
        bool operator()(const PGTO_Pair& a, const PGTO_Pair& b) const {
            return (std::fabs(a.U_P) > std::fabs(b.U_P));
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

        unsigned int index() const {
            const unsigned long v = ((a_prime * ERI_B_MAX) + b_prime) * ERI_P_MAX + p_prime;
            return v;
        }

    public:
        int a_prime : 6; // -31 ~ +32 まで OK
        int b_prime : 6;
        int p_prime : 6;
    };

    struct ContractScale_cmp {
        bool operator()(const ContractScale& rhs1, const ContractScale& rhs2) const {
            return (rhs1.index() < rhs2.index());
        }
    };
    
    typedef std::set<ContractScale, ContractScale_cmp> ContractScalesType;

    
    struct ContractState {
    public:
        explicit ContractState(const TlAngularMomentumVector& R = TlAngularMomentumVector(),
                               const unsigned int ap =0, const unsigned int bp =0, const unsigned int pp =0,
                               const unsigned int cp =0, const unsigned int dp =0, const unsigned int qp =0)
            : r_index(R.index()), r(R.angularMomentum()),
              a_prime(ap), b_prime(bp), p_prime(pp),
              c_prime(cp), d_prime(dp), q_prime(qp),
              isCached_index_abp_(false), cached_index_abp_(0),
              isCached_index_cdq_(false), cached_index_cdq_(0),
              isCached_index_r_(false), cached_index_r_(0),
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
              isCached_index_abp_(rhs.isCached_index_abp_), cached_index_abp_(rhs.cached_index_abp_),
              isCached_index_cdq_(rhs.isCached_index_cdq_), cached_index_cdq_(rhs.cached_index_cdq_),
              isCached_index_r_(rhs.isCached_index_r_), cached_index_r_(rhs.cached_index_r_),
              isCached_index_(rhs.isCached_index_), cached_index_(rhs.cached_index_) {
        }

    public:
        void setABP(const unsigned int ap, const unsigned int bp, const unsigned int pp) {
            this->a_prime = ap;
            this->b_prime = bp;
            this->p_prime = pp;
            this->isCached_index_abp_ = false;
        }

        void setCDQ(const unsigned int cp, const unsigned int dp, const unsigned int qp) {
            this->c_prime = cp;
            this->d_prime = dp;
            this->q_prime = qp;
            this->isCached_index_cdq_ = false;
        }

        void setR(const TlAngularMomentumVector& R) {
            this->r_index = R.index();
            this->r = R.angularMomentum();
            this->isCached_index_r_ = false;
        }

        unsigned int getCprime() const {
            return this->c_prime;
        }
        unsigned int getDprime() const {
            return this->d_prime;
        }
        unsigned int getQprime() const {
            return this->q_prime;
        }
        
        std::string debugOut() const {
            return TlUtils::format("[%u(%u); (%u %u %u|%u %u %u)]",
                                   r, r_index,
                                   a_prime, b_prime, p_prime,
                                   c_prime, d_prime, q_prime);
        }

        // TODO: hotspot
        std::size_t index() const {
            const bool check = (this->isCached_index_abp_ &
                                this->isCached_index_cdq_ &
                                this->isCached_index_r_ &
                                this->isCached_index_);
            if (check != true) {
                if (this->isCached_index_abp_ != true) {
                    this->cached_index_abp_ =
                        (((((( this->a_prime)*(ERI_B_PRIME_MAX)
                             + this->b_prime)*(ERI_P_PRIME_MAX)
                             + this->p_prime)*(ERI_A_PRIME_MAX)
                                            )*(ERI_B_PRIME_MAX)
                                            )*(ERI_P_PRIME_MAX)
                                            )*ERI_NUM_OF_AMVS;
                    this->isCached_index_abp_ = true;
                }
                if (this->isCached_index_cdq_ != true) {
                    this->cached_index_cdq_ =
                        ((( this->c_prime)*(ERI_B_PRIME_MAX)
                          + this->d_prime)*(ERI_P_PRIME_MAX)
                          + this->q_prime)*ERI_NUM_OF_AMVS;
                    this->isCached_index_cdq_ = true;
                }
                if (this->isCached_index_r_ != true) {
                    const std::size_t r = this->r;
                    this->cached_index_r_ = (r)*(r+1)*(r+2)/6 + this->r_index;
                    this->isCached_index_r_ = true;
                }

                this->cached_index_ =
                      this->cached_index_abp_
                    + this->cached_index_cdq_
                    + this->cached_index_r_;
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
        unsigned int r_index;
        unsigned int r;
        unsigned int a_prime;
        unsigned int b_prime;
        unsigned int p_prime;
        unsigned int c_prime;
        unsigned int d_prime;
        unsigned int q_prime;

    private:
        mutable bool isCached_index_abp_;
        mutable std::size_t cached_index_abp_;
        mutable bool isCached_index_cdq_;
        mutable std::size_t cached_index_cdq_;
        mutable bool isCached_index_r_;
        mutable std::size_t cached_index_r_;
        mutable bool isCached_index_;
        mutable std::size_t cached_index_;
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

        std::size_t index() const {
            const int index =
                ((((((( this->a_bar)  *(ERI_B_BAR_MAX)
                      + this->b_bar)  *(ERI_A_MAX)
                      + this->a)      *(ERI_B_MAX)
                      + this->b)      *(ERI_P_MAX)
                      + this->p)      *(ERI_A_PRIME_MAX)
                      + this->a_prime)*(ERI_B_PRIME_MAX)
                      + this->b_prime)*(ERI_P_PRIME_MAX)
                      + this->p_prime;
            return index;
        }
        
    public:
        int a_bar : 8;
        int b_bar : 8;
        int a : 8;
        int b : 8;
        int p : 8;
        int a_prime : 8;
        int b_prime : 8;
        int p_prime : 8;
    };

    struct ERI_State_cmp {
        bool operator()(const ERI_State& rhs1, const ERI_State& rhs2) const {
            return (rhs1.index() < rhs2.index());
        }
    };

    typedef std::vector<std::vector<std::vector<double> > > EriDataType;

        
    // construct & destruct ----------------------------------------------------
public:
    DfEriEngine();
    ~DfEriEngine();
    
    static CGTO_Pair getCGTO_pair(const TlOrbitalInfoObject& orbInfoconst,
                                  const index_type shellIndexP,
                                  const index_type shellIndexQ = -1,
                                  const double threshold = 1.0E-20);
    
    void calc(const Query& qAB, const Query& qCD,
              const CGTO_Pair& IJ, const CGTO_Pair& KL);

    void calcGrad(const Query& qAB, const Query& qCD,
                  const CGTO_Pair& IJ, const CGTO_Pair& KL);

    void setPrimitiveLevelThreshold(const double threshold) {
        this->primitiveLevelThreshold_ = threshold;
    }
    double getPrimitiveLevelThreshold() const {
        return this->primitiveLevelThreshold_;
    }
    
private:
    void initialize();
    
    void calc(const Query& qAB, const Query& qCD);

    void calcGrad(const Query& qAB, const Query& qCD);

    void calcGrad_sub(const Query& qAB, const Query& qCD);
    
    void copyResultsToOutputBuffer(const Query& qAB,
                                   const Query& qCD,
                                   double* pOutput);
    
    void calcE4CQ();

    // [0]^(m)テーブルを作成する ===============================================
    void calc0m();

    // [r]^(0)テーブルを作成する ===============================================
    void calcR0();
    void calcRM(const TlAngularMomentumVector& r, const int m);
    unsigned int indexRM(const TlAngularMomentumVector& amv, const int m) const;
    int initiativeRM(const TlAngularMomentumVector& amv) const;

    // contract ================================================================
    void choice(const int a_bar, const int b_bar,
                const int a, const int b, const int p,
                const int a_prime, const int b_prime, const int p_prime,
                ContractScalesType* pContractList);
    unsigned int index_contract(const int a_prime, const int b_prime, const int p_prime) const;

    void contract(const DfEriEngine::Query& qAB,
                  const DfEriEngine::Query& qCD);
    void contract_bra(const DfEriEngine::Query& qAB,
                      const TlAngularMomentumVector& r,
                      const int a_prime, const int b_prime, const int p_prime,
                      const std::size_t nR_dash_index);
    void contract_ket(const DfEriEngine::Query& qCD,
                      const ContractState& cs, const std::vector<double>& KQ_values);
    
    // calc PQ
    void calcPQ(const DfEriEngine::Query& qAB,
                const DfEriEngine::Query& qCD);
    void transpose(const DfEriEngine::Query& qAB,
                   const DfEriEngine::Query& qCD);
    
    // ERI
    void calcERI(const int a_bar, const int b_bar,
                 const int a, const int b, const int p,
                 const int a_prime, const int b_prime, const int p_prime,
                 EriDataType* pERI);

    int index(const TlAngularMomentumVector& a_bar, const TlAngularMomentumVector& b_bar,
              const TlAngularMomentumVector& a, const TlAngularMomentumVector& b,
              const TlAngularMomentumVector& p) const;
    
    void ERI_EQ44(const ERI_State eriState, EriDataType* pERI);
    void ERI_EQ45(const ERI_State eriState, EriDataType* pERI);
    void ERI_EQ46(const ERI_State eriState, EriDataType* pERI); // 非常用
    void ERI_EQ43(const ERI_State eriState, EriDataType* pERI);
    void ERI_EQ47(const ERI_State eriState, EriDataType* pERI);

    // transform 6D to 5D
    void transform6Dto5D(const DfEriEngine::Query& qAB,
                         const DfEriEngine::Query& qCD,
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

    void compD(const DfEriEngine::Query& qAB,
               const DfEriEngine::Query& qCD);
    
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
    ContractScalesType bra_contractScales_;
    ContractScalesType ket_contractScales_;

    double* pContractBraCoef_;
    
    //n_R_dash_Type n_R_dash_;
    nR_dash_Type nR_dash_;

    /// abp(r)cdq を格納
    /// indexはContractState::index()で求める。
    double* p_abpRcdq_; 
    
    // for recursive relations
    EriDataType ERI_bra_;
    EriDataType ERI_ket_;
    
    std::map<ERI_State, bool, ERI_State_cmp> isCalcdERI_;
    int ERI_batch_;

    // see J. Chem. Phys., 105, 2726 (1996), eq33
    double primitiveLevelThreshold_;

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
    std::size_t maxSizeOf_nR_dash_;
    int maxNumOfAMVs_;
    int maxERI_batch_;
#endif //
};


#endif // DFERIENGINE
