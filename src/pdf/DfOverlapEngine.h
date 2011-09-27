#ifndef DFOVERLAPENGINE_H
#define DFOVERLAPENGINE_H

#include <vector>
#include <bitset>
#include "TlAngularMomentumVector.h"
#include "TlPosition.h"

/// 結果出力用のサイズ
/// d軌道の6Dから5Dへの変換領域にも利用するため、
/// d軌道のサイズを6Dベースとして確保している。
//#define OVP_OUTPUT_BUFFER_SIZE (5832) // = 3 * 3 * 3 * 6(d) * 6 * 6

/// a^の最大数
#define OVP_A_BAR_MAX (2 +1)
/// b^の最大数
#define OVP_B_BAR_MAX (2 +1)
/// c^の最大数
#define OVP_C_BAR_MAX (2 +1)
/// aの最大数 (d) +1次微分(1) +"0スタート"
#define OVP_A_MAX (6 +1)
/// bの最大数
#define OVP_B_MAX (6 +1)
/// cの最大数
#define OVP_C_MAX (6 +1)
///
#define OVP_STATE_MAX (OVP_A_BAR_MAX * OVP_B_BAR_MAX * OVP_C_BAR_MAX * OVP_A_MAX * OVP_B_MAX * OVP_C_MAX)

/// see. S.Obara and A.Saika, J. Chem. Phys., 84, 3963 (1986)
class DfOverlapEngine {
public:
    struct Query {
    public:
        Query(int a_bar, int b_bar, int c_bar, int a, int b, int c)
            : shellTypeA_bar(a_bar), shellTypeB_bar(b_bar), shellTypeC_bar(c_bar),
              shellTypeA(a), shellTypeB(b), shellTypeC(c){
        }
        
    public:
        int shellTypeA_bar : 8;
        int shellTypeB_bar : 8;
        int shellTypeC_bar : 8;
        int shellTypeA : 8; // -31 ~ +32 まで OK
        int shellTypeB : 8;
        int shellTypeC : 8;
    };

    struct PGTO {
    public:
        explicit PGTO(double c = 0.0, double z = 0.0) : coef(c), zeta(z) {
        };

    public:
        double coef;
        double zeta;
    };

    typedef std::vector<PGTO> PGTOs;

    struct OvpState {
        OvpState(int inAbar, int inBbar, int inCbar,
                 int inA, int inB, int inC)
            : a_bar(inAbar), b_bar(inBbar), c_bar(inCbar), a(inA), b(inB), c(inC) {
            assert((0 <= a_bar) && (a_bar < OVP_A_BAR_MAX));
            assert((0 <= b_bar) && (b_bar < OVP_B_BAR_MAX));
            assert((0 <= c_bar) && (c_bar < OVP_C_BAR_MAX));
            assert((0 <= a) && (a < OVP_A_MAX));
            assert((0 <= b) && (b < OVP_B_MAX));
            assert((0 <= c) && (c < OVP_C_MAX));
        }

        int index() const {
            const int answer = ((((a_bar *OVP_B_BAR_MAX
                                  +b_bar)*OVP_C_BAR_MAX
                                  +c_bar)*OVP_A_MAX
                                  +a)    *OVP_B_MAX
                                  +b)    *OVP_C_MAX +c;
            return answer;
        }
        
    public:
        int a_bar : 8;
        int b_bar : 8;
        int c_bar : 8;
        int a : 8;
        int b : 8;
        int c : 8;
    };

    typedef std::vector<std::vector<std::vector<double> > > OvpDataType;
    
public:
    DfOverlapEngine();
    ~DfOverlapEngine();

    void calc(const Query& query,
              const TlPosition& A,
              const TlPosition& B,
              const TlPosition& C,
              const PGTOs& PGTOs_A,
              const PGTOs& PGTOs_B,
              const PGTOs& PGTOs_C);

private:
    void calc(const int aBar, const int bBar, const int cBar,
              const int a, const int b, const int c);

    int index(const TlAngularMomentumVector& amvAbar,
              const TlAngularMomentumVector& amvBbar,
              const TlAngularMomentumVector& amvCbar,
              const TlAngularMomentumVector& amvA,
              const TlAngularMomentumVector& amvB,
              const TlAngularMomentumVector& amvC) const;

    void copyResultsToOutputBuffer();

    void transform6Dto5D();
    void transform6Dto5D_i(const int I_, const int J_, const int K_,
                           const int J, const int K,
                           const double* pInput, double* pOutput);
    void transform6Dto5D_j(const int I_, const int J_, const int K_,
                           const int I, const int K,
                           const double* pInput, double* pOutput);
    void transform6Dto5D_k(const int I_, const int J_, const int K_,
                           const int I, const int J,
                           const double* pInput, double* pOutput);

    int initiative(const TlAngularMomentumVector& amv) const;

    void EQ20A(const OvpState& state);
    void EQ20B(const OvpState& state);
    void EQ20C(const OvpState& state);

    void gradA(const OvpState& state);
    void gradB(const OvpState& state);
    void gradC(const OvpState& state);
    
public:
    double* WORK; 
    
private:
    static const double MAX_EXPONENT;

    static const int OUTPUT_BUFFER_SIZE;
//     static const int A_BAR_MAX;
//     static const int B_BAR_MAX;
//     static const int C_BAR_MAX;
//     static const int A_MAX;
//     static const int B_MAX;
//     static const int C_MAX;
//     static const int OVP_STATE_MAX;
    
    // E1 = (1, 0, 0), (0, 1, 0), (0, 0, 1);
    static const TlAngularMomentumVector E1_[3];

    // E2 = (2, 0, 0), (0, 2, 0), (0, 0, 2);
    static const TlAngularMomentumVector E2_[3];

    static const double INV_SQRT3;

    int a_bar_;
    int b_bar_;
    int c_bar_;
    int a_;
    int b_;
    int c_;
    TlPosition A_;
    TlPosition B_;
    TlPosition C_;

    OvpDataType OVP_;
    std::vector<TlPosition> G_;
    std::vector<double> coef_;
    std::bitset<OVP_STATE_MAX> isCalcd_;

    std::vector<double> zetaA_;
    std::vector<double> zetaB_;
    std::vector<double> zetaC_;
    int numOfBatches_;
};


#endif // DFOVERLAPENGINE_H
