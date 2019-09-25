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

#ifndef DFHPQENGINE_H
#define DFHPQENGINE_H

#include <bitset>
#include <vector>
#include "TlAngularMomentumVector.h"
#include "TlAngularMomentumVectorSet.h"
#include "TlAtom.h"
#include "TlPosition.h"

/// 結果出力用のサイズ
/// d軌道の6Dから5Dへの変換領域にも利用するため、
/// d軌道のサイズを6Dベースとして確保している。
//#define HPQ_OUTPUT_BUFFER_SIZE (108) // = 6(d) * 6 * 3

/// a^の最大数
#define HPQ_A_BAR_MAX (1 + 1)
/// b^の最大数
#define HPQ_B_BAR_MAX (1 + 1)
/// aの最大数
#define HPQ_A_MAX (5 + 1)
/// bの最大数
#define HPQ_B_MAX (5 + 1)
/// Lの最大値
#define HPQ_L_MAX (10 + 1)
///
#define HPQ_OVP_STATE_MAX \
    (HPQ_A_BAR_MAX * HPQ_B_BAR_MAX * HPQ_A_MAX * HPQ_B_MAX)
#define HPQ_KIN_STATE_MAX \
    (HPQ_A_BAR_MAX * HPQ_B_BAR_MAX * HPQ_A_MAX * HPQ_B_MAX)
#define HPQ_NUC_STATE_MAX \
    (HPQ_A_BAR_MAX * HPQ_B_BAR_MAX * HPQ_A_MAX * HPQ_B_MAX * HPQ_L_MAX)

/// see. S.Obara and A.Saika, J. Chem. Phys., 84, 1963 (1986)
class DfHpqEngine {
   public:
    struct Query {
       public:
        Query(int a_bar, int b_bar, int a, int b)
            : shellTypeA_bar(a_bar),
              shellTypeB_bar(b_bar),
              shellTypeA(a),
              shellTypeB(b) {}

       public:
        int shellTypeA_bar : 6;
        int shellTypeB_bar : 6;
        int shellTypeA : 6;  // -31 ~ +32 まで OK
        int shellTypeB : 6;
    };

    struct PGTO {
       public:
        explicit PGTO(double c = 0.0, double z = 0.0) : coef(c), zeta(z){};

       public:
        double coef;
        double zeta;
    };

    typedef std::vector<PGTO> PGTOs;

    // OvpState ----------------------------------------------------------------
    struct OvpState {
        OvpState(int inA, int inB) : a(inA), b(inB) {}

        int index() const {
            const int answer = a * HPQ_B_MAX + b;
            return answer;
        }

       public:
        int a : 8;
        int b : 8;
    };
    typedef std::vector<std::vector<std::vector<double> > > OvpDataType;

    // KinState ----------------------------------------------------------------
    struct KinState {
        KinState(int inAbar, int inBbar, int inA, int inB)
            : a_bar(inAbar), b_bar(inBbar), a(inA), b(inB) {}

        int index() const {
            const int answer =
                ((a_bar * HPQ_B_BAR_MAX + b_bar) * HPQ_A_MAX + a) * HPQ_B_MAX +
                b;
            return answer;
        }

       public:
        int a_bar : 8;
        int b_bar : 8;
        int a : 8;
        int b : 8;
    };
    typedef std::vector<std::vector<std::vector<double> > > KinDataType;

    // NucState ----------------------------------------------------------------
    struct NucState {
        NucState(int inAbar, int inBbar, int inA, int inB, int inM = 0)
            : a_bar(inAbar), b_bar(inBbar), a(inA), b(inB), m(inM) {}

        int index() const {
            const int answer =
                (((a_bar * HPQ_B_BAR_MAX + b_bar) * HPQ_A_MAX + a) * HPQ_B_MAX +
                 b) *
                    HPQ_L_MAX +
                m;
            return answer;
        }

       public:
        int a_bar : 8;
        int b_bar : 8;
        int a : 8;
        int b : 8;
        int m : 8;
    };
    typedef std::vector<std::vector<std::vector<double> > > NucDataType;

   public:
    DfHpqEngine();
    ~DfHpqEngine();

    void calc(const Query& query, const TlPosition& A, const TlPosition& B,
              const PGTOs& PGTOsA, const PGTOs& PGTOsB,
              const std::vector<TlAtom>& Cs, const std::vector<TlAtom>& Xs);

    // for ESP
    void calc(const Query& query, const TlPosition& A, const TlPosition& B,
              const PGTOs& PGTOsA, const PGTOs& PGTOsB, const TlPosition& C);

    void calcKineticPart(const Query& query, const TlPosition& A,
                         const TlPosition& B, const PGTOs& PGTOsA,
                         const PGTOs& PGTOsB);
    void calcNuclearAttractionPart(const Query& query, const TlPosition& A,
                                   const TlPosition& B, const PGTOs& PGTOsA,
                                   const PGTOs& PGTOsB,
                                   const std::vector<TlAtom>& Cs);

   private:
    void prepare(const PGTOs& PGTOsA, const PGTOs& PGTOsB);
    void calcOvpSS();
    void calcKinSS();
    void calcNucSS();

    void calcOvp(const int a, const int b);
    void calcKin(const int a_bar, const int b_bar, const int a, const int b);
    void calcNuc(const int a_bar, const int b_bar, const int a, const int b,
                 const int m);

    void Ovp_EqA2_A(const OvpState& state);
    void Ovp_EqA2_B(const OvpState& state);

    void Kin_EqA12_A(const KinState& state);
    void Kin_EqA12_B(const KinState& state);
    void gradKinA(const KinState& state);
    void gradKinB(const KinState& state);

    void Nuc_EqA19_A(const NucState& state);
    void Nuc_EqA19_B(const NucState& state);
    void gradNucA(const NucState& state);
    void gradNucB(const NucState& state);

    int index(const TlAngularMomentumVector& amvA,
              const TlAngularMomentumVector& amvB) const;
    int index(const TlAngularMomentumVector& amvAbar,
              const TlAngularMomentumVector& amvBbar,
              const TlAngularMomentumVector& amvA,
              const TlAngularMomentumVector& amvB) const;

    void transform6Dto5D(double* pBuf);
    void transform6Dto5D_i(const int I_, const int J_, const int J,
                           const double* pInput, double* pOutput);
    void transform6Dto5D_j(const int I_, const int J_, const int I,
                           const double* pInput, double* pOutput);
    void transform10Fto7F_i(const int I_, const int J_, const int J,
                            const double* pInput, double* pOutput);
    void transform10Fto7F_j(const int I_, const int J_, const int I,
                            const double* pInput, double* pOutput);

    int initiative(const TlAngularMomentumVector& amv) const;
    void EQ20A(const OvpState& state);
    void EQ20B(const OvpState& state);
    void EQ20C(const OvpState& state);

    template <typename StateType, typename DataType>
    void addResultsToOutputBuffer(const DataType& fromBuf, double* toBuf);

   public:
    double* WORK_KIN;
    double* WORK_NUC;
    double* WORK_NUCX;

   private:
    static const int OUTPUT_BUFFER_SIZE;

    // E1 = (1, 0, 0), (0, 1, 0), (0, 0, 1);
    static const TlAngularMomentumVector E1_[3];

    // E2 = (2, 0, 0), (0, 2, 0), (0, 0, 2);
    static const TlAngularMomentumVector E2_[3];

    static const double PI3_2;
    static const double INV_SQRT3;

    int a_bar_;
    int b_bar_;
    int a_;
    int b_;
    TlPosition A_;
    TlPosition B_;

    TlPosition C_;
    double chargeC_;

    std::vector<double> zetaA_;
    std::vector<double> zetaB_;
    std::vector<double> zeta_;
    std::vector<double> coefAB_;
    std::vector<TlPosition> P_;

    OvpDataType OVP_;
    KinDataType KIN_;
    NucDataType NUC_;
    std::bitset<HPQ_OVP_STATE_MAX> isCalcdOvp_;
    std::bitset<HPQ_KIN_STATE_MAX> isCalcdKin_;
    std::bitset<HPQ_NUC_STATE_MAX> isCalcdNuc_;

    int numOfBatches_;
};

template <typename StateType, typename DataType>
void DfHpqEngine::addResultsToOutputBuffer(const DataType& fromBuf,
                                           double* toBuf) {
    const int a_bar = this->a_bar_;
    const int b_bar = this->b_bar_;
    const int a = this->a_;
    const int b = this->b_;

    const StateType state(a_bar, b_bar, a, b);
    const int stateIndex = state.index();

    const TlAngularMomentumVectorSet amvsAbar(a_bar);
    const TlAngularMomentumVectorSet amvsBbar(b_bar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int batch = this->numOfBatches_;
    int index = 0;
    for (int index_a_bar = 0; index_a_bar < amvsAbar.size(); ++index_a_bar) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(index_a_bar);

        for (int index_b_bar = 0; index_b_bar < amvsBbar.size();
             ++index_b_bar) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(index_b_bar);

            for (int index_a = 0; index_a < amvsA.size(); ++index_a) {
                const TlAngularMomentumVector amvA = amvsA.get(index_a);

                for (int index_b = 0; index_b < amvsB.size(); ++index_b) {
                    const TlAngularMomentumVector amvB = amvsB.get(index_b);

                    const int stateAmvsIndex =
                        this->index(amvAbar, amvBbar, amvA, amvB);

                    double value = 0.0;
                    for (int i = 0; i < batch; ++i) {
                        value += fromBuf[stateIndex][stateAmvsIndex][i];
                    }

                    toBuf[index] += value;
                    ++index;
                }  // b
            }      // a
        }
    }
}

#endif  // DFHPQENGINE_H
