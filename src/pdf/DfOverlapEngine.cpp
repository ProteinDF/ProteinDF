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

#include <cmath>
#include <limits>
#include "DfOverlapEngine.h"
#include "TlAngularMomentumVectorSet.h"

const double DfOverlapEngine::MAX_EXPONENT = std::log(std::numeric_limits<double>::max());

/// 結果出力用のサイズ
/// d軌道の6Dから5Dへの変換領域にも利用するため、
/// d軌道のサイズを6Dベースとして確保している。
/// = 3 * 3 * 3 * 3 * 6(d) * 6 * 6 * 6
const int DfOverlapEngine::OUTPUT_BUFFER_SIZE = 3 * 3 * 3 * 3 * 6 * 6 * 6 * 6;

const TlAngularMomentumVector DfOverlapEngine::E1_[3] = {
    TlAngularMomentumVector(1, 0, 0),
    TlAngularMomentumVector(0, 1, 0),
    TlAngularMomentumVector(0, 0, 1)
};

const TlAngularMomentumVector DfOverlapEngine::E2_[3] = {
    TlAngularMomentumVector(2, 0, 0),
    TlAngularMomentumVector(0, 2, 0),
    TlAngularMomentumVector(0, 0, 2)
};

const double DfOverlapEngine::INV_SQRT3 = 1.0 / sqrt(3.0);

DfOverlapEngine::DfOverlapEngine()
{
    this->OVP_.resize(OVP_STATE_MAX);

    this->WORK = new double[OUTPUT_BUFFER_SIZE];
    this->pTransformBuf_ = new double[OUTPUT_BUFFER_SIZE];
}


DfOverlapEngine::~DfOverlapEngine()
{
    delete[] this->WORK;
    this->WORK = NULL;

    delete[] this->pTransformBuf_;
    this->pTransformBuf_ = NULL;
}


void DfOverlapEngine::calc(const int diff1, const TlOrbitalInfoObject& orbInfo1, const index_type shell1,
                           const int diff2, const TlOrbitalInfoObject& orbInfo2, const index_type shell2,
                           const int diff3, const TlOrbitalInfoObject& orbInfo3, const index_type shell3,
                           const int diff4, const TlOrbitalInfoObject& orbInfo4, const index_type shell4)
{
    static const PGTO pgto0(1.0, 0.0);

    int shellType1 = 0;
    TlPosition pos1;
    PGTOs pgtos1;
    if (shell1 >= 0) {
        shellType1 = orbInfo1.getShellType(shell1);
        pos1 = orbInfo1.getPosition(shell1);
        pgtos1 = this->getPGTOs(orbInfo1, shell1);
    } else {
        pgtos1.resize(1);
        pgtos1[0] = pgto0;
    }

    int shellType2 = 0;
    TlPosition pos2;
    PGTOs pgtos2;
    if (shell2 >= 0) {
        shellType2 = orbInfo2.getShellType(shell2);
        pos2 = orbInfo2.getPosition(shell2);
        pgtos2 = this->getPGTOs(orbInfo2, shell2);
    } else {
        pgtos2.resize(1);
        pgtos2[0] = pgto0;
    }

    int shellType3 = 0;
    TlPosition pos3;
    PGTOs pgtos3;
    if (shell3 >= 0) {
        shellType3 = orbInfo3.getShellType(shell3);
        pos3 = orbInfo3.getPosition(shell3);
        pgtos3 = this->getPGTOs(orbInfo3, shell3);
    } else {
        pgtos3.resize(1);
        pgtos3[0] = pgto0;
    }

    int shellType4 = 0;
    TlPosition pos4;
    PGTOs pgtos4;
    if (shell4 >= 0) {
        shellType4 = orbInfo4.getShellType(shell4);
        pos4 = orbInfo4.getPosition(shell4);
        pgtos4 = this->getPGTOs(orbInfo4, shell4);
    } else {
        pgtos4.resize(1);
        pgtos4[0] = pgto0;
    }

    const Query query(diff1, diff2, diff3, diff4,
                      shellType1, shellType2, shellType3, shellType4);
    this->calc0(query, pos1, pos2, pos3, pos4,
                pgtos1, pgtos2, pgtos3, pgtos4);
}

void DfOverlapEngine::calc0(const Query& query,
                            const TlPosition& A,
                            const TlPosition& B,
                            const TlPosition& C,
                            const TlPosition& D,
                            const PGTOs& PGTOsA,
                            const PGTOs& PGTOsB,
                            const PGTOs& PGTOsC,
                            const PGTOs& PGTOsD)
{
    this->a_bar_ = query.shellTypeA_bar;
    this->b_bar_ = query.shellTypeB_bar;
    this->c_bar_ = query.shellTypeC_bar;
    this->d_bar_ =  query.shellTypeD_bar;
    this->a_ = query.shellTypeA;
    this->b_ = query.shellTypeB;
    this->c_ = query.shellTypeC;
    this->d_ = query.shellTypeD;
    
    this->A_ = A;
    this->B_ = B;
    this->C_ = C;
    this->D_ = D;
    const double A2 = A.squareDistanceFrom();
    const double B2 = B.squareDistanceFrom();
    const double C2 = C.squareDistanceFrom();
    const double D2 = D.squareDistanceFrom();

    const OvpState state(0, 0, 0, 0, 0, 0, 0, 0);
    const int stateIndex = state.index();
    
    //const TlAngularMomentumVector S(0, 0, 0);
    const int index = 0; // = this->index(S, S, S);
    this->OVP_[stateIndex].resize(1);
   
    const int numOfPGTOsA = PGTOsA.size();
    const int numOfPGTOsB = PGTOsB.size();
    const int numOfPGTOsC = PGTOsC.size();
    const int numOfPGTOsD = PGTOsD.size();
    const int numOfBatches = numOfPGTOsA * numOfPGTOsB * numOfPGTOsC * numOfPGTOsD;
    this->OVP_[stateIndex][index].resize(numOfBatches);
    this->G_.resize(numOfBatches);
    this->coef_.resize(numOfBatches);
    this->zetaA_.resize(numOfBatches);
    this->zetaB_.resize(numOfBatches);
    this->zetaC_.resize(numOfBatches);
    this->zetaD_.resize(numOfBatches);
    this->numOfBatches_ = numOfBatches;

    int batch = 0;
    for (int indexA = 0; indexA < numOfPGTOsA; ++indexA) {
        const double coefA = PGTOsA[indexA].coef;
        const double zetaA = PGTOsA[indexA].zeta;

        for (int indexB = 0; indexB < numOfPGTOsB; ++indexB) {
            const double coefB = PGTOsB[indexB].coef;
            const double zetaB = PGTOsB[indexB].zeta;
            const double zetaAB = zetaA + zetaB;
            
            for (int indexC = 0; indexC < numOfPGTOsC; ++indexC) {
                const double coefC = PGTOsC[indexC].coef;
                const double zetaC = PGTOsC[indexC].zeta;
                const double zetaABC = zetaAB + zetaC;

                for (int indexD = 0; indexD < numOfPGTOsD; ++indexD) {
                    const double coefD = PGTOsD[indexD].coef;
                    const double zetaD = PGTOsD[indexD].zeta;
                    
                    const double zeta = zetaABC + zetaD;
                    const double invZeta = 1.0 / zeta;
                    const TlPosition G = (zetaA * A + zetaB * B + zetaC * C + zetaD * D) * invZeta;
                    const double G2 = G.squareDistanceFrom();
                    const double kappa_exp = zeta * (G2 - (zetaA * A2 + zetaB * B2 + zetaC * C2 + zetaD * D2) * invZeta);
                    if (kappa_exp < DfOverlapEngine::MAX_EXPONENT) {
                        const double kappa = std::exp(kappa_exp);
                        const double base = M_PI * invZeta;
                        const double sss = coefA * coefB * coefC * coefD * base * std::sqrt(base) * kappa;
                        
                        // copy for batch operation
                        this->OVP_[stateIndex][index][batch] = sss;
                        this->G_[batch] = G;
                        this->coef_[batch] = 0.5 * invZeta;
                        this->zetaA_[batch] = zetaA;
                        this->zetaB_[batch] = zetaB;
                        this->zetaC_[batch] = zetaC;
                        this->zetaD_[batch] = zetaD;
                        ++batch;
                    }
                }
            }
        }
    }

    this->isCalcd_.reset();
    this->isCalcd_[stateIndex] = true;

    this->calc(this->a_bar_, this->b_bar_, this->c_bar_, this->d_bar_,
               this->a_, this->b_, this->c_, this->d_);

    this->copyResultsToOutputBuffer();
    this->transform6Dto5D();
}


DfOverlapEngine::PGTOs DfOverlapEngine::getPGTOs(const TlOrbitalInfoObject& orbitalInfo,
                                                 const int shellIndex)
{
    const int numOfContractions = orbitalInfo.getCgtoContraction(shellIndex);
    PGTOs pgtos(numOfContractions);

    for (int i = 0; i < numOfContractions; ++i) {
        const PGTO pgto(orbitalInfo.getCoefficient(shellIndex, i),
                        orbitalInfo.getExponent(shellIndex, i));
        pgtos[i] = pgto;
    }
    
    return pgtos;
}


void DfOverlapEngine::calc(const int a_bar, const int b_bar, const int c_bar, const int d_bar,
                           const int a, const int b, const int c, const int d)
{
    // std::cerr << TlUtils::format("DfOverlapEngine::calc() (%d %d %d, %d, %d, %d)",
    //                              a_bar, b_bar, c_bar,
    //                              a, b, c)
    //           << std::endl;
    if ((a_bar < 0) || (b_bar < 0) || (c_bar < 0) || (d_bar < 0) ||
        (a < 0) || (b < 0) || (c < 0) || (d < 0)) {
        return;
    }
    
    const OvpState state(a_bar, b_bar, c_bar, d_bar, a, b, c, d);
    const int stateIndex = state.index();
    if (this->isCalcd_[stateIndex] == true) {
        return;
    }

    if (a_bar > 0) {
        this->calc(a_bar-1, b_bar, c_bar, d_bar, a+1, b, c, d);
        this->calc(a_bar-1, b_bar, c_bar, d_bar, a-1, b, c, d);
        this->gradA(state);
    } else if (b_bar > 0) {
        this->calc(a_bar, b_bar-1, c_bar, d_bar, a, b+1, c, d);
        this->calc(a_bar, b_bar-1, c_bar, d_bar, a, b-1, c, d);
        this->gradB(state);
    } else if (c_bar > 0) {
        this->calc(a_bar, b_bar, c_bar-1, d_bar, a, b, c+1, d);
        this->calc(a_bar, b_bar, c_bar-1, d_bar, a, b, c-1, d);
        this->gradC(state);
    } else if (d_bar > 0) {
        this->calc(a_bar, b_bar, c_bar, d_bar-1, a, b, c, d+1);
        this->calc(a_bar, b_bar, c_bar, d_bar-1, a, b, c, d-1);
        this->gradD(state);
    } else if (a > 0) {
        this->calc(a_bar, b_bar, c_bar, d_bar, a-1, b,   c,   d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a-2, b,   c,   d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a-1, b-1, c,   d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a-1, b,   c-1, d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a-1, b,   c,   d-1);
        this->EQ20A(state);
    } else if (b > 0) {
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b-1, c,   d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a-1, b-1, c,   d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b-2, c,   d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b-1, c-1, d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b-1, c,   d-1);
        this->EQ20B(state);
    } else if (c > 0) {
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b,   c-1, d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a-1, b,   c-1, d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b-1, c-1, d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b,   c-2, d  );
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b,   c-1, d-1);
        this->EQ20C(state);
    } else if (d > 0) {
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b,   c,   d-1);
        this->calc(a_bar, b_bar, c_bar, d_bar, a-1, b,   c,   d-1);
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b-1, c,   d-1);
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b,   c-1, d-1);
        this->calc(a_bar, b_bar, c_bar, d_bar, a,   b,   c,   d-2);
        this->EQ20D(state);
    } else {
        std::cerr << TlUtils::format("Program Error: @DfOverlapEngine::calc() a=%d, b=%d, c=%d",
                                     a, b, c)
                  << std::endl;
        abort();
    }

    this->isCalcd_[stateIndex] = true;
}

int DfOverlapEngine::index(const TlAngularMomentumVector& amvAbar,
                           const TlAngularMomentumVector& amvBbar,
                           const TlAngularMomentumVector& amvCbar,
                           const TlAngularMomentumVector& amvDbar,
                           const TlAngularMomentumVector& amvA,
                           const TlAngularMomentumVector& amvB,
                           const TlAngularMomentumVector& amvC,
                           const TlAngularMomentumVector& amvD) const {
    //const int amAbar = amvAbar.angularMomentum();
    const int amBbar = amvBbar.angularMomentum();
    const int amCbar = amvCbar.angularMomentum();
    const int amDbar = amvDbar.angularMomentum();
    const int amA = amvA.angularMomentum();
    const int amB = amvB.angularMomentum();
    const int amC = amvC.angularMomentum();
    const int amD = amvD.angularMomentum();

    //const int aBars = TlAngularMomentumVectorSet(amAbar).size();
    const int bBars = TlAngularMomentumVectorSet(amBbar).size();
    const int cBars = TlAngularMomentumVectorSet(amCbar).size();
    const int dBars = TlAngularMomentumVectorSet(amDbar).size();
    const int as = TlAngularMomentumVectorSet(amA).size();
    const int bs = TlAngularMomentumVectorSet(amB).size();
    const int cs = TlAngularMomentumVectorSet(amC).size();
    const int ds = TlAngularMomentumVectorSet(amD).size();

    const int index = ((((((amvAbar.index() *bBars
                           +amvBbar.index())*cBars
                           +amvCbar.index())*dBars
                           +amvDbar.index())*as
                           +amvA.index())*bs
                           +amvB.index())*cs
                           +amvC.index())*ds
                           +amvD.index();
    return index;
}


void DfOverlapEngine::copyResultsToOutputBuffer()
{
    const int a_bar = this->a_bar_;
    const int b_bar = this->b_bar_;
    const int c_bar = this->c_bar_;
    const int d_bar = this->d_bar_;
    const int a = this->a_;
    const int b = this->b_;
    const int c = this->c_;
    const int d = this->d_;

    const OvpState state(a_bar, b_bar, c_bar, d_bar,
                         a, b, c, d);
    const int stateIndex = state.index();
    assert(0 <= stateIndex);
    assert(stateIndex < OVP_STATE_MAX);

    const TlAngularMomentumVectorSet amvsAbar(a_bar);
    const TlAngularMomentumVectorSet amvsBbar(b_bar);
    const TlAngularMomentumVectorSet amvsCbar(c_bar);
    const TlAngularMomentumVectorSet amvsDbar(d_bar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int batch = this->numOfBatches_;
    int index = 0;
    for (int index_a_bar = 0; index_a_bar < amvsAbar.size(); ++index_a_bar) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(index_a_bar);
        
        for (int index_b_bar = 0; index_b_bar < amvsBbar.size(); ++index_b_bar) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(index_b_bar);
            
            for (int index_c_bar = 0; index_c_bar < amvsCbar.size(); ++index_c_bar) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(index_c_bar);
                
                for (int index_d_bar = 0; index_d_bar < amvsDbar.size(); ++index_d_bar) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(index_d_bar);

                    for (int index_a = 0; index_a < amvsA.size(); ++index_a) {
                        const TlAngularMomentumVector amvA = amvsA.get(index_a);
                        
                        for (int index_b = 0; index_b < amvsB.size(); ++index_b) {
                            const TlAngularMomentumVector amvB = amvsB.get(index_b);
                            
                            for (int index_c = 0; index_c < amvsC.size(); ++index_c) {
                                const TlAngularMomentumVector amvC = amvsC.get(index_c);

                                for (int index_d = 0; index_d < amvsD.size(); ++index_d) {
                                    const TlAngularMomentumVector amvD = amvsD.get(index_d);

                                    const int stateAmvsIndex = this->index(amvAbar, amvBbar, amvCbar, amvDbar,
                                                                           amvA, amvB, amvC, amvD);
                                    assert(0 <= stateAmvsIndex);
                                    assert(stateAmvsIndex < static_cast<int>(this->OVP_[stateIndex].size()));
                                    assert(static_cast<int>(this->OVP_[stateIndex][stateAmvsIndex].size()) <= batch);

                                    double value = 0.0;
                                    for (int i = 0; i < batch; ++i) {
                                        value += this->OVP_[stateIndex][stateAmvsIndex][i];
                                    }
                                    assert(index < OUTPUT_BUFFER_SIZE);
                                    this->WORK[index] = value;
                                    ++index;
                                } // d
                            } // c
                        } // b
                    } // a
                }
            }
        }
    }
}    


void DfOverlapEngine::transform6Dto5D()
{
    const int a_bar = this->a_bar_;
    const int b_bar = this->b_bar_;
    const int c_bar = this->c_bar_;
    const int d_bar = this->d_bar_;
    const int a = this->a_;
    const int b = this->b_;
    const int c = this->c_;
    const int d = this->d_;
    assert(0 <= a_bar);
    assert(0 <= b_bar);
    assert(0 <= c_bar);
    assert(0 <= d_bar);
    assert(0 <= a);
    assert(0 <= b);
    assert(0 <= c);
    assert(0 <= d);
    
    if ((a < 2) && (b < 2) && (c < 2) && (d < 2)) {
        // nothing to do
        return;
    }

    // 6D ベースの要素数(1, 3, 6, 10, ...)
    // 5D ベースの場合は(2n +1: 1, 3, 5, 7, ...)
    int I_ = a_bar * (a_bar + 3) / 2 + 1;
    int J_ = b_bar * (b_bar + 3) / 2 + 1;
    int K_ = c_bar * (c_bar + 3) / 2 + 1;
    int L_ = d_bar * (d_bar + 3) / 2 + 1;
    int I = a * (a + 3) / 2 + 1;
    int J = b * (b + 3) / 2 + 1;
    int K = c * (c + 3) / 2 + 1;
    int L = d * (d + 3) / 2 + 1;

    assert(this->pTransformBuf_ != NULL);
    if (a == 2) {
        this->transform6Dto5D_i(I_, J_, K_, L_, J, K, L, this->WORK, this->pTransformBuf_);
        I = 5;
        const std::size_t end = I_ * J_ * K_ * L_ * I * J * K * L;
        assert(end <= static_cast<std::size_t>(OUTPUT_BUFFER_SIZE));
        std::copy(this->pTransformBuf_,
                  this->pTransformBuf_ + end,
                  this->WORK);
    }
    if (b == 2) {
        this->transform6Dto5D_j(I_, J_, K_, L_, I, K, L, this->WORK, this->pTransformBuf_);
        J = 5;
        const std::size_t end = I_ * J_ * K_ * L_ * I * J * K * L;
        assert(end <= static_cast<std::size_t>(OUTPUT_BUFFER_SIZE));
        std::copy(this->pTransformBuf_,
                  this->pTransformBuf_ + end,
                  this->WORK);
    }
    if (c == 2) {
        this->transform6Dto5D_k(I_, J_, K_, L_, I, J, L, this->WORK, this->pTransformBuf_);
        K = 5;
        const std::size_t end = I_ * J_ * K_ * L_ * I * J * K * L;
        assert(end <= static_cast<std::size_t>(OUTPUT_BUFFER_SIZE));
        std::copy(this->pTransformBuf_,
                  this->pTransformBuf_ + end,
                  this->WORK);
    }
    if (d == 2) {
        this->transform6Dto5D_l(I_, J_, K_, L_, I, J, K, this->WORK, this->pTransformBuf_);
        L = 5;
        const std::size_t end = I_ * J_ * K_ * L_ * I * J * K * L;
        assert(end <= static_cast<std::size_t>(OUTPUT_BUFFER_SIZE));
        std::copy(this->pTransformBuf_,
                  this->pTransformBuf_ + end,
                  this->WORK);
    }
}


void DfOverlapEngine::transform6Dto5D_i(const int I_, const int J_, const int K_, const int L_,
                                        const int J, const int K, const int L,
                                        const double* pInput, double* pOutput)
{
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int k_ = 0; k_ < K_; ++k_) {
                for (int l_ = 0; l_ < L_; ++l_) {
                
                    for (int j = 0; j < J; ++j) {
                        for (int k = 0; k < K; ++k) {
                            for (int l = 0; l < L; ++l) {
                                const double xx = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*6 +0)*J +j)*K +k)*L + l];
                                const double xy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ *l_)*6 +1)*J +j)*K +k)*L + l];
                                const double xz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ *l_)*6 +2)*J +j)*K +k)*L + l];
                                const double yy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ *l_)*6 +3)*J +j)*K +k)*L + l];
                                const double yz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ *l_)*6 +4)*J +j)*K +k)*L + l];
                                const double zz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ *l_)*6 +5)*J +j)*K +k)*L + l];
                                
                                const int xy_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*5 +0)*J +j)*K +k)*L +l;
                                const int xz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*5 +1)*J +j)*K +k)*L +l;
                                const int yz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*5 +2)*J +j)*K +k)*L +l;
                                const int xxyy_5d = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*5 +3)*J +j)*K +k)*L +l;
                                const int rr_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*5 +4)*J +j)*K +k)*L +l;
                                pOutput[xy_5d] = xy;
                                pOutput[xz_5d] = xz;
                                pOutput[yz_5d] = yz;
                                pOutput[xxyy_5d] = 0.5*(xx-yy);
                                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfOverlapEngine::transform6Dto5D_j(const int I_, const int J_, const int K_, const int L_,
                                        const int I, const int K, const int L,
                                        const double* pInput, double* pOutput)
{
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int k_ = 0; k_ < K_; ++k_) {
                for (int l_ = 0; l_ < L_; ++l_) {
                    
                    for (int i = 0; i < I; ++i) {
                        for (int k = 0; k < K; ++k) {
                            for (int l = 0; l < L; ++l) {
                                
                                const double xx = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*6 +0)*K +k)*L +l];
                                const double xy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*6 +1)*K +k)*L +l];
                                const double xz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*6 +2)*K +k)*L +l];
                                const double yy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*6 +3)*K +k)*L +l];
                                const double yz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*6 +4)*K +k)*L +l];
                                const double zz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*6 +5)*K +k)*L +l];
                                
                                const int xy_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*5 +0)*K +k)*L +l;
                                const int xz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*5 +1)*K +k)*L +l;
                                const int yz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*5 +2)*K +k)*L +l;
                                const int xxyy_5d = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*5 +3)*K +k)*L +l;
                                const int rr_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*5 +4)*K +k)*L +l;
                                pOutput[xy_5d] = xy;
                                pOutput[xz_5d] = xz;
                                pOutput[yz_5d] = yz;
                                pOutput[xxyy_5d] = 0.5*(xx-yy);
                                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfOverlapEngine::transform6Dto5D_k(const int I_, const int J_, const int K_, const int L_,
                                        const int I, const int J, const int L,
                                        const double* pInput, double* pOutput)
{
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int k_ = 0; k_ < K_; ++k_) {
                for (int l_ = 0; l_ < L_; ++l_) {

                    for (int i = 0; i < I; ++i) {
                        for (int j = 0; j < J; ++j) {
                            for (int l = 0; l < L; ++l) {
                                const double xx = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*6 +0)*L +l];
                                const double xy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*6 +1)*L +l];
                                const double xz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*6 +2)*L +l];
                                const double yy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*6 +3)*L +l];
                                const double yz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*6 +4)*L +l];
                                const double zz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*6 +5)*L +l];
                                
                                const int xy_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*5 +0)*L +l;
                                const int xz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*5 +1)*L +l;
                                const int yz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*5 +2)*L +l;
                                const int xxyy_5d = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*5 +3)*L +l;
                                const int rr_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*5 +4)*L +l;
                                pOutput[xy_5d] = xy;
                                pOutput[xz_5d] = xz;
                                pOutput[yz_5d] = yz;
                                pOutput[xxyy_5d] = 0.5*(xx-yy);
                                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfOverlapEngine::transform6Dto5D_l(const int I_, const int J_, const int K_, const int L_,
                                        const int I, const int J, const int K,
                                        const double* pInput, double* pOutput)
{
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int k_ = 0; k_ < K_; ++k_) {
                for (int l_ = 0; l_ < L_; ++l_) {

                    for (int i = 0; i < I; ++i) {
                        for (int j = 0; j < J; ++j) {
                            for (int k = 0; k < K; ++k) {
                                const double xx = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*6 +0];
                                const double xy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*6 +1];
                                const double xz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*6 +2];
                                const double yy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*6 +3];
                                const double yz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*6 +4];
                                const double zz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*6 +5];
                                
                                const int xy_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*5 +0;
                                const int xz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*5 +1;
                                const int yz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*5 +2;
                                const int xxyy_5d = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*5 +3;
                                const int rr_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ +l_)*I +i)*J +j)*K +k)*5 +4;
                                pOutput[xy_5d] = xy;
                                pOutput[xz_5d] = xz;
                                pOutput[yz_5d] = yz;
                                pOutput[xxyy_5d] = 0.5*(xx-yy);
                                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
                            }
                        }
                    }
                }
            }
        }
    }
}


int DfOverlapEngine::initiative(const TlAngularMomentumVector& amv) const {
    //const char x = amv.get(0);
    const char y = amv.get(1);
    const char z = amv.get(2);
    
    if (z > 0) {
        return 2;
    } else if (y > 0) {
        return 1;
    } else {
        return 0;
    }
}


void DfOverlapEngine::EQ20A(const OvpState& state)
{
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int cBar = state.c_bar;
    const int dBar = state.d_bar;
    const int a = state.a;
    const int b = state.b;
    const int c = state.c;
    const int d = state.d;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsCbar(cBar);
    const TlAngularMomentumVectorSet amvsDbar(dBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsD.size();
    this->OVP_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsCbar * numOfAmvsDbar
                                  * numOfAmvsA * numOfAmvsB * numOfAmvsC * numOfAmvsD);
    
    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar; ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar; ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvCbar_index = 0; amvCbar_index < numOfAmvsCbar; ++amvCbar_index) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(amvCbar_index);

                for (int amvDbar_index = 0; amvDbar_index < numOfAmvsDbar; ++amvDbar_index) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(amvDbar_index);

                    for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                        const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                        
                        const int i = this->initiative(amvA);
                        const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
                        const TlAngularMomentumVector amvA2 = amvA - this->E2_[i];
                        assert(amvA1.isExist() == true);
                        
                        for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                            const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                            const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];
                            
                            for (int amvC_index = 0; amvC_index < numOfAmvsC; ++amvC_index) {
                                const TlAngularMomentumVector amvC = amvsC.get(amvC_index);
                                const TlAngularMomentumVector amvC1 = amvC - this->E1_[i];
                                
                                for (int amvD_index = 0; amvD_index < numOfAmvsD; ++amvD_index) {
                                    const TlAngularMomentumVector amvD = amvsD.get(amvD_index);
                                    const TlAngularMomentumVector amvD1 = amvD - this->E1_[i];
                                    
                                    // batch -------------------------------------------------------
                                    const OvpState state1(aBar, bBar, cBar, dBar, a-1, b, c, d);
                                    const std::size_t state1_index = state1.index();
                                    const int index1 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA1, amvB, amvC, amvD);
                                    
                                    std::size_t state2_index = 0;
                                    int index2 = 0;
                                    if (amvA2.isExist() == true) {
                                        const OvpState state2(aBar, bBar, cBar, dBar, a-2, b, c, d);
                                        state2_index = state2.index();
                                        index2 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA2, amvB, amvC, amvD);
                                    }
                                    
                                    std::size_t state3_index = 0;
                                    int index3 = 0;
                                    if (amvB1.isExist() == true) {
                                        const OvpState state3(aBar, bBar, cBar, dBar, a-1, b-1, c, d);
                                        state3_index = state3.index();
                                        index3 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA1, amvB1, amvC, amvD);
                                    }
                                    
                                    std::size_t state4_index = 0;
                                    int index4 = 0;
                                    if (amvC1.isExist() == true) {
                                        const OvpState state4(aBar, bBar, cBar, dBar, a-1, b, c-1, d);
                                        state4_index = state4.index();
                                        index4 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA1, amvB, amvC1, amvD);
                                    }
                                    
                                    std::size_t state5_index = 0;
                                    int index5 = 0;
                                    if (amvD1.isExist() == true) {
                                        const OvpState state5(aBar, bBar, cBar, dBar, a-1, b, c, d-1);
                                        state5_index = state5.index();
                                        index5 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA1, amvB, amvC, amvD1);
                                    }

                                    const int index = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD);
                                    this->OVP_[stateIndex][index].resize(numOfBatches);
                                    for (int batch = 0; batch < numOfBatches; ++batch) {
                                        double answer = 0.0;
                                        const double coef = this->coef_[batch];
                                        
                                        // 1st term
                                        {
                                            answer += (this->G_[batch][i] - this->A_[i]) * this->OVP_[state1_index][index1][batch];
                                        }
                                        
                                        // 2nd term
                                        if (amvA2.isExist() == true) {
                                            answer += coef * amvA1.get(i) * this->OVP_[state2_index][index2][batch];
                                        }
                                        
                                        // 3rd term
                                        if (amvB1.isExist() == true) {
                                            answer += coef * amvB.get(i) * this->OVP_[state3_index][index3][batch];
                                        }
                                        
                                        // 4th term
                                        if (amvC1.isExist() == true) {
                                            answer += coef * amvC.get(i) * this->OVP_[state4_index][index4][batch];
                                        }
                                        
                                        // 5th term
                                        if (amvD1.isExist() == true) {
                                            answer += coef * amvD.get(i) * this->OVP_[state5_index][index5][batch];
                                        }

                                        this->OVP_[stateIndex][index][batch] = answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void DfOverlapEngine::EQ20B(const OvpState& state)
{
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int cBar = state.c_bar;
    const int dBar = state.d_bar;
    const int a = state.a;
    const int b = state.b;
    const int c = state.c;
    const int d = state.d;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsCbar(cBar);
    const TlAngularMomentumVectorSet amvsDbar(dBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsD.size();
    this->OVP_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsCbar * numOfAmvsDbar
                                  * numOfAmvsA * numOfAmvsB * numOfAmvsC * numOfAmvsD);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar; ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar; ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvCbar_index = 0; amvCbar_index < numOfAmvsCbar; ++amvCbar_index) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(amvCbar_index);

                for (int amvDbar_index = 0; amvDbar_index < numOfAmvsDbar; ++amvDbar_index) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(amvDbar_index);

                    for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                        const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                        const int i = this->initiative(amvB);
                        const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];
                        const TlAngularMomentumVector amvB2 = amvB - this->E2_[i];
                        
                        for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                            const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                            const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
                            
                            for (int amvC_index = 0; amvC_index < numOfAmvsC; ++amvC_index) {
                                const TlAngularMomentumVector amvC = amvsC.get(amvC_index);
                                const TlAngularMomentumVector amvC1 = amvC - this->E1_[i];

                                for (int amvD_index = 0; amvD_index < numOfAmvsD; ++amvD_index) {
                                    const TlAngularMomentumVector amvD = amvsD.get(amvD_index);
                                    const TlAngularMomentumVector amvD1 = amvD - this->E1_[i];
                            
                                    // batch -------------------------------------------------------
                                    const OvpState state1(aBar, bBar, cBar, dBar, a, b-1, c, d);
                                    const std::size_t state1_index = state1.index();
                                    const int index1 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB1, amvC, amvD);
                            
                                    std::size_t state2_index = 0;
                                    int index2 = 0;
                                    if (amvA1.isExist() == true) {
                                        const OvpState state2(aBar, bBar, cBar, dBar, a-1, b-1, c, d);
                                        state2_index = state2.index();
                                        index2 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA1, amvB1, amvC, amvD);
                                    }
                            
                                    std::size_t state3_index = 0;
                                    int index3 = 0;
                                    if (amvB2.isExist() == true) {
                                        const OvpState state3(aBar, bBar, cBar, dBar, a, b-2, c, d);
                                        state3_index = state3.index();
                                        index3 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB2, amvC, amvD);
                                    }
                            
                                    std::size_t state4_index = 0;
                                    int index4 = 0;
                                    if (amvC1.isExist() == true) {
                                        const OvpState state4(aBar, bBar, cBar, dBar, a, b-1, c -1, d);
                                        state4_index = state4.index();
                                        index4 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB1, amvC1, amvD);
                                    }
                            
                                    std::size_t state5_index = 0;
                                    int index5 = 0;
                                    if (amvD1.isExist() == true) {
                                        const OvpState state5(aBar, bBar, cBar, dBar, a, b-1, c, d -1);
                                        state5_index = state5.index();
                                        index5 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB1, amvC, amvD1);
                                    }

                                    const int index = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD);
                                    this->OVP_[stateIndex][index].resize(numOfBatches);
                                    for (int batch = 0; batch < numOfBatches; ++batch) {
                                        double answer = 0.0;
                                        const double coef = this->coef_[batch];
                                        
                                        // 1st term
                                        {
                                            answer += (this->G_[batch][i] - this->B_[i]) * this->OVP_[state1_index][index1][batch];
                                        }
                                        
                                        // 2nd term
                                        if (amvA1.isExist() == true) {
                                            answer += coef * amvA.get(i) * this->OVP_[state2_index][index2][batch];
                                        }
                                
                                        // 3rd term
                                        if (amvB2.isExist() == true) {
                                            answer += coef * amvB1.get(i) * this->OVP_[state3_index][index3][batch];
                                        }
                                        
                                        // 4th term
                                        if (amvC1.isExist() == true) {
                                            answer += coef * amvC.get(i) * this->OVP_[state4_index][index4][batch];
                                        }
                                
                                        // 5th term
                                        if (amvD1.isExist() == true) {
                                            answer += coef * amvD.get(i) * this->OVP_[state5_index][index5][batch];
                                        }

                                        this->OVP_[stateIndex][index][batch] = answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

            
void DfOverlapEngine::EQ20C(const OvpState& state)
{
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int cBar = state.c_bar;
    const int dBar = state.d_bar;
    const int a = state.a;
    const int b = state.b;
    const int c = state.c;
    const int d = state.d;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsCbar(cBar);
    const TlAngularMomentumVectorSet amvsDbar(dBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsD.size();
    this->OVP_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsCbar * numOfAmvsDbar
                                  * numOfAmvsA * numOfAmvsB * numOfAmvsC * numOfAmvsD);
    
    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar; ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar; ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvCbar_index = 0; amvCbar_index < numOfAmvsCbar; ++amvCbar_index) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(amvCbar_index);

                for (int amvDbar_index = 0; amvDbar_index < numOfAmvsDbar; ++amvDbar_index) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(amvDbar_index);

                    for (int amvC_index = 0; amvC_index < numOfAmvsC; ++amvC_index) {
                        const TlAngularMomentumVector amvC = amvsC.get(amvC_index);
                        const int i = this->initiative(amvC);
                        const TlAngularMomentumVector amvC1 = amvC - this->E1_[i];
                        const TlAngularMomentumVector amvC2 = amvC - this->E2_[i];
                        
                        for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                            const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                            const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
                            
                            for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                                const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                                const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];

                                for (int amvD_index = 0; amvD_index < numOfAmvsD; ++amvD_index) {
                                    const TlAngularMomentumVector amvD = amvsD.get(amvD_index);
                                    const TlAngularMomentumVector amvD1 = amvD - this->E1_[i];
                            
                                    // batch -------------------------------------------------------
                                    const OvpState state1(aBar, bBar, cBar, dBar, a, b, c-1, d);
                                    std::size_t state1_index = state1.index();
                                    const int index1 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC1, amvD);
                            
                                    std::size_t state2_index = 0;
                                    int index2 = 0;
                                    if (amvA1.isExist() == true) {
                                        const OvpState state2(aBar, bBar, cBar, dBar, a-1, b, c-1, d);
                                        state2_index = state2.index();
                                        index2 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA1, amvB, amvC1, amvD);
                                    }
                            
                                    std::size_t state3_index = 0;
                                    int index3 = 0;
                                    if (amvB1.isExist() == true) {
                                        const OvpState state3(aBar, bBar, cBar, dBar, a, b-1, c-1, d);
                                        state3_index = state3.index();
                                        index3 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB1, amvC1, amvD);
                                    }
                            
                                    std::size_t state4_index = 0;
                                    int index4 = 0;
                                    if (amvC2.isExist() == true) {
                                        const OvpState state4(aBar, bBar, cBar, dBar, a, b, c -2, d);
                                        state4_index = state4.index();
                                        index4 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC2, amvD);
                                    }
                                    
                                    std::size_t state5_index = 0;
                                    int index5 = 0;
                                    if (amvD1.isExist() == true) {
                                        const OvpState state5(aBar, bBar, cBar, dBar, a, b, c-1, d -1);
                                        state5_index = state5.index();
                                        index5 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC1, amvD1);
                                    }
                            
                                    const int index = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD);
                                    this->OVP_[stateIndex][index].resize(numOfBatches);
                                    for (int batch = 0; batch < numOfBatches; ++batch) {
                                        double answer = 0.0;
                                        const double coef = this->coef_[batch];
                                        
                                        // 1st term
                                        {
                                            answer += (this->G_[batch][i] - this->C_[i]) * this->OVP_[state1_index][index1][batch];
                                        }
                                        
                                        // 2nd term
                                        if (amvA1.isExist() == true) {
                                            answer += coef * amvA.get(i) * this->OVP_[state2_index][index2][batch];
                                        }
                                        
                                        // 3rd term
                                        if (amvB1.isExist() == true) {
                                            answer += coef * amvB.get(i) * this->OVP_[state3_index][index3][batch];
                                        }
                                        
                                        // 4th term
                                        if (amvC2.isExist() == true) {
                                            answer += coef * amvC1.get(i) * this->OVP_[state4_index][index4][batch];
                                        }
                                
                                        // 5th term
                                        if (amvD1.isExist() == true) {
                                            answer += coef * amvD.get(i) * this->OVP_[state5_index][index5][batch];
                                        }
                                        
                                        this->OVP_[stateIndex][index][batch] = answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfOverlapEngine::EQ20D(const OvpState& state)
{
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int cBar = state.c_bar;
    const int dBar = state.d_bar;
    const int a = state.a;
    const int b = state.b;
    const int c = state.c;
    const int d = state.d;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsCbar(cBar);
    const TlAngularMomentumVectorSet amvsDbar(dBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsD.size();
    this->OVP_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsCbar * numOfAmvsDbar
                                  * numOfAmvsA * numOfAmvsB * numOfAmvsC * numOfAmvsD);
    
    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar; ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar; ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvCbar_index = 0; amvCbar_index < numOfAmvsCbar; ++amvCbar_index) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(amvCbar_index);

                for (int amvDbar_index = 0; amvDbar_index < numOfAmvsDbar; ++amvDbar_index) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(amvDbar_index);

                    for (int amvD_index = 0; amvD_index < numOfAmvsD; ++amvD_index) {
                        const TlAngularMomentumVector amvD = amvsD.get(amvD_index);
                        const int i = this->initiative(amvD);
                        const TlAngularMomentumVector amvD1 = amvD - this->E1_[i];
                        const TlAngularMomentumVector amvD2 = amvD - this->E2_[i];
                        
                        for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                            const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                            const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
                            
                            for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                                const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                                const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];

                                for (int amvC_index = 0; amvC_index < numOfAmvsC; ++amvC_index) {
                                    const TlAngularMomentumVector amvC = amvsC.get(amvC_index);
                                    const TlAngularMomentumVector amvC1 = amvC - this->E1_[i];
                            
                                    // batch -------------------------------------------------------
                                    const OvpState state1(aBar, bBar, cBar, dBar, a, b, c, d -1);
                                    std::size_t state1_index = state1.index();
                                    const int index1 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD1);
                            
                                    std::size_t state2_index = 0;
                                    int index2 = 0;
                                    if (amvA1.isExist() == true) {
                                        const OvpState state2(aBar, bBar, cBar, dBar, a-1, b, c, d-1);
                                        state2_index = state2.index();
                                        index2 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA1, amvB, amvC, amvD1);
                                    }
                            
                                    std::size_t state3_index = 0;
                                    int index3 = 0;
                                    if (amvB1.isExist() == true) {
                                        const OvpState state3(aBar, bBar, cBar, dBar, a, b-1, c, d-1);
                                        state3_index = state3.index();
                                        index3 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB1, amvC, amvD1);
                                    }
                            
                                    std::size_t state4_index = 0;
                                    int index4 = 0;
                                    if (amvC1.isExist() == true) {
                                        const OvpState state4(aBar, bBar, cBar, dBar, a, b, c-1, d-1);
                                        state4_index = state4.index();
                                        index4 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC1, amvD1);
                                    }
                                    
                                    std::size_t state5_index = 0;
                                    int index5 = 0;
                                    if (amvD2.isExist() == true) {
                                        const OvpState state5(aBar, bBar, cBar, dBar, a, b, c, d-2);
                                        state5_index = state5.index();
                                        index5 = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD2);
                                    }
                            
                                    const int index = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD);
                                    this->OVP_[stateIndex][index].resize(numOfBatches);
                                    for (int batch = 0; batch < numOfBatches; ++batch) {
                                        double answer = 0.0;
                                        const double coef = this->coef_[batch];
                                        
                                        // 1st term
                                        {
                                            answer += (this->G_[batch][i] - this->D_[i]) * this->OVP_[state1_index][index1][batch];
                                        }
                                        
                                        // 2nd term
                                        if (amvA1.isExist() == true) {
                                            answer += coef * amvA.get(i) * this->OVP_[state2_index][index2][batch];
                                        }
                                        
                                        // 3rd term
                                        if (amvB1.isExist() == true) {
                                            answer += coef * amvB.get(i) * this->OVP_[state3_index][index3][batch];
                                        }
                                        
                                        // 4th term
                                        if (amvC1.isExist() == true) {
                                            answer += coef * amvC.get(i) * this->OVP_[state4_index][index4][batch];
                                        }
                                
                                        // 5th term
                                        if (amvD2.isExist() == true) {
                                            answer += coef * amvD1.get(i) * this->OVP_[state5_index][index5][batch];
                                        }
                                        
                                        this->OVP_[stateIndex][index][batch] = answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfOverlapEngine::gradA(const OvpState& state)
{
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int cBar = state.c_bar;
    const int dBar = state.d_bar;
    const int a = state.a;
    const int b = state.b;
    const int c = state.c;
    const int d = state.d;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsCbar(cBar);
    const TlAngularMomentumVectorSet amvsDbar(dBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsC.size();
    this->OVP_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsCbar * numOfAmvsDbar
                                  * numOfAmvsA * numOfAmvsB * numOfAmvsC * numOfAmvsD);
    
    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar; ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);
        const int i = this->initiative(amvAbar);
        const TlAngularMomentumVector amvAbar1 = amvAbar - this->E1_[i];
        assert(amvAbar1.isExist() == true);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar; ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvCbar_index = 0; amvCbar_index < numOfAmvsCbar; ++amvCbar_index) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(amvCbar_index);

                for (int amvDbar_index = 0; amvDbar_index < numOfAmvsDbar; ++amvDbar_index) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(amvDbar_index);

                    for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                        const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                        const TlAngularMomentumVector amvA1  = amvA - this->E1_[i];
                        const TlAngularMomentumVector amvA1p = amvA + this->E1_[i];

                        for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                            const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                        
                            for (int amvC_index = 0; amvC_index < numOfAmvsC; ++amvC_index) {
                                const TlAngularMomentumVector amvC = amvsC.get(amvC_index);
                            
                                for (int amvD_index = 0; amvD_index < numOfAmvsD; ++amvD_index) {
                                    const TlAngularMomentumVector amvD = amvsD.get(amvD_index);

                                    // batch -------------------------------------------------------
                                    const OvpState state1(aBar-1, bBar, cBar, dBar, a+1, b, c, d);
                                    const std::size_t state1_index = state1.index();
                                    const int index1 = this->index(amvAbar1, amvBbar, amvCbar, amvDbar, amvA1p, amvB, amvC, amvD);
                                    
                                    std::size_t state2_index = 0;
                                    int index2 = 0;
                                    if (amvA1.isExist() == true) {
                                        const OvpState state2(aBar-1, bBar, cBar, dBar, a-1, b, c, d);
                                        state2_index = state2.index();
                                        index2 = this->index(amvAbar1, amvBbar, amvCbar, amvDbar, amvA1, amvB, amvC, amvD);
                                    }
                            
                                    const int index = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD);
                                    this->OVP_[stateIndex][index].resize(numOfBatches);
                                    for (int batch = 0; batch < numOfBatches; ++batch) {
                                        double answer = 0.0;
                                
                                        // 1st term
                                        answer += 2.0 * this->zetaA_[batch] * this->OVP_[state1_index][index1][batch];
                                
                                        // 2nd term
                                        if (amvA1.isExist() == true) {
                                            answer -= amvA.get(i) * this->OVP_[state2_index][index2][batch];
                                        }
                                
                                        this->OVP_[stateIndex][index][batch] = answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfOverlapEngine::gradB(const OvpState& state)
{
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int cBar = state.c_bar;
    const int dBar = state.d_bar;
    const int a = state.a;
    const int b = state.b;
    const int c = state.c;
    const int d = state.d;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsCbar(cBar);
    const TlAngularMomentumVectorSet amvsDbar(dBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsD.size();
    this->OVP_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsCbar * numOfAmvsDbar
                                  * numOfAmvsA * numOfAmvsB * numOfAmvsC * numOfAmvsD);
    
    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar; ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar; ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);
            const int i = this->initiative(amvBbar);
            const TlAngularMomentumVector amvBbar1 = amvBbar - this->E1_[i];
            assert(amvBbar1.isExist() == true);

            for (int amvCbar_index = 0; amvCbar_index < numOfAmvsCbar; ++amvCbar_index) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(amvCbar_index);

                for (int amvDbar_index = 0; amvDbar_index < numOfAmvsDbar; ++amvDbar_index) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(amvDbar_index);

                    for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                        const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                        
                        for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                            const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                            const TlAngularMomentumVector amvB1  = amvB - this->E1_[i];
                            const TlAngularMomentumVector amvB1p = amvB + this->E1_[i];
                            
                            for (int amvC_index = 0; amvC_index < numOfAmvsC; ++amvC_index) {
                                const TlAngularMomentumVector amvC = amvsC.get(amvC_index);
                            
                                for (int amvD_index = 0; amvD_index < numOfAmvsD; ++amvD_index) {
                                    const TlAngularMomentumVector amvD = amvsD.get(amvD_index);

                                    // batch -------------------------------------------------------
                                    const OvpState state1(aBar, bBar-1, cBar, dBar, a, b+1, c, d);
                                    const std::size_t state1_index = state1.index();
                                    const int index1 = this->index(amvAbar, amvBbar1, amvCbar, amvDbar, amvA, amvB1p, amvC, amvD);
                                    
                                    std::size_t state2_index = 0;
                                    int index2 = 0;
                                    if (amvB1.isExist() == true) {
                                        const OvpState state2(aBar, bBar-1, cBar, dBar, a, b-1, c, d);
                                        state2_index = state2.index();
                                        index2 = this->index(amvAbar, amvBbar1, amvCbar, amvDbar, amvA, amvB1, amvC, amvD);
                                    }
                            
                                    const int index = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD);
                                    this->OVP_[stateIndex][index].resize(numOfBatches);
                                    for (int batch = 0; batch < numOfBatches; ++batch) {
                                        double answer = 0.0;
                                        
                                        // 1st term
                                        answer += 2.0 * this->zetaB_[batch] * this->OVP_[state1_index][index1][batch];
                                        
                                        // 2nd term
                                        if (amvB1.isExist() == true) {
                                            answer -= amvB.get(i) * this->OVP_[state2_index][index2][batch];
                                        }
                                
                                        this->OVP_[stateIndex][index][batch] = answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfOverlapEngine::gradC(const OvpState& state)
{
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int cBar = state.c_bar;
    const int dBar = state.d_bar;
    const int a = state.a;
    const int b = state.b;
    const int c = state.c;
    const int d = state.d;
    
    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsCbar(cBar);
    const TlAngularMomentumVectorSet amvsDbar(dBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsD.size();
    this->OVP_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsCbar * numOfAmvsDbar
                                  * numOfAmvsA * numOfAmvsB * numOfAmvsC * numOfAmvsD);
    
    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar; ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar; ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvCbar_index = 0; amvCbar_index < numOfAmvsCbar; ++amvCbar_index) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(amvCbar_index);
                const int i = this->initiative(amvCbar);
                const TlAngularMomentumVector amvCbar1 = amvCbar - this->E1_[i];
                assert(amvCbar1.isExist() == true);

                for (int amvDbar_index = 0; amvDbar_index < numOfAmvsDbar; ++amvDbar_index) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(amvDbar_index);

                    for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                        const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                        
                        for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                            const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                        
                            for (int amvC_index = 0; amvC_index < numOfAmvsC; ++amvC_index) {
                                const TlAngularMomentumVector amvC = amvsC.get(amvC_index);
                                const TlAngularMomentumVector amvC1  = amvC - this->E1_[i];
                                const TlAngularMomentumVector amvC1p = amvC + this->E1_[i];
                                
                                for (int amvD_index = 0; amvD_index < numOfAmvsD; ++amvD_index) {
                                    const TlAngularMomentumVector amvD = amvsD.get(amvD_index);
                                    
                                    // batch -------------------------------------------------------
                                    const OvpState state1(aBar, bBar, cBar-1, dBar, a, b, c+1, d);
                                    const std::size_t state1_index = state1.index();
                                    const int index1 = this->index(amvAbar, amvBbar, amvCbar1, amvDbar, amvA, amvB, amvC1p, amvD);
                                    
                                    std::size_t state2_index = 0;
                                    int index2 = 0;
                                    if (amvC1.isExist() == true) {
                                        const OvpState state2(aBar, bBar, cBar-1, dBar, a, b, c-1, d);
                                        state2_index = state2.index();
                                        index2 = this->index(amvAbar, amvBbar, amvCbar1, amvDbar, amvA, amvB, amvC1, amvD);
                                    }
                                    
                                    const int index = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD);
                                    this->OVP_[stateIndex][index].resize(numOfBatches);
                                    for (int batch = 0; batch < numOfBatches; ++batch) {
                                        double answer = 0.0;
                                        
                                        // 1st term
                                        answer += 2.0 * this->zetaC_[batch] * this->OVP_[state1_index][index1][batch];
                                        
                                        // 2nd term
                                        if (amvC1.isExist() == true) {
                                            answer -= amvC.get(i) * this->OVP_[state2_index][index2][batch];
                                        }
                                        
                                        this->OVP_[stateIndex][index][batch] = answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfOverlapEngine::gradD(const OvpState& state)
{
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int cBar = state.c_bar;
    const int dBar = state.d_bar;
    const int a = state.a;
    const int b = state.b;
    const int c = state.c;
    const int d = state.d;
    
    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsCbar(cBar);
    const TlAngularMomentumVectorSet amvsDbar(dBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsD.size();
    this->OVP_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsCbar * numOfAmvsDbar
                                  * numOfAmvsA * numOfAmvsB * numOfAmvsC * numOfAmvsD);
    
    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar; ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar; ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvCbar_index = 0; amvCbar_index < numOfAmvsCbar; ++amvCbar_index) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(amvCbar_index);

                for (int amvDbar_index = 0; amvDbar_index < numOfAmvsDbar; ++amvDbar_index) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(amvDbar_index);
                    const int i = this->initiative(amvDbar);
                    const TlAngularMomentumVector amvDbar1 = amvDbar - this->E1_[i];
                    assert(amvDbar1.isExist() == true);

                    for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                        const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                        
                        for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                            const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                        
                            for (int amvC_index = 0; amvC_index < numOfAmvsC; ++amvC_index) {
                                const TlAngularMomentumVector amvC = amvsC.get(amvC_index);
                                    
                                for (int amvD_index = 0; amvD_index < numOfAmvsD; ++amvD_index) {
                                    const TlAngularMomentumVector amvD = amvsD.get(amvD_index);
                                    const TlAngularMomentumVector amvD1  = amvD - this->E1_[i];
                                    const TlAngularMomentumVector amvD1p = amvD + this->E1_[i];
                                
                                    // batch -------------------------------------------------------
                                    const OvpState state1(aBar, bBar, cBar, dBar-1, a, b, c, d+1);
                                    const std::size_t state1_index = state1.index();
                                    const int index1 = this->index(amvAbar, amvBbar, amvCbar, amvDbar1, amvA, amvB, amvC, amvD1p);
                                    
                                    std::size_t state2_index = 0;
                                    int index2 = 0;
                                    if (amvD1.isExist() == true) {
                                        const OvpState state2(aBar, bBar, cBar, dBar-1, a, b, c, d-1);
                                        state2_index = state2.index();
                                        index2 = this->index(amvAbar, amvBbar, amvCbar, amvDbar1, amvA, amvB, amvC, amvD1);
                                    }
                                    
                                    const int index = this->index(amvAbar, amvBbar, amvCbar, amvDbar, amvA, amvB, amvC, amvD);
                                    this->OVP_[stateIndex][index].resize(numOfBatches);
                                    for (int batch = 0; batch < numOfBatches; ++batch) {
                                        double answer = 0.0;
                                        
                                        // 1st term
                                        answer += 2.0 * this->zetaD_[batch] * this->OVP_[state1_index][index1][batch];
                                        
                                        // 2nd term
                                        if (amvD1.isExist() == true) {
                                            answer -= amvD.get(i) * this->OVP_[state2_index][index2][batch];
                                        }
                                        
                                        this->OVP_[stateIndex][index][batch] = answer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


