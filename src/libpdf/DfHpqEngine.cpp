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

#include "DfHpqEngine.h"
#include <cmath>
#include "TlFmt.h"

/// 結果出力用のサイズ
/// F軌道の10Fから7Fへの変換領域にも利用するため、
/// F軌道のサイズを10Fベースとして確保している。
const int DfHpqEngine::OUTPUT_BUFFER_SIZE = 10 * 10 * 3;

const TlAngularMomentumVector DfHpqEngine::E1_[3] = {
    TlAngularMomentumVector(1, 0, 0), TlAngularMomentumVector(0, 1, 0),
    TlAngularMomentumVector(0, 0, 1)};

const TlAngularMomentumVector DfHpqEngine::E2_[3] = {
    TlAngularMomentumVector(2, 0, 0), TlAngularMomentumVector(0, 2, 0),
    TlAngularMomentumVector(0, 0, 2)};

const double DfHpqEngine::PI3_2 = M_PI * std::sqrt(M_PI);
const double DfHpqEngine::INV_SQRT3 = 1.0 / sqrt(3.0);

DfHpqEngine::DfHpqEngine() {
    this->OVP_.resize(HPQ_OVP_STATE_MAX);
    this->KIN_.resize(HPQ_KIN_STATE_MAX);
    this->NUC_.resize(HPQ_NUC_STATE_MAX);

    this->WORK_KIN = new double[DfHpqEngine::OUTPUT_BUFFER_SIZE];
    this->WORK_NUC = new double[DfHpqEngine::OUTPUT_BUFFER_SIZE];
    this->WORK_NUCX = new double[DfHpqEngine::OUTPUT_BUFFER_SIZE];
}

DfHpqEngine::~DfHpqEngine() {
    delete[] this->WORK_KIN;
    this->WORK_KIN = NULL;

    delete[] this->WORK_NUC;
    this->WORK_NUC = NULL;

    delete[] this->WORK_NUCX;
    this->WORK_NUCX = NULL;
}

void DfHpqEngine::calc(const Query& query, const TlPosition& A,
                       const TlPosition& B, const PGTOs& PGTOsA,
                       const PGTOs& PGTOsB, const std::vector<TlAtom>& Cs,
                       const std::vector<TlAtom>& Xs) {
    this->a_bar_ = query.shellTypeA_bar;
    this->b_bar_ = query.shellTypeB_bar;
    this->a_ = query.shellTypeA;
    this->b_ = query.shellTypeB;

    this->A_ = A;
    this->B_ = B;

    this->prepare(PGTOsA, PGTOsB);

    // clean up
    {
        const int dimAbar = this->a_bar_ * (this->a_bar_ + 3) / 2 + 1;
        const int dimBbar = this->b_bar_ * (this->b_bar_ + 3) / 2 + 1;
        const int dimA = this->a_ * (this->a_ + 3) / 2 + 1;
        const int dimB = this->b_ * (this->b_ + 3) / 2 + 1;
        const int length = dimAbar * dimBbar * dimA * dimB;
        std::fill(this->WORK_KIN, this->WORK_KIN + length, 0.0);
        std::fill(this->WORK_NUC, this->WORK_NUC + length, 0.0);
        std::fill(this->WORK_NUCX, this->WORK_NUCX + length, 0.0);
    }

    this->isCalcdOvp_.reset();
    this->isCalcdKin_.reset();
    this->calcOvpSS();
    this->calcKinSS();

    this->calcKin(this->a_bar_, this->b_bar_, this->a_, this->b_);
    this->addResultsToOutputBuffer<KinState, KinDataType>(this->KIN_,
                                                          this->WORK_KIN);

    const int numOfRealAtoms = Cs.size();
    for (int i = 0; i < numOfRealAtoms; ++i) {
        this->C_ = Cs[i].getPosition();
        this->chargeC_ = Cs[i].getCharge();
        this->isCalcdNuc_.reset();
        this->calcNucSS();
        this->calcNuc(this->a_bar_, this->b_bar_, this->a_, this->b_, 0);
        this->addResultsToOutputBuffer<NucState, NucDataType>(this->NUC_,
                                                              this->WORK_NUC);
    }

    const int numOfDummyAtoms = Xs.size();
    for (int i = 0; i < numOfDummyAtoms; ++i) {
        this->C_ = Xs[i].getPosition();
        this->chargeC_ = Xs[i].getCharge();
        this->isCalcdNuc_.reset();
        this->calcNucSS();
        this->calcNuc(this->a_bar_, this->b_bar_, this->a_, this->b_, 0);
        this->addResultsToOutputBuffer<NucState, NucDataType>(this->NUC_,
                                                              this->WORK_NUCX);
    }

    this->transform6Dto5D(this->WORK_KIN);
    this->transform6Dto5D(this->WORK_NUC);
    this->transform6Dto5D(this->WORK_NUCX);
}

void DfHpqEngine::calc(const Query& query, const TlPosition& A,
                       const TlPosition& B, const PGTOs& PGTOsA,
                       const PGTOs& PGTOsB, const TlPosition& C) {
    this->a_bar_ = query.shellTypeA_bar;
    this->b_bar_ = query.shellTypeB_bar;
    this->a_ = query.shellTypeA;
    this->b_ = query.shellTypeB;

    this->A_ = A;
    this->B_ = B;

    this->prepare(PGTOsA, PGTOsB);

    // clean up
    {
        const int dimAbar = this->a_bar_ * (this->a_bar_ + 3) / 2 + 1;
        const int dimBbar = this->b_bar_ * (this->b_bar_ + 3) / 2 + 1;
        const int dimA = this->a_ * (this->a_ + 3) / 2 + 1;
        const int dimB = this->b_ * (this->b_ + 3) / 2 + 1;
        const int length = dimAbar * dimBbar * dimA * dimB;
        std::fill(this->WORK_NUC, this->WORK_NUC + length, 0.0);
    }

    this->isCalcdOvp_.reset();
    this->calcOvpSS();

    this->C_ = C;
    this->chargeC_ = 1.0;
    this->isCalcdNuc_.reset();
    this->calcNucSS();
    this->calcNuc(this->a_bar_, this->b_bar_, this->a_, this->b_, 0);
    this->addResultsToOutputBuffer<NucState, NucDataType>(this->NUC_,
                                                          this->WORK_NUC);
    this->transform6Dto5D(this->WORK_NUC);
}

void DfHpqEngine::calcKineticPart(const Query& query, const TlPosition& A,
                                  const TlPosition& B, const PGTOs& PGTOsA,
                                  const PGTOs& PGTOsB) {
    this->a_bar_ = query.shellTypeA_bar;
    this->b_bar_ = query.shellTypeB_bar;
    this->a_ = query.shellTypeA;
    this->b_ = query.shellTypeB;

    this->A_ = A;
    this->B_ = B;

    this->prepare(PGTOsA, PGTOsB);

    // clean up
    {
        const int dimAbar = this->a_bar_ * (this->a_bar_ + 3) / 2 + 1;
        const int dimBbar = this->b_bar_ * (this->b_bar_ + 3) / 2 + 1;
        const int dimA = this->a_ * (this->a_ + 3) / 2 + 1;
        const int dimB = this->b_ * (this->b_ + 3) / 2 + 1;
        const int length = dimAbar * dimBbar * dimA * dimB;
        std::fill(this->WORK_KIN, this->WORK_KIN + length, 0.0);
    }

    this->isCalcdOvp_.reset();
    this->isCalcdKin_.reset();
    this->calcOvpSS();
    this->calcKinSS();

    this->calcKin(this->a_bar_, this->b_bar_, this->a_, this->b_);
    this->addResultsToOutputBuffer<KinState, KinDataType>(this->KIN_,
                                                          this->WORK_KIN);

    this->transform6Dto5D(this->WORK_KIN);
}

void DfHpqEngine::calcNuclearAttractionPart(
    const Query& query, const TlPosition& A, const TlPosition& B,
    const PGTOs& PGTOsA, const PGTOs& PGTOsB, const std::vector<TlAtom>& Cs) {
    this->a_bar_ = query.shellTypeA_bar;
    this->b_bar_ = query.shellTypeB_bar;
    this->a_ = query.shellTypeA;
    this->b_ = query.shellTypeB;

    this->A_ = A;
    this->B_ = B;

    this->prepare(PGTOsA, PGTOsB);

    // clean up
    {
        const int dimAbar = this->a_bar_ * (this->a_bar_ + 3) / 2 + 1;
        const int dimBbar = this->b_bar_ * (this->b_bar_ + 3) / 2 + 1;
        const int dimA = this->a_ * (this->a_ + 3) / 2 + 1;
        const int dimB = this->b_ * (this->b_ + 3) / 2 + 1;
        const int length = dimAbar * dimBbar * dimA * dimB;
        std::fill(this->WORK_NUC, this->WORK_NUC + length, 0.0);
    }

    this->isCalcdOvp_.reset();
    this->isCalcdKin_.reset();
    this->calcOvpSS();
    this->calcKinSS();

    const int numOfRealAtoms = Cs.size();
    for (int i = 0; i < numOfRealAtoms; ++i) {
        this->C_ = Cs[i].getPosition();
        this->chargeC_ = Cs[i].getCharge();
        this->isCalcdNuc_.reset();
        this->calcNucSS();
        this->calcNuc(this->a_bar_, this->b_bar_, this->a_, this->b_, 0);
        this->addResultsToOutputBuffer<NucState, NucDataType>(this->NUC_,
                                                              this->WORK_NUC);
    }

    this->transform6Dto5D(this->WORK_NUC);
}

void DfHpqEngine::prepare(const PGTOs& PGTOsA, const PGTOs& PGTOsB) {
    const int numOfPGTOsA = PGTOsA.size();
    const int numOfPGTOsB = PGTOsB.size();
    const int numOfBatches = numOfPGTOsA * numOfPGTOsB;
    this->zetaA_.resize(numOfBatches);
    this->zetaB_.resize(numOfBatches);
    this->zeta_.resize(numOfBatches);
    this->coefAB_.resize(numOfBatches);
    this->P_.resize(numOfBatches);
    this->numOfBatches_ = numOfBatches;

    int batch = 0;
    for (int indexA = 0; indexA < numOfPGTOsA; ++indexA) {
        const double coefA = PGTOsA[indexA].coef;
        const double zetaA = PGTOsA[indexA].zeta;

        for (int indexB = 0; indexB < numOfPGTOsB; ++indexB) {
            const double coefB = PGTOsB[indexB].coef;
            const double zetaB = PGTOsB[indexB].zeta;

            // copy for batch operation
            this->zetaA_[batch] = zetaA;
            this->zetaB_[batch] = zetaB;
            const double zeta = zetaA + zetaB;
            this->zeta_[batch] = zeta;
            this->coefAB_[batch] = coefA * coefB;
            this->P_[batch] = (zetaA * this->A_ + zetaB * this->B_) / zeta;
            ++batch;
        }
    }
}

int DfHpqEngine::initiative(const TlAngularMomentumVector& amv) const {
    // const char x = amv.get(0);
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

int DfHpqEngine::index(const TlAngularMomentumVector& amvA,
                       const TlAngularMomentumVector& amvB) const {
    // const int amA = amvA.angularMomentum();
    const int amB = amvB.angularMomentum();

    // const int as = TlAngularMomentumVectorSet(amA).size();
    const int bs = TlAngularMomentumVectorSet(amB).size();

    const int index = amvA.index() * bs + amvB.index();
    return index;
}

int DfHpqEngine::index(const TlAngularMomentumVector& amvAbar,
                       const TlAngularMomentumVector& amvBbar,
                       const TlAngularMomentumVector& amvA,
                       const TlAngularMomentumVector& amvB) const {
    // const int amAbar = amvAbar.angularMomentum();
    const int amBbar = amvBbar.angularMomentum();
    const int amA = amvA.angularMomentum();
    const int amB = amvB.angularMomentum();

    // const int aBars = TlAngularMomentumVectorSet(amAbar).size();
    const int bBars = TlAngularMomentumVectorSet(amBbar).size();
    const int as = TlAngularMomentumVectorSet(amA).size();
    const int bs = TlAngularMomentumVectorSet(amB).size();

    const int index =
        ((amvAbar.index() * bBars + amvBbar.index()) * as + amvA.index()) * bs +
        amvB.index();
    return index;
}

void DfHpqEngine::transform6Dto5D(double* pBuf) {
    const int a_bar = this->a_bar_;
    const int b_bar = this->b_bar_;
    const int a = this->a_;
    const int b = this->b_;

    if ((a < 2) && (b < 2)) {
        // nothing to do
        return;
    }

    double pTempBuf[OUTPUT_BUFFER_SIZE];

    // 6D ベースの要素数(1, 3, 6, 10, ...)
    // 5D ベースの場合は(2n +1: 1, 3, 5, 7, ...)
    int I_ = a_bar * (a_bar + 3) / 2 + 1;
    int J_ = b_bar * (b_bar + 3) / 2 + 1;
    int I = a * (a + 3) / 2 + 1;
    int J = b * (b + 3) / 2 + 1;

    if (a == 2) {
        this->transform6Dto5D_i(I_, J_, J, pBuf, pTempBuf);
        I = 5;
        const std::size_t end = I_ * J_ * I * J;
        std::copy(pTempBuf, pTempBuf + end, pBuf);
    }
    if (b == 2) {
        this->transform6Dto5D_j(I_, J_, I, pBuf, pTempBuf);
        J = 5;
        const std::size_t end = I_ * J_ * I * J;
        std::copy(pTempBuf, pTempBuf + end, pBuf);
    }

    if (a == 3) {
        this->transform10Fto7F_i(I_, J_, J, pBuf, pTempBuf);
        I = 7;
        const std::size_t end = I_ * J_ * I * J;
        std::copy(pTempBuf, pTempBuf + end, pBuf);
    }
    if (b == 3) {
        this->transform10Fto7F_j(I_, J_, I, pBuf, pTempBuf);
        J = 7;
        const std::size_t end = I_ * J_ * I * J;
        std::copy(pTempBuf, pTempBuf + end, pBuf);
    }
}

void DfHpqEngine::transform6Dto5D_i(const int I_, const int J_, const int J,
                                    const double* pInput, double* pOutput) {
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int j = 0; j < J; ++j) {
                const double xx = pInput[((i_ * J_ + j_) * 6 + 0) * J + j];
                const double xy = pInput[((i_ * J_ + j_) * 6 + 1) * J + j];
                const double xz = pInput[((i_ * J_ + j_) * 6 + 2) * J + j];
                const double yy = pInput[((i_ * J_ + j_) * 6 + 3) * J + j];
                const double yz = pInput[((i_ * J_ + j_) * 6 + 4) * J + j];
                const double zz = pInput[((i_ * J_ + j_) * 6 + 5) * J + j];

                const int xy_5d = ((i_ * J_ + j_) * 5 + 0) * J + j;
                const int xz_5d = ((i_ * J_ + j_) * 5 + 1) * J + j;
                const int yz_5d = ((i_ * J_ + j_) * 5 + 2) * J + j;
                const int xxyy_5d = ((i_ * J_ + j_) * 5 + 3) * J + j;
                const int rr_5d = ((i_ * J_ + j_) * 5 + 4) * J + j;
                pOutput[xy_5d] = xy;
                pOutput[xz_5d] = xz;
                pOutput[yz_5d] = yz;
                pOutput[xxyy_5d] = 0.5 * (xx - yy);
                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
            }
        }
    }
}

void DfHpqEngine::transform6Dto5D_j(const int I_, const int J_, const int I,
                                    const double* pInput, double* pOutput) {
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int i = 0; i < I; ++i) {
                const double xx = pInput[((i_ * J_ + j_) * I + i) * 6 + 0];
                const double xy = pInput[((i_ * J_ + j_) * I + i) * 6 + 1];
                const double xz = pInput[((i_ * J_ + j_) * I + i) * 6 + 2];
                const double yy = pInput[((i_ * J_ + j_) * I + i) * 6 + 3];
                const double yz = pInput[((i_ * J_ + j_) * I + i) * 6 + 4];
                const double zz = pInput[((i_ * J_ + j_) * I + i) * 6 + 5];

                const int xy_5d = ((i_ * J_ + j_) * I + i) * 5 + 0;
                const int xz_5d = ((i_ * J_ + j_) * I + i) * 5 + 1;
                const int yz_5d = ((i_ * J_ + j_) * I + i) * 5 + 2;
                const int xxyy_5d = ((i_ * J_ + j_) * I + i) * 5 + 3;
                const int rr_5d = ((i_ * J_ + j_) * I + i) * 5 + 4;
                pOutput[xy_5d] = xy;
                pOutput[xz_5d] = xz;
                pOutput[yz_5d] = yz;
                pOutput[xxyy_5d] = 0.5 * (xx - yy);
                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
            }
        }
    }
}

void DfHpqEngine::transform10Fto7F_i(const int I_, const int J_, const int J,
                                     const double* pInput, double* pOutput) {
    const double route_23 = std::sqrt(2.0 / 3.0);
    const double route_25 = std::sqrt(2.0 / 5.0);
    const double route_1_15 = std::sqrt(1.0 / 15.0);

    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int j = 0; j < J; ++j) {
                const double xxx = pInput[((i_ * J_ + j_) * 10 + 0) * J + j];
                const double xxy = pInput[((i_ * J_ + j_) * 10 + 1) * J + j];
                const double xxz = pInput[((i_ * J_ + j_) * 10 + 2) * J + j];
                const double xyy = pInput[((i_ * J_ + j_) * 10 + 3) * J + j];
                const double xyz = pInput[((i_ * J_ + j_) * 10 + 4) * J + j];
                const double xzz = pInput[((i_ * J_ + j_) * 10 + 5) * J + j];
                const double yyy = pInput[((i_ * J_ + j_) * 10 + 6) * J + j];
                const double yyz = pInput[((i_ * J_ + j_) * 10 + 7) * J + j];
                const double yzz = pInput[((i_ * J_ + j_) * 10 + 8) * J + j];
                const double zzz = pInput[((i_ * J_ + j_) * 10 + 9) * J + j];

                const int z3_7f = ((i_ * J_ + j_) * 7 + 0) * J + j;
                const int xz2_7f = ((i_ * J_ + j_) * 7 + 1) * J + j;
                const int yz2_7f = ((i_ * J_ + j_) * 7 + 2) * J + j;
                const int x2y_y3_7f = ((i_ * J_ + j_) * 7 + 3) * J + j;
                const int x3_xy2_7f = ((i_ * J_ + j_) * 7 + 4) * J + j;
                const int xyz_7f = ((i_ * J_ + j_) * 7 + 5) * J + j;
                const int x2z_y2z_7f = ((i_ * J_ + j_) * 7 + 6) * J + j;

                pOutput[z3_7f] =
                    0.5 * route_1_15 * (5 * zzz - 3 * (xxz + yyz + zzz));
                pOutput[xz2_7f] =
                    0.25 * route_25 * (5 * xzz - (xxx + xyy + xzz));
                pOutput[yz2_7f] =
                    0.25 * route_25 * (5 * yzz - (xxy + yyy + yzz));
                pOutput[x2y_y3_7f] = 0.25 * route_23 * (3 * xxy - yyy);
                pOutput[x3_xy2_7f] = 0.25 * route_23 * (xxx - 3 * xyy);
                pOutput[xyz_7f] = xyz;
                pOutput[x2z_y2z_7f] = 0.5 * (xxz - yyz);
            }
        }
    }
}

void DfHpqEngine::transform10Fto7F_j(const int I_, const int J_, const int I,
                                     const double* pInput, double* pOutput) {
    const double route_23 = std::sqrt(2.0 / 3.0);
    const double route_25 = std::sqrt(2.0 / 5.0);
    const double route_1_15 = std::sqrt(1.0 / 15.0);

    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int i = 0; i < I; ++i) {
                const double xxx = pInput[((i_ * J_ + j_) * I + i) * 10 + 0];
                const double xxy = pInput[((i_ * J_ + j_) * I + i) * 10 + 1];
                const double xxz = pInput[((i_ * J_ + j_) * I + i) * 10 + 2];
                const double xyy = pInput[((i_ * J_ + j_) * I + i) * 10 + 3];
                const double xyz = pInput[((i_ * J_ + j_) * I + i) * 10 + 4];
                const double xzz = pInput[((i_ * J_ + j_) * I + i) * 10 + 5];
                const double yyy = pInput[((i_ * J_ + j_) * I + i) * 10 + 6];
                const double yyz = pInput[((i_ * J_ + j_) * I + i) * 10 + 7];
                const double yzz = pInput[((i_ * J_ + j_) * I + i) * 10 + 8];
                const double zzz = pInput[((i_ * J_ + j_) * I + i) * 10 + 9];

                const int z3_7f = ((i_ * J_ + j_) * I + i) * 7 + 0;
                const int xz2_7f = ((i_ * J_ + j_) * I + i) * 7 + 1;
                const int yz2_7f = ((i_ * J_ + j_) * I + i) * 7 + 2;
                const int x2y_y3_7f = ((i_ * J_ + j_) * I + i) * 7 + 3;
                const int x3_xy2_7f = ((i_ * J_ + j_) * I + i) * 7 + 4;
                const int xyz_7f = ((i_ * J_ + j_) * I + i) * 7 + 5;
                const int x2z_y2z_7f = ((i_ * J_ + j_) * I + i) * 7 + 6;

                pOutput[z3_7f] =
                    0.5 * route_1_15 * (5 * zzz - 3 * (xxz + yyz + zzz));
                pOutput[xz2_7f] =
                    0.25 * route_25 * (5 * xzz - (xxx + xyy + xzz));
                pOutput[yz2_7f] =
                    0.25 * route_25 * (5 * yzz - (xxy + yyy + yzz));
                pOutput[x2y_y3_7f] = 0.25 * route_23 * (3 * xxy - yyy);
                pOutput[x3_xy2_7f] = 0.25 * route_23 * (xxx - 3 * xyy);
                pOutput[xyz_7f] = xyz;
                pOutput[x2z_y2z_7f] = 0.5 * (xxz - yyz);
            }
        }
    }
}

void DfHpqEngine::calcOvp(const int a, const int b) {
    if ((a < 0) || (b < 0)) {
        return;
    }

    const OvpState state(a, b);
    const int stateIndex = state.index();
    if (this->isCalcdOvp_[stateIndex] == true) {
        return;
    }

    if (a > 0) {
        this->calcOvp(a - 1, b);
        this->calcOvp(a - 2, b);
        this->calcOvp(a - 1, b - 1);
        this->Ovp_EqA2_A(state);
    } else if (b > 0) {
        this->calcOvp(a, b - 1);
        this->calcOvp(a - 1, b - 1);
        this->calcOvp(a, b - 2);
        this->Ovp_EqA2_B(state);
    } else {
        std::cerr << TlUtils::format(
                         "Program Error: @DfHpqEngine::calcOvp() a=%d, b=%d", a,
                         b)
                  << std::endl;
        abort();
    }

    this->isCalcdOvp_[stateIndex] = true;
}

void DfHpqEngine::calcKin(const int a_bar, const int b_bar, const int a,
                          const int b) {
    if ((a_bar < 0) || (b_bar < 0) || (a < 0) || (b < 0)) {
        return;
    }

    const KinState state(a_bar, b_bar, a, b);
    const int stateIndex = state.index();
    if (this->isCalcdKin_[stateIndex] == true) {
        return;
    }

    if (a_bar > 0) {
        this->calcKin(a_bar - 1, b_bar, a + 1, b);
        this->calcKin(a_bar - 1, b_bar, a - 1, b);
        this->gradKinA(state);
    } else if (b_bar > 0) {
        this->calcKin(a_bar, b_bar - 1, a, b + 1);
        this->calcKin(a_bar, b_bar - 1, a, b - 1);
        this->gradKinB(state);
    } else if (a > 0) {
        this->calcKin(a_bar, b_bar, a - 1, b);
        this->calcKin(a_bar, b_bar, a - 2, b);
        this->calcKin(a_bar, b_bar, a - 1, b - 1);
        this->calcOvp(a, b);
        this->calcOvp(a - 2, b);
        this->Kin_EqA12_A(state);
    } else if (b > 0) {
        this->calcKin(a_bar, b_bar, a, b - 1);
        this->calcKin(a_bar, b_bar, a - 1, b - 1);
        this->calcKin(a_bar, b_bar, a, b - 2);
        this->calcOvp(a, b);
        this->calcOvp(a, b - 2);
        this->Kin_EqA12_B(state);
    } else {
        std::cerr << TlUtils::format(
                         "Program Error: @DfHpqEngine::calcKin() a=%d, b=%d", a,
                         b)
                  << std::endl;
        abort();
    }

    this->isCalcdNuc_[stateIndex] = true;
}

void DfHpqEngine::calcNuc(const int a_bar, const int b_bar, const int a,
                          const int b, const int m) {
    if ((a_bar < 0) || (b_bar < 0) || (a < 0) || (b < 0) || (m < 0)) {
        return;
    }

    const NucState state(a_bar, b_bar, a, b, m);
    const int stateIndex = state.index();
    if (this->isCalcdNuc_[stateIndex] == true) {
        return;
    }

    if (a_bar > 0) {
        this->calcNuc(a_bar - 1, b_bar, a + 1, b, m);
        this->calcNuc(a_bar - 1, b_bar, a - 1, b, m);
        this->gradNucA(state);
    } else if (b_bar > 0) {
        this->calcNuc(a_bar, b_bar - 1, a, b + 1, m);
        this->calcNuc(a_bar, b_bar - 1, a, b - 1, m);
        this->gradNucB(state);
    } else if (a > 0) {
        this->calcNuc(a_bar, b_bar, a - 1, b, m);
        this->calcNuc(a_bar, b_bar, a - 1, b, m + 1);
        this->calcNuc(a_bar, b_bar, a - 2, b, m);
        this->calcNuc(a_bar, b_bar, a - 2, b, m + 1);
        this->calcNuc(a_bar, b_bar, a - 1, b - 1, m);
        this->calcNuc(a_bar, b_bar, a - 1, b - 1, m + 1);
        this->Nuc_EqA19_A(state);
    } else if (b > 0) {
        this->calcNuc(a_bar, b_bar, a, b - 1, m);
        this->calcNuc(a_bar, b_bar, a, b - 1, m + 1);
        this->calcNuc(a_bar, b_bar, a - 1, b - 1, m);
        this->calcNuc(a_bar, b_bar, a - 1, b - 1, m + 1);
        this->calcNuc(a_bar, b_bar, a, b - 2, m);
        this->calcNuc(a_bar, b_bar, a, b - 2, m + 1);
        this->Nuc_EqA19_B(state);
    } else {
        std::cerr
            << TlUtils::format(
                   "Program Error: @DfHpqEngine::calcNuc() a=%d, b=%d, m=%d", a,
                   b, m)
            << std::endl;
        abort();
    }

    this->isCalcdNuc_[stateIndex] = true;
}

// ovp(0, 0)
void DfHpqEngine::calcOvpSS() {
    const OvpState state(0, 0);
    const std::size_t stateIndex = state.index();
    const int index = 0;  // = this->index((0, 0, 0), (0, 0, 0));
    this->OVP_[stateIndex].resize(1);

    const double AB2 = this->B_.squareDistanceFrom(this->A_);

    const int numOfBatches = this->numOfBatches_;
    this->OVP_[stateIndex][index].resize(numOfBatches);
    for (int batch = 0; batch < numOfBatches; ++batch) {
        const double coefAB = this->coefAB_[batch];
        const double zeta = this->zeta_[batch];
        const double invZeta = 1.0 / zeta;
        const double invZeta3_2 = invZeta * std::sqrt(invZeta);
        const double eta = this->zetaA_[batch] * this->zetaB_[batch] / zeta;

        const double value = coefAB * PI3_2 * invZeta3_2 * std::exp(-eta * AB2);
        this->OVP_[stateIndex][index][batch] = value;
        //         std::cerr << TlUtils::format("OvpSS = % f", value) <<
        //         std::endl;
    }
    this->isCalcdOvp_[stateIndex] = true;
}

void DfHpqEngine::calcKinSS() {
    OvpState ovpState(0, 0);
    const std::size_t ovpStateIndex = ovpState.index();
    const int ovpIndex = 0;  // = this->index((0, 0, 0), (0, 0, 0));

    const KinState state(0, 0, 0, 0);
    const std::size_t stateIndex = state.index();
    const int index = 0;  // = this->index((0, 0), (0, 0));
    this->KIN_[stateIndex].resize(1);

    const double AB2 = this->B_.squareDistanceFrom(this->A_);
    const int numOfBatches = this->numOfBatches_;
    this->KIN_[stateIndex][index].resize(numOfBatches);
    for (int batch = 0; batch < numOfBatches; ++batch) {
        const double zeta = this->zeta_[batch];
        const double eta = this->zetaA_[batch] * this->zetaB_[batch] / zeta;
        const double coef = eta * (3.0 - 2.0 * eta * AB2);

        this->KIN_[stateIndex][index][batch] =
            coef * this->OVP_[ovpStateIndex][ovpIndex][batch];
    }
    this->isCalcdKin_[stateIndex] = true;
}

void DfHpqEngine::calcNucSS() {
    const OvpState ovpState(0, 0);
    const std::size_t ovpStateIndex = ovpState.index();
    const int ovpIndex = 0;  // = this->index((0, 0, 0), (0, 0, 0));

    const double PI_1_2 = 2.0 / std::sqrt(M_PI);
    double fmtBuf[HPQ_L_MAX];
    TlFmt& fmt = TlFmt::getInstance();
    const int sumOfAngularMomentums =
        this->a_bar_ + this->b_bar_ + this->a_ + this->b_;
    assert(sumOfAngularMomentums < HPQ_L_MAX);
    const int numOfBatches = this->numOfBatches_;

    std::size_t stateIndexPool[HPQ_L_MAX];  // 後のループで使う。再計算防止。
    int index = 0;                          // = this->index((0,0), (0,0))
    for (int m = 0; m <= sumOfAngularMomentums; ++m) {
        NucState state(0, 0, 0, 0, m);
        const std::size_t stateIndex = state.index();
        stateIndexPool[m] = stateIndex;
        this->NUC_[stateIndex].resize(1);
        this->NUC_[stateIndex][index].resize(numOfBatches);
    }

    for (int batch = 0; batch < numOfBatches; ++batch) {
        const double zeta = this->zeta_[batch];
        const double sqrt_zeta = std::sqrt(zeta);

        const double charge = -this->chargeC_;  // 符号は負
        const double ovpSS = this->OVP_[ovpStateIndex][ovpIndex]
                                       [batch];  // 軌道係数はここに含まれる
        const double coef = charge * PI_1_2 * sqrt_zeta * ovpSS;

        const TlPosition P = this->P_[batch];
        const double U = zeta * P.squareDistanceFrom(this->C_);
        fmt.getFmT(sumOfAngularMomentums, U, fmtBuf);

        //         std::cerr << TlUtils::format("calcNucSS() coef=% f, charge=%
        //         f, sqrt_zeta=% f, ovp=% f",
        //                                      coef, charge, sqrt_zeta, ovpSS)
        //                   << std::endl;
        //         std::cerr << TlUtils::format("calcNucSS() zeta=% f, P(% f,
        //         %f, %f) C(% f, % f, % f)",
        //                                      zeta, P.x(), P.y(), P.z(),
        //                                      this->C_.x(), this->C_.y(),
        //                                      this->C_.z())
        //                   << std::endl;
        for (int m = 0; m <= sumOfAngularMomentums; ++m) {
            const std::size_t stateIndex = stateIndexPool[m];
            const double value = coef * fmtBuf[m];
            //             std::cerr << TlUtils::format("calcNucSS() value=% f,
            //             m=%d, coef=% f, fmt=% f",
            //                                          value, m, coef,
            //                                          fmtBuf[m])
            //                       << std::endl;
            this->NUC_[stateIndex][index][batch] = value;
        }
    }

    for (int m = 0; m <= sumOfAngularMomentums; ++m) {
        const NucState state(0, 0, 0, 0, m);
        const std::size_t stateIndex = state.index();
        this->isCalcdNuc_[stateIndex] = true;
    }
}

// eq.A2 for A
void DfHpqEngine::Ovp_EqA2_A(const OvpState& state) {
    const int stateIndex = state.index();
    const int a = state.a;
    const int b = state.b;
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->OVP_[stateIndex].resize(numOfAmvsA * numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
        const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
        const int i = this->initiative(amvA);
        const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
        const TlAngularMomentumVector amvA2 = amvA - this->E2_[i];

        for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
            const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
            const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];

            // batch -----------------------------------------------------------
            const OvpState state1(a - 1, b);
            const std::size_t state1_index = state1.index();
            const int index1 = this->index(amvA1, amvB);

            std::size_t state2_index = 0;
            int index2 = 0;
            if (amvA2.isExist() == true) {
                const OvpState state2(a - 2, b);
                state2_index = state2.index();
                index2 = this->index(amvA2, amvB);
            }

            std::size_t state3_index = 0;
            int index3 = 0;
            if (amvB1.isExist() == true) {
                const OvpState state3(a - 1, b - 1);
                state3_index = state3.index();
                index3 = this->index(amvA1, amvB1);
            }

            const int index = this->index(amvA, amvB);
            this->OVP_[stateIndex][index].resize(numOfBatches);
            for (int batch = 0; batch < numOfBatches; ++batch) {
                double answer = 0.0;
                const double zeta = this->zeta_[batch];
                const double invZeta2 = 1.0 / (2.0 * zeta);

                // 1st term
                {
                    answer += (this->P_[batch][i] - this->A_[i]) *
                              this->OVP_[state1_index][index1][batch];
                }

                // 2nd term
                if (amvA2.isExist() == true) {
                    answer += invZeta2 * amvA1.get(i) *
                              this->OVP_[state2_index][index2][batch];
                }

                // 3rd term
                if (amvB1.isExist() == true) {
                    answer += invZeta2 * amvB.get(i) *
                              this->OVP_[state3_index][index3][batch];
                }

                this->OVP_[stateIndex][index][batch] = answer;
            }
        }
    }
}

// eq.A2 for B
void DfHpqEngine::Ovp_EqA2_B(const OvpState& state) {
    const int stateIndex = state.index();
    const int a = state.a;
    const int b = state.b;
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->OVP_[stateIndex].resize(numOfAmvsA * numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
        const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
        const int i = this->initiative(amvB);
        const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];
        const TlAngularMomentumVector amvB2 = amvB - this->E2_[i];

        for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
            const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
            const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];

            // batch -----------------------------------------------------------
            const OvpState state1(a, b - 1);
            const std::size_t state1_index = state1.index();
            const int index1 = this->index(amvA, amvB1);

            std::size_t state2_index = 0;
            int index2 = 0;
            if (amvA1.isExist() == true) {
                const OvpState state2(a - 1, b - 1);
                state2_index = state2.index();
                index2 = this->index(amvA1, amvB1);
            }

            std::size_t state3_index = 0;
            int index3 = 0;
            if (amvB2.isExist() == true) {
                const OvpState state3(a, b - 2);
                state3_index = state3.index();
                index3 = this->index(amvA, amvB2);
            }

            const int index = this->index(amvA, amvB);
            this->OVP_[stateIndex][index].resize(numOfBatches);
            for (int batch = 0; batch < numOfBatches; ++batch) {
                double answer = 0.0;
                const double zeta = this->zeta_[batch];
                const double invZeta2 = 1.0 / (2.0 * zeta);

                // 1st term
                {
                    answer += (this->P_[batch][i] - this->B_[i]) *
                              this->OVP_[state1_index][index1][batch];
                }

                // 2nd term
                if (amvA1.isExist() == true) {
                    answer += invZeta2 * amvA.get(i) *
                              this->OVP_[state2_index][index2][batch];
                }

                // 3rd term
                if (amvB2.isExist() == true) {
                    answer += invZeta2 * amvB1.get(i) *
                              this->OVP_[state3_index][index3][batch];
                }

                this->OVP_[stateIndex][index][batch] = answer;
            }
        }
    }
}

void DfHpqEngine::Kin_EqA12_A(const KinState& state) {
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int a = state.a;
    const int b = state.b;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->KIN_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsA *
                                  numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar;
         ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar;
             ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                const int i = this->initiative(amvA);
                const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
                const TlAngularMomentumVector amvA2 = amvA - this->E2_[i];

                for (int amvB_index = 0; amvB_index < numOfAmvsB;
                     ++amvB_index) {
                    const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                    const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];

                    // batch
                    // -----------------------------------------------------------
                    const KinState state1(aBar, bBar, a - 1, b);
                    const std::size_t state1_index = state1.index();
                    const int index1 =
                        this->index(amvAbar, amvBbar, amvA1, amvB);

                    std::size_t state2_index = 0;
                    int index2 = 0;
                    if (amvA2.isExist() == true) {
                        const KinState state2(aBar, bBar, a - 2, b);
                        state2_index = state2.index();
                        index2 = this->index(amvAbar, amvBbar, amvA2, amvB);
                    }

                    std::size_t state3_index = 0;
                    int index3 = 0;
                    if (amvB1.isExist() == true) {
                        const KinState state3(aBar, bBar, a - 1, b - 1);
                        state3_index = state3.index();
                        index3 = this->index(amvAbar, amvBbar, amvA1, amvB1);
                    }

                    const OvpState state4(a, b);
                    const std::size_t state4_index = state4.index();
                    const int index4 =
                        this->index(amvAbar, amvBbar, amvA, amvB);

                    std::size_t state5_index = 0;
                    int index5 = 0;
                    if (amvA2.isExist() == true) {
                        const OvpState state5(a - 2, b);
                        state5_index = state5.index();
                        index5 = this->index(amvAbar, amvBbar, amvA2, amvB);
                    }

                    const int index = this->index(amvAbar, amvBbar, amvA, amvB);
                    this->KIN_[stateIndex][index].resize(numOfBatches);
                    for (int batch = 0; batch < numOfBatches; ++batch) {
                        double answer = 0.0;
                        const double zeta = this->zeta_[batch];
                        const double eta =
                            this->zetaA_[batch] * this->zetaB_[batch] / zeta;
                        const double zetaA = this->zetaA_[batch];
                        const double invZeta2 = 1.0 / (2.0 * zeta);

                        // 1st term
                        {
                            answer += (this->P_[batch][i] - this->A_[i]) *
                                      this->KIN_[state1_index][index1][batch];
                        }

                        // 2nd term
                        if (amvA2.isExist() == true) {
                            answer += invZeta2 * amvA1.get(i) *
                                      this->KIN_[state2_index][index2][batch];
                        }

                        // 3rd term
                        if (amvB1.isExist() == true) {
                            answer += invZeta2 * amvB.get(i) *
                                      this->KIN_[state3_index][index3][batch];
                        }

                        // 4th term
                        {
                            answer += 2.0 * eta *
                                      this->OVP_[state4_index][index4][batch];
                        }

                        // 5th term
                        if (amvA2.isExist() == true) {
                            answer -= (eta / zetaA) * amvA1.get(i) *
                                      this->OVP_[state5_index][index5][batch];
                        }

                        this->KIN_[stateIndex][index][batch] = answer;
                    }
                }
            }
        }
    }
}

void DfHpqEngine::Kin_EqA12_B(const KinState& state) {
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int a = state.a;
    const int b = state.b;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->KIN_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsA *
                                  numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar;
         ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar;
             ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                const int i = this->initiative(amvB);
                const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];
                const TlAngularMomentumVector amvB2 = amvB - this->E2_[i];

                for (int amvA_index = 0; amvA_index < numOfAmvsA;
                     ++amvA_index) {
                    const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                    const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];

                    // batch
                    // -----------------------------------------------------------
                    const KinState state1(aBar, bBar, a, b - 1);
                    const std::size_t state1_index = state1.index();
                    const int index1 =
                        this->index(amvAbar, amvBbar, amvA, amvB1);

                    std::size_t state2_index = 0;
                    int index2 = 0;
                    if (amvA1.isExist() == true) {
                        const KinState state2(aBar, bBar, a - 1, b - 1);
                        state2_index = state2.index();
                        index2 = this->index(amvAbar, amvBbar, amvA1, amvB1);
                    }

                    std::size_t state3_index = 0;
                    int index3 = 0;
                    if (amvB2.isExist() == true) {
                        const KinState state3(aBar, bBar, a, b - 2);
                        state3_index = state3.index();
                        index3 = this->index(amvAbar, amvBbar, amvA, amvB2);
                    }

                    const OvpState state4(a, b);
                    const std::size_t state4_index = state4.index();
                    const int index4 =
                        this->index(amvAbar, amvBbar, amvA, amvB);

                    std::size_t state5_index = 0;
                    int index5 = 0;
                    if (amvB2.isExist() == true) {
                        const OvpState state5(a, b - 2);
                        state5_index = state5.index();
                        index5 = this->index(amvAbar, amvBbar, amvA, amvB2);
                    }

                    const int index = this->index(amvAbar, amvBbar, amvA, amvB);
                    this->KIN_[stateIndex][index].resize(numOfBatches);
                    for (int batch = 0; batch < numOfBatches; ++batch) {
                        double answer = 0.0;
                        const double zetaB = this->zetaB_[batch];
                        const double zeta = this->zeta_[batch];
                        const double eta =
                            this->zetaA_[batch] * this->zetaB_[batch] / zeta;
                        const double invZeta2 = 1.0 / (2.0 * zeta);

                        // 1st term
                        {
                            answer += (this->P_[batch][i] - this->B_[i]) *
                                      this->KIN_[state1_index][index1][batch];
                        }

                        // 2nd term
                        if (amvA1.isExist() == true) {
                            answer += invZeta2 * amvA.get(i) *
                                      this->KIN_[state2_index][index2][batch];
                        }

                        // 3rd term
                        if (amvB2.isExist() == true) {
                            answer += invZeta2 * amvB1.get(i) *
                                      this->KIN_[state3_index][index3][batch];
                        }

                        // 4th term
                        {
                            answer += 2.0 * eta *
                                      this->OVP_[state4_index][index4][batch];
                        }

                        // 5th term
                        if (amvB2.isExist() == true) {
                            answer -= (eta / zetaB) * amvB1.get(i) *
                                      this->OVP_[state5_index][index5][batch];
                        }

                        this->KIN_[stateIndex][index][batch] = answer;
                    }
                }
            }
        }
    }
}

void DfHpqEngine::Nuc_EqA19_A(const NucState& state) {
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int a = state.a;
    const int b = state.b;
    const int m = state.m;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->NUC_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsA *
                                  numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar;
         ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar;
             ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                const int i = this->initiative(amvA);
                const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
                const TlAngularMomentumVector amvA2 = amvA - this->E2_[i];

                for (int amvB_index = 0; amvB_index < numOfAmvsB;
                     ++amvB_index) {
                    const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                    const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];

                    // batch
                    // -----------------------------------------------------------
                    const NucState state1(aBar, bBar, a - 1, b, m);
                    const std::size_t state1_index = state1.index();
                    const int index1 =
                        this->index(amvAbar, amvBbar, amvA1, amvB);

                    const NucState state2(aBar, bBar, a - 1, b, m + 1);
                    const std::size_t state2_index = state2.index();
                    const int index2 = index1;

                    std::size_t state31_index = 0;
                    std::size_t state32_index = 0;
                    int index3 = 0;
                    if (amvA2.isExist() == true) {
                        const NucState state31(aBar, bBar, a - 2, b, m);
                        const NucState state32(aBar, bBar, a - 2, b, m + 1);
                        state31_index = state31.index();
                        state32_index = state32.index();
                        index3 = this->index(amvAbar, amvBbar, amvA2, amvB);
                    }

                    std::size_t state41_index = 0;
                    std::size_t state42_index = 0;
                    int index4 = 0;
                    if (amvB1.isExist() == true) {
                        const NucState state41(aBar, bBar, a - 1, b - 1, m);
                        const NucState state42(aBar, bBar, a - 1, b - 1, m + 1);
                        state41_index = state41.index();
                        state42_index = state42.index();
                        index4 = this->index(amvAbar, amvBbar, amvA1, amvB1);
                    }

                    const int index = this->index(amvAbar, amvBbar, amvA, amvB);
                    this->NUC_[stateIndex][index].resize(numOfBatches);
                    for (int batch = 0; batch < numOfBatches; ++batch) {
                        double answer = 0.0;
                        const double zetaA = this->zetaA_[batch];
                        const double zetaB = this->zetaB_[batch];
                        const double zeta = zetaA + zetaB;
                        const double invZeta2 = 1.0 / (2.0 * zeta);

                        // 1st term
                        {
                            const double term1 =
                                (this->P_[batch][i] - this->A_[i]) *
                                this->NUC_[state1_index][index1][batch];
                            answer += term1;
                        }

                        // 2nd term
                        {
                            const double term2 =
                                (this->P_[batch][i] - this->C_[i]) *
                                this->NUC_[state2_index][index2][batch];
                            answer -= term2;
                        }

                        // 3rd term
                        if (amvA2.isExist() == true) {
                            const double coef = invZeta2 * amvA1.get(i);
                            const double term3 =
                                coef *
                                (this->NUC_[state31_index][index3][batch] -
                                 this->NUC_[state32_index][index3][batch]);
                            answer += term3;
                        }

                        // 4th term
                        if (amvB1.isExist() == true) {
                            const double coef = invZeta2 * amvB.get(i);
                            answer +=
                                coef *
                                (this->NUC_[state41_index][index4][batch] -
                                 this->NUC_[state42_index][index4][batch]);
                        }

                        this->NUC_[stateIndex][index][batch] = answer;
                    }
                }
            }
        }
    }
}

void DfHpqEngine::Nuc_EqA19_B(const NucState& state) {
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int a = state.a;
    const int b = state.b;
    const int m = state.m;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->NUC_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsA *
                                  numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar;
         ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar;
             ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvB_index = 0; amvB_index < numOfAmvsB; ++amvB_index) {
                const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                const int i = this->initiative(amvB);
                const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];
                const TlAngularMomentumVector amvB2 = amvB - this->E2_[i];

                for (int amvA_index = 0; amvA_index < numOfAmvsA;
                     ++amvA_index) {
                    const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                    const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];

                    // batch
                    // -----------------------------------------------------------
                    const NucState state1(aBar, bBar, a, b - 1, m);
                    const std::size_t state1_index = state1.index();
                    const int index1 =
                        this->index(amvAbar, amvBbar, amvA, amvB1);

                    const NucState state2(aBar, bBar, a, b - 1, m + 1);
                    const std::size_t state2_index = state2.index();
                    const int index2 = index1;

                    std::size_t state31_index = 0;
                    std::size_t state32_index = 0;
                    int index3 = 0;
                    if (amvA1.isExist() == true) {
                        const NucState state31(aBar, bBar, a - 1, b - 1, m);
                        const NucState state32(aBar, bBar, a - 1, b - 1, m + 1);
                        state31_index = state31.index();
                        state32_index = state32.index();
                        index3 = this->index(amvAbar, amvBbar, amvA1, amvB1);
                    }

                    std::size_t state41_index = 0;
                    std::size_t state42_index = 0;
                    int index4 = 0;
                    if (amvB2.isExist() == true) {
                        const NucState state41(aBar, bBar, a, b - 2, m);
                        const NucState state42(aBar, bBar, a, b - 2, m + 1);
                        state41_index = state41.index();
                        state42_index = state42.index();
                        index4 = this->index(amvAbar, amvBbar, amvA, amvB2);
                    }

                    const int index = this->index(amvAbar, amvBbar, amvA, amvB);
                    this->NUC_[stateIndex][index].resize(numOfBatches);
                    for (int batch = 0; batch < numOfBatches; ++batch) {
                        double answer = 0.0;
                        const double zetaA = this->zetaA_[batch];
                        const double zetaB = this->zetaB_[batch];
                        const double zeta = zetaA + zetaB;
                        const double invZeta2 = 1.0 / (2.0 * zeta);

                        // 1st term
                        {
                            answer += (this->P_[batch][i] - this->B_[i]) *
                                      this->NUC_[state1_index][index1][batch];
                        }

                        // 2nd term
                        {
                            answer -= (this->P_[batch][i] - this->C_[i]) *
                                      this->NUC_[state2_index][index2][batch];
                        }

                        // 3rd term
                        if (amvA1.isExist() == true) {
                            const double coef = invZeta2 * amvA.get(i);
                            answer +=
                                coef *
                                (this->NUC_[state31_index][index3][batch] -
                                 this->NUC_[state32_index][index3][batch]);
                        }

                        // 4th term
                        if (amvB2.isExist() == true) {
                            const double coef = invZeta2 * amvB1.get(i);
                            answer +=
                                coef *
                                (this->NUC_[state41_index][index4][batch] -
                                 this->NUC_[state42_index][index4][batch]);
                        }

                        this->NUC_[stateIndex][index][batch] = answer;
                    }
                }
            }
        }
    }
}

void DfHpqEngine::gradKinA(const KinState& state) {
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int a = state.a;
    const int b = state.b;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->KIN_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsA *
                                  numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar;
         ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);
        const int i = this->initiative(amvAbar);
        const TlAngularMomentumVector amvAbar1 = amvAbar - this->E1_[i];

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar;
             ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
                const TlAngularMomentumVector amvA1p = amvA + this->E1_[i];

                for (int amvB_index = 0; amvB_index < numOfAmvsB;
                     ++amvB_index) {
                    const TlAngularMomentumVector amvB = amvsB.get(amvB_index);

                    // batch
                    // -----------------------------------------------------------
                    const KinState state1(aBar - 1, bBar, a + 1, b);
                    const std::size_t state1_index = state1.index();
                    const int index1 =
                        this->index(amvAbar1, amvBbar, amvA1p, amvB);

                    std::size_t state2_index = 0;
                    int index2 = 0;
                    if (amvA1.isExist() == true) {
                        const KinState state2(aBar - 1, bBar, a - 1, b);
                        state2_index = state2.index();
                        index2 = this->index(amvAbar1, amvBbar, amvA1, amvB);
                    }

                    const int index = this->index(amvAbar, amvBbar, amvA, amvB);
                    this->KIN_[stateIndex][index].resize(numOfBatches);
                    for (int batch = 0; batch < numOfBatches; ++batch) {
                        double answer = 0.0;

                        // 1st term
                        answer += 2.0 * this->zetaA_[batch] *
                                  this->KIN_[state1_index][index1][batch];

                        // 2nd term
                        if (amvA1.isExist() == true) {
                            answer -= amvA.get(i) *
                                      this->KIN_[state2_index][index2][batch];
                        }

                        this->KIN_[stateIndex][index][batch] = answer;
                    }
                }
            }
        }
    }
}

void DfHpqEngine::gradKinB(const KinState& state) {
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int a = state.a;
    const int b = state.b;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->KIN_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsA *
                                  numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar;
         ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar;
             ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);
            const int i = this->initiative(amvBbar);
            const TlAngularMomentumVector amvBbar1 = amvBbar - this->E1_[i];

            for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                const TlAngularMomentumVector amvA = amvsA.get(amvA_index);

                for (int amvB_index = 0; amvB_index < numOfAmvsB;
                     ++amvB_index) {
                    const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                    const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];
                    const TlAngularMomentumVector amvB1p = amvB + this->E1_[i];

                    // batch
                    // -----------------------------------------------------------
                    const KinState state1(aBar, bBar - 1, a, b + 1);
                    const std::size_t state1_index = state1.index();
                    const int index1 =
                        this->index(amvAbar, amvBbar1, amvA, amvB1p);

                    std::size_t state2_index = 0;
                    int index2 = 0;
                    if (amvB1.isExist() == true) {
                        const KinState state2(aBar, bBar - 1, a, b - 1);
                        state2_index = state2.index();
                        index2 = this->index(amvAbar, amvBbar1, amvA, amvB1);
                    }

                    const int index = this->index(amvAbar, amvBbar, amvA, amvB);
                    this->KIN_[stateIndex][index].resize(numOfBatches);
                    for (int batch = 0; batch < numOfBatches; ++batch) {
                        double answer = 0.0;

                        // 1st term
                        answer += 2.0 * this->zetaB_[batch] *
                                  this->KIN_[state1_index][index1][batch];

                        // 2nd term
                        if (amvB1.isExist() == true) {
                            answer -= amvB.get(i) *
                                      this->KIN_[state2_index][index2][batch];
                        }

                        this->KIN_[stateIndex][index][batch] = answer;
                    }
                }
            }
        }
    }
}

void DfHpqEngine::gradNucA(const NucState& state) {
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int a = state.a;
    const int b = state.b;
    const int m = state.m;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->NUC_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsA *
                                  numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar;
         ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);
        const int i = this->initiative(amvAbar);
        const TlAngularMomentumVector amvAbar1 = amvAbar - this->E1_[i];

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar;
             ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);

            for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                const TlAngularMomentumVector amvA = amvsA.get(amvA_index);
                const TlAngularMomentumVector amvA1 = amvA - this->E1_[i];
                const TlAngularMomentumVector amvA1p = amvA + this->E1_[i];

                for (int amvB_index = 0; amvB_index < numOfAmvsB;
                     ++amvB_index) {
                    const TlAngularMomentumVector amvB = amvsB.get(amvB_index);

                    // batch
                    // -----------------------------------------------------------
                    const NucState state1(aBar - 1, bBar, a + 1, b, m);
                    const std::size_t state1_index = state1.index();
                    const int index1 =
                        this->index(amvAbar1, amvBbar, amvA1p, amvB);

                    std::size_t state2_index = 0;
                    int index2 = 0;
                    if (amvA1.isExist() == true) {
                        const NucState state2(aBar - 1, bBar, a - 1, b, m);
                        state2_index = state2.index();
                        index2 = this->index(amvAbar1, amvBbar, amvA1, amvB);
                    }

                    const int index = this->index(amvAbar, amvBbar, amvA, amvB);
                    this->NUC_[stateIndex][index].resize(numOfBatches);
                    for (int batch = 0; batch < numOfBatches; ++batch) {
                        double answer = 0.0;

                        // 1st term
                        answer += 2.0 * this->zetaA_[batch] *
                                  this->NUC_[state1_index][index1][batch];

                        // 2nd term
                        if (amvA1.isExist() == true) {
                            answer -= amvA.get(i) *
                                      this->NUC_[state2_index][index2][batch];
                        }

                        this->NUC_[stateIndex][index][batch] = answer;
                    }
                }
            }
        }
    }
}

void DfHpqEngine::gradNucB(const NucState& state) {
    const int stateIndex = state.index();
    const int aBar = state.a_bar;
    const int bBar = state.b_bar;
    const int a = state.a;
    const int b = state.b;
    const int m = state.m;

    const TlAngularMomentumVectorSet amvsAbar(aBar);
    const TlAngularMomentumVectorSet amvsBbar(bBar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    this->NUC_[stateIndex].resize(numOfAmvsAbar * numOfAmvsBbar * numOfAmvsA *
                                  numOfAmvsB);

    const int numOfBatches = this->numOfBatches_;
    for (int amvAbar_index = 0; amvAbar_index < numOfAmvsAbar;
         ++amvAbar_index) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(amvAbar_index);

        for (int amvBbar_index = 0; amvBbar_index < numOfAmvsBbar;
             ++amvBbar_index) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(amvBbar_index);
            const int i = this->initiative(amvBbar);
            const TlAngularMomentumVector amvBbar1 = amvBbar - this->E1_[i];

            for (int amvA_index = 0; amvA_index < numOfAmvsA; ++amvA_index) {
                const TlAngularMomentumVector amvA = amvsA.get(amvA_index);

                for (int amvB_index = 0; amvB_index < numOfAmvsB;
                     ++amvB_index) {
                    const TlAngularMomentumVector amvB = amvsB.get(amvB_index);
                    const TlAngularMomentumVector amvB1 = amvB - this->E1_[i];
                    const TlAngularMomentumVector amvB1p = amvB + this->E1_[i];

                    // batch
                    // -----------------------------------------------------------
                    const NucState state1(aBar, bBar - 1, a, b + 1, m);
                    const std::size_t state1_index = state1.index();
                    const int index1 =
                        this->index(amvAbar, amvBbar1, amvA, amvB1p);

                    std::size_t state2_index = 0;
                    int index2 = 0;
                    if (amvB1.isExist() == true) {
                        const NucState state2(aBar, bBar - 1, a, b - 1, m);
                        state2_index = state2.index();
                        index2 = this->index(amvAbar, amvBbar1, amvA, amvB1);
                    }

                    const int index = this->index(amvAbar, amvBbar, amvA, amvB);
                    this->NUC_[stateIndex][index].resize(numOfBatches);
                    for (int batch = 0; batch < numOfBatches; ++batch) {
                        double answer = 0.0;

                        // 1st term
                        answer += 2.0 * this->zetaB_[batch] *
                                  this->NUC_[state1_index][index1][batch];

                        // 2nd term
                        if (amvB1.isExist() == true) {
                            answer -= amvB.get(i) *
                                      this->NUC_[state2_index][index2][batch];
                        }

                        this->NUC_[stateIndex][index][batch] = answer;
                    }
                }
            }
        }
    }
}
