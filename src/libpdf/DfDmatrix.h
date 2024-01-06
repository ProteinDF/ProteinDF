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

#ifndef DFDMATRIX_H
#define DFDMATRIX_H

#include <cmath>
#include <string>
#include <vector>

#include "CnError.h"
#include "DfObject.h"
#include "TlTime.h"
#include "common.h"

/// 分子軌道の占有数の決定と密度行列の作成を行うクラス
///
/// 収束させるための処理として、軌道の重なりの方法、射影演算子法を用いる
class DfDmatrix : public DfObject {
public:
    DfDmatrix(TlSerializeData* pPdfParam);
    virtual ~DfDmatrix();

public:
    virtual void run();

protected:
    template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
    void run_impl();

    // occupation
    // ---------------------------------------------------------------
    template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
    void makeOccupation(const DfObject::RUN_TYPE runType);

    template <typename GeneralMatrix, typename Vector>
    Vector getOccupationUsingOverlap(DfObject::RUN_TYPE runType);

    template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
    Vector getOccupationUsingProjection(DfObject::RUN_TYPE runType);

    // density matrix
    // -----------------------------------------------------------
    template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
    void generateDensityMatrix(DfObject::RUN_TYPE runType);

    template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
    SymmetricMatrix calcDensMatrix(DfObject::RUN_TYPE runType, const GeneralMatrix& C, double base);

    // others
    // -------------------------------------------------------------------
    virtual void checkOccupation(const TlDenseVectorObject& prevOcc, const TlDenseVectorObject& currOcc);
    virtual void printOccupation(const TlDenseVectorObject& occ);

    void printTwoVectors(const std::vector<double>& a, const std::vector<double>& b, const std::string& title,
                         int pnumcol);

    // --------------------------------------------------------------------------
protected:
    enum ORBITAL_CORRESPONDENCE_METHOD { OCM_NONE,
                                         OCM_OVERLAP,
                                         OCM_PROJECTION };

protected:
    ORBITAL_CORRESPONDENCE_METHOD orbitalCorrespondenceMethod_;
};

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
void DfDmatrix::run_impl() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->makeOccupation<GeneralMatrix, SymmetricMatrix, Vector>(RUN_RKS);
            this->generateDensityMatrix<GeneralMatrix, SymmetricMatrix, Vector>(RUN_RKS);
            break;

        case METHOD_UKS:
            this->makeOccupation<GeneralMatrix, SymmetricMatrix, Vector>(RUN_UKS_ALPHA);
            this->generateDensityMatrix<GeneralMatrix, SymmetricMatrix, Vector>(RUN_UKS_ALPHA);

            this->makeOccupation<GeneralMatrix, SymmetricMatrix, Vector>(RUN_UKS_BETA);
            this->generateDensityMatrix<GeneralMatrix, SymmetricMatrix, Vector>(RUN_UKS_BETA);
            break;

        case METHOD_ROKS:
            this->makeOccupation<GeneralMatrix, SymmetricMatrix, Vector>(RUN_ROKS_CLOSED);
            this->generateDensityMatrix<GeneralMatrix, SymmetricMatrix, Vector>(RUN_ROKS_CLOSED);

            this->makeOccupation<GeneralMatrix, SymmetricMatrix, Vector>(RUN_ROKS_OPEN);
            this->generateDensityMatrix<GeneralMatrix, SymmetricMatrix, Vector>(RUN_ROKS_OPEN);

            // ROKS_alpha,beta
            {
                SymmetricMatrix PC =
                    DfObject::getSpinDensityMatrix<SymmetricMatrix>(RUN_ROKS_CLOSED, this->m_nIteration);
                SymmetricMatrix PO = DfObject::getSpinDensityMatrix<SymmetricMatrix>(RUN_ROKS_OPEN, this->m_nIteration);

                SymmetricMatrix PA = PC + PO;
                DfObject::saveSpinDensityMatrix(RUN_ROKS_ALPHA, this->m_nIteration, PA);

                SymmetricMatrix PB = PC;
                DfObject::saveSpinDensityMatrix(RUN_ROKS_BETA, this->m_nIteration, PB);
            }
            break;

        default:
            CnErr.abort();
            break;
    }
}

// ----------------------------------------------------------------------------
// occ
// ----------------------------------------------------------------------------
template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
void DfDmatrix::makeOccupation(const DfObject::RUN_TYPE runType) {
    Vector currOcc;
    switch (this->orbitalCorrespondenceMethod_) {
        case OCM_OVERLAP:
            this->log_.info(" orbital correspondence method: MO-overlap");
            currOcc = this->getOccupationUsingOverlap<GeneralMatrix, Vector>(runType);
            currOcc.save(this->getOccupationPath(runType));
            break;

        case OCM_PROJECTION:
            this->log_.info(" orbital correspondence method: MO-projection");
            currOcc = this->getOccupationUsingProjection<GeneralMatrix, SymmetricMatrix, Vector>(runType);
            currOcc.save(this->getOccupationPath(runType));
            break;

        default:
            this->log_.info(" orbital correspondence method: none");
            this->log_.info(TlUtils::format(" use occ: %s", this->getOccupationPath(runType).c_str()));
            break;
    }
}

// 軌道の重なりの対応のルーチン ====================================
template <typename GeneralMatrix, typename Vector>
Vector DfDmatrix::getOccupationUsingOverlap(DfObject::RUN_TYPE runType) {
    this->log_.info(" MO overlap method is started.");
    const int nNumOfMOs = this->m_nNumOfMOs;

    this->log_.info(" load previous C' matrix");
    GeneralMatrix prevCprime(nNumOfMOs, nNumOfMOs);
    {
        // read current orbital in orthonormal basis
        GeneralMatrix Cprime;
        Cprime = DfObject::getCprimeMatrix<GeneralMatrix>(runType, this->m_nIteration);

        // read previous orbital in orthonormal basis
        if ((this->m_nIteration != 1) ||
            (this->initialGuessType_ == GUESS_LCAO || this->initialGuessType_ == GUESS_HUCKEL)) {
            prevCprime = DfObject::getCprimeMatrix<GeneralMatrix>(runType, this->m_nIteration - 1);
        } else if (this->m_nIteration == 1) {
            prevCprime = Cprime;
        }

        // calculate and print MO overlap between previous and current orbitals
        // C <-- C_pre^dagger * C_cur, in program  D <-- B^dagger * A
        prevCprime.transposeInPlace();
        prevCprime *= Cprime;
    }

    Vector prevOcc = this->getOccVtr<Vector>(runType);

    // construct occupation number of current orbital with MO overlap matrix
    // 旧MO(pre)がどの新MO(crr)との重なりが一番大きいかを探す
    this->log_.info(" check overlap");
    Vector currOcc(nNumOfMOs);
    {
        std::vector<bool> g(nNumOfMOs, false);
        bool bListHeaderOutput = false;
        for (int pre = 0; pre < nNumOfMOs; ++pre) {
            const double dPrevOcc = prevOcc.get(pre);
            if ((std::fabs(dPrevOcc - 1.00) < 1.0e-10) ||  // 占有数が 1 or 2 の時
                (std::fabs(dPrevOcc - 2.00) < 1.0e-10)) {
                const TlDenseVector_Lapack prevCprimeRowVector =
                    prevCprime.template getRowVector_tmpl<TlDenseVector_Lapack>(pre);
                double xval = 0.0;
                int xord = -1;
                for (int crr = 0; crr < nNumOfMOs; ++crr) {
                    const double val = std::fabs(prevCprimeRowVector.get(crr));
                    if ((g[crr] == false) && (val > xval)) {
                        xval = val;  // fabs(prevCprime(pre, crr));
                        xord = crr;
                    }
                }

                if (xord == -1) {
                    this->log_.info(TlUtils::format(" crr MO %d th is not corresponded!\n", pre));
                    // CnErr.abort("DfDmatrix", "", "", "MO Overlap is not
                    // corresponded");
                }

                if (xval < 0.3) {
                    if (bListHeaderOutput == false) {
                        this->log_.info(" ##### MO Overlap is less than 0.3 #####\n");
                        bListHeaderOutput = true;
                    }
                    this->log_.info(
                        TlUtils::format("pre MO %6d th -> crr MO %6d th %12.8lf\n", (pre + 1), (xord + 1), xval));
                }
                currOcc.set(xord, dPrevOcc);
                g[xord] = true;
            }
        }
    }

    this->log_.info(" check occupation vectors");
    this->checkOccupation(prevOcc, currOcc);

    this->log_.info(" finish");
    return currOcc;
}

// 射影演算子法 ====================================================
template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
Vector DfDmatrix::getOccupationUsingProjection(const DfObject::RUN_TYPE runType) {
    this->log_.info("orbital_overlap_method is mo-projection");
    const index_type numOfMOs = this->m_nNumOfMOs;

    // read molecular orbital ^(n)
    GeneralMatrix C = DfObject::getCMatrix<GeneralMatrix>(runType, this->m_nIteration);
    GeneralMatrix Ct = C.transpose();

    SymmetricMatrix S = DfObject::getSpqMatrix<SymmetricMatrix>();

    // calculation of projection diagonal
    std::vector<double> pd(numOfMOs);

    // read input density matrix
    const SymmetricMatrix D = DfObject::getPInMatrix<SymmetricMatrix>(runType, this->m_nIteration);
    const GeneralMatrix SDS = S * D * S;

    // diagonal
    const GeneralMatrix judge = Ct * SDS * C;
    for (DfObject::index_type i = 0; i < numOfMOs; ++i) {
        pd[i] = judge.get(i, i);
    }

    for (DfObject::index_type i = 0; i < numOfMOs; ++i) {
        if (pd[i] < 0.0) {
            pd[i] *= -1.0;
        }
    }

    // set order of MO (0から始まる軌道の番号)
    std::vector<double> pd_ord(numOfMOs);
    for (DfObject::index_type k = 0; k < numOfMOs; ++k) {
        pd_ord[k] = k;
    }

    // sort pd
    for (DfObject::index_type i = 0; i < numOfMOs - 1; ++i) {
        for (DfObject::index_type j = i; j < numOfMOs; ++j) {
            if (pd[i] < pd[j]) {
                std::swap(pd[i], pd[j]);
                std::swap(pd_ord[i], pd_ord[j]);
            }
        }
    }

    // output
    double occupation = 0.0;
    std::string title = "";
    switch (runType) {
        case RUN_RKS:
            occupation = 2.0;
            title = "projection diagonal of MO";
            break;

        case RUN_UKS_ALPHA:
            occupation = 1.0;
            title = "projection diagonal of MO(alpha)";
            break;

        case RUN_UKS_BETA:
            occupation = 1.0;
            title = "projection diagonal of MO(beta)";
            break;

        case RUN_ROKS_CLOSED:
            occupation = 2.0;
            title = "projection diagonal of MO(closed)";
            break;

        case RUN_ROKS_OPEN:
            occupation = 1.0;
            title = "projection diagonal of MO(open)";
            break;

        default:
            this->log_.critical("program error.");
            break;
    }
    this->printTwoVectors(pd_ord, pd, title, 10);

    // store electrons
    Vector currOcc(numOfMOs);
    for (DfObject::index_type k = 0; k < numOfMOs; ++k) {
        currOcc.set(static_cast<index_type>(pd_ord[k]), occupation);
    }

    // check
    {
        Vector prevOcc = this->getOccVtr<Vector>(runType);
        const double sumOfPrevOcc = prevOcc.sum();
        const double sumOfCurrOcc = currOcc.sum();

        this->log_.info(TlUtils::format("prev occ = %10.2f", sumOfPrevOcc));
        this->log_.info(TlUtils::format("prev occ = %10.2f", sumOfCurrOcc));
        if (std::fabs(sumOfPrevOcc - sumOfCurrOcc) > 1.0E-5) {
            this->log_.info("projection occupation was failed.");
            CnErr.abort();
        }
    }

    return currOcc;
}

// ----------------------------------------------------------------------------
// density matrix
// ----------------------------------------------------------------------------
template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
inline void DfDmatrix::generateDensityMatrix(const DfObject::RUN_TYPE runType) {
    switch (runType) {
        case RUN_RKS: {
            GeneralMatrix C = DfObject::getCMatrix<GeneralMatrix>(runType, this->m_nIteration);
            SymmetricMatrix P = this->calcDensMatrix<GeneralMatrix, SymmetricMatrix, Vector>(runType, C, 2.0);
            this->saveSpinDensityMatrix(runType, this->m_nIteration, P);

            P *= 2.0;
            this->savePOutMatrix(runType, this->m_nIteration, P);
        } break;

        case RUN_UKS_ALPHA:
        case RUN_UKS_BETA: {
            GeneralMatrix C = DfObject::getCMatrix<GeneralMatrix>(runType, this->m_nIteration);
            SymmetricMatrix P = this->calcDensMatrix<GeneralMatrix, SymmetricMatrix, Vector>(runType, C, 1.0);
            this->saveSpinDensityMatrix(runType, this->m_nIteration, P);
            this->savePOutMatrix(runType, this->m_nIteration, P);
        } break;

        case RUN_ROKS_CLOSED: {
            GeneralMatrix C = DfObject::getCMatrix<GeneralMatrix>(RUN_ROKS, this->m_nIteration);
            SymmetricMatrix P = this->calcDensMatrix<GeneralMatrix, SymmetricMatrix, Vector>(RUN_ROKS_CLOSED, C, 2.0);
            this->saveSpinDensityMatrix(runType, this->m_nIteration, P);
            P *= 2.0;
            this->savePOutMatrix(runType, this->m_nIteration, P);

        } break;

        case RUN_ROKS_OPEN: {
            GeneralMatrix C = DfObject::getCMatrix<GeneralMatrix>(RUN_ROKS, this->m_nIteration);
            SymmetricMatrix P = this->calcDensMatrix<GeneralMatrix, SymmetricMatrix, Vector>(RUN_ROKS_OPEN, C, 1.0);
            this->savePOutMatrix(runType, this->m_nIteration, P);
        } break;

        default:
            CnErr.abort("program error.");
            break;
    }
}

template <typename GeneralMatrix, typename SymmetricMatrix, typename Vector>
inline SymmetricMatrix DfDmatrix::calcDensMatrix(DfObject::RUN_TYPE runType, const GeneralMatrix& C, double base) {
    const TlMatrixObject::index_type nNumOfMOs = C.getNumOfCols();

    SymmetricMatrix E(nNumOfMOs);
    {
        Vector occ;
        occ.load(this->getOccupationPath(runType));

        TlMatrixObject::index_type max_k = std::min(nNumOfMOs, occ.getSize());
        for (TlMatrixObject::index_type k = 0; k < max_k; ++k) {
            if (std::fabs(occ.get(k) - base) < 1.0e-10) {
                E.set(k, k, 1.0);
            }
        }
    }
    const SymmetricMatrix P = C * E * C.transpose();
    return P;
}

#ifdef HAVE_VIENNACL
template <>
inline TlDenseSymmetricMatrix_ViennaCL
DfDmatrix::calcDensMatrix<TlDenseGeneralMatrix_ViennaCL, TlDenseSymmetricMatrix_ViennaCL, TlDenseVector_ViennaCL>(
    DfObject::RUN_TYPE runType, const TlDenseGeneralMatrix_ViennaCL& C, double base) {
    const TlMatrixObject::index_type nNumOfMOs = C.getNumOfCols();

    TlDenseSymmetricMatrix_ViennaCL E(nNumOfMOs);
    {
        TlDenseVector_Eigen occ;
        occ.load(this->getOccupationPath(runType));

        TlMatrixObject::index_type max_k = std::min(nNumOfMOs, occ.getSize());
        TlDenseSymmetricMatrix_Eigen tmpE(nNumOfMOs);
        for (TlMatrixObject::index_type k = 0; k < max_k; ++k) {
            if (std::fabs(occ.get(k) - base) < 1.0e-10) {
                tmpE.set(k, k, 1.0);
            }
        }
        E = tmpE;
    }
    const TlDenseSymmetricMatrix_ViennaCL P = C * E * C.transpose();
    return P;
}
#endif  // HAVE_VIENNACL

#endif  // DFDMATRIX_H
