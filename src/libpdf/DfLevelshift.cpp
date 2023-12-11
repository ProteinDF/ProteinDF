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

#include "DfLevelshift.h"

#include "CnError.h"
#include "Fl_Tbl_Fragment.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"

DfLevelshift::DfLevelshift(TlSerializeData* pPdfParam, int num_iter)
    : DfObject(pPdfParam) {}

DfLevelshift::~DfLevelshift() {}

void DfLevelshift::DfLshiftMain() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->main(RUN_RKS, this->m_nIteration);
            break;

        case METHOD_UKS:
            this->main(RUN_UKS_ALPHA, this->m_nIteration);
            this->main(RUN_UKS_BETA, this->m_nIteration);
            break;

        case METHOD_ROKS:
            this->main(RUN_ROKS, this->m_nIteration);
            break;

        default:
            CnErr.abort();
            break;
    }
}

// for extended QCLO method
void DfLevelshift::DfLshiftQclo(const std::string& fragname, int norbcut) {
    this->m_nNumOfMOs = norbcut;

    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->main(RUN_RKS, this->m_nIteration, fragname);
            break;

        case METHOD_UKS:
            this->main(RUN_UKS_ALPHA, this->m_nIteration, fragname);
            this->main(RUN_UKS_BETA, this->m_nIteration, fragname);
            break;

        case METHOD_ROKS:
            this->main(RUN_ROKS, this->m_nIteration, fragname);
            break;

        default:
            CnErr.abort();
            break;
    }
}

void DfLevelshift::main(const RUN_TYPE runType, int iteration,
                        const std::string& fragname, bool bPdfQcloMode) {
    TlSerializeData& pdfParam = *(this->pPdfParam_);

    // construct level shift matrix for "F' matrix"
    const double ls_closed_mo = pdfParam["level_shift/closed_mo"].getDouble();
    const double ls_open_mo = pdfParam["level_shift/open_mo"].getDouble();
    const double ls_virtual_mo = pdfParam["level_shift/virtual_mo"].getDouble();

    const double delta_group_closed = pdfParam["level_shift/delta_group_closed"].getDouble();
    const double delta_group_open = pdfParam["level_shift/delta_group_open"].getDouble();
    const double delta_group_virtual = pdfParam["level_shift/delta_group_virtual"].getDouble();

    this->log_.info("construct Level Shift Operator.");

    this->log_.info("RKS,UKS calculation");
    this->log_.info(TlUtils::format(" level shift for closed  MO = %8.2lf", ls_closed_mo));
    this->log_.info(TlUtils::format(" level shift for open    MO = %8.2lf", ls_open_mo));
    this->log_.info(TlUtils::format(" level shift for virtual MO = %8.2lf\n", ls_virtual_mo));
    this->log_.info(TlUtils::format(" level shift between closed  MO = %8.2lf", delta_group_closed));
    this->log_.info(TlUtils::format(" level shift between open    MO = %8.2lf", delta_group_open));
    this->log_.info(TlUtils::format(" level shift between virtual MO = %8.2lf", delta_group_virtual));

    // check numOfMOs
    if (this->m_nNumOfMOs <= 0) {
        this->log_.info(TlUtils::format("#MOs = %d", this->m_nNumOfMOs));
        this->log_.info("read #MOs from X matrix");
        TlDenseGeneralMatrix_Lapack X = DfObject::getXMatrix<TlDenseGeneralMatrix_Lapack>();
        const index_type newNumOfMOs = X.getNumOfCols();
        this->log_.info(TlUtils::format("update MOs: %d", newNumOfMOs));
        this->m_nNumOfMOs = newNumOfMOs;
    }

    const index_type numOfMOs = this->m_nNumOfMOs;
    this->log_.info(TlUtils::format("# MOs: %d", numOfMOs));

    // prepare "beta" vector (level-shift values)
    TlDenseSymmetricMatrix_Lapack beta(numOfMOs);
    {
        this->log_.info("load occupation vector");
        TlDenseVector_Lapack vOcc;
        if (runType == RUN_ROKS) {
            vOcc.load(DfObject::getOccupationPath(RUN_ROKS_CLOSED));
            TlDenseVector_Lapack occ_open;
            occ_open.load(DfObject::getOccupationPath(RUN_ROKS_OPEN));
            vOcc += occ_open;
        } else {
            vOcc.load(DfObject::getOccupationPath(runType));
        }
        if (vOcc.getSize() < numOfMOs) {
            this->log_.info(TlUtils::format("occ resize to %d", numOfMOs));
            vOcc.resize(numOfMOs);
        }

        double shift_closed = 0.0;
        double shift_open = 0.0;
        double shift_virtual = 0.0;

        this->log_.info("construct beta vector");
        int orb_id = 0;
        for (int k = 0; k < numOfMOs; k++) {
            if (bPdfQcloMode == true) {
                int frag_id = -1;
                Fl_Tbl_Fragment Tfrag(Fl_Geometry((*this->pPdfParam_)["coordinates"]));
                if (Tfrag.getFragment(k) != frag_id) {
                    continue;
                }
            }

            if (std::fabs(vOcc.get(k)) < 1.0E-10) {
                beta.add(orb_id, orb_id, ls_virtual_mo + shift_virtual);
                shift_virtual += delta_group_virtual;
            } else if (fabs(vOcc.get(k) - 2.0) < 1.0E-10) {
                beta.add(orb_id, orb_id, ls_closed_mo + shift_closed);
                shift_closed += delta_group_closed;
            } else if (fabs(vOcc.get(k) - 1.0) < 1.0E-10) {
                beta.add(orb_id, orb_id, ls_open_mo + shift_open);
                shift_open += delta_group_open;
            } else {
                CnErr.abort("DfLevelshift", "", "", "occ is illegal");
            }

            orb_id++;
        }
    }

    // prepare "C'" matrix
    this->log_.info("prepare C' matrix");
    TlDenseGeneralMatrix_Lapack Cprime(numOfMOs, numOfMOs);
    {
        if (iteration == 1 && ((this->initialGuessType_ != GUESS_LCAO) &&
                               (this->initialGuessType_ != GUESS_HUCKEL))) {
            this->log_.info("construct level shift operator with fukue's Rou at iteration==1");
            this->log_.info("level shift value is simply add to F' matrix (guess Rou is solution of DFT)");

            // construct "C'" matrix for initial Rou
            const int numOfMOs = this->m_nNumOfMOs;
            for (int k = 0; k < numOfMOs; ++k) {
                Cprime.set(k, k, 1.0);
            }
        } else {
            // "read previous C' matrix"
            if (TlFile::isExistFile(DfObject::getCprimeMatrixPath(runType, iteration - 1)) == true) {
                this->log_.info(TlUtils::format("load C' matrix: %d", iteration - 1));
                Cprime = DfObject::getCprimeMatrix<TlDenseGeneralMatrix_Lapack>(runType, iteration - 1);
            } else {
                this->log_.info(TlUtils::format("load C matrix: %d", iteration - 1));
                const TlDenseGeneralMatrix_Lapack C = DfObject::getCMatrix<TlDenseGeneralMatrix_Lapack>(runType, iteration - 1);

                this->log_.info("load Xinv matrix");
                const TlDenseGeneralMatrix_Lapack Xinv = DfObject::getXInvMatrix<TlDenseGeneralMatrix_Lapack>();

                this->log_.info("make C' matrix");
                Cprime = Xinv * C;
            }

            if (Cprime.getNumOfCols() != numOfMOs) {
                this->log_.warn(TlUtils::format("C' matrix size: (%d, %d)", Cprime.getNumOfRows(), Cprime.getNumOfCols()));
                this->log_.warn(TlUtils::format("num of MOs = %d", numOfMOs));
                this->log_.warn("the dimension is not consistency, but continue...");
            }
        }
    }

    // calculate "C' * beta * C'^dagger" = "C' * (C' * beta)^dagger"
    this->log_.info("calc L matrix");
    {
        // calc. "C' * beta"
        TlDenseGeneralMatrix_Lapack Cprime_beta = Cprime * beta;
        // TlDenseGeneralMatrix_Lapack Cprime_beta(numOfMOs, numOfMOs);
        // for (int i = 0; i < numOfAOs; ++i) {
        //     for (int j = 0; j < numOfMOs; ++j) {
        //         Cprime_beta.set(i, j, Cprime.get(i, j) * beta.get(j));
        //     }
        // }

        // (C' * beta)^dagger
        Cprime_beta.transposeInPlace();

        // calc. "C' * (C'*beta)^dagger"
        Cprime *= Cprime_beta;
    }

    // "read F' matrix"
    this->log_.info("load F' matrix");
    TlDenseSymmetricMatrix_Lapack Fprime;
    {
        Fprime = DfObject::getFprimeMatrix<TlDenseSymmetricMatrix_Lapack>(runType, iteration);

        if ((Fprime.getNumOfRows() != numOfMOs) || (Fprime.getNumOfCols() != numOfMOs)) {
            this->log_.warn(TlUtils::format("rowDim of F' matrix = %d", Fprime.getNumOfRows()));
            this->log_.warn(TlUtils::format("colDIm of F' matrix = %d", Fprime.getNumOfCols()));
            this->log_.warn(TlUtils::format("number_mo_basis = %d", this->m_nNumOfMOs));
            this->log_.warn("DfLevelshift dimension is not consistency, but continue");
        }
    }

    // construct "F' + Level shift matrix "
    this->log_.info("construct F' matrix");
    Fprime += Cprime;

    // write "shifted F' matrix"
    DfObject::saveFprimeMatrix(runType, iteration, Fprime);
}
