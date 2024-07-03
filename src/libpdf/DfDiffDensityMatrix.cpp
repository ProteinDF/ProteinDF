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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#ifdef HAVE_LAPACK
#include "tl_dense_symmetric_matrix_lapack.h"
#endif  // HAVA_LAPACK

#ifdef HAVE_EIGEN
#include "tl_dense_symmetric_matrix_eigen.h"
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
#include "tl_dense_symmetric_matrix_viennacl.h"
#endif  // HAVE_VIENNACL

#include "DfDiffDensityMatrix.h"
#include "TlFile.h"

DfDiffDensityMatrix::DfDiffDensityMatrix(TlSerializeData* pPdfData)
    : DfObject(pPdfData) {
}

DfDiffDensityMatrix::~DfDiffDensityMatrix() {}

void DfDiffDensityMatrix::exec() {
    switch (this->linearAlgebraPackage_) {
#ifdef HAVE_LAPACK
        case LAP_LAPACK:
            this->exec_lapack();
            break;
#endif  // HAVE_LAPACK

#ifdef HAVE_EIGEN
        case LAP_EIGEN:
            this->exec_eigen();
            break;
#endif  // HAVE_EIGEN

        default:
            this->log_.critical(TlUtils::format("program error: @%s,%d", __FILE__, __LINE__));
            this->log_.critical(TlUtils::format("linear algebra package: %d", this->linearAlgebraPackage_));
            CnErr.abort();
            break;
    }
}

#ifdef HAVE_LAPACK
void DfDiffDensityMatrix::exec_lapack() {
    this->log_.info("make diff matrix by using LAPACK");
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_RKS, this->m_nIteration);
            break;

        case METHOD_UKS:
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_UKS_ALPHA, this->m_nIteration);
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_UKS_BETA, this->m_nIteration);
            break;

        case METHOD_ROKS:
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_ROKS_CLOSED, this->m_nIteration);
            this->calc<TlDenseSymmetricMatrix_Lapack>(RUN_ROKS_OPEN, this->m_nIteration);
            break;

        default:
            std::cerr << "program error. " << __FILE__ << ":" << __LINE__
                      << std::endl;
            abort();
            break;
    }
}
#endif // HAVE_LAPACK

#ifdef HAVE_EIGEN
void DfDiffDensityMatrix::exec_eigen() {
    this->log_.info("make diff matrix by using LAPACK");
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->calc<TlDenseSymmetricMatrix_Eigen>(RUN_RKS, this->m_nIteration);
            break;

        case METHOD_UKS:
            this->calc<TlDenseSymmetricMatrix_Eigen>(RUN_UKS_ALPHA, this->m_nIteration);
            this->calc<TlDenseSymmetricMatrix_Eigen>(RUN_UKS_BETA, this->m_nIteration);
            break;

        case METHOD_ROKS:
            this->calc<TlDenseSymmetricMatrix_Eigen>(RUN_ROKS_CLOSED, this->m_nIteration);
            this->calc<TlDenseSymmetricMatrix_Eigen>(RUN_ROKS_OPEN, this->m_nIteration);
            break;

        default:
            std::cerr << "program error. " << __FILE__ << ":" << __LINE__
                      << std::endl;
            abort();
            break;
    }
}
#endif // HAVE_EIGEN
