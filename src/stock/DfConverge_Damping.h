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

#ifndef DFCONVERGE_DAMPING
#define DFCONVERGE_DAMPING

#include <cassert>

#include "DfConverge.h"

class DfConverge_Damping : public DfConverge {
   public:
    DfConverge_Damping(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping();

   protected:
    virtual void convergeRhoTilde();
    virtual void convergeKSMatrix();
    virtual void convergePMatrix();

   protected:
    template <class VectorType>
    void convergeRhoTilde(DfObject::RUN_TYPE runType);

    template <class SymmetricMatrixType>
    void convergeKSMatrix(DfObject::RUN_TYPE runType);

    template <class SymmetricMatrix>
    void convergePMatrix(DfObject::RUN_TYPE runType);

   protected:
    virtual double getDampingFactor() const;

   protected:
    int m_nStartIteration;
    double m_dDampingFactor;
};

template <class VectorType>
void DfConverge_Damping::convergeRhoTilde(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->m_nStartIteration) {
        this->log_.info("damping to cd coefficient");

        // read previous rho
        VectorType prevRho;
        prevRho = DfObject::getRho<VectorType>(runType, iteration - 1);
        assert(prevRho.getSize() == this->m_nNumOfAux);

        // read current rho
        VectorType currRho;
        currRho = DfObject::getRho<VectorType>(runType, iteration);
        assert(currRho.getSize() == this->m_nNumOfAux);

        // get damped rho
        const double dampingFactor = this->getDampingFactor();
        this->log_.info(TlUtils::format(" damping factor = %f", dampingFactor));
        currRho *= (1.0 - dampingFactor);
        prevRho *= dampingFactor;
        currRho += prevRho;

        // write damped current rho
        DfObject::saveRho(runType, iteration, currRho);
    }
}

template <class SymmetricMatrixType>
void DfConverge_Damping::convergeKSMatrix(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->m_nStartIteration) {
        this->log_.info("damping to kohn-sham matrix");

        // Fpq damping
        SymmetricMatrixType currFpq;
        currFpq =
            DfObject::getFpqMatrix<SymmetricMatrixType>(runType, iteration);
        SymmetricMatrixType prevFpq;
        prevFpq =
            DfObject::getFpqMatrix<SymmetricMatrixType>(runType, iteration - 1);

        // get damped Fpq
        const double dampingFactor = this->getDampingFactor();
        this->log_.info(TlUtils::format(" damping factor = %f", dampingFactor));
        currFpq *= (1.0 - dampingFactor);
        prevFpq *= dampingFactor;
        currFpq += prevFpq;

        // write damped current Fpq
        DfObject::saveFpqMatrix(runType, iteration, currFpq);
    }
}

template <class SymmetricMatrix>
void DfConverge_Damping::convergePMatrix(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->m_nStartIteration) {
        this->log_.info(" damping to density matrix");

        SymmetricMatrix currPpq =
            DfObject::getPpqMatrix<SymmetricMatrix>(runType, iteration - 1);

        SymmetricMatrix prevPpq =
            DfObject::getPpqMatrix<SymmetricMatrix>(runType, iteration - 2);

        SymmetricMatrix currP_spin =
            this->getSpinDensityMatrix<SymmetricMatrix>(runType,
                                                        this->m_nIteration - 1);
        SymmetricMatrix prevP_spin =
            this->getSpinDensityMatrix<SymmetricMatrix>(runType,
                                                        this->m_nIteration - 2);

        const double dampingFactor = this->getDampingFactor();
        this->log_.info(TlUtils::format(" damping factor = %f", dampingFactor));
        currPpq *= (1.0 - dampingFactor);
        prevPpq *= dampingFactor;
        currPpq += prevPpq;

        currP_spin *= (1.0 - dampingFactor);
        prevP_spin *= dampingFactor;
        currP_spin += prevP_spin;

        DfObject::savePpqMatrix(runType, iteration - 1, currPpq);
        DfObject::saveSpinDensityMatrix<SymmetricMatrix>(runType, iteration - 1, currP_spin);
    }
}

#endif  // DFCONVERGE_DAMPING
