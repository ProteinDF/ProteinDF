#ifndef DFCONVERGE_DAMPING
#define DFCONVERGE_DAMPING

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

    template <class SymmetricMatrixType>
    void convergePMatrix(DfObject::RUN_TYPE runType);

protected:
    int m_nStartIteration;
    double m_dDampingFactor;
};


template <class VectorType>
void DfConverge_Damping::convergeRhoTilde(const DfObject::RUN_TYPE runType)
{
    const int iteration = this->m_nIteration;

    if (iteration >= this->m_nStartIteration) {
        this->log_.info("damping to cd coefficient");

        // read previous rho
        VectorType prevRho;
        prevRho = DfObject::getRho<VectorType>(runType, iteration -1);
        assert(prevRho.getSize() == this->m_nNumOfAux);

        // read current rho
        VectorType currRho;
        currRho = DfObject::getRho<VectorType>(runType, iteration);
        assert(currRho.getSize() == this->m_nNumOfAux);

        // get damped rho
        this->log_.info(TlUtils::format(" damping factor = %f", this->m_dDampingFactor));
        currRho *= (1.0 - this->m_dDampingFactor);
        prevRho *= this->m_dDampingFactor;
        currRho += prevRho;

        // write damped current rho
        DfObject::saveRho(runType, iteration, currRho);
    }
}


template <class SymmetricMatrixType>
void DfConverge_Damping::convergeKSMatrix(const DfObject::RUN_TYPE runType)
{
    const int iteration = this->m_nIteration;

    if (iteration >= this->m_nStartIteration) {
        this->log_.info("damping to kohn-sham matrix");

        // Fpq damping
        SymmetricMatrixType currFpq;
        currFpq = DfObject::getFpqMatrix<SymmetricMatrixType>(runType, iteration);
        SymmetricMatrixType prevFpq;
        prevFpq = DfObject::getFpqMatrix<SymmetricMatrixType>(runType, iteration -1);

        // get damped Fpq
        this->log_.info(TlUtils::format(" damping factor = %f", this->m_dDampingFactor));
        currFpq *= (1.0 - this->m_dDampingFactor);
        prevFpq *= this->m_dDampingFactor;
        currFpq += prevFpq;

        // write damped current Fpq
        DfObject::saveFpqMatrix(runType, iteration, currFpq);
    }
}


template <class SymmetricMatrixType>
void DfConverge_Damping::convergePMatrix(const DfObject::RUN_TYPE runType)
{
    const int iteration = this->m_nIteration;

    if (iteration >= (this->m_nStartIteration +1)) {
        this->log_.info(" damping to density matrix");

        // Fpq damping
        SymmetricMatrixType currPpq;
        currPpq = DfObject::getPpqMatrix<SymmetricMatrixType>(runType, iteration -1);
        SymmetricMatrixType prevPpq;
        prevPpq = DfObject::getPpqMatrix<SymmetricMatrixType>(runType, iteration -2);

        // get damped Ppq
        this->log_.info(TlUtils::format(" damping factor = %f", this->m_dDampingFactor));
        currPpq *= (1.0 - this->m_dDampingFactor);
        prevPpq *= this->m_dDampingFactor;
        currPpq += prevPpq;

        // write damped current Ppq
        DfObject::savePpqMatrix(runType, iteration -1, currPpq);
    }
}

#endif // DFCONVERGE_DAMPING
