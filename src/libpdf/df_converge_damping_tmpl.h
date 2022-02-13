#ifndef DF_CONVERGE_DAMPING_TMPL_H
#define DF_CONVERGE_DAMPING_TMPL_H

#include <cassert>

#include "DfConverge.h"

template <class SymmetricMatrix, class Vector>
class DfConverge_Damping_tmpl : public DfConverge {
   public:
    DfConverge_Damping_tmpl(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_tmpl();

   protected:
    virtual void convergeRhoTilde();
    virtual void convergeKSMatrix();
    virtual void convergePMatrix();

   protected:
    void convergeRhoTilde(const DfObject::RUN_TYPE runType);
    void convergeKSMatrix(const DfObject::RUN_TYPE runType);
    void convergePMatrix(const DfObject::RUN_TYPE runType);

   protected:
    virtual void getDampingFactor();

   protected:
    int startIteration_;
    double dampingFactor_;
};

///////////////////////////////////////////////////////////////////////////////
// implementation
///////////////////////////////////////////////////////////////////////////////

template <class SymmetricMatrix, class Vector>
DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::DfConverge_Damping_tmpl(TlSerializeData* pPdfParam)
    : DfConverge(pPdfParam) {
    this->startIteration_ = std::max((*this->pPdfParam_)["scf_acceleration/damping/number_of_damping"].getInt(), 2);

    this->getDampingFactor();
}

template <class SymmetricMatrix, class Vector>
DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::~DfConverge_Damping_tmpl() {}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::getDampingFactor() {
    this->dampingFactor_ = (*this->pPdfParam_)["scf_acceleration/damping/damping_factor"].getDouble();
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::convergeRhoTilde() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->convergeRhoTilde(DfObject::RUN_RKS);
            break;
        case METHOD_UKS:
            this->convergeRhoTilde(DfObject::RUN_UKS_ALPHA);
            this->convergeRhoTilde(DfObject::RUN_UKS_BETA);
            break;
        case METHOD_ROKS:
            this->convergeRhoTilde(DfObject::RUN_ROKS_CLOSED);
            this->convergeRhoTilde(DfObject::RUN_ROKS_OPEN);
            break;
        default:
            std::cerr << "program error. @DfConverge_Damping::convergeRhoTilde()"
                      << std::endl;
            break;
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::convergeKSMatrix() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->convergeKSMatrix(DfObject::RUN_RKS);
            break;
        case METHOD_UKS:
            this->convergeKSMatrix(DfObject::RUN_UKS_ALPHA);
            this->convergeKSMatrix(DfObject::RUN_UKS_BETA);
            break;
        case METHOD_ROKS:
            this->convergeKSMatrix(DfObject::RUN_ROKS_CLOSED);
            this->convergeKSMatrix(DfObject::RUN_ROKS_OPEN);
            break;
        default:
            std::cerr << "program error. @DfConverge_Damping::convergeKSMatrix()"
                      << std::endl;
            break;
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::convergePMatrix() {
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->convergePMatrix(DfObject::RUN_RKS);
            break;
        case METHOD_UKS:
            this->convergePMatrix(DfObject::RUN_UKS_ALPHA);
            this->convergePMatrix(DfObject::RUN_UKS_BETA);
            break;
        case METHOD_ROKS:
            this->convergePMatrix(DfObject::RUN_ROKS_CLOSED);
            this->convergePMatrix(DfObject::RUN_ROKS_OPEN);
            break;
        default:
            std::cerr << "program error. @DfConverge_Damping::convergePMatrix()"
                      << std::endl;
            break;
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::convergeRhoTilde(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->startIteration_) {
        this->log_.info("damping to cd coefficient");

        // read previous rho
        Vector prevRho;
        prevRho = DfObject::getRho<Vector>(runType, iteration - 1);
        assert(prevRho.getSize() == this->m_nNumOfAux);

        // read current rho
        Vector currRho;
        currRho = DfObject::getRho<Vector>(runType, iteration);
        assert(currRho.getSize() == this->m_nNumOfAux);

        // get damped rho
        const double dampingFactor = this->dampingFactor_;
        this->log_.info(TlUtils::format(" damping factor = %f", dampingFactor));
        currRho *= (1.0 - dampingFactor);
        prevRho *= dampingFactor;
        currRho += prevRho;

        // write damped current rho
        DfObject::saveRho(runType, iteration, currRho);
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::convergeKSMatrix(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->startIteration_) {
        this->log_.info("damping to kohn-sham matrix");

        // Fpq damping
        SymmetricMatrix currFpq;
        currFpq = DfObject::getFpqMatrix<SymmetricMatrix>(runType, iteration);
        SymmetricMatrix prevFpq;
        prevFpq = DfObject::getFpqMatrix<SymmetricMatrix>(runType, iteration - 1);

        // get damped Fpq
        const double dampingFactor = this->dampingFactor_;
        this->log_.info(TlUtils::format(" damping factor = %f", dampingFactor));
        currFpq *= (1.0 - dampingFactor);
        prevFpq *= dampingFactor;
        currFpq += prevFpq;

        // write damped current Fpq
        DfObject::saveFpqMatrix(runType, iteration, currFpq);
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_tmpl<SymmetricMatrix, Vector>::convergePMatrix(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->startIteration_) {
        this->log_.info(" damping to density matrix");

        SymmetricMatrix currPpq = DfObject::getPpqMatrix<SymmetricMatrix>(runType, iteration - 1);
        SymmetricMatrix prevPpq = DfObject::getPpqMatrix<SymmetricMatrix>(runType, iteration - 2);

        SymmetricMatrix currP_spin = this->getSpinDensityMatrix<SymmetricMatrix>(runType, this->m_nIteration - 1);
        SymmetricMatrix prevP_spin = this->getSpinDensityMatrix<SymmetricMatrix>(runType, this->m_nIteration - 2);

        const double dampingFactor = this->dampingFactor_;
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

#endif  // DF_CONVERGE_DAMPING_TMPL_H
