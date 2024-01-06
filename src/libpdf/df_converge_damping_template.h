#ifndef DF_CONVERGE_DAMPING_TEMPLATE_H
#define DF_CONVERGE_DAMPING_TEMPLATE_H

#include <cassert>

#include "df_converge.h"

template <class SymmetricMatrix, class Vector>
class DfConverge_Damping_Template : public DfConverge {
public:
    DfConverge_Damping_Template(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Template();

protected:
    virtual void convergeRhoTilde();
    virtual void convergeKSMatrix();
    virtual void convergePMatrix();

protected:
    void convergeRhoTilde(DfObject::RUN_TYPE runType);
    void convergeKSMatrix(DfObject::RUN_TYPE runType);
    void convergePMatrix(DfObject::RUN_TYPE runType);

protected:
    virtual double getDampingFactor();

protected:
    int startIteration_;
    double dampingFactor_;
};

template <class SymmetricMatrix, class Vector>
DfConverge_Damping_Template<SymmetricMatrix, Vector>::DfConverge_Damping_Template(TlSerializeData* pPdfParam)
    : DfConverge(pPdfParam) {
    const TlSerializeData& pdfParam = *pPdfParam;
    this->startIteration_ = std::max(pdfParam["scf_acceleration/damping/number_of_damping"].getInt(), 2);
    this->dampingFactor_ = pdfParam["scf_acceleration/damping/damping_factor"].getDouble();
}

template <class SymmetricMatrix, class Vector>
DfConverge_Damping_Template<SymmetricMatrix, Vector>::~DfConverge_Damping_Template() {}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeRhoTilde() {
    this->log_.info("damping to the rho~");
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
            std::cerr << "program error. @DfConverge_Damping_Template::convergeRhoTilde()"
                      << std::endl;
            break;
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeKSMatrix() {
    this->log_.info("damping to the Fock(KS) matrix");
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
            std::cerr
                << "program error. @DfConverge_Damping_Template::convergeKSMatrix()"
                << std::endl;
            break;
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergePMatrix() {
    this->log_.info("damping to the density matrix");
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
            std::cerr << "program error. @DfConverge_Damping_Template::convergePMatrix()"
                      << std::endl;
            break;
    }
}

template <class SymmetricMatrix, class Vector>
double DfConverge_Damping_Template<SymmetricMatrix, Vector>::getDampingFactor() {
    const double dampingFactor = (*this->pPdfParam_)["scf_acceleration/damping/damping_factor"].getDouble();

    return dampingFactor;
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeRhoTilde(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->startIteration_) {
        this->log_.info("damping to cd coefficient");

        // read previous rho
        Vector prevRho = DfObject::getRho<Vector>(runType, iteration - 1);
        assert(prevRho.getSize() == this->m_nNumOfAux);

        // read current rho
        Vector currRho = DfObject::getRho<Vector>(runType, iteration);
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

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeKSMatrix(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->startIteration_) {
        this->log_.info("damping to kohn-sham matrix");

        // Fpq damping
        SymmetricMatrix currFpq = DfObject::getFpqMatrix<SymmetricMatrix>(runType, iteration);
        SymmetricMatrix prevFpq = DfObject::getFpqMatrix<SymmetricMatrix>(runType, iteration - 1);

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

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergePMatrix(const DfObject::RUN_TYPE runType) {
    const int iteration = this->m_nIteration;

    if (iteration >= this->startIteration_) {
        this->log_.info(" damping to density matrix");

        SymmetricMatrix P0_in = DfObject::getPInMatrix<SymmetricMatrix>(runType, iteration - 1);
        SymmetricMatrix P0_out = DfObject::getPOutMatrix<SymmetricMatrix>(runType, iteration - 1);

        // SymmetricMatrix currP_spin = this->getSpinDensityMatrix<SymmetricMatrix>(runType, this->m_nIteration - 1);
        // SymmetricMatrix prevP_spin = this->getSpinDensityMatrix<SymmetricMatrix>(runType, this->m_nIteration - 2);

        const double dampingFactor = this->getDampingFactor();
        this->log_.info(TlUtils::format(" damping factor = %f", dampingFactor));
        P0_in *= dampingFactor;
        P0_out *= (1.0 - dampingFactor);
        P0_in += P0_out;

        DfObject::savePInMatrix(runType, iteration, P0_in);

        // currP_spin *= (1.0 - dampingFactor);
        // prevP_spin *= dampingFactor;
        // currP_spin += prevP_spin;

        // DfObject::saveSpinDensityMatrix<SymmetricMatrix>(runType, iteration - 1, currP_spin);
    }
}

#endif  // DF_CONVERGE_DAMPING_TEMPLATE_H
