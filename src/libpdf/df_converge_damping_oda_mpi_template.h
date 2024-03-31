#ifndef DF_CONVERGE_DAMPING_ODA_MPI_TEMPLATE_H
#define DF_CONVERGE_DAMPING_ODA_MPI_TEMPLATE_H

#include "TlCommunicate.h"
#include "df_converge_damping_oda_template.h"

// ----------------------------------------------------------------------------
// template class
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector>
class DfConverge_Damping_Oda_MpiTemplate : public DfConverge_Damping_OdaTemplate<SymmetricMatrix, Vector> {
public:
    DfConverge_Damping_Oda_MpiTemplate(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Oda_MpiTemplate();

protected:
    virtual void convergeRhoTilde();
    virtual void convergeKSMatrix();
    virtual void convergePMatrix();

protected:
    virtual double getDampingFactor();
};

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector>
DfConverge_Damping_Oda_MpiTemplate<SymmetricMatrix, Vector>::DfConverge_Damping_Oda_MpiTemplate(TlSerializeData* pPdfParam)
    : DfConverge_Damping_OdaTemplate<SymmetricMatrix, Vector>(pPdfParam) {
}

template <class SymmetricMatrix, class Vector>
DfConverge_Damping_Oda_MpiTemplate<SymmetricMatrix, Vector>::~DfConverge_Damping_Oda_MpiTemplate() {}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Oda_MpiTemplate<SymmetricMatrix, Vector>::convergeRhoTilde() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeRhoTilde();
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Oda_MpiTemplate<SymmetricMatrix, Vector>::convergeKSMatrix() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeKSMatrix();
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Oda_MpiTemplate<SymmetricMatrix, Vector>::convergePMatrix() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergePMatrix();
    }
}

template <class SymmetricMatrix, class Vector>
double DfConverge_Damping_Oda_MpiTemplate<SymmetricMatrix, Vector>::getDampingFactor() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    double dampingFactor = 0.0;
    if (rComm.isMaster()) {
        dampingFactor = DfConverge_Damping_OdaTemplate<SymmetricMatrix, Vector>::getDampingFactor();
    }
    rComm.broadcast(dampingFactor);

    return dampingFactor;
}

#endif  // DF_CONVERGE_DAMPING_ODA_MPI_TEMPLATE_H
