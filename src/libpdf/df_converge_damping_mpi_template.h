#ifndef DF_CONVERGE_DAMPING_MPI_TEMPLATE_H
#define DF_CONVERGE_DAMPING_MPI_TEMPLATE_H

#include "TlCommunicate.h"
#include "df_converge_damping_template.h"

// ----------------------------------------------------------------------------
// template class
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector>
class DfConverge_Damping_Mpi_Template : public DfConverge_Damping_Template<SymmetricMatrix, Vector> {
public:
    DfConverge_Damping_Mpi_Template(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Mpi_Template();

protected:
    virtual void convergeRhoTilde();
    virtual void convergeKSMatrix();
    virtual void convergePMatrix();

protected:
    // virtual double getDampingFactor();
};

// ----------------------------------------------------------------------------
// implementation
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector>
DfConverge_Damping_Mpi_Template<SymmetricMatrix, Vector>::DfConverge_Damping_Mpi_Template(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Template<SymmetricMatrix, Vector>(pPdfParam) {
    // this->log_.info("converge: damping method (MPI)");
}

template <class SymmetricMatrix, class Vector>
DfConverge_Damping_Mpi_Template<SymmetricMatrix, Vector>::~DfConverge_Damping_Mpi_Template() {}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Mpi_Template<SymmetricMatrix, Vector>::convergeRhoTilde() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeRhoTilde();
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Mpi_Template<SymmetricMatrix, Vector>::convergeKSMatrix() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeKSMatrix();
    }
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_Mpi_Template<SymmetricMatrix, Vector>::convergePMatrix() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergePMatrix();
    }
}

// template <class SymmetricMatrix, class Vector>
// double DfConverge_Damping_Mpi_Template<SymmetricMatrix, Vector>::getDampingFactor() {
//     TlCommunicate& rComm = TlCommunicate::getInstance();

//     double dampingFactor = 0.0;
//     if (rComm.isMaster()) {
//         dampingFactor = DfConverge_Damping_Template<SymmetricMatrix, Vector>::getDampingFactor();
//     }
//     rComm.broadcast(dampingFactor);

//     return dampingFactor;
// }

#endif  // DF_CONVERGE_DAMPING_MPI_TEMPLATE_H
