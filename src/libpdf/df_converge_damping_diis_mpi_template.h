#ifndef DF_CONVERGE_DAMPING_DIIS_MPI_TEMPLATE_H
#define DF_CONVERGE_DAMPING_DIIS_MPI_TEMPLATE_H

#include "TlCommunicate.h"
#include "df_converge_damping_diis_template.h"

template <class Matrix, class SymmetricMatrix, class Vector>
class DfConverge_Damping_Diis_Mpi_Template : public DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector> {
public:
    DfConverge_Damping_Diis_Mpi_Template(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Diis_Mpi_Template();

protected:
    virtual void convergeRhoTilde(const DfObject::RUN_TYPE runType);
    virtual void convergeKSMatrix(const DfObject::RUN_TYPE runType);
    virtual void convergePMatrix(const DfObject::RUN_TYPE runType);
};

template <class Matrix, class SymmetricMatrix, class Vector>
DfConverge_Damping_Diis_Mpi_Template<Matrix, SymmetricMatrix, Vector>::DfConverge_Damping_Diis_Mpi_Template(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>(pPdfParam) {
}

template <class Matrix, class SymmetricMatrix, class Vector>
DfConverge_Damping_Diis_Mpi_Template<Matrix, SymmetricMatrix, Vector>::~DfConverge_Damping_Diis_Mpi_Template() {
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Diis_Mpi_Template<Matrix, SymmetricMatrix, Vector>::convergeRhoTilde(const DfObject::RUN_TYPE runType) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergeRhoTilde(runType);
    }
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Diis_Mpi_Template<Matrix, SymmetricMatrix, Vector>::convergeKSMatrix(const DfObject::RUN_TYPE runType) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergeKSMatrix(runType);
    }
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Diis_Mpi_Template<Matrix, SymmetricMatrix, Vector>::convergePMatrix(const DfObject::RUN_TYPE runType) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergePMatrix(runType);
    }
}

#endif  // DF_CONVERGE_DAMPING_DIIS_MPI_TEMPLATE_H
