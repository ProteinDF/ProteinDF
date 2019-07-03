#ifndef DF_TOTAL_ENERGY_PARALLEL_TMPL_H
#define DF_TOTAL_ENERGY_PARALLEL_TMPL_H

#include "df_total_energy_tmpl.h"
#include "TlCommunicate.h"

// ----------------------------------------------------------------------------
// template class
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector, class DfOverlapType>
class DfTotalEnergy_Parallel_tmpl : public DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType> {
   public:
    DfTotalEnergy_Parallel_tmpl(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy_Parallel_tmpl();

   protected:
    virtual SymmetricMatrix getSpinDensityMatrix(const DfObject::RUN_TYPE runType, const int iteration);
    virtual Vector getRho();
    virtual Vector getEps(const DfObject::RUN_TYPE runType);

    virtual void saveParams();
};

// ----------------------------------------------------------------------------
// constructor & destructor
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector, class DfOverlapType>
DfTotalEnergy_Parallel_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::DfTotalEnergy_Parallel_tmpl(TlSerializeData* pPdfParam) 
: DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::DfTotalEnergy_tmpl(pPdfParam) {

}

template <class GeneralMatrix, class SymmetricMatrix, class Vector, class DfOverlapType>
DfTotalEnergy_Parallel_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::~DfTotalEnergy_Parallel_tmpl() 
{
}

// ----------------------------------------------------------------------------
// density matrix
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector, class DfOverlapType>
SymmetricMatrix DfTotalEnergy_Parallel_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::getSpinDensityMatrix(const DfObject::RUN_TYPE runType, const int iteration) {
    SymmetricMatrix matrix;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        matrix = DfObject::getSpinDensityMatrix<SymmetricMatrix>(runType, iteration);
    }
    rComm.broadcast(&matrix);

    return matrix;
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector, class DfOverlapType>
Vector DfTotalEnergy_Parallel_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::getRho() {
    Vector rho;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        rho = DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::getRho();
    }
    rComm.broadcast(&rho);

    return rho;
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector, class DfOverlapType>
Vector DfTotalEnergy_Parallel_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::getEps(const DfObject::RUN_TYPE runType) {
    Vector eps;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        eps = DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::getEps(runType);
    }
    rComm.broadcast(&eps);

    return eps;
}

// ----------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector, class DfOverlapType>
void DfTotalEnergy_Parallel_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::saveParams() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::saveParams();
    }
}

#endif // DF_TOTAL_ENERGY_PARALLEL_TMPL_H