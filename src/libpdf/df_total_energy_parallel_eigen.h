#ifndef DF_TOTAL_ENERGY_PARALLEL_EIGEN_H
#define DF_TOTAL_ENERGY_PARALLEL_EIGEN_H

#include "df_total_energy_parallel_tmpl.h"

#include "DfOverlapX.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"

class DfTotalEnergy_Parallel_Eigen
    : public DfTotalEnergy_Parallel_tmpl<TlDenseGeneralMatrix_Eigen,
                                         TlDenseSymmetricMatrix_Eigen,
                                         TlDenseVector_Eigen, DfOverlapX> {
   public:
    DfTotalEnergy_Parallel_Eigen(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy_Parallel_Eigen();
};

#endif  // DF_TOTAL_ENERGY_PARALLEL_EIGEN_H