#ifndef DF_TOTAL_ENERGY_EIGEN_H
#define DF_TOTAL_ENERGY_EIGEN_H

#include "DfOverlapX.h"
#include "df_total_energy_tmpl.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"

class DfTotalEnergy_Eigen
    : public DfTotalEnergy_tmpl<
          TlDenseGeneralMatrix_Eigen, TlDenseSymmetricMatrix_Eigen,
          TlDenseVector_Eigen, DfOverlapX> {
   public:
    DfTotalEnergy_Eigen(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy_Eigen();
};

#endif  // DF_TOTAL_ENERGY_EIGEN_H
