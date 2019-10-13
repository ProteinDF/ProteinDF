#ifndef DF_POPULATION_EIGEN_H
#define DF_POPULATION_EIGEN_H

#include "df_population_tmpl.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"

class DfPopulation_Eigen
    : public DfPopulation_tmpl<TlDenseSymmetricMatrix_Eigen,
                               TlDenseVector_Eigen> {
   public:
    DfPopulation_Eigen(TlSerializeData* pPdfParam);
    virtual ~DfPopulation_Eigen();
};

#endif  // DF_POPULATION_EIGEN_H
