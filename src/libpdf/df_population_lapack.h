#ifndef DF_POPULATION_LAPACK_H
#define DF_POPULATION_LAPACK_H

#include "df_population_tmpl.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfPopulation_Lapack
    : public DfPopulation_tmpl<TlDenseSymmetricMatrix_Lapack,
                               TlDenseVector_Lapack> {
   public:
    DfPopulation_Lapack(TlSerializeData* pPdfParam);
    virtual ~DfPopulation_Lapack();
};

#endif  // DF_POPULATION_LAPACK_H
