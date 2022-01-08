#ifndef DF_POPULATION_SCALAPACK_H
#define DF_POPULATION_SCALAPACK_H

#include "df_population_tmpl.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_vector_scalapack.h"

typedef DfPopulation_tmpl<TlDenseSymmetricMatrix_Scalapack, TlDenseVector_Scalapack> DfPopulationTmplScalapack;

class DfPopulation_Scalapack : public DfPopulationTmplScalapack {
public:
    DfPopulation_Scalapack(TlSerializeData* pPdfParam);
    virtual ~DfPopulation_Scalapack();
};

#endif  // DF_POPULATION_SCALAPACK_H
