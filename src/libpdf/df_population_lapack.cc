#include "df_population_lapack.h"

DfPopulation_Lapack::DfPopulation_Lapack(TlSerializeData* pPdfParam)
    : DfPopulation_tmpl<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(
          pPdfParam) {}

DfPopulation_Lapack::~DfPopulation_Lapack() {}
