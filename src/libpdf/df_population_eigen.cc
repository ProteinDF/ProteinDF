#include "df_population_eigen.h"

DfPopulation_Eigen::DfPopulation_Eigen(TlSerializeData* pPdfParam)
    : DfPopulation_tmpl<TlDenseSymmetricMatrix_Eigen, TlDenseVector_Eigen>(
          pPdfParam) {}

DfPopulation_Eigen::~DfPopulation_Eigen() {}
