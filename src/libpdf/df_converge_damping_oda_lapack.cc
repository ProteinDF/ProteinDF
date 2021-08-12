#include "df_converge_damping_oda_lapack.h"

DfConverge_Damping_Oda_Lapack::DfConverge_Damping_Oda_Lapack(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Oda_tmpl<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(pPdfParam) {
}

DfConverge_Damping_Oda_Lapack::~DfConverge_Damping_Oda_Lapack() {}
