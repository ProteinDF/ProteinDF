#include "df_converge_damping_lapack.h"

DfConverge_Damping_Lapack::DfConverge_Damping_Lapack(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Template<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(pPdfParam) {
}

DfConverge_Damping_Lapack::~DfConverge_Damping_Lapack() {}
