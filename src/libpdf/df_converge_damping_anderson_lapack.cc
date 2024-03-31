#include "df_converge_damping_anderson_lapack.h"

DfConverge_Damping_Anderson_Lapack::DfConverge_Damping_Anderson_Lapack(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Anderson_Template<TlDenseGeneralMatrix_Lapack, TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(pPdfParam) {
}

DfConverge_Damping_Anderson_Lapack::~DfConverge_Damping_Anderson_Lapack() {
}
