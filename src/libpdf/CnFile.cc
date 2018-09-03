#include "CnFile.h"

#include "tl_matrix_object.h"

CnFile::CnFile() {}

CnFile::~CnFile() {}

// ----------------------------------------------------------------------------
// template<class SymmetricMatrix>
// void CnFile::getDensityMatrix(const DfObject::RUN_TYPE, int itr, TlMatrixObject* pP) {
//     SymmetricMatrixType P;

//     switch (runType) {
//     case RUN_RKS:
//     {
//         P = 0.5 * this->locdMatrix<SymmetricMatrix>(DfObject::getPpqMatrixPath(RUN_RKS, itr));
//     }
//     break;

//       case RUN_UKS_ALPHA:
//       case RUN_UKS_BETA:
//         P = this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(runType, itr);
//         break;

//       case RUN_ROKS_ALPHA: {
//         P = 0.5 * this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
//                       RUN_ROKS_CLOSED, itr);
//         P += this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(RUN_ROKS_OPEN,
//                                                                itr);
//       } break;

//       case RUN_ROKS_BETA: {
//         P = 0.5 * this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
//                       RUN_ROKS_CLOSED, itr);
//       } break;

//       default:
//         this->log_.critical(
//             TlUtils::format("Program Error: %s:%d", __FILE__, __LINE__));
//         CnErr.abort();
//     }
// }

// template<class MatrixType>
// void CnFile::loadMatrix(const std::string& path) {
//     MatrixType 
// }

void CnFile::saveMatrix(const std::string& path, const TlMatrixObject& matrix) {
    matrix.save(path);
}
