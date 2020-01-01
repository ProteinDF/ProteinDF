#ifndef DF_CD_MATRIX_OBJECT_H
#define DF_CD_MATRIX_OBJECT_H

#include "DfObject.h"
class DfCDMatrixObject : public DfObject {};

template <typename SymmetricMatrix, typename Vector>
Vector DfCdkMatrix::getScreenedDensityMatrix(const RUN_TYPE runType,
                                             const PQ_PairArray& I2PR) {
  SymmetricMatrix P;
  switch (runType) {
    case RUN_RKS:
      P = 0.5 *
          this->getPpqMatrix<SymmetricMatrix>(RUN_RKS, this->m_nIteration - 1);
      break;

    case RUN_UKS_ALPHA:
    case RUN_UKS_BETA:
      P = this->getPpqMatrix<SymmetricMatrix>(runType, this->m_nIteration - 1);
      break;

    case RUN_ROKS_ALPHA: {
      P = 0.5 *
          this->getPpqMatrix<SymmetricMatrix>(RUN_ROKS_CLOSED,
                                              this->m_nIteration - 1);
      P += this->getPpqMatrix<TlDenseSymmetricMatrix_Lapack>(
          RUN_ROKS_OPEN, this->m_nIteration - 1);
    } break;

    case RUN_ROKS_BETA: {
      P = 0.5 *
          this->getPpqMatrix<SymmetricMatrix>(RUN_ROKS_CLOSED,
                                              this->m_nIteration - 1);
    } break;

    default:
      this->log_.critical(
          TlUtils::format("Program Error: %s:%d", __FILE__, __LINE__));
      CnErr.abort();
  }

  const std::size_t numOfI = I2PR.size();
  Vector answer(numOfI);

  for (std::size_t i = 0; i < numOfI; ++i) {
    const Index2& pair = I2PR[i];
    const index_type r = pair.index1();
    const index_type c = pair.index2();
    answer.set(i, P.get(r, c));
  }

  return answer;
}

template <typename SymmetricMatrix typename Vector>
Vector DfCDMatrixObject::getScreenedDensityMatrix(const PQ_PairArray& I2PQ) {
  const SymmetricMatrix P = this->getPMatrix();
  const std::size_t numOfI = I2PQ.size();
  Vector answer(numOfI);

  for (std::size_t i = 0; i < numOfI; ++i) {
    const Index2& pair = I2PQ[i];
    const index_type r = pair.index1();
    const index_type c = pair.index2();
    const double coef = (r != c) ? 2.0 : 1.0;
    answer.set(i, coef * P.get(r, c));
  }

  return answer;
}

#endif  // DF_CD_MATRIX_OBJECT_H
