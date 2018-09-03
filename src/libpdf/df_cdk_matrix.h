#ifndef DF_CDK_MATRIX_H
#define DF_CDK_MATRIX_H

#include "DfObject.h"
#include "common.h"
#include "datatype.h"

class DfTaskCtrl;

class DfCdkMatrix : public DfObject {
 public:
  DfCdkMatrix(TlSerializeData* pPdfParam);
  virtual ~DfCdkMatrix();

 public:
  void getK();

 protected:
  template <typename SymmetricMatrix, typename Vector,
            typename SparseGeneralMatrix, typename SparseSymmetricMatrix>
  void getK_method();

  template <typename SymmetricMatrix, typename Vector,
            typename SparseGeneralMatrix, typename SparseSymmetricMatrix>
  void getK_L(DfObject::RUN_TYPE runType);

  template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector>
  void getK_byLjk(const RUN_TYPE runType);

  template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector,
            typename SparseGeneralMatrix>
  void getK_byLjk_useTransMatrix(const RUN_TYPE runType);

  template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector,
            typename SparseGeneralMatrix, typename SparseSymmetricMatrix>
  void getK_byLjk_useSparseMatrix(const RUN_TYPE runType);

 protected:
  DfTaskCtrl* getDfTaskCtrlObject() const;
  PQ_PairArray getI2PQ(const std::string& filepath);

  template <typename SymmetricMatrix, typename Vector>
  SymmetricMatrix getCholeskyVector(const Vector& L_col,
                                    const PQ_PairArray& I2PQ);

  template <typename SparseGeneralMatrix>
  SparseGeneralMatrix getTrans_I2PQ_Matrix(const PQ_PairArray& I2PQ);

  template <typename SymmetricMatrix, typename Vector,
            typename SparseGeneralMatrix>
  SymmetricMatrix convert_I2PQ(const SparseGeneralMatrix& I2PQ_mat,
                               const Vector& L);

 protected:
  bool useMmapMatrix_;
  int sparseMatrixLevel_;
};

#endif  // DF_CDK_MATRIX_H
