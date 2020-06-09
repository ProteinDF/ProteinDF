#ifndef DF_CDK_MATRIX_PARALLEL_H
#define DF_CDK_MATRIX_PARALLEL_H

#include "df_cdk_matrix.h"

class DfCdkMatrix_Parallel : public DfCdkMatrix {
 public:
  DfCdkMatrix_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfCdkMatrix_Parallel();

 public:
  virtual void getK();

  // ---------------------------------------------------------------------------
  // task control
  // ---------------------------------------------------------------------------
 protected:
  DfTaskCtrl* getDfTaskCtrlObject() const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 protected:
  // virtual void finalize(TlDenseSymmetricMatrixObject* pK);
};

#endif  // DF_CDK_MATRIX_PARALLEL_H
