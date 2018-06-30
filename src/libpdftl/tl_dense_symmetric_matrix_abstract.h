#ifndef TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H
#define TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H

#include "tl_matrix_object.h"

class TlDenseSymmetricMatrixAbstract : public TlMatrixObject {
 public:
  TlDenseSymmetricMatrixAbstract() : TlMatrixObject(TlMatrixObject::RLHD) {}
  virtual ~TlDenseSymmetricMatrixAbstract() {}

 public:
  virtual bool load(const std::string& filePath, double* pBuf,
                    const TlMatrixObject::size_type length);
  virtual bool save(const std::string& filePath, const double* pBuf,
                    const TlMatrixObject::size_type length) const;
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H
