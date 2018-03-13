#ifndef TL_DENSE_GENERAL_MATRIX_ABSTRACT_H
#define TL_DENSE_GENERAL_MATRIX_ABSTRACT_H

#include "tl_matrix_object.h"
#include "tl_matrix_utils.h"

class TlDenseGeneralMatrixAbstract : public TlMatrixObject {
 public:
  TlDenseGeneralMatrixAbstract() : TlMatrixObject(TlMatrixObject::CSFD) {}
  virtual ~TlDenseGeneralMatrixAbstract() {}

public:
    virtual void resize(const TlMatrixObject::index_type row,
                        const TlMatrixObject::index_type col) =0;

public:
    bool load(const std::string& filePath);
    bool save(const std::string& filePath) const;

 public:
  virtual bool load(const std::string& filePath, double* pBuf,
                    const TlMatrixObject::size_type length);
  virtual bool save(const std::string& filePath, const double* pBuf,
                    const TlMatrixObject::size_type length) const;
};

#endif  // TL_DENSE_GENERAL_MATRIX_ABSTRACT_H
