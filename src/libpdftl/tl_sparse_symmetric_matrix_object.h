#ifndef TL_SPARSE_SYMMETRIC_MATRIX_OBJECT_H
#define TL_SPARSE_SYMMETRIC_MATRIX_OBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "TlSerializeData.h"
#include "tl_matrix_object.h"

class TlSparseMatrix_ImplObject;

class TlSparseSymmetricMatrixObject : public TlMatrixObject {
   public:
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
    TlSparseSymmetricMatrixObject(TlSparseMatrix_ImplObject* pImpl = NULL);
    virtual ~TlSparseSymmetricMatrixObject();

   public:
    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
    virtual TlMatrixObject::index_type getNumOfRows() const;
    virtual TlMatrixObject::index_type getNumOfCols() const;
    virtual void resize(const TlMatrixObject::index_type dim);

    virtual double get(const TlMatrixObject::index_type row,
                       const TlMatrixObject::index_type col) const;
    virtual void set(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);
    virtual void add(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
    virtual bool load(const std::string& filePath);
    virtual bool save(const std::string& filePath) const;

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
   protected:
    TlSparseMatrix_ImplObject* pImpl_;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlSparseSymmetricMatrixObject& mat);

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_OBJECT_H
