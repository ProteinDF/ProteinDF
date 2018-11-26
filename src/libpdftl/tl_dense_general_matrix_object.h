#ifndef TL_DENSE_GENERAL_MATRIX_OBJECT_H
#define TL_DENSE_GENERAL_MATRIX_OBJECT_H

#include "tl_dense_matrix_impl_object.h"
#include "tl_dense_vector_object.h"
#include "tl_matrix_object.h"

class TlDenseGeneralMatrixObject : public TlMatrixObject {
 public:
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
  TlDenseGeneralMatrixObject(TlDenseMatrix_ImplObject* pImpl = NULL);
  virtual ~TlDenseGeneralMatrixObject();

 public:
  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
  virtual TlMatrixObject::index_type getNumOfRows() const;
  virtual TlMatrixObject::index_type getNumOfCols() const;
  virtual void resize(const TlMatrixObject::index_type row,
                      const TlMatrixObject::index_type col);

  virtual double get(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col) const;
  virtual void set(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);
  virtual void add(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);

  /// ブロック行列を返す
  ///
  /// @param[in] row 始点となる行
  /// @param[in] col 始点となる列
  /// @param[in] rowDistance 取得する行数
  /// @param[in] colDistance 取得する列数
  /// @return row_distance × col_distance 次元のブロック行列
  void block(const TlMatrixObject::index_type row,
             const TlMatrixObject::index_type col,
             const TlMatrixObject::index_type rowDistance,
             const TlMatrixObject::index_type colDistance,
             TlDenseGeneralMatrixObject* pOut) const;

  /// 行列要素を指定された位置に上書きする
  ///
  /// @param[in] row 始点となる行
  /// @param[in] col 始点となる列
  /// @param[in] ref 参照行列
  void block(const TlMatrixObject::index_type row,
             const TlMatrixObject::index_type col,
             const TlDenseGeneralMatrixObject& ref);

  template <class VectorType>
  VectorType getRowVector(const TlMatrixObject::index_type row) const;

  template <class VectorType>
  VectorType getColVector(const TlMatrixObject::index_type col) const;

  // ---------------------------------------------------------------------------
  // Operations
  // ---------------------------------------------------------------------------
  virtual std::vector<double> diagonals() const;
  virtual double sum() const;
  virtual double trace() const;
  virtual double getRMS() const;
  virtual double getMaxAbsoluteElement(
      TlMatrixObject::index_type* outRow = NULL,
      TlMatrixObject::index_type* outCol = NULL) const;
  virtual void transposeInPlace();

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
  virtual bool load(const std::string& filePath);
  virtual bool save(const std::string& filePath) const;
  virtual bool saveText(const std::string& filePath) const;
  virtual void saveText(std::ostream& os) const;

#ifdef HAVE_HDF5
  virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);
  virtual bool saveHdf5(const std::string& filepath,
                        const std::string& h5path) const;
#endif  // HAVE_HDF5

  virtual void dump(double* buf, const std::size_t size) const;
  virtual void restore(const double* buf, const std::size_t size);

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlDenseMatrix_ImplObject* pImpl_;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlDenseGeneralMatrixObject& mat);

// -----------------------------------------------------------------------------
template <class VectorType>
VectorType TlDenseGeneralMatrixObject::getRowVector(
    const TlMatrixObject::index_type row) const {
  const TlMatrixObject::index_type size = this->getNumOfCols();
  VectorType v(size);
  for (TlMatrixObject::index_type i = 0; i < size; ++i) {
    v.set(i, this->get(row, i));
  }

  return v;
}

template <class VectorType>
VectorType TlDenseGeneralMatrixObject::getColVector(
    const TlMatrixObject::index_type col) const {
  const TlMatrixObject::index_type size = this->getNumOfRows();
  VectorType v(size);
  for (TlMatrixObject::index_type i = 0; i < size; ++i) {
    v.set(i, this->get(i, col));
  }

  return v;
}

#endif  // TL_DENSE_GENERAL_MATRIX_OBJECT_H
