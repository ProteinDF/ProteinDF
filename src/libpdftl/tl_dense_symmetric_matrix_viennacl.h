#ifndef TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_H
#define TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_H

#include "tl_dense_symmetric_matrix_object.h"

class TlDenseGeneralMatrix_ViennaCL;
class TlDenseVector_ViennaCL;
class TlDenseSymmetricMatrix_Eigen;
class TlSparseSymmetricMatrix_ViennaCL;

class TlDenseSymmetricMatrix_ViennaCL : public TlDenseSymmetricMatrixObject {
  // ---------------------------------------------------------------------------
  //
  // ---------------------------------------------------------------------------
 public:
  enum EIG_METHOD { EIG_QR, EIG_POWERITERATION };

  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseSymmetricMatrix_ViennaCL(
      const TlMatrixObject::index_type dim = 1);
  TlDenseSymmetricMatrix_ViennaCL(const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL(const TlDenseGeneralMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL(const TlSparseSymmetricMatrix_ViennaCL& rhs);
#ifdef HAVE_EIGEN
  TlDenseSymmetricMatrix_ViennaCL(const TlDenseSymmetricMatrix_Eigen& rhs);
#endif  // HAVE_EIGEN

  virtual ~TlDenseSymmetricMatrix_ViennaCL();

  virtual void vtr2mat(const std::vector<double>& vtr);

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseSymmetricMatrix_ViennaCL& operator=(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL& operator=(
      const TlDenseSymmetricMatrix_Eigen& rhs);

  const TlDenseSymmetricMatrix_ViennaCL operator+(
      const TlDenseSymmetricMatrix_ViennaCL& rhs) const;
  const TlDenseSymmetricMatrix_ViennaCL operator-(
      const TlDenseSymmetricMatrix_ViennaCL& rhs) const;
  const TlDenseSymmetricMatrix_ViennaCL operator*(
      const TlDenseSymmetricMatrix_ViennaCL& rhs) const;

  TlDenseSymmetricMatrix_ViennaCL& operator+=(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL& operator-=(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL& operator*=(const double coef);
  TlDenseSymmetricMatrix_ViennaCL& operator/=(const double coef);
  TlDenseSymmetricMatrix_ViennaCL& operator*=(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  const TlDenseSymmetricMatrix_ViennaCL& dotInPlace(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);

  bool eig(TlDenseVector_ViennaCL* pEigVal,
           TlDenseGeneralMatrix_ViennaCL* pEigVec,
           EIG_METHOD eigMethod = EIG_QR) const;
  TlDenseSymmetricMatrix_ViennaCL inverse() const;

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------
 protected:
  bool eig_powerIteration(TlDenseVector_ViennaCL* pEigVal,
                          TlDenseGeneralMatrix_ViennaCL* pEigVec) const;

  bool eig_QR(TlDenseVector_ViennaCL* pEigVal,
              TlDenseGeneralMatrix_ViennaCL* pEigVec) const;

  // ---------------------------------------------------------------------------
  // Others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_ViennaCL;
  friend class TlDenseSymmetricMatrix_Eigen;
  friend class TlSparseSymmetricMatrix_ViennaCL;

  friend TlDenseGeneralMatrix_ViennaCL operator*(
      const TlDenseGeneralMatrix_ViennaCL& mat1,
      const TlDenseSymmetricMatrix_ViennaCL& mat2);
  friend TlDenseGeneralMatrix_ViennaCL operator*(
      const TlDenseSymmetricMatrix_ViennaCL& mat1,
      const TlDenseGeneralMatrix_ViennaCL& mat2);
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_H
