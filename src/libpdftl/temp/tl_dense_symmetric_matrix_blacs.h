// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef TL_DENSE_SYMMETRIC_MATRIX_BLACS_H
#define TL_DENSE_SYMMETRIC_MATRIX_BLACS_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <list>
#include "tl_dense_general_matrix_blacs.h"
#include "tl_dense_vector_blacs.h"
#include "tl_dense_vector_lapack.h"
#include "tl_sparse_symmetric_matrix.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

class TlDenseSymmetricMatrix_blacs : public TlDenseGeneralMatrix_blacs {
 public:
  enum DIAGONAL_METHOD { QR, DIVIDE_AND_CONQUER };

 public:
  explicit TlDenseSymmetricMatrix_blacs(int dim = 1);
  TlDenseSymmetricMatrix_blacs(const TlDenseSymmetricMatrix_blacs& rhs);
  TlDenseSymmetricMatrix_blacs(const TlDenseGeneralMatrix_blacs& rhs);
  explicit TlDenseSymmetricMatrix_blacs(const TlDistributedVector& rhs,
                                        int dim);

  virtual ~TlDenseSymmetricMatrix_blacs();

 protected:
  virtual size_type getNumOfElements() const;

 public:
  void resize(int nSize);

 public:
  // need to call by all processes.
  TlDenseSymmetricMatrix_blacs& operator=(
      const TlDenseSymmetricMatrix_blacs& rhs);
  TlDenseSymmetricMatrix_blacs& operator=(
      const TlDenseGeneralMatrix_blacs& rhs);

  /// 要素の絶対値の最大値を返す
  ///
  /// @param[out] outRow 該当要素の行
  /// @param[out] outCol 該当要素の列
  /// @return 要素の絶対値の最大値
  virtual double getMaxAbsoluteElement(int* pOutRow = NULL,
                                       int* pOutCol = NULL) const;

  virtual double getRMS() const;

  // need to call by all processes.
  virtual double get(index_type row, index_type col) const;
  virtual void set(index_type row, index_type col, double value);
  virtual void add(index_type row, index_type col, double value);

  /** ローカルから該当する要素があれば値を返す
   *
   *  範囲外であれば0.0を返す
   *  すべてのプロセスが同時に呼ぶ必要はない
   */
  virtual double getLocal(index_type row, index_type col) const;

  /// 指定した行の要素から構成されるベクトルを返す
  ///
  /// @param[in] nRow 指定する行
  virtual TlDenseVector_Lapack getRowVector(index_type row) const;

  /// 指定した列の要素から構成されるベクトルを返す
  ///
  /// @param[in] nCol 指定する列
  virtual TlDenseVector_Lapack getColumnVector(index_type col) const;

  /// 全行列を各プロセスに均等分割された疎行列を返す
  // TlSparseSymmetricMatrix getPartialMatrix(double threshold = 1.0E-16) const;

  /// 指定された要素を持つ疎行列を返す
  // void getPartialMatrix(TlSparseSymmetricMatrix& ioMatrix) const;
  bool getSparseMatrixX(TlSparseSymmetricMatrix* pMatrix,
                        bool isFinalize = false) const;

  TlDenseVector_Lapack getPartialMatrix(int* pStartRow, int* pEndRow,
                                        int* pStartCol, int* pEndCol) const;

  /// 各ノードが与えた疎行列を大域行列に加算する。
  ///
  /// 全ノードがこの関数を呼び出す必要がある。
  void mergeSparseMatrix(const TlSparseSymmetricMatrix& M);

  virtual void mergeSparseMatrixAsync(const TlSparseMatrix* pMatrix,
                                      bool isFinalize = false);

  /// 固有値を求める
  ///
  /// @param[out] pEigVal 固有値が格納されたベクトル
  /// @param[out] pEigVec 固有値ベクトルが格納された行列
  /// @retval true 固有値が求められた
  /// @retval false エラーが発生した
  virtual bool eig(TlDenseVector_Lapack* pEigVal,
                   TlDenseGeneralMatrix_blacs* pEigVec,
                   DIAGONAL_METHOD method = DIVIDE_AND_CONQUER) const;

  TlDenseSymmetricMatrix_blacs inverse() const;

  // const TlDenseSymmetricMatrix_blacs& dot(const TlDenseSymmetricMatrix_blacs&
  // X);
  virtual double sum() const;

  /// 指定された要素(グローバル)がどのプロセスが所有しているかを返す
  virtual int getProcIdForIndex(index_type globalRow,
                                index_type globalCol) const;

  virtual size_type getIndex(index_type globalRow, index_type globalCol) const;

 private:
  // TlDenseGeneralMatrix_blacs choleskyFactorization(
  //     const double threshold = 1.0E-16) const;
  // TlDenseGeneralMatrix_blacs choleskyFactorization_mod(
  //     const double threshold = 1.0E-16) const;
 public:
  TlDenseGeneralMatrix_blacs choleskyFactorization_mod2(
      const double threshold = 1.0E-16) const;

 public:
  /// 指定された入力ストリームが読み込み可能かどうかを返す
  ///
  /// @param[in,out] ifs 入力ファイルストリーム
  /// @retval true 読み取り可能
  /// @retval false 読み取り不可能
  // static bool isLoadable(std::ifstream& ifs);
  // static bool isLoadable(const std::string& rFilePath);

  virtual bool load(const std::string& sFilePath);
  virtual bool load(std::fstream& fs);
  virtual bool save(const std::string& sFilePath) const;

// protected:
//     virtual void saveElements(TlMatrixFile* pFileMatrix, const
//     std::vector<TlMatrixObject::MatrixElement>& elements) const;

#ifdef HAVE_HDF5
 public:
  virtual bool saveHdf5(const std::string& filepath,
                        const std::string& h5path) const;
  virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);

 protected:
  // bool saveHdf5(const std::string& filepath, const std::string& h5path, int
  // saveMatType) const;
  virtual size_type getArrayIndex(const index_type row,
                                  const index_type col) const;
// void saveElements(TlHdf5Utils* pH5, const std::string& path,
//                   const std::vector<TlMatrixObject::MatrixElement>&
//                   elements) const;
#endif  // HAVE_HDF5

 protected:
  virtual std::vector<TlMatrixObject::MatrixElement> getMatrixElementsInLocal()
      const;

 public:
  // need to call by all processes.
  template <typename T>
  void print(T& out) const;

 protected:
  bool load_RLHD(std::ifstream& ifs);
  bool load_CLHD(std::ifstream& ifs);

 protected:
  virtual bool loadLocal(const std::string& filePath);
  virtual bool saveLocal(const std::string& filePath) const;

  // friend
 protected:
#ifdef HAVE_SCALAPACK
  friend TlDenseGeneralMatrix_blacs operator*(
      const TlDenseSymmetricMatrix_blacs& X,
      const TlDenseSymmetricMatrix_blacs& Y);
  friend TlDenseGeneralMatrix_blacs operator*(
      const TlDenseSymmetricMatrix_blacs& X,
      const TlDenseGeneralMatrix_blacs& Y);
  friend TlDenseGeneralMatrix_blacs operator*(
      const TlDenseGeneralMatrix_blacs& X,
      const TlDenseSymmetricMatrix_blacs& Y);
  friend TlDistributedVector operator*(const TlDenseSymmetricMatrix_blacs& A,
                                       const TlDistributedVector& X);

  friend TlDenseSymmetricMatrix_blacs operator*(
      double X, const TlDenseSymmetricMatrix_blacs& Y);
  friend TlDenseSymmetricMatrix_blacs operator*(
      const TlDenseSymmetricMatrix_blacs& X, double Y);

  friend TlDenseGeneralMatrix_blacs multiplicationByScaLapack(
      const TlDenseGeneralMatrix_blacs& X,
      const TlDenseSymmetricMatrix_blacs& Y);
  friend TlDenseGeneralMatrix_blacs multiplicationByScaLapack(
      const TlDenseSymmetricMatrix_blacs& X,
      const TlDenseGeneralMatrix_blacs& Y);

  /// 対称行列の積を求める
  // friend TlMatrix multiplicationByLapack(const TlMatrix_Symmetric& X, const
  // TlMatrix& Y);  friend TlMatrix multiplicationByLapack(const TlMatrix& X,
  // const TlMatrix_Symmetric& Y);

  /// 対称行列の固有値を求める(QR法)
  ///
  /// @param[in] inMatrix 対称行列
  /// @param[out] outEigVal 固有値が格納されたベクトル
  /// @param[out] outEigVec 固有値ベクトルが格納された行列
  /// @retval true 固有値が求められた
  /// @retval false エラーが発生した
  friend bool diagonalByScaLapack_QR(
      const TlDenseSymmetricMatrix_blacs& inMatrix,
      TlDenseVector_Lapack* outEigVal, TlDenseGeneralMatrix_blacs* outEigVec);

  /// 対称行列の固有値を求める(Divide&Conquer)
  ///
  /// @param[in] inMatrix 対称行列
  /// @param[out] outEigVal 固有値が格納されたベクトル
  /// @param[out] outEigVec 固有値ベクトルが格納された行列
  /// @retval true 固有値が求められた
  /// @retval false エラーが発生した
  friend bool diagonalByScaLapack_DC(
      const TlDenseSymmetricMatrix_blacs& inMatrix,
      TlDenseVector_Lapack* outEigVal, TlDenseGeneralMatrix_blacs* outEigVec);

  friend TlDenseSymmetricMatrix_blacs inverseByScaLapack(const TlDenseSymmetricMatrix_blacs& inoutMatrix);
#endif  // HAVE_SCALAPACK
};

template <typename T>
void TlDenseSymmetricMatrix_blacs::print(T& out) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  const int nNumOfDim = this->getNumOfRows();  // == this->getNumOfCols()

  if (rComm.isMaster() == true) {
    out << "\n\n";
  }
  for (int ord = 0; ord < nNumOfDim; ord += 10) {
    if (rComm.isMaster() == true) {
      out << "       ";
      for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); ++j) {
        out << TlUtils::format("   %5d th", j + 1);
      }
      out << "\n"
          << " ----";

      for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); ++j) {
        out << "-----------";
      }
      out << "----\n";
    }

    for (int i = 0; i < nNumOfDim; ++i) {
      if (rComm.isMaster() == true) {
        out << TlUtils::format(" %5d  ", i + 1);
      }

      for (int j = ord; ((j < ord + 10) && (j < nNumOfDim)); ++j) {
        if (j > i) {
          if (rComm.isMaster() == true) {
            out << "    ----   ";
          }
        } else {
          const double dValue = this->get(i, j);
          if (rComm.isMaster() == true) {
            out << TlUtils::format(" %10.6lf", dValue);
          }
        }
      }
      if (rComm.isMaster() == true) {
        out << "\n";
      }
    }
    if (rComm.isMaster() == true) {
      out << "\n\n";
    }
  }

  if (rComm.isMaster() == true) {
    out.flush();
  }
}

#endif  // TL_DENSE_SYMMETRIC_MATRIX_BLACS_H
