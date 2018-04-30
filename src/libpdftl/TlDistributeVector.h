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

#ifndef TLDISTRIBUTEVECTOR_H
#define TLDISTRIBUTEVECTOR_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlVectorObject.h"

class TlDistributeVector : public TlVectorObject {
 protected:
  typedef std::valarray<double> DataType;

 public:
  explicit TlDistributeVector(size_type nSize = 1);
  TlDistributeVector(const TlDistributeVector& rhs);
  explicit TlDistributeVector(const std::vector<double>& rhs, int globalSize);
  explicit TlDistributeVector(const TlVector& rhs);

  virtual ~TlDistributeVector();

  static void setSystemBlockSize(int blockSize);

 public:
  virtual void resize(size_type nSize);

  TlDistributeVector& operator=(const TlDistributeVector& rhs);
  TlDistributeVector& operator=(const TlVector& rhs);

  size_type getSize() const;

  void set(size_type index, double value);
  double get(size_type index) const;
  virtual void add(const size_type index, const double value);

  double operator[](size_type index) const;
  double& operator[](size_type index);

  TlDistributeVector& operator+=(const TlDistributeVector& rhs);
  TlDistributeVector& operator-=(const TlDistributeVector& rhs);
  TlDistributeVector& operator*=(double dCoef);
  TlDistributeVector& operator/=(double dCoef);

 public:
  virtual bool load(const std::string& sFilePath);
  virtual bool save(const std::string& sFilePath) const;

  /// STLのvector<double>を返す
  /// 全ノードで呼び出される必要がある
  /// @ret 全ノードで同一のベクトルが返される
  std::vector<double> getVector() const;

 protected:
  virtual void initialize();
  virtual int getIndex(int nGlobalRow, int nGlobalCol) const;
  int getProcIdForIndex(const index_type globalIndex) const;

  virtual bool load(std::ifstream& ifs);
  virtual bool save(std::ofstream& ofs) const;

 public:
  // friend
  friend TlDistributeVector operator*(const TlDistributeMatrix& A,
                                      const TlDistributeVector& X);
  friend TlDistributeVector operator*(const TlDistributeVector& X,
                                      const TlDistributeMatrix& A);
  friend TlDistributeVector operator*(const TlDistributeSymmetricMatrix& A,
                                      const TlDistributeVector& X);
  // friend TlDistributeVector operator*(const TlDistributeVector& X, const
  // TlDistributeSymmetricMatrix& A);
  friend double operator*(const TlDistributeVector& X,
                          const TlDistributeVector& Y);
  friend TlDistributeVector operator*(const TlDistributeVector& X,
                                      const double& Y);
  friend TlDistributeVector operator*(const double& X,
                                      const TlDistributeVector& Y);
  friend TlDistributeVector operator+(const TlDistributeVector& X,
                                      const TlDistributeVector& Y);
  friend TlDistributeVector operator-(const TlDistributeVector& X,
                                      const TlDistributeVector& Y);

 protected:
  /// MPI通信タグ
  enum {
    // for copy constructor
    // TAG_CONSTRUCTOR_SIZE   = 10001,
    // TAG_CONSTRUCTOR_INDEX  = 10002,
    // TAG_CONSTRUCTOR_VALUE  = 10003,
    // for LOAD
    TAG_LOAD_SIZE = 11011,
    TAG_LOAD_ROWCOLS = 11012,
    TAG_LOAD_VALUES = 11013,
    TAG_LOAD_END = 11014
    // for save
    // TAG_SAVE_HANDSHAKE     = 10021,
    // TAG_SAVE_HANDSHAKE_OK  = 10022,
    // TAG_SAVE_DATA_ROWS     = 10023,
    // TAG_SAVE_DATA_COLS     = 10024,
    // TAG_SAVE_DATA_ROWINDEX = 10025,
    // TAG_SAVE_DATA_COLINDEX = 10026,
    // TAG_SAVE_DATA          = 10027,
  };

 protected:
  static int systemBlockSize_;
  static const std::size_t FILE_BUFFER_SIZE;

 protected:
  int m_nContext;
  int m_pDESC[9];

  int m_nRows;  // 大域行列の行数
  int m_nCols;  // 大域行列の列数

  int m_nMyRows;  // ローカル行列の行数
  int m_nMyCols;  // ローカル行列の列数

  DataType data_;

  std::vector<int>
      m_RowIndexTable;  // m_RowIndexTable[local_index] = global_index
  std::vector<int> m_ColIndexTable;

  double m_dTempVar;

 protected:
  int m_nRank;         // プロセスのランク
  int m_nProc;         // 総プロセス数
  int m_nProcGridRow;  // プロセスグリッドの行数
  int m_nProcGridCol;  // プロセスグリッドの列数
  int m_nMyProcRow;    // プロセスグリッドにおける自分の行数
  int m_nMyProcCol;    // プロセスグリッドにおける自分の列数
  int m_nBlockSize;    // ブロックサイズ

 protected:
  std::string filePath_;

  friend class TlDistributeMatrix;
  friend class TlDistributeSymmetricMatrix;
};

#endif  // TLDISTRIBUTEVECTOR_H
