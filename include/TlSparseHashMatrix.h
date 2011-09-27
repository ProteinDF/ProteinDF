#ifndef TLSPARSEHASHMATRIX_H
#define TLSPARSEHASHMATRIX_H

#include <cassert>
#include <ext/hash_map>
//#include "TlHashMap.h"

class TlSparseHashMatrix {
public:
  typedef __gnu_cxx::hash_map<unsigned long, double> ContainerType;
  typedef ContainerType::const_iterator const_iterator;
  typedef ContainerType::iterator iterator;

public:	 
  TlSparseHashMatrix(int numOfRows =1, int numOfCols =1);
  TlSparseHashMatrix(const TlSparseHashMatrix& rhs);
  ~TlSparseHashMatrix();

public:
  /// 最初のイテレータを返す
  const_iterator begin() const {
    return this->container.begin();
  }

  /// 最初のイテレータを返す
  iterator begin() {
    return this->container.begin();
  }

  /// 最後のイテレータを返す
  const_iterator end() const {
    return this->container.end();
  }

  /// 最後のイテレータを返す
  iterator end() {
    return this->container.end();
  }

public:
  TlSparseHashMatrix& operator=(const TlSparseHashMatrix& rhs);

  int getNumOfRows() const;
  int getNumOfCols() const;
  int getSize() const;
  void resize(int row, int col);

  double& operator()(int row, int col);

public:
  unsigned long index(const int row, const int col) const {
    const unsigned int r = static_cast<unsigned int>(row);
    const unsigned int c = static_cast<unsigned int>(col);

    return ((static_cast<unsigned long>(r) << INT_BITS) + c);
  }

  void index(const unsigned long i, int* pRow, int* pCol) const {
    assert(pRow != NULL);
    assert(pCol != NULL);
    
    const unsigned int r = static_cast<unsigned int>(i >> INT_BITS);
    const unsigned int c = static_cast<unsigned int>(i & TlSparseHashMatrix::MAX_INT);
    
    *pRow = static_cast<int>(r);
    *pCol = static_cast<int>(c);
  }

protected:
  int numOfRows;
  int numOfCols;

  //TlHashMap container;
  //__gnu_cxx::hash_map<unsigned long, int, std::tr1::hash<unsigned long> > container;
  ContainerType container;

  static const int INT_BITS;
  static const int MAX_INT;

  friend class TlCommunicate;
};

#endif // TLSPARSEHASHMATRIX_H
