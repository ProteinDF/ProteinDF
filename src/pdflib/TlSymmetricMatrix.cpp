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

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "TlSymmetricMatrix.h"
#include "TlUtils.h"
#include "TlLogging.h"

////////////////////////////////////////////////////////////////////////
//
TlSymmetricMatrix::TlSymmetricMatrix(int dim)
    : TlMatrix(dim, dim, NULL)
{
    assert(dim >= 0);

    this->initialize();
}


TlSymmetricMatrix::TlSymmetricMatrix(const TlSymmetricMatrix& rhs)
        : TlMatrix(rhs.getNumOfRows(), rhs.getNumOfCols(), NULL)
{
    this->initialize(false);
    std::copy(rhs.data_, rhs.data_ + rhs.getNumOfElements(), this->data_);
}


TlSymmetricMatrix::TlSymmetricMatrix(const TlMatrix& rhs)
    : TlMatrix(rhs.getNumOfRows(), rhs.getNumOfCols(), NULL)
{
    assert(rhs.getNumOfRows() == rhs.getNumOfCols());

    this->initialize(false);

    //std::cerr << "Full matrix is converted to Symmetric matrix..." << std::endl;
    const int dim = this->getNumOfRows();
#pragma omp parallel for
    for (int row = 0; row < dim; ++row) {
        for (int col = 0; col <= row; ++col) {
            this->set(row, col, rhs.get(row, col));
        }
    }
}


TlSymmetricMatrix::TlSymmetricMatrix(const TlSerializeData& data)
    : TlMatrix(1, 1, NULL)
{
    this->m_nRows = std::max(data["row"].getInt(), 1);
    this->m_nCols = std::max(data["col"].getInt(), 1);
    assert(this->getNumOfRows() == this->getNumOfCols());
    this->initialize(false);

    size_type index = 0;
    const index_type maxRow = this->getNumOfRows();
    for (index_type row = 0; row < maxRow; ++row) {
        for (index_type col = 0; col <= row; ++col) {
            this->set(row, col, data["data"].getAt(index).getDouble());
            ++index;
        }
    }
}


TlSymmetricMatrix::TlSymmetricMatrix(const TlVector& vct, index_type dim)
    : TlMatrix(dim, dim, NULL)
{
    assert(0 <= dim);

    const size_type size = this->getNumOfElements();
    assert(vct.getSize() == size);

    this->initialize(false);

#ifndef __FUJITSU // cannot use OpenMP in constructor
#pragma omp parallel for
#endif // __FUJITSU 
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] = vct[index];
    }
}


TlSymmetricMatrix::~TlSymmetricMatrix()
{
    // バッファの削除は親クラスで行う。
    //this->clear();
}


void TlSymmetricMatrix::resize(const index_type dim)
{
    assert(dim > 0);

    TlSymmetricMatrix oldMatrix(*this);

    this->clear();
    this->m_nRows = dim;
    this->m_nCols = dim;
    this->initialize(true);

    const index_type nMaxDimForCopy = std::min<index_type>(oldMatrix.getNumOfRows(), dim);

#pragma omp parallel for
    for (index_type i = 0; i < nMaxDimForCopy; ++i) {
        for (index_type j = 0; j <= i; ++j) {
            this->set(i, j, oldMatrix.get(i, j));
        }
    }
}


std::size_t TlSymmetricMatrix::index(index_type row,
                                     index_type col) const
{
    assert((0 <= row) && (row < this->m_nRows));
    assert((0 <= col) && (col < this->m_nCols));
    
    // This class treats 'U' type matrix.
    // if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; in Fortran
    // index = row + col * (col +1) /2; // in C/C++
    
    if (row > col) {
        std::swap(row, col);
    }
    unsigned int t = col * (col +1);
    t = t >> 1; // == 't /= 2'
    return row + t;
}

std::size_t TlSymmetricMatrix::vtr_index(index_type row,
                                         index_type col,
                                         index_type dim) 
{
    assert((0 <= row) && (row < dim));
    assert((0 <= col) && (col < dim));

    // This class treats 'U' type matrix.
    // if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; in Fortran
    // index = row + col * (col +1) /2; // in C/C++
    
    if (row > col) {
        std::swap(row, col);
    }
    unsigned int t = col * (col +1);
    t = t >> 1; // == 't /= 2'
    return row + t;
}


double TlSymmetricMatrix::sum() const
{
    double answer = std::accumulate(this->data_,
                                    this->data_ + this->getNumOfElements(),
                                    0.0);
    answer *= 2.0;
    answer -= this->trace();

    return answer;
}


TlSymmetricMatrix& TlSymmetricMatrix::operator=(const TlSymmetricMatrix& rhs)
{
    if (&rhs != this) {
        this->clear();
        this->m_nRows = rhs.getNumOfRows();
        this->m_nCols = rhs.getNumOfCols();
        this->initialize(false);

        const std::size_t size = this->getNumOfElements();
        std::copy(rhs.data_, rhs.data_ + size, this->data_);
    }

    return *this;
}


TlSymmetricMatrix& TlSymmetricMatrix::operator+=(const TlSymmetricMatrix& rhs)
{
    assert(this->m_nRows == rhs.m_nRows);
    assert(this->m_nCols == rhs.m_nCols);

    const size_type size = this->getNumOfElements();

#ifdef _OPENMP
    const size_type quot = size / MAX_LOOP;
    const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
    {
        for (size_type block = 0; block < quot; ++block) {
            const size_type index_base = block * MAX_LOOP;
#pragma omp for
            for (int i = 0; i < MAX_LOOP; ++i) {
                const size_type index = index_base + i;
                this->data_[index] += rhs.data_[index];
            }
        }

        const size_type index_base = quot * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < rem; ++i) {
            const size_type index = index_base + i;
            this->data_[index] += rhs.data_[index];
        }
    }
#else
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] += rhs.data_[index];
    }
#endif // _OPENMP
    return (*this);
}


// TlSymmetricMatrix& TlSymmetricMatrix::operator+=(const TlSparseSymmetricMatrix& rhs)
// {
//     std::cerr << "call TlSymmetricMatrix::operator+=(const TlSparseSymmetricMatrix&)" << std::endl;
//     assert(this->getNumOfRows() == rhs.getNumOfRows());
//     assert(this->getNumOfCols() == rhs.getNumOfCols());

//     index_type row = 0;
//     index_type col = 0;
//     TlSparseSymmetricMatrix::const_iterator itEnd = rhs.end();
//     for (TlSparseSymmetricMatrix::const_iterator it = rhs.begin(); it != itEnd; ++it) {
//         rhs.index(it->first, &row, &col);
//         this->add(row, col, it->second);
//     }

//     return *this;
// }


TlSymmetricMatrix& TlSymmetricMatrix::operator-=(const TlSymmetricMatrix& rhs)
{
    assert(this->m_nRows == rhs.m_nRows);
    assert(this->m_nCols == rhs.m_nCols);

    const size_type size = this->getNumOfElements();

#ifdef _OPENMP
    const size_type quot = size / MAX_LOOP;
    const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
    {
        for (size_type block = 0; block < quot; ++block) {
            const size_type index_base = block * MAX_LOOP;
#pragma omp for
            for (int i = 0; i < MAX_LOOP; ++i) {
                const size_type index = index_base + i;
                this->data_[index] -= rhs.data_[index];
            }
        }

        const size_type index_base = quot * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < rem; ++i) {
            const size_type index = index_base + i;
            this->data_[index] -= rhs.data_[index];
        }
    }
#else
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] -= rhs.data_[index];
    }
#endif // _OPENMP
    
    return (*this);
}


TlSymmetricMatrix& TlSymmetricMatrix::operator*=(const double& rhs)
{
    const size_type size = this->getNumOfElements();

#ifdef _OPENMP
    const size_type quot = size / MAX_LOOP;
    const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
    {
        for (size_type block = 0; block < quot; ++block) {
            const size_type index_base = block * MAX_LOOP;
#pragma omp for
            for (int i = 0; i < MAX_LOOP; ++i) {
                const size_type index = index_base + i;
                this->data_[index] *= rhs;
            }
        }

        const size_type index_base = quot * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < rem; ++i) {
            const size_type index = index_base + i;
            this->data_[index] *= rhs;
        }
    }
#else
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] *= rhs;
    }
#endif // _OPENMP
    
    return (*this);
}


TlSymmetricMatrix& TlSymmetricMatrix::operator/=(const double& rhs)
{
    return (this->operator*=(1.0 / rhs));
}


TlSymmetricMatrix operator+(const TlSymmetricMatrix& X, const TlSymmetricMatrix& Y)
{
    assert(X.m_nRows == Y.m_nCols);
    assert(X.m_nCols == Y.m_nCols);

    TlSymmetricMatrix answer(X);
    answer += Y;

    return answer;
}


TlMatrix operator+(const TlMatrix& X, const TlSymmetricMatrix& Y)
{
    assert(X.getNumOfRows() == Y.getNumOfRows());
    assert(X.getNumOfCols() == Y.getNumOfCols());

    TlMatrix answer(X);
    const int numOfRows = Y.getNumOfRows();
#pragma omp parallel for
    for (int row = 0; row < numOfRows; ++row) {
        // col != row
        for (int col = 0; col < row; ++col) {
            const double tmp = Y.get(row, col);
            answer.add(row, col, tmp);
            answer.add(col, row, tmp);
        }
        // col == row
        answer.add(row, row, Y.get(row, row));
    }

    return answer;
}


TlSymmetricMatrix operator-(const TlSymmetricMatrix& X, const TlSymmetricMatrix& Y)
{
    assert(X.getNumOfRows() == Y.getNumOfRows());
    assert(X.getNumOfCols() == Y.getNumOfCols());

    TlSymmetricMatrix answer(X);
    answer -= Y;

    return answer;
}


TlMatrix operator-(const TlMatrix& X, const TlSymmetricMatrix& Y)
{
    return (X + (-1.0 * Y));
}


TlMatrix operator-(const TlSymmetricMatrix& X, const TlMatrix& Y)
{
    return ((-1.0 * Y) + X);
}


// X x Y
TlMatrix operator*(const TlSymmetricMatrix& X, const TlSymmetricMatrix& Y)
{
    assert(X.getNumOfCols() == Y.getNumOfRows());

    TlMatrix answer(X.m_nRows, Y.m_nCols);

#ifdef HAVE_LAPACK
    // use LAPACK
    TlMatrix tmpY = Y;
    answer = multiplicationByLapack(TlMatrix(X), tmpY);
#else
    // not use LAPACK
    {
        for (int row = 0; row < X.m_nRows; row++) {
            for (int col = 0; col < Y.m_nCols; col++) {
                for (int t = 0; t < X.m_nCols; t++) {
                    answer(row, col) += X(row, t) * Y(t, col);
                }
            }
        }
    }
#endif // HAVE_LAPACK

    return answer;
}


TlMatrix operator*(const TlMatrix& X, const TlSymmetricMatrix& Y)
{
    assert(X.getNumOfCols() == Y.getNumOfRows());
    TlMatrix answer(X.getNumOfRows(), Y.getNumOfCols());

#ifdef HAVE_LAPACK
    answer = multiplicationByLapack(X, Y);
#else
    // not use LAPACK
    const int nMaxNumOfRows = X.getNumOfRows();
    const int nMaxNumOfCols = Y.getNumOfCols();
    const int nMaxNumOfTemp = X.getNumOfCols();
    for (int row = 0; row < nMaxNumOfRows; ++row) {
        for (int col = 0; col < nMaxNumOfCols; ++col) {
            for (int t = 0; t < nMaxNumOfTemp; ++t) {
                answer(row, col) += X(row, t) * Y(t, col);
            }
        }
    }
#endif // HAVE_LAPACK

    return answer;
}


TlMatrix operator*(const TlSymmetricMatrix& X, const TlMatrix& Y)
{
    assert(X.getNumOfCols() == Y.getNumOfRows());
    TlMatrix answer(X.getNumOfRows(), Y.getNumOfCols());

#ifdef HAVE_LAPACK
    answer = multiplicationByLapack(X, Y);
#else
    // not use LAPACK
    for (int row = 0; row < X.getNumOfRows(); row++) {
        for (int col = 0; col < Y.getNumOfCols(); col++) {
            for (int t = 0; t < X.getNumOfCols(); t++) {
                answer(row, col) += X(row, t) * Y(t, col);
            }
        }
    }
#endif // HAVE_LAPACK

    return answer;
}


TlVector operator*(const TlSymmetricMatrix& A, const TlVector& X)
{
    assert(A.getNumOfCols() == X.getSize());
    TlVector answer(A.getNumOfRows());

#ifdef HAVE_LAPACK
    answer = multiplicationByLapack(A, X);
#else
    // not use LAPACK
    for (int row = 0; row < X.getNumOfRows(); ++row) {
        double tmp = 0.0;
        for (int col = 0; col < Y.getNumOfCols(); ++col) {
            tmp += A.get(row, col) * X.get(col);
        }
        answer.set(row, tmp);
    }
#endif // HAVE_LAPACK

    return answer;
}


TlVector operator*(const TlVector& X, const TlSymmetricMatrix& A)
{
    return (A * X);
}


TlSymmetricMatrix operator*(const TlSymmetricMatrix& X, double Y)
{
    TlSymmetricMatrix answer(X);
    answer *= Y;

    return answer;
}


TlSymmetricMatrix operator*(double X, const TlSymmetricMatrix& Y)
{
    return (Y * X);
}


const TlSymmetricMatrix& TlSymmetricMatrix::dot(const TlSymmetricMatrix& X)
{
    assert(X.getNumOfRows() == this->getNumOfRows());
    assert(X.getNumOfCols() == this->getNumOfCols());

    const std::size_t size = this->getNumOfElements();

#ifdef _OPENMP
    // use OpenMP
    const size_type quot = size / MAX_LOOP;
    const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
    {
        for (size_type block = 0; block < quot; ++block) {
            const size_type index_base = block * MAX_LOOP;
#pragma omp for
            for (int i = 0; i < MAX_LOOP; ++i) {
                const size_type index = index_base + i;
                this->data_[index] *= X.data_[index];
            }
        }

        const size_type index_base = quot * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < rem; ++i) {
            const size_type index = index_base + i;
            this->data_[index] *= X.data_[index];
        }
    }
#else
    // not use OpenMP
    for (std::size_t i = 0; i < size; ++i) {
        this->data_[i] *= X.data_[i];
    }
#endif // _OPENMP
    
    return (*this);
}


TlSymmetricMatrix dot(const TlSymmetricMatrix& X, const TlSymmetricMatrix& Y)
{
    assert(X.getNumOfRows() == Y.getNumOfRows());
    assert(X.getNumOfCols() == Y.getNumOfCols());

    TlSymmetricMatrix Z = X;
    Z.dot(Y);

    return Z;
}

/**
   入力されたstd::ifstream が読み込み可能なフォーマットかチェックする。
   チェック後、渡されたifstream の読み込み位置は先頭に戻される。
 */
bool TlSymmetricMatrix::isLoadable(std::ifstream& ifs)
{
    if (ifs.is_open() != true) {
        return false;
    }

    enum {RSFD=0, CSFD, RLHD, CLHD, RUHD, CUHD, RSFS, CSFS, RLHS, CLHS, RUHS, CUHS, LHSC};
    int nType = 0;

    ifs.seekg(0, std::ios::beg);
    const bool bHeader = TlSymmetricMatrix::getHeaderInfo(ifs, &nType);
    ifs.seekg(0, std::ios::beg);

    bool bJudge = false;
    if (bHeader == true) {
        if (nType == RLHD) {
            bJudge =  true;
        }
    }

    return bJudge;
}


bool TlSymmetricMatrix::isLoadable(const std::string& rFilePath)
{
    std::ifstream ifs;
    ifs.open(rFilePath.c_str());
    if (ifs.fail()) {
#ifdef DEBUG
        std::cerr << "[error] TlMatrix::load(): could not open file. " << rFilePath << std::endl;
#endif // DEBUG
        return false;
    }

    const bool bAnswer = TlSymmetricMatrix::isLoadable(ifs);
    ifs.close();
    return bAnswer;
}


bool TlSymmetricMatrix::getHeaderInfo(std::ifstream& ifs, int* pType,
                                      int* pNumOfRows, int* pNumOfCols)
{
    bool bAnswer = true;

    // get file size
    std::ifstream::pos_type nFileSize = 0;
    {
        ifs.seekg(0, std::ios_base::beg);
        std::ifstream::pos_type begin = ifs.tellg();
        ifs.seekg(0, std::ios_base::end);
        std::ifstream::pos_type end = ifs.tellg();

        nFileSize = end - begin;
    }

    int  nType = 0;
    bool bCheckVariableType = false;
    std::ifstream::pos_type nStartContentPos = 0;

    // int case:
    {
        int nRows = 0;
        int nCols = 0;
        ifs.seekg(0, std::ios_base::beg);
        ifs.read((char*)&nType, sizeof(int));
        ifs.read((char*)&nRows, sizeof(int));
        ifs.read((char*)&nCols, sizeof(int));

        const std::ifstream::pos_type nEstimatedFileSize =
            std::ifstream::pos_type(sizeof(int) *3)
            + (std::ifstream::pos_type(nRows) * std::ifstream::pos_type(nCols +1)
               * std::ifstream::pos_type(sizeof(double) / 2));
        if ((nRows == nCols) && (nEstimatedFileSize == nFileSize)) {
            if (pType != NULL) {
                *pType = nType;
            }
            if (pNumOfRows != NULL) {
                *pNumOfRows = nRows;
            }
            if (pNumOfCols != NULL) {
                *pNumOfCols = nCols;
            }
            bCheckVariableType = true;
            nStartContentPos = sizeof(int) * 3;
        }
    }

    // long case:
    {
        long nRows = 0;
        long nCols = 0;
        ifs.seekg(0, std::ios_base::beg);
        ifs.read((char*)&nType, sizeof(int));
        ifs.read((char*)&nRows, sizeof(long));
        ifs.read((char*)&nCols, sizeof(long));

        const std::ifstream::pos_type nEstimatedFileSize = (sizeof(int) + sizeof(long) * 2) + (nRows * (nCols +1) * sizeof(double) / 2);
        if ((nRows == nCols) && (nEstimatedFileSize == nFileSize)) {
            if (pType != NULL) {
                *pType = nType;
            }
            if (pNumOfRows != NULL) {
                *pNumOfRows = static_cast<int>(nRows);
            }
            if (pNumOfCols != NULL) {
                *pNumOfCols = static_cast<int>(nCols);
            }
            bCheckVariableType = true;
            nStartContentPos = sizeof(int) + sizeof(long) * 2;
        }
    }

    if (bCheckVariableType == true) {
        ifs.seekg(nStartContentPos, std::ios_base::beg);
    } else {
        //std::cerr << "file size mismatch." << std::endl;
        bAnswer = false;
    }

    return bAnswer;
}

bool TlSymmetricMatrix::load(const std::string& sFilePath)
{
    std::ifstream ifs;

    ifs.open(sFilePath.c_str(), std::ifstream::binary | std::ifstream::in);
    if (ifs.fail()) {
#ifdef DEBUG
        std::cerr << "[error] TlSymmetricMatrix::load(): could not open file. " << sFilePath << std::endl;
#endif // DEBUG
        return false;
    }

    bool bAnswer = this->load(ifs);
    ifs.close();

    if (bAnswer == false) {
        TlLogging& log = TlLogging::getInstance();
        log.critical(TlUtils::format("TlSymmetricMatrix::load() is not supported: %s",
                                     sFilePath.c_str()));
        std::abort();

        return false;
    }

    return true;
}

bool TlSymmetricMatrix::load(std::ifstream& ifs)
{
    bool bAnswer = true;

    // binary mode
    enum {RSFD, CSFD, RLHD, CLHD, RUHD, CUHD, RSFS, CSFS, RLHS, CLHS, RUHS, CUHS, LHSC};

    // read header
    int nType = 0;
    int nRows = 0;
    int nCols = 0;

    bAnswer = TlSymmetricMatrix::getHeaderInfo(ifs, &nType, &nRows, &nCols);

    switch (nType) {
    case RLHD:
        //std::cerr << "load RLHD" << std::endl;
        break;

    default:
        std::cerr << "unknown matrix type(" << nType << ")." << std::endl;
        bAnswer = false;
        break;
    }

    // const int nMatrixSize = this->m_nRows + (this->m_nRows * (this->m_nRows -1)) / 2;
    // ifs.read(reinterpret_cast<char*>(&(this->m_aMatrix[0])), nMatrixSize);
    if (bAnswer == true) {
        assert(nRows == nCols);
        this->clear();
        this->m_nRows = nRows;
        this->m_nCols = nCols;
        this->initialize();

        const int numOfRows = this->getNumOfRows();
        for (int row = 0; row < numOfRows; ++row) {
            for (int col = 0; col <= row; ++col) {
                double tmp;
                ifs.read((char*)&tmp, sizeof(double));
                (*this)(row, col) = tmp;
            }
        }
    }

    return bAnswer;
}

bool TlSymmetricMatrix::save(const std::string& sFilePath) const
{
    bool bAnswer = true;

    std::ofstream ofs;
    ofs.open(sFilePath.c_str(), std::ofstream::out | std::ofstream::binary);
    bAnswer = this->save(ofs);
    ofs.close();

//   std::cerr << "TlSymmetricMatrix::save() sFilePath = " << sFilePath << std::endl;

    return bAnswer;
}

bool TlSymmetricMatrix::save(std::ofstream& ofs) const
{
    bool bAnswer = true;

    const int nType = 2; // means RLHD
    ofs.write(reinterpret_cast<const char*>(&nType), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&this->m_nRows), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&this->m_nCols), sizeof(int));

    assert(this->m_nRows == this->m_nCols);
// const int nMatrixSize = this->m_nRows + (this->m_nRows * (this->m_nRows -1)) / 2;
// ofs.write(reinterpret_cast<const char*>(&(this->m_aMatrix[0])), nMatrixSize);
    const std::size_t maxRow = this->m_nRows;
    for (std::size_t row = 0; row < maxRow; ++row) {
        for (std::size_t col = 0; col <= row; ++col) {
            const double tmp = (*this)(row, col);
            ofs.write(reinterpret_cast<const char*>(&tmp), sizeof(double)) ;
        }
    }

    return bAnswer;
}


TlSerializeData TlSymmetricMatrix::getSerialize() const
{
    TlSerializeData data;
    data["row"] = this->getNumOfRows();
    data["col"] = this->getNumOfCols();
    data["type"] = "RLHD";

    TlSerializeData tmp;
    const size_type size = this->getNumOfElements();
    for (size_type index = 0; index < size; ++index) {
        tmp.pushBack(this->data_[index]);
    }
    data["data"] = tmp;

    return data;
}


/**
 *  出力ストリーム
 */
std::ostream& operator <<(std::ostream& out, const TlSymmetricMatrix& rhs)
{
    rhs.print(out);
    return out;
}

bool TlSymmetricMatrix::diagonal(TlVector* pEigVal, TlMatrix* pEigVec) const
{
#ifdef HAVE_LAPACK
    return diagonalByLapack(*this, pEigVal, pEigVec);
#else
    std::cerr << "please WRITE diagnal argorithm. stop." << std::endl;
    abort();
#endif // HAVE_LAPACK
}

bool TlSymmetricMatrix::inverse()
{
#ifdef HAVE_LAPACK
    // using LAPACK
    return inverseByLapack(*this);
#else
    // without LAPACK
    std::cerr << "sorry. this code is not implemented." << std::endl;
    abort();
    return false;
#endif // HAVE_LAPACK  
}

TlSymmetricMatrix& TlSymmetricMatrix::transpose()
{
    return (*this);
}

TlVector TlSymmetricMatrix::getMaxAbsoluteVectorOnEachRow() const
{
    const int nSize = this->getNumOfRows();

    TlVector v(nSize);

    for (int row = 0; row < nSize; ++row) {
        // たいていの場合、対角要素が大きいので、これを最初にもってくる
        double dMax = std::fabs((*this)(row, row));

        for (int col = (row -1); col >= 0; --col) {
            dMax = std::max(dMax, std::fabs((*this)(row, col)));
        }

        v[row] = dMax;
    }

    return v;
}

// 要素の絶対値の最大値を返す
// int* outRow, outCol にはその要素がある行と列を返す
double TlSymmetricMatrix::getMaxAbsoluteElement(int* outRow, int* outCol) const
{
    //std::cout << "TlSymmetricMatrix::getMaxAbsoluteElement()" << std::endl;
    double dAnswer = -1.0;
    int nMaxRow = 0;
    int nMaxCol = 0;

    for (int row = 0; row < this->m_nRows; row++) {
        for (int col = 0; col <= row; col++) {
            double val = fabs((*this)(row, col));
            if (dAnswer < val) {
                dAnswer = val;
                nMaxRow = row;
                nMaxCol = col;
            }
        }
    }

    if (outRow != NULL) {
        *outRow = nMaxRow;
    }
    if (outCol != NULL) {
        *outCol = nMaxCol;
    }

    return dAnswer;
}

double TlSymmetricMatrix::getMaxAbsoluteElementByIndex(int index) const
{
    double dAnswer = 0.0;

    for (int row = 0; row < this->m_nRows; ++row) {
        double val = fabs((*this)(row, index));
        dAnswer = std::max(dAnswer, val);
    }

    return dAnswer;
}

////////////////////////////////////////////////////////////////////////
// NOTICE:
// friend functions
#ifdef HAVE_LAPACK
// for LAPACK
extern "C" {
    // dspmv
    // DSPMV  performs the matrix-vector operation
    //  y := alpha*A*x + beta*y,
    // where alpha and beta are scalars, x and y are n element vectors and
    // A is an n by n symmetric matrix, supplied in packed form.
    void dspmv_(const char* UPLO, const int* N, const double* ALPHA,
                const double* AP, const double* X,
                const int* INCX, double* BETA, double* Y, const int* INCY);

    // performs one of the matrix-matrix operations
    //    C := alpha*A*B + beta*C,
    // or
    //    C := alpha*B*A + beta*C,
    // where alpha and beta are scalars,  A is a symmetric matrix and  B and
    // C are  m by n matrices.
    void dsymm_(const char* SIDE, const char* UPLO, const int* M, const int* N,
                const double* ALPHA, const double* A, const int* LDA,
                const double* B, const int* LDB, const double* BETA ,
                double* C, const int* LDC);

    // DSPEV
    void dspev_(const char* JOBZ, const char* UPLO,
                const int* N, double* AP, double* W,
                double* Z, const int* LDZ, double* WORK,
                int* INFO);
    // DSPEVD
    // void dspevd_(const char* JOBZ, const char* UPLO,
    //              const int* N, double* AP, double* W,
    //              double* Z, const int* LDZ, double* WORK, int* LWORK,
    //              int* IWORK, int* LIWORK, int* INFO);

    // DSPTRF computes the factorization of a real symmetric matrix A stored
    // in packed format using the Bunch-Kaufman diagonal pivoting method:
    //    A = U*D*U**T  or  A = L*D*L**T
    // where U (or L) is a product of permutation and unit upper (lower)
    // triangular matrices, and D is symmetric and block diagonal with
    // 1-by-1 and 2-by-2 diagonal blocks.
    void dsptrf_(const char* UPLO, const int* N, double* PA, int* IPIV, int* INFO);

    // DSPTRI computes the inverse of a real symmetric indefinite matrix
    // A in packed storage using the factorization A = U*D*U**T or
    // A = L*D*L**T computed by DSPTRF.
    void dsptri_(const char* UPLO, const int* N, double* PA, int* IPIV, double* WORK, int* INFO);
}


TlMatrix multiplicationByLapack(const TlSymmetricMatrix& A, const TlMatrix& B)
{
    assert(A.getNumOfCols() == B.getNumOfRows());

    TlMatrix C(A.getNumOfRows(), B.getNumOfCols());

    // use LAPACK
    const char SIDE = 'L';                 // L means "C := alpha*A*B + beta*C",
                                           // R means "C := alpha*B*A + beta*C"
    const char UPLO = 'U';                 // L means the lower triangular part of the symmetric matrix
                                           // U means the upper triangular part of the symmetric matrix
    const int M = C.getNumOfRows();        // the number of rows of the matrix  C
    const int N = C.getNumOfCols();        // the number of columns of the matrix C
    const double ALPHA = 1.0;              // ALPHA specifies the scalar alpha
    const int LDA = M;                     // When  SIDE = 'L' or 'l'  then LDA must be at least max(1, M),
                                           // otherwise  LDA must be at least  max(1, N).

    double* pA = new double[M * M];
    for (int c = 0; c < M; ++c) {
        const int base = c * M;
        for (int r = 0; r <= c;  ++r) {
            pA[base + r] = A.get(r, c);
        }
    }

    const double* pB = const_cast<TlMatrix&>(B).data_; // DIMENSION (LDB, N)
    const int LDB = M;                     // max(1, M)
    const double BETA = 0.0;               // BETA  specifies the scalar  beta
    double* pC = C.data_;                  // DIMENSION (LDC, N)
    const int LDC = M;                     // max(1, M)

    dsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, pA, &LDA, pB, &LDB, &BETA, pC, &LDC);

    delete[] pA;
    pA = NULL;

    return C;
}


TlMatrix multiplicationByLapack(const TlMatrix& B, const TlSymmetricMatrix& A)
{
    assert(B.getNumOfCols() == A.getNumOfRows());

    TlMatrix C(B.getNumOfRows(), A.getNumOfCols());

    // use LAPACK
    const char SIDE = 'R';                 // L means "C := alpha*A*B + beta*C",
                                           // R means "C := alpha*B*A + beta*C"
    const char UPLO = 'U';                 // L means the lower triangular part of the symmetric matrix
                                           // U means the upper triangular part of the symmetric matrix
    const int M = C.getNumOfRows();        // the number of rows of the matrix  C
    const int N = C.getNumOfCols();        // the number of columns of the matrix C
    const double ALPHA = 1.0;              // ALPHA specifies the scalar alpha
    const int LDA = N;                     // When  SIDE = 'L' or 'l'  then LDA must be at least max(1, M),
                                           // otherwise  LDA must be at least  max(1, N).
    double* pA = new double[N * N];
    for (int c = 0; c < N; ++c) {
        const int base = c * N;
        for (int r = 0; r <= c;  ++r) {
            pA[base + r] = A.get(r, c);
        }
    }

    const double* pB = const_cast<TlMatrix&>(B).data_; // DIMENSION (LDB, N)
    const int LDB = M;                     // max(1, M)
    const double BETA = 0.0;               // BETA  specifies the scalar  beta
    double* pC = C.data_;                  // DIMENSION (LDC, N)
    const int LDC = M;                     // max(1, M)

    dsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, pA, &LDA, pB, &LDB, &BETA, pC, &LDC);

    delete[] pA;
    pA = NULL;

    return C;
}


TlVector multiplicationByLapack(const TlSymmetricMatrix& A, const TlVector& X)
{
    assert(A.getNumOfCols() == X.getSize());

    // use LAPACK
    const char UPLO = 'U';                 // L means the lower triangular part of the symmetric matrix
                                           // U means the upper triangular part of the symmetric matrix
    const int N = A.getNumOfRows();
    const double ALPHA = 1.0;              // ALPHA specifies the scalar alpha
    const double* AP = A.data_;
    const double* pX = X.data_;
    const int INCX = 1;
    double BETA = 1.0;
    TlVector answer(N);
    double* Y = answer.data_;
    const int INCY = 1;

    dspmv_(&UPLO, &N, &ALPHA, AP, pX, &INCX,
           &BETA, Y, &INCY);
    
    return answer;
}


bool diagonalByLapack(const TlSymmetricMatrix& inMatrix,
                      TlVector* outEigVal, TlMatrix* outEigVec)
{
    assert(outEigVal != NULL);
    assert(outEigVec != NULL);
    assert(inMatrix.getNumOfRows() == inMatrix.getNumOfCols());
    assert(inMatrix.getNumOfRows() >= 2);

    bool bAnswer = false;

    const char JOBZ = 'V';                         // 固有値と固有ベクトルを計算する。
    const char UPLO = 'U';
    const int N = inMatrix.getNumOfRows();     // 行列Aの次数(N>=0)
    double* AP = new double[N * (N+1) / 2];
    std::copy(inMatrix.data_, inMatrix.data_ + (N * (N+1) / 2), AP);

    outEigVal->resize(N);
    double* W = outEigVal->data_;                  // (出力用) INFO=0のとき, Wに固有値が昇順で入る。大きさN

    const int LDZ = N;
    outEigVec->resize(LDZ, N);
    double* Z = outEigVec->data_;

    double* WORK = new double[3 * N];              // (作業/出力用)
    int INFO =0;                                   // (出力用) =0: 正常終了, >0: 収束しなかった

    dspev_(&JOBZ, &UPLO, &N, AP, W, Z, &LDZ, WORK, &INFO);

    delete[] WORK;
    WORK = NULL;
    delete[] AP;
    AP = NULL;

    if (INFO == 0) {
        bAnswer = true;
    } else {
#ifdef DEBUG
        std::cerr << "[error] INFO != 0 at TlSymmetricMatrix::diagonalByLapack(). abort." << std::endl;
#endif //DEBUG
        abort();
    }

    return bAnswer;
}


bool inverseByLapack(TlSymmetricMatrix& X)
{
    // (input)
    // 'U':  Upper triangle of A is stored
    // 'L':  Lower triangle of A is stored.
    char UPLO = 'U';

    //(input) The order of the matrix A.  N >= 0.
    const int N = X.getNumOfRows();

    // (input/output) The order of the matrix A.  N >= 0.
    // if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
    // if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
    // On exit, if INFO = 0, the triangular factor U or L from the
    // Cholesky factorization A = U**T*U or A = L*L**T, in the same
    // storage format as A.
    double* AP = X.data_;

    // (output) for dsptrf_
    // INTEGER array, dimension (N)
    // Details of the interchanges and the block structure of D.
    // If IPIV(k) > 0, then rows and columns k and IPIV(k) were
    // interchanged and D(k,k) is a 1-by-1 diagonal block.
    // If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
    // columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
    // is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
    // IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
    // interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
    //
    // (input)  for dsptri_
    // INTEGER array, dimension (N)
    // Details of the interchanges and the block structure of D
    // as determined by DSPTRF.
    int* IPIV = new int[N];

    // (workspace) DOUBLE PRECISION array, dimension (N)
    double* WORK = new double[N];

    // (output)
    // = 0:  successful exit
    // < 0:  if INFO = -i, the i-th argument had an illegal value
    // > 0:  if INFO = i, the leading minor of order i is not
    // positive definite, and the factorization could not be completed.
    int INFO = 0;

    // execute LAPACK
    bool bAnswer = false;
    dsptrf_(&UPLO, &N, AP, IPIV, &INFO);
    if (INFO == 0) {
        dsptri_(&UPLO, &N, AP, IPIV, WORK, &INFO);
        if (INFO == 0) {
            bAnswer = true;
        } else {
            std::cerr << "inverseByLapack() failed.: dsptri() return code = " << INFO << std::endl;
            abort();
        }
    } else {
        std::cerr << "inverseByLapack() failed.: dsptrf() return code = " << INFO << std::endl;
        abort();
    }

    // finalize
    delete[] IPIV;
    IPIV = NULL;
    delete[] WORK;
    WORK = NULL;

    return bAnswer;
}

// int choleskyFactorization(TlSymmetricMatrix* A,
//                           std::vector<int>* pPivot)
// {
//     // *  UPLO    (input) CHARACTER*1
//     // *          = 'U':  Upper triangle of A is stored;
//     // *          = 'L':  Lower triangle of A is stored.
//     char UPLO = 'L';

//     // *  N       (input) INTEGER
//     // *          The order of the matrix A.  N >= 0.
//     const int N = A->getNumOfRows();

//     // (input/output) The order of the matrix A.  N >= 0.
//     // if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//     // if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
//     // On exit, if INFO = 0, the triangular factor U or L from the
//     // Cholesky factorization A = U**T*U or A = L*L**T, in the same
//     // storage format as A.
//     double* AP = A->data_;

//     // INTEGER array, dimension (N)
//     // Details of the interchanges and the block structure of D
//     // as determined by DSPTRF.
//     //int* IPIV = new int[N];
//     pPivot->resize(N);
//     int* IPIV = &((*pPivot)[0]);

//     int INFO = 0;

//     dsptrf_(&UPLO, &N, AP, IPIV, &INFO);

//     //bool answer = (INFO == 0);
//     // if (INFO < 0) {
//     //     std::cerr << (-INFO) << "th argument had an illegal value." << std::endl;
//     // } else if (INFO > 0) {
//     //     std::cerr << "the leading mirror of order " << INFO << " is not positive definite."
//     //               << std::endl;
//     // }

//     return INFO;
// }


TlMatrix TlSymmetricMatrix::choleskyFactorization()
{
    const index_type dim = this->getNumOfRows();
    TlMatrix L(dim, dim);
    for (index_type j = 0; j < dim; ++j) {
        // calc L_jj
        double s = this->get(j, j);
        for (index_type k = 0; k < j; ++k) {
            s -= L.get(j, k) * L.get(j, k);
        }
        if (s < 0.0) {
            std::cerr << "CholeskyFactorization() s < 0" << std::endl;
            abort();
        }
        L.set(j, j, std::sqrt(s));
        
        // calc L_ij (i > j)
        const double L_jj = L.get(j, j);
        for (index_type i = j +1; i < dim; ++i) {
            double s = this->get(i, j);
            for (index_type k = 0; k < j; ++k) {
                s -= L.get(i, k) * L.get(j, k);
            }
            L.set(i, j, s / L_jj);
        }
    }

    return L;
}

// Harbrecht, Peter, Schneider, 2011
TlMatrix TlSymmetricMatrix::choleskyFactorization2(const double threshold) const
{
    const index_type N = this->getNumOfRows();
    TlVector d = this->getDiagonalElements();
    double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    TlMatrix L(N, N);
    index_type m = 0;
    double sum_ll = 0.0;

    while (error > threshold) {
        std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                       pivot.end());
        const index_type i = it - pivot.begin();
        std::swap(pivot[m], pivot[i]);
        
        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(m, pivot[m], l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;
        
        for (index_type i = m +1; i < N; ++i) {
            double sum_ll = 0.0;
            for (index_type j = 0; j < m; ++j) {
                sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
            }
            const double l_m_pi = (this->get(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
            L.set(m, pivot[i], l_m_pi);
            
            d[pivot[i]] -= l_m_pi * l_m_pi;
        }
            
        error = 0.0;
        for (index_type i = m +1; i < N; ++i) {
            error += d[pivot[i]];
        }
        ++m;
    }

    L.transpose();
    L.resize(N, m);

    return L;
}


TlMatrix TlSymmetricMatrix::choleskyFactorization2omp(const double threshold) const
{
    const index_type N = this->getNumOfRows();
    TlVector d = this->getDiagonalElements();
    double error = d.sum();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    TlMatrix L(N, N);
    index_type m = 0;
    double sum_ll = 0.0;
    double inv_l_m_pm = 0.0;

#pragma omp parallel
    {
        while (error > threshold) {
#pragma omp single
            {
                std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                               pivot.end());
                const index_type i = it - pivot.begin();
                std::swap(pivot[m], pivot[i]);
                
                const double l_m_pm = std::sqrt(d[pivot[m]]);
                L.set(m, pivot[m], l_m_pm);
                
                //const double inv_l_m_pm = 1.0 / l_m_pm;
                inv_l_m_pm = 1.0 / l_m_pm;
            }

#pragma omp single // この箇所をOpenMP並列にすると遅くなる。要調査
            {
                for (index_type i = m +1; i < N; ++i) {
                    sum_ll = 0.0;
                    for (index_type j = 0; j < m; ++j) {
                        sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
                    }
                    const double l_m_pi = (this->get(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
                    L.set(m, pivot[i], l_m_pi);
                    
                    d[pivot[i]] -= l_m_pi * l_m_pi;
                }
            }

#pragma omp single
            {
                ++m;
                error = 0.0;
            }
#pragma omp for schedule(runtime) reduction(+:error)
            for (index_type i = m +1; i < N; ++i) {
                error += d[pivot[i]];
            }
        }
    }
        
    L.transpose();
    L.resize(N, m);

    return L;
}


// TlMatrix TlSymmetricMatrix::choleskyFactorization2omp(const double threshold) const
// {
//     const index_type N = this->getNumOfRows();
//     TlVector d = this->getDiagonalElements();
//     double error = d.sum();
//     std::vector<TlVector::size_type> pivot(N);
//     for (index_type i = 0; i < N; ++i) {
//         pivot[i] = i;
//     }

//     TlMatrix L(N, N);
//     index_type m = 0;
//     double sum_ll = 0.0;
//     double inv_l_m_pm = 0.0;

// #pragma omp parallel
//     {
//         while (error > threshold) {
// #pragma omp single
//             {
//                 std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
//                                                                                pivot.end());
//                 const index_type i = it - pivot.begin();
//                 std::swap(pivot[m], pivot[i]);
                
//                 const double l_m_pm = std::sqrt(d[pivot[m]]);
//                 L.set(m, pivot[m], l_m_pm);
                
//                 //const double inv_l_m_pm = 1.0 / l_m_pm;
//                 inv_l_m_pm = 1.0 / l_m_pm;
//             }

// #pragma omp single
//             {
//                 for (index_type i = m +1; i < N; ++i) {
//                     sum_ll = 0.0;
//                     for (index_type j = 0; j < m; ++j) {
//                         sum_ll += L.get(j, pivot[m]) * L.get(j, pivot[i]);
//                     }   
//                     const double l_m_pi = (this->get(pivot[m], pivot[i]) - sum_ll) * inv_l_m_pm;
//                     L.set(m, pivot[i], l_m_pi);
                        
//                     d[pivot[i]] -= l_m_pi * l_m_pi;
//                 }
//             }

// #pragma omp single
//             {
//                 error = 0.0;
//                 for (index_type i = m +1; i < N; ++i) {
//                     error += d[pivot[i]];
//                 }
//                 ++m;
//             }
//         }
//     }
        
//     L.transpose();
//     L.resize(N, m);

//     return L;
// }

#endif // HAVE_LAPACK


