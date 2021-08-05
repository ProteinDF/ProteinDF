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

#ifndef TL_DENSE_SYMMETRIC_MATRIX_IMPL_SCALAPACK_H
#define TL_DENSE_SYMMETRIC_MATRIX_IMPL_SCALAPACK_H

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "tl_dense_general_matrix_impl_scalapack.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

class TlDenseVector_ImplLapack;

class TlDenseSymmetricMatrix_ImplScalapack
    : public TlDenseGeneralMatrix_ImplScalapack {
public:
    enum DIAGONAL_METHOD { QR,
                           DIVIDE_AND_CONQUER };

    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlDenseSymmetricMatrix_ImplScalapack(
        const TlMatrixObject::index_type dim = 0);
    TlDenseSymmetricMatrix_ImplScalapack(
        const TlDenseSymmetricMatrix_ImplScalapack& rhs);
    TlDenseSymmetricMatrix_ImplScalapack(
        const TlDenseGeneralMatrix_ImplScalapack& rhs);
    //   TlDenseSymmetricMatrix_ImplScalapack(const TlDenseVector_ImplScalapack&
    //   v,
    //                                        const TlMatrixObject::index_type
    //                                        dim);
    virtual ~TlDenseSymmetricMatrix_ImplScalapack();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
public:
    // virtual TlMatrixObject::size_type getNumOfElements() const;
    virtual double get(TlMatrixObject::index_type row,
                       TlMatrixObject::index_type col) const;
    virtual void set(TlMatrixObject::index_type row,
                     TlMatrixObject::index_type col, const double value);
    virtual void add(TlMatrixObject::index_type row,
                     TlMatrixObject::index_type col, double value);

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    // TlDenseSymmetricMatrix_ImplScalapack& operator*=(const double coef);
    // const TlDenseGeneralMatrix_ImplLapack operator*=(
    //   const TlDenseSymmetricMatrix_ImplScalapack& rhs);

    // ---------------------------------------------------------------------------
    // operations
    // ---------------------------------------------------------------------------
public:
    // virtual double sum() const;
    // virtual double getRMS() const;
    // virtual double getMaxAbsoluteElement(
    //     TlMatrixObject::index_type* outRow,
    //     TlMatrixObject::index_type* outCol) const;

    // const TlDenseGeneralMatrix_ImplLapack& dotInPlace(
    //     const TlDenseGeneralMatrix_ImplLapack& rhs);
    // TlDenseSymmetricMatrix_ImplScalapack transpose() const;
    // TlDenseSymmetricMatrix_ImplScalapack inverse() const;

    bool eig(TlDenseVector_ImplLapack* pEigVal,
             TlDenseGeneralMatrix_ImplScalapack* pEigVec,
             const DIAGONAL_METHOD method = DIVIDE_AND_CONQUER) const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
public:
    virtual bool load(const std::string& filePath);
    virtual bool save(const std::string& filePath) const;

public:
    // virtual void dump(double* buf, const std::size_t size) const;
    // virtual void restore(const double* buf, const std::size_t size);

    // ---------------------------------------------------------------------------
    // protected
    // ---------------------------------------------------------------------------
protected:
    virtual TlMatrixObject::size_type index(
        TlMatrixObject::index_type row, TlMatrixObject::index_type col) const;

    virtual bool load(std::ifstream& ifs);

    virtual std::vector<TlMatrixObject::MatrixElement>
    getMatrixElementsInLocal2() const;
    // ---------------------------------------------------------------------------
    // private
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // friend
    // ---------------------------------------------------------------------------
    friend bool diagonalByScaLapack_QR(
        const TlDenseSymmetricMatrix_ImplScalapack& inMatrix,
        TlDenseVector_ImplLapack* outEigVal,
        TlDenseGeneralMatrix_ImplScalapack* outEigVec);
    friend bool diagonalByScaLapack_DC(
        const TlDenseSymmetricMatrix_ImplScalapack& inMatrix,
        TlDenseVector_ImplLapack* outEigVal,
        TlDenseGeneralMatrix_ImplScalapack* outEigVec);

    // friend TlDenseGeneralMatrix_ImplLapack operator*(
    //     const TlDenseSymmetricMatrix_ImplScalapack& rhs1,
    //     const TlDenseGeneralMatrix_ImplLapack& rhs2);
    // friend TlDenseGeneralMatrix_ImplLapack operator*(
    //     const TlDenseGeneralMatrix_ImplLapack& rhs1,
    //     const TlDenseSymmetricMatrix_ImplScalapack& rhs2);
    //
    // friend TlDenseVector_ImplLapack operator*(
    //     const TlDenseSymmetricMatrix_ImplScalapack& mat,
    //     const TlDenseVector_ImplLapack& vec);
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IMPL_SCALAPACK_H
