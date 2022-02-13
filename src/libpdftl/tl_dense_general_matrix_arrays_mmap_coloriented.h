#ifndef TL_DENSE_GENERAL_MATRIX_ARRAYS_MMAP_COLORIENTED_H
#define TL_DENSE_GENERAL_MATRIX_ARRAYS_MMAP_COLORIENTED_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <valarray>

#include "tl_dense_matrix_arrays_mmap_object.h"

class TlDenseGeneralMatrix_arrays_mmap_ColOriented : public TlDenseMatrix_arrays_mmap_Object {
public:
    explicit TlDenseGeneralMatrix_arrays_mmap_ColOriented(const std::string& baseFilePath, index_type row = 1,
                                                          index_type col = 1, int numOfSubunits = 1, int subunitID = 0,
                                                          int reservedRows = 0);

    virtual ~TlDenseGeneralMatrix_arrays_mmap_ColOriented();

public:
    virtual void resize(index_type row, index_type col);
    void reserveRowSize(const index_type reserveRowSize);

    index_type getNumOfRows() const {
        return this->getSizeOfVector();
    };

    index_type getNumOfCols() const {
        return this->getNumOfVectors();
    };

    virtual void set(index_type row, index_type col, double value);
    virtual void add(index_type row, index_type col, double value);
    virtual double get(index_type row, index_type col) const;

    virtual std::vector<double> getColVector(const index_type col) const;
    virtual std::size_t getColVector(const index_type col, double* pBuf, std::size_t maxCount) const;

    virtual void setRowVector(const index_type row, const std::valarray<double>& v);

public:
    /// TlDenseGeneralMatrix_Lapackオブジェクトを返す(for debug)
    TlDenseGeneralMatrix_Lapack getTlMatrixObject() const;
};

std::ostream& operator<<(std::ostream& stream, const TlDenseGeneralMatrix_arrays_mmap_ColOriented& mat);

// bool RowVectorMatrix2CSFD_mmap(const std::string& rvmBasePath, const std::string& csfdPath, bool verbose = false,
//                                bool showProgress = false);

#endif  // TL_DENSE_GENERAL_MATRIX_ARRAYS_MMAP_COLORIENTED_H
