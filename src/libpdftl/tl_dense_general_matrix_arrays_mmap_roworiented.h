#ifndef TL_DENSE_GENERAL_MATRIX_ARRAYS_MMAP_ROWORIENTED_H
#define TL_DENSE_GENERAL_MATRIX_ARRAYS_MMAP_ROWORIENTED_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <valarray>

#include "tl_dense_general_matrix_mmap.h"
#include "tl_dense_matrix_arrays_mmap_object.h"

class TlDenseGeneralMatrix_arrays_mmap_RowOriented : public TlDenseMatrix_arrays_mmap_Object {
public:
    explicit TlDenseGeneralMatrix_arrays_mmap_RowOriented(const std::string& baseFilePath, index_type row,
                                                          index_type col, int numOfSubunits = 1, int subunitID = 0,
                                                          int reservedCols = 0);
    TlDenseGeneralMatrix_arrays_mmap_RowOriented(const std::string& filePath);

    virtual ~TlDenseGeneralMatrix_arrays_mmap_RowOriented();

public:
    virtual void resize(index_type row, index_type col);
    void reserveColSize(const index_type reserveColSize);

    index_type getNumOfRows() const {
        return this->getNumOfVectors();
    };

    index_type getNumOfCols() const {
        return this->getSizeOfVector();
    };

    virtual void set(index_type row, index_type col, double value);
    virtual void add(index_type row, index_type col, double value);
    virtual double get(index_type row, index_type col) const;

    virtual std::vector<double> getRowVector(const index_type row) const;
    virtual std::size_t getRowVector(const index_type row, double* pBuf, std::size_t maxCount) const;

    virtual void setColVector(const index_type col, const std::valarray<double>& v);

public:
    /// TlDenseGeneralMatrix_Lapackオブジェクトを返す(for debug)
    TlDenseGeneralMatrix_Lapack getTlMatrixObject() const;

    /// TlDenseGeneralMatrix_arrays_ColOriented形式で保存
    // void saveByTlDenseGeneralMatrix_arrays_ColOriented(const std::string& basename) const;

public:
    std::string tempDir() const {
        return this->tempDir_;
    };
    void tempDir(const std::string& tempDir) {
        this->tempDir_ = tempDir;
    };

    void convertMemoryLayout(const std::string& tempCsfdPath, const bool verbose = false, const bool showProgress = false) const;
    void set2csfd(TlDenseGeneralMatrix_mmap* pOutMat, const bool verbose = false,
                  const bool showProgress = false) const;

private:
    std::string tempDir_;
    mutable std::string tempCsfdMatPath_;  // used for convert to CSFD matrix
    mutable bool isRemoveTempCsfdMat_;
};

std::ostream& operator<<(std::ostream& stream, const TlDenseGeneralMatrix_arrays_mmap_RowOriented& mat);

bool convert2csfd(const std::string& rvmBasePath, const int unit, const std::string& outputPath, const bool verbose,
                  const bool showProgress);
void copy2csfd(const TlMatrixObject::index_type numOfRows, const TlMatrixObject::index_type numOfCols,
               const int numOfSubunits, const int sizeOfChunk, const std::string& inMatPath, const int unit,
               TlDenseGeneralMatrix_mmap* pOutMat, const bool verbose);
bool transpose2CSFD(const std::string& rvmBasePath, const std::string& outputMatrixPath, const bool verbose = false,
                    const bool showProgress = false);

#endif  // TL_DENSE_GENERAL_MATRIX_ARRAYS_MMAP_ROWORIENTED_H
