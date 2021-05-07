#include "tl_dense_general_matrix_arrays_mmap_roworiented.h"

#include <cassert>
#include <iostream>

#include "TlFile.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_lapack.h"
// #include "tl_dense_general_matrix_mmap.h"

TlDenseGeneralMatrix_arrays_mmap_RowOriented::TlDenseGeneralMatrix_arrays_mmap_RowOriented(
    const std::string& baseFilePath, const index_type row, const index_type col, const int numOfSubunits,
    const int subunitID, const int reservedCols)
    : TlDenseMatrix_arrays_mmap_Object(baseFilePath, row, col, numOfSubunits, subunitID, reservedCols) {
}

TlDenseGeneralMatrix_arrays_mmap_RowOriented::TlDenseGeneralMatrix_arrays_mmap_RowOriented(const std::string& filePath)
    : TlDenseMatrix_arrays_mmap_Object(filePath) {
}

// TlDenseGeneralMatrix_arrays_mmap_RowOriented::
//     TlDenseGeneralMatrix_arrays_mmap_RowOriented(
//         const TlDenseGeneralMatrix_Lapack& rhs, const int numOfSubunits,
//         const int subunitID, bool isUsingMemManager)
//     : TlDenseMatrix_arrays_mmap_Object(rhs.getNumOfRows(),
//     rhs.getNumOfCols(),
//                                        numOfSubunits, subunitID,
//                                        isUsingMemManager) {
//     const index_type numOfRows = rhs.getNumOfRows();
//     const index_type numOfCols = rhs.getNumOfCols();
//     for (index_type r = 0; r < numOfRows; ++r) {
//         for (index_type c = 0; c < numOfCols; ++c) {
//             this->set(r, c, rhs.get(r, c));
//         }
//     }
// }

// TlDenseGeneralMatrix_arrays_mmap_RowOriented::
//     TlDenseGeneralMatrix_arrays_mmap_RowOriented(
//         const TlDenseGeneralMatrix_arrays_mmap_RowOriented& rhs)
//     : TlDenseMatrix_arrays_Object(rhs) {}

TlDenseGeneralMatrix_arrays_mmap_RowOriented::~TlDenseGeneralMatrix_arrays_mmap_RowOriented() {
}

void TlDenseGeneralMatrix_arrays_mmap_RowOriented::resize(const index_type row, const index_type col) {
    TlDenseMatrix_arrays_mmap_Object::resize(row, col);
}

void TlDenseGeneralMatrix_arrays_mmap_RowOriented::reserveColSize(const index_type reserveColSize) {
    TlDenseMatrix_arrays_mmap_Object::reserveVectorSize(reserveColSize);
}

void TlDenseGeneralMatrix_arrays_mmap_RowOriented::set(const index_type row, const index_type col, const double value) {
    TlDenseMatrix_arrays_mmap_Object::set_to_vm(row, col, value);
}

void TlDenseGeneralMatrix_arrays_mmap_RowOriented::add(const index_type row, const index_type col, const double value) {
    TlDenseMatrix_arrays_mmap_Object::add_to_vm(row, col, value);
}

double TlDenseGeneralMatrix_arrays_mmap_RowOriented::get(const index_type row, const index_type col) const {
    return TlDenseMatrix_arrays_mmap_Object::get_from_vm(row, col);
}

std::vector<double> TlDenseGeneralMatrix_arrays_mmap_RowOriented::getRowVector(const index_type row) const {
    return TlDenseMatrix_arrays_mmap_Object::getVector(row);
}

std::size_t TlDenseGeneralMatrix_arrays_mmap_RowOriented::getRowVector(const index_type row, double* pBuf,
                                                                       std::size_t maxCount) const {
    return TlDenseMatrix_arrays_mmap_Object::getVector(row, pBuf, maxCount);
}

void TlDenseGeneralMatrix_arrays_mmap_RowOriented::setColVector(const index_type col, const std::valarray<double>& v) {
    TlDenseMatrix_arrays_mmap_Object::setAcrossMultipleVectors(col, v);
}

TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_arrays_mmap_RowOriented::getTlMatrixObject() const {
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlDenseGeneralMatrix_Lapack answer(numOfRows, numOfCols);

    for (index_type r = 0; r < numOfRows; ++r) {
        TlDenseVector_Lapack v = TlDenseMatrix_arrays_mmap_Object::getVector(r);
        assert(v.getSize() == numOfCols);
        for (index_type c = 0; c < numOfCols; ++c) {
            answer.set(r, c, v.get(c));
        }
    }

    return answer;
}

// void TlDenseGeneralMatrix_arrays_mmap_RowOriented::saveByTlDenseGeneralMatrix_arrays_ColOriented(
//     const std::string& basename) const {
//     TlDenseMatrix_arrays_mmap_Object::saveByTheOtherType(basename);
// }

// -----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream, const TlDenseGeneralMatrix_arrays_mmap_RowOriented& mat) {
    const TlMatrixObject::index_type numOfRows = mat.getNumOfRows();
    const TlMatrixObject::index_type numOfCols = mat.getNumOfCols();

    for (TlMatrixObject::index_type ord = 0; ord < numOfCols; ord += 10) {
        stream << "       ";
        for (TlMatrixObject::index_type j = ord; ((j < ord + 10) && (j < numOfCols)); ++j) {
            stream << TlUtils::format("   %5d th", j + 1);
        }
        stream << "\n ----";

        for (TlMatrixObject::index_type j = ord; ((j < ord + 10) && (j < numOfCols)); ++j) {
            stream << "-----------";
        }
        stream << "----\n";

        for (TlMatrixObject::index_type i = 0; i < numOfRows; ++i) {
            stream << TlUtils::format(" %5d  ", i + 1);

            for (TlMatrixObject::index_type j = ord; ((j < ord + 10) && (j < numOfCols)); ++j) {
                if (mat.getSubunitID(i) == mat.getSubunitID()) {
                    stream << TlUtils::format(" %10.6lf", mat.get(i, j));

                } else {
                    stream << " ----------";
                }
            }
            stream << "\n";
        }
        stream << "\n\n";
    }

    return stream;
}

// -----------------------------------------------------------------------------
bool convert2csfd(const std::string& rvmBasePath, const int unit, const std::string& outputPath, const bool verbose,
                  const bool showProgress) {
    bool answer = false;

    // check loadable
    bool isLoadable = false;
    TlMatrixObject::index_type numOfRows = 0;
    TlMatrixObject::index_type numOfCols = 0;
    int numOfSubunits = 0;
    int sizeOfChunk = 0;
    {
        int subunitID = 0;
        const std::string inputPath0 = TlDenseMatrix_arrays_mmap_Object::getFileName(rvmBasePath, unit);
        TlMatrixObject::index_type sizeOfVector, numOfVectors;
        isLoadable = TlDenseMatrix_arrays_mmap_Object::isLoadable(inputPath0, &numOfVectors, &sizeOfVector,
                                                                  &numOfSubunits, &subunitID, &sizeOfChunk);
        if (isLoadable != true) {
            std::cerr << "can not open file: " << inputPath0 << std::endl;
            answer = false;
        }

        // row-vector matrix
        numOfRows = numOfVectors;
        numOfCols = sizeOfVector;
    }

    {
        // load & set
        TlDenseGeneralMatrix_arrays_mmap_RowOriented partMat(rvmBasePath, 1, 1, numOfSubunits, unit);
        if (verbose) {
            std::cerr << "load partial matrix " << std::endl;
        }

        std::vector<double> chunkBuf(numOfCols * sizeOfChunk);
        std::vector<double> transBuf(numOfCols * sizeOfChunk);
        const int numOfLocalChunks =
            TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(numOfRows, numOfSubunits, sizeOfChunk);

        const TlMatrixObject::index_type localNumOfRows = numOfLocalChunks * sizeOfChunk;
        // std::cerr << TlUtils::format("localNumOfRows=%d, numOfLocalChunks=%d, sizeOfChunk=%d", localNumOfRows,
        //                              numOfLocalChunks, sizeOfChunk)
        //           << std::endl;
        TlDenseGeneralMatrix_mmap outMat(outputPath, localNumOfRows, numOfCols);
        if (verbose) {
            std::cerr << "output matrix has been prepared by mmap." << std::endl;
        }

        for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
            const TlMatrixObject::index_type chunkStartRow = sizeOfChunk * (numOfSubunits * chunk + unit);

            if (chunkStartRow < numOfRows) {
                partMat.getChunk(chunkStartRow, &(chunkBuf[0]), numOfCols * sizeOfChunk);

                // change memory layout
                const TlMatrixObject::index_type readRowChunks = std::min(sizeOfChunk, numOfRows - chunkStartRow);
                TlUtils::changeMemoryLayout(&(chunkBuf[0]), readRowChunks, numOfCols, &(transBuf[0]));

                TlDenseGeneralMatrix_Eigen tmpMat(readRowChunks, numOfCols, &(transBuf[0]));
                // std::cerr << TlUtils::format("chunk=%d; %d, %d", chunk, readRowChunks, tmpMat.getNumOfRows())
                //           << std::endl;
                outMat.block(chunk * sizeOfChunk, 0, tmpMat);
            }

            if (showProgress) {
                TlUtils::progressbar(float(chunk) / numOfLocalChunks);
            }
        }

        if (showProgress) {
            TlUtils::progressbar(1.0);
            std::cout << std::endl;
        }
    }

    return answer;
}

void copy2csfd(const TlMatrixObject::index_type numOfRows, const TlMatrixObject::index_type numOfCols,
               const int numOfSubunits, const int sizeOfChunk, const std::string& inMatPath, const int unit,
               TlDenseGeneralMatrix_mmap* pOutMat, const bool verbose) {
    TlDenseGeneralMatrix_mmap inMat(inMatPath);
    const TlMatrixObject::index_type numOfLocalRows = inMat.getNumOfRows();

    const int numOfLocalChunks =
        TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(numOfRows, numOfSubunits, sizeOfChunk);
    TlDenseGeneralMatrix_Eigen tmpMat;
    for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
        TlMatrixObject::index_type row = sizeOfChunk * chunk;
        const TlMatrixObject::index_type chunkStartRow = sizeOfChunk * (numOfSubunits * chunk + unit);
        TlMatrixObject::index_type rowDistance = std::min(sizeOfChunk, numOfRows - chunkStartRow);
        inMat.block(row, 0, rowDistance, numOfCols, &tmpMat);

        // std::cerr << TlUtils::format("chunk: %d/%d", chunk, numOfLocalChunks - 1) << std::endl;
        // std::cerr << TlUtils::format("chunkStartRow=%d, numOfLocalRows=%d, row=%d, rowDistance", chunkStartRow,
        //                              numOfLocalRows, row, rowDistance)
        //           << std::endl;
        pOutMat->block(chunkStartRow, 0, tmpMat);
    }
}

//
bool transpose2CSFD(const std::string& rvmBasePath, const std::string& outputMatrixPath, const bool verbose,
                    const bool showProgress) {
    bool answer = false;

    // 初期データ読み込み
    TlMatrixObject::index_type numOfRows = 0;
    TlMatrixObject::index_type numOfCols = 0;
    int numOfSubunits = 0;
    int sizeOfChunk = 0;
    bool isLoadable = false;
    {
        int subunitID = 0;
        const std::string inputPath0 = TlDenseMatrix_arrays_mmap_Object::getFileName(rvmBasePath, subunitID);
        TlMatrixObject::index_type sizeOfVector, numOfVectors;
        isLoadable = TlDenseMatrix_arrays_mmap_Object::isLoadable(inputPath0, &numOfVectors, &sizeOfVector,
                                                                  &numOfSubunits, &subunitID, &sizeOfChunk);
        if (isLoadable != true) {
            std::cerr << "can not open file: " << inputPath0 << std::endl;
            return false;
        }

        numOfRows = numOfVectors;
        numOfCols = sizeOfVector;

        if (verbose) {
            std::cerr << "rows: " << numOfRows << std::endl;
            std::cerr << "cols: " << numOfCols << std::endl;
            std::cerr << "units: " << numOfSubunits << std::endl;
            std::cerr << "chunk: " << sizeOfChunk << std::endl;
        }
    }

    // prepare CSFD file
    if (TlFile::isExistFile(outputMatrixPath)) {
        if (verbose) {
            std::cerr << "file overwrite: " << outputMatrixPath << std::endl;
        }
        TlFile::remove(outputMatrixPath);
    }
    TlDenseGeneralMatrix_mmap outMat(outputMatrixPath, numOfRows, numOfCols);
    if (verbose) {
        std::cerr << "output matrix has been prepared by mmap." << std::endl;
    }

    // 最終書き込み
    for (int unit = 0; unit < numOfSubunits; ++unit) {
        // std::cerr << "copy: " << unit << std::endl;
        const std::string tempMatPath = TlUtils::format("Ljk.temp.%d", unit);

        convert2csfd(rvmBasePath, unit, tempMatPath, verbose, showProgress);
        copy2csfd(numOfRows, numOfCols, numOfSubunits, sizeOfChunk, tempMatPath, unit, &outMat, verbose);
        TlFile::remove(tempMatPath);
    }

    return answer;
}
