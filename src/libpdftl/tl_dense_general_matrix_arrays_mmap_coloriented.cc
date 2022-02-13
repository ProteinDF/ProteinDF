#include "tl_dense_general_matrix_arrays_mmap_coloriented.h"

#include <cassert>
#include <iostream>

#include "TlFile.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_general_matrix_mmap.h"

TlDenseGeneralMatrix_arrays_mmap_ColOriented::TlDenseGeneralMatrix_arrays_mmap_ColOriented(
    const std::string& baseFilePath, const index_type row, const index_type col, const int numOfSubunits,
    const int subunitID, const int reservedRows)
    : TlDenseMatrix_arrays_mmap_Object(baseFilePath, col, row, numOfSubunits, subunitID, reservedRows) {
}

TlDenseGeneralMatrix_arrays_mmap_ColOriented::~TlDenseGeneralMatrix_arrays_mmap_ColOriented() {
}

void TlDenseGeneralMatrix_arrays_mmap_ColOriented::resize(const index_type row, const index_type col) {
    TlDenseMatrix_arrays_mmap_Object::resize(col, row);
}

void TlDenseGeneralMatrix_arrays_mmap_ColOriented::reserveRowSize(const index_type reserveRowSize) {
    TlDenseMatrix_arrays_mmap_Object::reserveVectorSize(reserveRowSize);
}

void TlDenseGeneralMatrix_arrays_mmap_ColOriented::set(const index_type row, const index_type col, const double value) {
    TlDenseMatrix_arrays_mmap_Object::set_to_vm(col, row, value);
}

void TlDenseGeneralMatrix_arrays_mmap_ColOriented::add(const index_type row, const index_type col, const double value) {
    TlDenseMatrix_arrays_mmap_Object::add_to_vm(col, row, value);
}

double TlDenseGeneralMatrix_arrays_mmap_ColOriented::get(const index_type row, const index_type col) const {
    return TlDenseMatrix_arrays_mmap_Object::get_from_vm(col, row);
}

std::vector<double> TlDenseGeneralMatrix_arrays_mmap_ColOriented::getColVector(const index_type col) const {
    return TlDenseMatrix_arrays_mmap_Object::getVector(col);
}

std::size_t TlDenseGeneralMatrix_arrays_mmap_ColOriented::getColVector(const index_type col, double* pBuf,
                                                                       std::size_t maxCount) const {
    return TlDenseMatrix_arrays_mmap_Object::getVector(col, pBuf, maxCount);
}

void TlDenseGeneralMatrix_arrays_mmap_ColOriented::setRowVector(const index_type row, const std::valarray<double>& v) {
    TlDenseMatrix_arrays_mmap_Object::setAcrossMultipleVectors(row, v);
}

TlDenseGeneralMatrix_Lapack TlDenseGeneralMatrix_arrays_mmap_ColOriented::getTlMatrixObject() const {
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlDenseGeneralMatrix_Lapack answer(numOfRows, numOfCols);

    for (index_type c = 0; c < numOfCols; ++c) {
        TlDenseVector_Lapack v = TlDenseMatrix_arrays_mmap_Object::getVector(c);
        assert(v.getSize() == numOfRows);
        for (index_type r = 0; r < numOfRows; ++r) {
            answer.set(r, c, v.get(r));
        }
    }

    return answer;
}

// -----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream, const TlDenseGeneralMatrix_arrays_mmap_ColOriented& mat) {
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
// bool RowVectorMatrix2CSFD_mmap(const std::string& rvmBasePath, const std::string& csfdPath, bool verbose,
//                                bool showProgress) {
//     // check
//     TlMatrixObject::index_type numOfRows = 0;
//     TlMatrixObject::index_type numOfCols = 0;
//     int numOfSubunits = 0;
//     int sizeOfChunk = 0;
//     {
//         int subunitID = 0;
//         const std::string inputPath0 = TlDenseMatrix_arrays_mmap_Object::getFileName(rvmBasePath, subunitID);
//         TlMatrixObject::index_type sizeOfVector, numOfVectors;
//         const bool isLoadable = TlDenseMatrix_arrays_mmap_Object::isLoadable(inputPath0, &numOfVectors,
//         &sizeOfVector,
//                                                                              &numOfSubunits, &subunitID,
//                                                                              &sizeOfChunk);
//         if (isLoadable != true) {
//             std::cerr << "can not open file: " << inputPath0 << std::endl;
//             return false;
//         }

//         // row-vector matrix
//         numOfRows = numOfVectors;
//         numOfCols = sizeOfVector;

//         if (verbose) {
//             std::cerr << "rows: " << numOfRows << std::endl;
//             std::cerr << "cols: " << numOfCols << std::endl;
//             std::cerr << "units: " << numOfSubunits << std::endl;
//             std::cerr << "chunk: " << sizeOfChunk << std::endl;
//         }
//     }

//     // prepare output
//     if (TlFile::isExistFile(csfdPath)) {
//         if (verbose) {
//             std::cerr << "file overwrite: " << csfdPath << std::endl;
//         }
//         TlFile::remove(csfdPath);
//     }
//     TlDenseGeneralMatrix_mmap fileMat(csfdPath, numOfRows, numOfCols);
//     if (verbose) {
//         std::cerr << "output matrix has been prepared by mmap." << std::endl;
//     }

//     // load & set
//     for (int unit = 0; unit < numOfSubunits; ++unit) {
//         if (verbose) {
//             std::cerr << TlUtils::format("%d / %d", unit + 1, numOfSubunits) << std::endl;
//         }

//         TlDenseGeneralMatrix_arrays_mmap_ColOriented partMat(rvmBasePath, 1, 1, numOfSubunits, unit);
//         if (verbose) {
//             std::cerr << "load partial matrix " << std::endl;
//         }

//         // std::vector<double> vtr(numOfCols);
//         // for (TlMatrixObject::index_type r = unit; r < numOfRows; ++r) {
//         //     if (i == m.getSubunitID(r)) {
//         //         m.getVector(r, &(vtr[0]), numOfCols);
//         //         fileMat.setRowVector(r, vtr);
//         //     }

//         //     if (showProgress) {
//         //         TlUtils::progressbar(float(r) / numOfRows);
//         //     }
//         // }
//         std::vector<double> chunkBuf(numOfCols * sizeOfChunk);
//         std::vector<double> transBuf(numOfCols * sizeOfChunk);
//         const int numOfLocalChunks =
//             TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(numOfRows, numOfSubunits, sizeOfChunk);
//         for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
//             const TlMatrixObject::index_type chunkStartRow = sizeOfChunk * (numOfSubunits * chunk + unit);

//             if (chunkStartRow < numOfRows) {
//                 partMat.getChunk(chunkStartRow, &(chunkBuf[0]), numOfCols * sizeOfChunk);

//                 // change memory layout
//                 const TlMatrixObject::index_type readRowChunks = std::min(sizeOfChunk, numOfRows - chunkStartRow);
//                 TlUtils::changeMemoryLayout(&(chunkBuf[0]), readRowChunks, numOfCols, &(transBuf[0]));

//                 TlDenseGeneralMatrix_Eigen tmpMat(readRowChunks, numOfCols, &(transBuf[0]));
//                 fileMat.block(chunkStartRow, 0, tmpMat);
//             }

//             TlUtils::progressbar(float(chunk) / numOfLocalChunks);
//         }

//         if (showProgress) {
//             TlUtils::progressbar(1.0);
//             std::cout << std::endl;
//         }
//     }

//     if (verbose) {
//         std::cerr << "end." << std::endl;
//     }

//     return true;
// }
