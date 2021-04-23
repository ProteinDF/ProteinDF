#include "tl_matrix_utils.h"

#include <iostream>

#include "TlUtils.h"
#include "tl_dense_matrix_arrays_object.h"
// #include "nullptr.h"

// -----------------------------------------------------------------------------
// public
// -----------------------------------------------------------------------------
bool TlMatrixUtils::isLoadable(const std::string& filePath, const TlMatrixObject::MatrixType matrixType) {
    // TlMatrixObject::MatrixType loadMatrixType;
    TlMatrixObject::HeaderInfo headerInfo;
    const TlMatrixUtils::FileSize headerSize = TlMatrixUtils::getHeaderInfo(filePath, &headerInfo);

    bool answer = false;
    if ((headerSize > 0) && (headerInfo.matrixType == matrixType)) {
        answer = true;
    }
    return answer;
}

TlMatrixUtils::FileSize TlMatrixUtils::getHeaderInfo(const std::string& filepath,
                                                     TlMatrixObject::HeaderInfo* pHeaderInfo) {
    FileSize headerSize = 0;

    std::ios_base::sync_with_stdio(false);
    std::fstream fs;
    fs.open(filepath.c_str(), std::ios::in | std::ios::binary);

    if (!fs.fail()) {
        headerSize = TlMatrixUtils::getHeaderInfo(fs, pHeaderInfo);
    } else {
        std::cerr << TlUtils::format("cannot open matrix file: %s @%s:%d", filepath.c_str(), __FILE__, __LINE__)
                  << std::endl;
    }

    fs.close();

    return headerSize;
}

// TlMatrixUtils::FileSize TlMatrixUtils::getHeaderInfo(const std::string& filepath,
//                                                      TlMatrixObject::MatrixType* pMatrixType,
//                                                      TlMatrixObject::index_type* pNumOfRows,
//                                                      TlMatrixObject::index_type* pNumOfCols, std::size_t*
//                                                      pNumOfItems, int* pNumOfSubunits, int* pSubunitId, int*
//                                                      pSizeOfChunk) {
//     FileSize headerSize = 0;

//     std::ios_base::sync_with_stdio(false);
//     std::fstream fs;
//     fs.open(filepath.c_str(), std::ios::in | std::ios::binary);

//     if (!fs.fail()) {
//         headerSize = TlMatrixUtils::getHeaderInfo(fs, pMatrixType, pNumOfRows, pNumOfCols, pNumOfItems,
//         pNumOfSubunits,
//                                                   pSubunitId, pSizeOfChunk);
//     } else {
//         std::cerr << TlUtils::format("cannot open matrix file: %s @%s:%d", filepath.c_str(), __FILE__, __LINE__)
//                   << std::endl;
//     }

//     fs.close();

//     return headerSize;
// }

// TlMatrixUtils::FileSize TlMatrixUtils::getHeaderInfo(std::fstream& fs, TlMatrixObject::MatrixType* pMatrixType,
//                                                      TlMatrixObject::index_type* pNumOfRows,
//                                                      TlMatrixObject::index_type* pNumOfCols, std::size_t*
//                                                      pNumOfItems, int* pNumOfSubunits, int* pSubunitId, int*
//                                                      pSizeOfChunk) {
//     const FileSize headerSize = TlMatrixUtils::getHeaderSize_templ1<std::fstream>(
//         fs, pMatrixType, pNumOfRows, pNumOfCols, pNumOfItems, pNumOfSubunits, pSubunitId, pSizeOfChunk);
//     fs.seekg(headerSize, std::ios_base::beg);

//     return headerSize;
// }

TlMatrixUtils::FileSize TlMatrixUtils::getHeaderInfo(std::fstream& fs, TlMatrixObject::HeaderInfo* pHeaderInfo) {
    const FileSize headerSize = TlMatrixUtils::getHeaderSize_templ1<std::fstream>(fs, pHeaderInfo);
    fs.seekg(headerSize, std::ios_base::beg);

    return headerSize;
}

bool TlMatrixUtils::saveMatrix(const std::string& filepath, const TlMatrixObject::MatrixType matrixType,
                               const TlMatrixObject::index_type rows, const TlMatrixObject::index_type cols,
                               const double* pBuf, const std::size_t numOfItems) {
    const char nMatrixType = static_cast<char>(matrixType);
    std::ios_base::sync_with_stdio(false);

    std::fstream fs(filepath.c_str(), std::ios::out | std::ios::binary);
    fs.write(&nMatrixType, sizeof(char));
    fs.write(reinterpret_cast<const char*>(&rows), sizeof(TlMatrixObject::index_type));
    fs.write(reinterpret_cast<const char*>(&cols), sizeof(TlMatrixObject::index_type));
    fs.write(reinterpret_cast<const char*>(&pBuf), sizeof(double) * numOfItems);

    fs.close();

    return true;
}

void TlMatrixUtils::CSFD2RSFD(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                              const double* pBufIn, double* pBufOut) {
    for (TlMatrixObject::index_type c = 0; c < col; ++c) {
        for (TlMatrixObject::index_type r = 0; r < row; ++r) {
            const TlMatrixObject::size_type rsfd_index = c + col * r;
            pBufOut[rsfd_index] = pBufIn[r + row * c];
        }
    }
}

void TlMatrixUtils::RSFD2CSFD(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                              const double* pBufIn, double* pBufOut) {
    for (TlMatrixObject::index_type c = 0; c < col; ++c) {
        for (TlMatrixObject::index_type r = 0; r < row; ++r) {
            const TlMatrixObject::size_type rsfd_index = c + col * r;
            pBufOut[r + row * c] = pBufIn[rsfd_index];
        }
    }
}

// -----------------------------------------------------------------------------
// protected
// -----------------------------------------------------------------------------
template <typename StreamType>
TlMatrixUtils::FileSize TlMatrixUtils::getHeaderSize_templ1(StreamType& s, TlMatrixObject::HeaderInfo* pHeaderInfo) {
    // TlMatrixObject::MatrixType matrixType;
    // TlMatrixObject::index_type rows = 0;
    // TlMatrixObject::index_type cols = 0;
    // std::size_t numOfItems = 0;
    // int numOfSubunits = 0;
    // int subunitId = 0;
    // int sizeOfChunk = 0;

    FileSize headerSize = 0;
    // char case:
    {
        // char type = 0;
        headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, char, int>(s, pHeaderInfo);
        // std::cout << TlUtils::format("%s: %d (%d, %d)", "ci", type, rows,
        // cols)
        //           << std::endl;

        if (headerSize == 0) {
            headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, char, long>(s, pHeaderInfo);
            // std::cout << TlUtils::format("%s: %d (%d, %d)", "cl", type, rows,
            // cols)
            //           << std::endl;
        }
        // matrixType = static_cast<TlMatrixObject::MatrixType>(type);
    }

    // int case:
    if (headerSize == 0) {
        // int type = 0;
        headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, int, int>(s, pHeaderInfo);
        // std::cout << TlUtils::format("%s: %d (%d, %d)", "ii", type, rows,
        // cols)
        //           << std::endl;

        if (headerSize == 0) {
            headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, int, long>(s, pHeaderInfo);
            // std::cout << TlUtils::format("%s: %d (%d, %d)", "il", type, rows,
            // cols)
            //           << std::endl;
        }
        // matrixType = static_cast<TlMatrixObject::MatrixType>(type);
    }

    // long case:
    if (headerSize == 0) {
        // long type = 0;
        headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, long, int>(s, pHeaderInfo);
        // std::cout << TlUtils::format("%s: %d (%d, %d)", "li", type, rows,
        // cols)
        //           << std::endl;

        if (headerSize == 0) {
            headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, long, long>(s, pHeaderInfo);
            // std::cout << TlUtils::format("%s: %d (%d, %d)", "ll", type, rows,
            // cols)
            //           << std::endl;
        }
        // matrixType = static_cast<TlMatrixObject::MatrixType>(type);
    }

    // if (headerSize > 0) {
    //     if (pMatrixType != NULL) {
    //         *pMatrixType = matrixType;
    //     }
    //     if (pNumOfRows != NULL) {
    //         *pNumOfRows = rows;
    //     }
    //     if (pNumOfCols != NULL) {
    //         *pNumOfCols = cols;
    //     }
    //     if (pNumOfItems != NULL) {
    //         *pNumOfItems = numOfItems;
    //     }
    //     if (pNumOfSubunits != NULL) {
    //         *pNumOfSubunits = numOfSubunits;
    //     }
    //     if (pSubunitId != NULL) {
    //         *pSubunitId = subunitId;
    //     }
    //     if (pSizeOfChunk != NULL) {
    //         *pSizeOfChunk = sizeOfChunk;
    //     }
    // }

    return headerSize;
}

template <typename StreamType, typename MatrixType, typename IndexType>
TlMatrixUtils::FileSize TlMatrixUtils::getHeaderSize_templ2(StreamType& s, TlMatrixObject::HeaderInfo* pHeaderInfo) {
    MatrixType matrixType = 0;
    TlMatrixObject::HeaderInfo headerInfo;

    // get file size
    FileSize fileSize = 0;
    {
        s.seekg(0, std::ios_base::end);
        fileSize = s.tellg();
    }

    s.seekg(0, std::ios_base::beg);
    s.read((char*)&(matrixType), sizeof(MatrixType));
    headerInfo.matrixType = static_cast<TlMatrixObject::MatrixType>(matrixType);
    // std::cerr << "matrixType (233): " << (int)matrixType << std::endl;
    FileSize headerSize = sizeof(MatrixType);

    switch (matrixType) {
        case TlMatrixObject::RSFD:
        case TlMatrixObject::CSFD:
        case TlMatrixObject::RLHD:
        case TlMatrixObject::RUHD:
        case TlMatrixObject::CLHD:
        case TlMatrixObject::CUHD: {
            IndexType rows = 0;
            IndexType cols = 0;
            s.read((char*)&(rows), sizeof(IndexType));
            s.read((char*)&(cols), sizeof(IndexType));

            headerInfo.numOfRows = rows;
            headerInfo.numOfCols = cols;

            headerSize += sizeof(IndexType) * 2;
        } break;

        case TlMatrixObject::COOF:
        case TlMatrixObject::COOS: {
            IndexType rows = 0;
            IndexType cols = 0;
            std::size_t numOfItems = 0;

            s.read((char*)&(rows), sizeof(IndexType));
            s.read((char*)&(cols), sizeof(IndexType));
            s.read(reinterpret_cast<char*>(&numOfItems), sizeof(std::size_t));

            headerInfo.numOfRows = rows;
            headerInfo.numOfCols = cols;
            headerInfo.numOfItems = numOfItems;

            headerSize += sizeof(IndexType) * 2 + sizeof(std::size_t);
        } break;

        case TlMatrixObject::ABGD: {
            IndexType numOfVectors = 0;
            IndexType sizeOfVector = 0;
            IndexType reservedSizeOfVector = 0;
            int numOfSubunits = 0;
            int subunitId = 0;
            int sizeOfChunk = 0;

            s.read(reinterpret_cast<char*>(&numOfVectors), sizeof(IndexType));
            s.read(reinterpret_cast<char*>(&sizeOfVector), sizeof(IndexType));
            s.read(reinterpret_cast<char*>(&reservedSizeOfVector), sizeof(IndexType));
            s.read(reinterpret_cast<char*>(&numOfSubunits), sizeof(int));
            s.read(reinterpret_cast<char*>(&subunitId), sizeof(int));
            s.read(reinterpret_cast<char*>(&sizeOfChunk), sizeof(int));

            headerInfo.numOfVectors = numOfVectors;
            headerInfo.sizeOfVector = sizeOfVector;
            headerInfo.reservedSizeOfVector = reservedSizeOfVector;
            headerInfo.numOfSubunits = numOfSubunits;
            headerInfo.subunitId = subunitId;
            headerInfo.sizeOfChunk = sizeOfChunk;

            headerSize += sizeof(IndexType) * 3 + sizeof(int) * 3;
        } break;

        default: {
            std::cerr << TlUtils::format("program error. %s@%d", __FILE__, __LINE__) << std::endl;
        } break;
    }

    if (pHeaderInfo != NULL) {
        *pHeaderInfo = headerInfo;
    }
    const FileSize estimatedFileSize = headerSize + TlMatrixUtils::estimateFileSize(headerInfo);
    // std::cerr << TlUtils::format("%s@%d [%d, %d] type: %d, fs: %ld/%ld %d/%d/%d", __FILE__, __LINE__,
    //                              headerInfo.numOfRows, headerInfo.numOfCols, headerInfo.matrixType, fileSize,
    //                              estimatedFileSize, headerInfo.numOfSubunits, headerInfo.subunitId,
    //                              headerInfo.sizeOfChunk)
    //           << std::endl;

    FileSize answer = 0;
    if (estimatedFileSize == fileSize) {
        answer = headerSize;
    }

    return answer;
}

std::size_t TlMatrixUtils::estimateFileSize(const TlMatrixObject::HeaderInfo& headerInfo) {
    std::size_t answer = 0;

    const std::size_t size_index = sizeof(TlMatrixObject::index_type);
    const std::size_t size_double = sizeof(double);

    switch (headerInfo.matrixType) {
        case TlMatrixObject::RSFD:
        case TlMatrixObject::CSFD:
            answer = size_double * headerInfo.numOfRows * headerInfo.numOfCols;
            break;

        case TlMatrixObject::RLHD:
        case TlMatrixObject::RUHD:
        case TlMatrixObject::CLHD:
        case TlMatrixObject::CUHD: {
            TlMatrixObject::index_type dim = headerInfo.numOfRows;
            if (dim != headerInfo.numOfCols) {
                answer = 0;
                break;
            }
            answer = size_double * dim * (dim + 1) / 2;
        } break;

        case TlMatrixObject::COOF:
        case TlMatrixObject::COOS:
            answer = (size_index * 2 + size_double) * headerInfo.numOfItems;
            break;

        case TlMatrixObject::ABGD: {
            const std::size_t numOfChunks = TlDenseMatrix_arrays_Object::getNumOfLocalChunks(
                headerInfo.numOfVectors, headerInfo.numOfSubunits, headerInfo.sizeOfChunk);
            // std::cerr << TlUtils::format("chunks: %d", numOfChunks) << std::endl;
            answer = size_double * (numOfChunks * headerInfo.sizeOfChunk) * headerInfo.reservedSizeOfVector;
        } break;

        default:
            std::cerr << TlUtils::format("program error: matrix_type=%d %s@%d", static_cast<int>(headerInfo.matrixType),
                                         __FILE__, __LINE__)
                      << std::endl;
            break;
    }

    return answer;
}
