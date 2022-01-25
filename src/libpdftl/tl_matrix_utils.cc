#include "tl_matrix_utils.h"

#include <iostream>

#include "TlLogging.h"
#include "TlUtils.h"
#include "tl_dense_matrix_arrays_object.h"
// #include "nullptr.h"

// -----------------------------------------------------------------------------
// public
// -----------------------------------------------------------------------------
bool TlMatrixUtils::isLoadable(const std::string& filePath, const TlMatrixObject::MatrixType matrixType) {
    TlMatrixObject::HeaderInfo headerInfo;
    const bool isLoadable = TlMatrixUtils::getHeaderInfo(filePath, &headerInfo);

    bool answer = false;
    if ((isLoadable == true) && (headerInfo.matrixType == matrixType)) {
        answer = true;
    }

    return answer;
}

bool TlMatrixUtils::getHeaderInfo(const std::string& filepath, TlMatrixObject::HeaderInfo* pHeaderInfo) {
    bool answer = false;

    std::ios_base::sync_with_stdio(false);
    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ios::in | std::ios::binary);

    if (!ifs.fail()) {
        answer = TlMatrixUtils::getHeaderInfo(ifs, pHeaderInfo);
    } else {
        std::cerr << TlUtils::format("cannot open matrix file: %s @%s:%d", filepath.c_str(), __FILE__, __LINE__)
                  << std::endl;
    }

    ifs.close();

    return answer;
}

bool TlMatrixUtils::getHeaderInfo(std::ifstream& ifs, TlMatrixObject::HeaderInfo* pHeaderInfo) {
    TlMatrixObject::HeaderInfo headerInfo;
    const bool answer = TlMatrixUtils::getHeaderSize_template<std::ifstream>(ifs, &headerInfo);
    if (answer == true) {
        const std::size_t headerSize = headerInfo.headerSize;
        ifs.seekg(headerSize, std::ios_base::beg);
    }

    if (pHeaderInfo != NULL) {
        *pHeaderInfo = headerInfo;
    }

    return answer;
}

bool TlMatrixUtils::getHeaderInfo(std::fstream& fs, TlMatrixObject::HeaderInfo* pHeaderInfo) {
    TlMatrixObject::HeaderInfo headerInfo;
    const bool answer = TlMatrixUtils::getHeaderSize_template<std::fstream>(fs, &headerInfo);
    if (answer == true) {
        const std::size_t headerSize = headerInfo.headerSize;
        fs.seekg(headerSize, std::ios_base::beg);
    }

    if (pHeaderInfo != NULL) {
        *pHeaderInfo = headerInfo;
    }

    return answer;
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
bool TlMatrixUtils::getHeaderSize_template(StreamType& s, TlMatrixObject::HeaderInfo* pHeaderInfo) {
    TlLogging& log = TlLogging::getInstance();

    bool answer = false;

    // char case:
    {
        // char type = 0;
        answer = TlMatrixUtils::getHeaderSizeByType_template<StreamType, char, int>(s, pHeaderInfo);
        log.debug(TlUtils::format("ci: %s(%d) (%d, %d)",
                                  TlMatrixObject::matrixTypeStr(pHeaderInfo->matrixType).c_str(), pHeaderInfo->matrixType, pHeaderInfo->numOfRows, pHeaderInfo->numOfCols));

        if (answer == false) {
            answer = TlMatrixUtils::getHeaderSizeByType_template<StreamType, char, long>(s, pHeaderInfo);
            log.debug(TlUtils::format("cl: %s(%d) (%d, %d)",
                                      TlMatrixObject::matrixTypeStr(pHeaderInfo->matrixType).c_str(), pHeaderInfo->matrixType, pHeaderInfo->numOfRows, pHeaderInfo->numOfCols));
        }
    }

    // int case:
    if (answer == false) {
        answer = TlMatrixUtils::getHeaderSizeByType_template<StreamType, int, int>(s, pHeaderInfo);
        log.debug(TlUtils::format("ii: %s(%d) (%d, %d)",
                                  TlMatrixObject::matrixTypeStr(pHeaderInfo->matrixType).c_str(), pHeaderInfo->matrixType, pHeaderInfo->numOfRows, pHeaderInfo->numOfCols));

        if (answer == false) {
            answer = TlMatrixUtils::getHeaderSizeByType_template<StreamType, int, long>(s, pHeaderInfo);
            log.debug(TlUtils::format("il: %s(%d) (%d, %d)",
                                      TlMatrixObject::matrixTypeStr(pHeaderInfo->matrixType).c_str(), pHeaderInfo->matrixType, pHeaderInfo->numOfRows, pHeaderInfo->numOfCols));
        }
    }

    // long case:
    if (answer == false) {
        answer = TlMatrixUtils::getHeaderSizeByType_template<StreamType, long, int>(s, pHeaderInfo);
        log.debug(TlUtils::format("li: %s(%d) (%d, %d)",
                                  TlMatrixObject::matrixTypeStr(pHeaderInfo->matrixType).c_str(), pHeaderInfo->matrixType, pHeaderInfo->numOfRows, pHeaderInfo->numOfCols));

        if (answer == false) {
            answer = TlMatrixUtils::getHeaderSizeByType_template<StreamType, long, long>(s, pHeaderInfo);
            log.debug(TlUtils::format("ll: %s(%d) (%d, %d)",
                                      TlMatrixObject::matrixTypeStr(pHeaderInfo->matrixType).c_str(), pHeaderInfo->matrixType, pHeaderInfo->numOfRows, pHeaderInfo->numOfCols));
        }
    }

    return answer;
}

template <typename StreamType, typename MatrixType, typename IndexType>
bool TlMatrixUtils::getHeaderSizeByType_template(StreamType& s, TlMatrixObject::HeaderInfo* pHeaderInfo) {
    TlLogging& log = TlLogging::getInstance();

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
    std::ifstream::pos_type pos_type = s.tellg();
    headerInfo.matrixType = static_cast<TlMatrixObject::MatrixType>(matrixType);
    // std::cerr << "matrixType: " << TlMatrixObject::matrixTypeStr(headerInfo.matrixType) << std::endl;

    FileSize headerSize = sizeof(char);
    bool answer = false;

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

            headerSize = sizeof(MatrixType) + sizeof(IndexType) * 2;
            const FileSize estimatedFileSize = headerSize + TlMatrixUtils::estimateFileSize(headerInfo);
            if (estimatedFileSize == fileSize) {
                headerInfo.headerSize = headerSize;
                headerInfo.version = 0;

                answer = true;
            }
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

            headerSize = sizeof(MatrixType) + sizeof(IndexType) * 2 + sizeof(std::size_t);
            const FileSize estimatedFileSize = headerSize + TlMatrixUtils::estimateFileSize(headerInfo);
            if (estimatedFileSize == fileSize) {
                headerInfo.headerSize = headerSize;
                headerInfo.version = 0;

                answer = true;
            }
        } break;

        case TlMatrixObject::ABGD: {
            // old format
            {
                IndexType numOfVectors = 0;
                IndexType sizeOfVector = 0;
                IndexType reservedSizeOfVector = 0;
                int numOfSubunits = 0;
                int subunitId = 0;
                int sizeOfChunk = 0;

                // s.seekg(pos_type);
                s.read(reinterpret_cast<char*>(&numOfVectors), sizeof(IndexType));
                s.read(reinterpret_cast<char*>(&sizeOfVector), sizeof(IndexType));
                // skip reading "reservedSizeOfVector" for old version
                reservedSizeOfVector = sizeOfVector;
                s.read(reinterpret_cast<char*>(&numOfSubunits), sizeof(int));
                s.read(reinterpret_cast<char*>(&subunitId), sizeof(int));
                s.read(reinterpret_cast<char*>(&sizeOfChunk), sizeof(int));
                log.debug(TlUtils::format("ABGD(rev.1) (%d, %d, %d) [%d/%d] @%d",
                                          numOfVectors, sizeOfVector, reservedSizeOfVector, subunitId, numOfSubunits, sizeOfChunk));

                headerInfo.numOfVectors = numOfVectors;
                headerInfo.sizeOfVector = sizeOfVector;
                headerInfo.reservedSizeOfVector = reservedSizeOfVector;
                headerInfo.numOfSubunits = numOfSubunits;
                headerInfo.subunitId = subunitId;
                headerInfo.sizeOfChunk = sizeOfChunk;

                headerSize = sizeof(MatrixType) + sizeof(IndexType) * 2 + sizeof(int) * 3;
                const FileSize estimatedFileSize = headerSize + TlMatrixUtils::estimateFileSize(headerInfo);
                if (estimatedFileSize == fileSize) {
                    headerInfo.headerSize = headerSize;
                    headerInfo.version = 0;

                    answer = true;

                    log.debug(TlUtils::format("ABGD(rev.1) check filesize: (%d, %d) [%d/%d] @%d",
                                              numOfVectors, sizeOfVector, subunitId, numOfSubunits, sizeOfChunk));
                    std::cerr << "This file format is obsolete. Please update the format." << std::endl;
                }
            }

            //
            if (answer == false) {
                IndexType numOfVectors = 0;
                IndexType sizeOfVector = 0;
                IndexType reservedSizeOfVector = 0;
                int numOfSubunits = 0;
                int subunitId = 0;
                int sizeOfChunk = 0;

                s.seekg(pos_type);
                s.read(reinterpret_cast<char*>(&numOfVectors), sizeof(IndexType));
                s.read(reinterpret_cast<char*>(&sizeOfVector), sizeof(IndexType));
                s.read(reinterpret_cast<char*>(&reservedSizeOfVector), sizeof(IndexType));

                s.read(reinterpret_cast<char*>(&numOfSubunits), sizeof(int));
                s.read(reinterpret_cast<char*>(&subunitId), sizeof(int));
                s.read(reinterpret_cast<char*>(&sizeOfChunk), sizeof(int));
                log.debug(TlUtils::format("ABGD(rev.2) (%d, %d, %d) [%d/%d] @%d",
                                          numOfVectors, sizeOfVector, reservedSizeOfVector, subunitId, numOfSubunits, sizeOfChunk));

                headerInfo.numOfVectors = numOfVectors;
                headerInfo.sizeOfVector = sizeOfVector;
                headerInfo.reservedSizeOfVector = reservedSizeOfVector;
                headerInfo.numOfSubunits = numOfSubunits;
                headerInfo.subunitId = subunitId;
                headerInfo.sizeOfChunk = sizeOfChunk;

                headerSize = sizeof(MatrixType) + sizeof(IndexType) * 3 + sizeof(int) * 3;
                const FileSize estimatedFileSize = headerSize + TlMatrixUtils::estimateFileSize(headerInfo);
                log.debug(TlUtils::format("ABGD(rev.2) check filesize: %ld (%ld) = %dx%d", fileSize, estimatedFileSize, numOfVectors, sizeOfVector));
                if (estimatedFileSize == fileSize) {
                    headerInfo.headerSize = headerSize;
                    headerInfo.version = 1;

                    answer = true;
                }
            }
        } break;

        default: {
            log.debug(TlUtils::format("unknown matrix file format. %s@%d", __FILE__, __LINE__));
        } break;
    }

    if (pHeaderInfo != NULL) {
        *pHeaderInfo = headerInfo;
    }

    // std::cerr << TlUtils::format("%s@%d [%d, %d] type: %d, fs: %ld/%ld %d/%d/%d", __FILE__, __LINE__,
    //                              headerInfo.numOfRows, headerInfo.numOfCols, headerInfo.matrixType, fileSize,
    //                              estimatedFileSize, headerInfo.numOfSubunits, headerInfo.subunitId,
    //                              headerInfo.sizeOfChunk)
    //           << std::endl;

    return answer;
}

std::size_t TlMatrixUtils::estimateFileSize(const TlMatrixObject::HeaderInfo& headerInfo) {
    TlLogging& log = TlLogging::getInstance();

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
            log.debug(TlUtils::format("TlMatrixUtils::estimateFileSize(): numOfVectors: %d, subunits: %d, sizeOfChunk: %d",
                                      headerInfo.numOfVectors, headerInfo.numOfSubunits, headerInfo.sizeOfChunk));
            const std::size_t numOfLocalChunks =
                TlDenseMatrix_arrays_Object::getNumOfLocalChunks(headerInfo.numOfVectors, headerInfo.numOfSubunits, headerInfo.sizeOfChunk);
            log.debug(TlUtils::format("TlMatrixUtils::estimateFileSize: chunks: %d", numOfLocalChunks));
            answer = size_double * (numOfLocalChunks * headerInfo.sizeOfChunk) * headerInfo.reservedSizeOfVector;
        } break;

        default:
            log.critical(TlUtils::format("program error: matrix_type=%d %s@%d", static_cast<int>(headerInfo.matrixType),
                                         __FILE__, __LINE__));
            break;
    }

    return answer;
}
