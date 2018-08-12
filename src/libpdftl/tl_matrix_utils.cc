#include "tl_matrix_utils.h"
#include <iostream>
#include "TlUtils.h"
// #include "nullptr.h"

// -----------------------------------------------------------------------------
// public
// -----------------------------------------------------------------------------
bool TlMatrixUtils::isLoadable(const std::string& filePath,
                               const TlMatrixObject::MatrixType matrixType) {
  TlMatrixObject::MatrixType loadMatrixType;
  const TlMatrixUtils::FileSize headerSize =
      TlMatrixUtils::getHeaderInfo(filePath, &loadMatrixType);

  bool answer = false;
  if ((headerSize > 0) && (loadMatrixType == matrixType)) {
    answer = true;
  }
  return answer;
}

TlMatrixUtils::FileSize TlMatrixUtils::getHeaderInfo(
    const std::string& filepath, TlMatrixObject::MatrixType* pMatrixType,
    TlMatrixObject::index_type* pNumOfRows,
    TlMatrixObject::index_type* pNumOfCols,
    std::size_t* pNumOfItems) {
  FileSize headerSize = 0;

  std::ios_base::sync_with_stdio(false);
  std::fstream fs;
  fs.open(filepath.c_str(), std::ios::in | std::ios::binary);

  if (!fs.fail()) {
    headerSize =
        TlMatrixUtils::getHeaderInfo(fs, pMatrixType, pNumOfRows, pNumOfCols, pNumOfItems);
  } else {
    std::cerr << TlUtils::format("could not open file. %s @%s,%d",
                                 filepath.c_str(), __FILE__, __LINE__)
              << std::endl;
  }

  fs.close();

  return headerSize;
}

TlMatrixUtils::FileSize TlMatrixUtils::getHeaderInfo(
    std::fstream& fs, TlMatrixObject::MatrixType* pMatrixType,
    TlMatrixObject::index_type* pNumOfRows,
    TlMatrixObject::index_type* pNumOfCols,
    std::size_t* pNumOfItems) {
  const FileSize headerSize = TlMatrixUtils::getHeaderSize_templ1<std::fstream>(
      fs, pMatrixType, pNumOfRows, pNumOfCols, pNumOfItems);
  fs.seekg(headerSize, std::ios_base::beg);

  return headerSize;
}

bool TlMatrixUtils::saveMatrix(const std::string& filepath,
                               const TlMatrixObject::MatrixType matrixType,
                               const TlMatrixObject::index_type rows,
                               const TlMatrixObject::index_type cols,
                               const double* pBuf,
                               const std::size_t numOfItems) {
  const char nMatrixType = static_cast<char>(matrixType);
  std::ios_base::sync_with_stdio(false);

  std::fstream fs(filepath.c_str(), std::ios::out | std::ios::binary);
  fs.write(&nMatrixType, sizeof(char));
  fs.write(reinterpret_cast<const char*>(&rows),
           sizeof(TlMatrixObject::index_type));
  fs.write(reinterpret_cast<const char*>(&cols),
           sizeof(TlMatrixObject::index_type));
  fs.write(reinterpret_cast<const char*>(&pBuf), sizeof(double) * numOfItems);

  fs.close();

  return true;
}

void TlMatrixUtils::CSFD2RSFD(const TlMatrixObject::index_type row,
                              const TlMatrixObject::index_type col,
                              const double* pBufIn, double* pBufOut) {
  for (TlMatrixObject::index_type c = 0; c < col; ++c) {
    for (TlMatrixObject::index_type r = 0; r < row; ++r) {
      const TlMatrixObject::size_type rsfd_index = c + col * r;
      pBufOut[rsfd_index] = pBufIn[r + row * c];
    }
  }
}

void TlMatrixUtils::RSFD2CSFD(const TlMatrixObject::index_type row,
                              const TlMatrixObject::index_type col,
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
TlMatrixUtils::FileSize TlMatrixUtils::getHeaderSize_templ1(
    StreamType& s, TlMatrixObject::MatrixType* pMatrixType,
    TlMatrixObject::index_type* pNumOfRows,
    TlMatrixObject::index_type* pNumOfCols,
    std::size_t* pNumOfItems) {
  TlMatrixObject::MatrixType matrixType;
  TlMatrixObject::index_type rows = 0;
  TlMatrixObject::index_type cols = 0;
  std::size_t numOfItems = 0;

  FileSize headerSize = 0;
  // char case:
  {
    char type = 0;
    headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, char, int>(
        s, &type, &rows, &cols, &numOfItems);
    // std::cout << TlUtils::format("%s: %d (%d, %d)", "ci", type, rows, cols)
    //           << std::endl;

    if (headerSize == 0) {
      headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, char, long>(
          s, &type, &rows, &cols, &numOfItems);
      // std::cout << TlUtils::format("%s: %d (%d, %d)", "cl", type, rows, cols)
      //           << std::endl;
    }
    matrixType = static_cast<TlMatrixObject::MatrixType>(type);
  }

  // int case:
  if (headerSize == 0) {
    int type = 0;
    headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, int, int>(
        s, &type, &rows, &cols, &numOfItems);
    // std::cout << TlUtils::format("%s: %d (%d, %d)", "ii", type, rows, cols)
    //           << std::endl;

    if (headerSize == 0) {
      headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, int, long>(
          s, &type, &rows, &cols, &numOfItems);
      // std::cout << TlUtils::format("%s: %d (%d, %d)", "il", type, rows, cols)
      //           << std::endl;
    }
    matrixType = static_cast<TlMatrixObject::MatrixType>(type);
  }

  // long case:
  if (headerSize == 0) {
    long type = 0;
    headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, long, int>(
        s, &type, &rows, &cols, &numOfItems);
    // std::cout << TlUtils::format("%s: %d (%d, %d)", "li", type, rows, cols)
    //           << std::endl;

    if (headerSize == 0) {
      headerSize = TlMatrixUtils::getHeaderSize_templ2<StreamType, long, long>(
          s, &type, &rows, &cols, &numOfItems);
      // std::cout << TlUtils::format("%s: %d (%d, %d)", "ll", type, rows, cols)
      //           << std::endl;
    }
    matrixType = static_cast<TlMatrixObject::MatrixType>(type);
  }

  if (headerSize > 0) {
    if (pMatrixType != NULL) {
      *pMatrixType = matrixType;
    }
    if (pNumOfRows != NULL) {
      *pNumOfRows = rows;
    }
    if (pNumOfCols != NULL) {
      *pNumOfCols = cols;
    }
  }

  return headerSize;
}

template <typename StreamType, typename MatrixType, typename IndexType>
TlMatrixUtils::FileSize TlMatrixUtils::getHeaderSize_templ2(
    StreamType& s, MatrixType* pMatrixType,
    TlMatrixObject::index_type* pNumOfRows,
    TlMatrixObject::index_type* pNumOfCols,
    std::size_t* pNumOfItems) {
  MatrixType matrixType = 0;
  IndexType rows = 0;
  IndexType cols = 0;
  std::size_t numOfItems = 0;

  // get file size
  FileSize fileSize = 0;
  {
    s.seekg(0, std::ios_base::end);
    fileSize = s.tellg();
  }

  s.seekg(0, std::ios_base::beg);
  s.read((char*)&(matrixType), sizeof(MatrixType));
  s.read((char*)&(rows), sizeof(IndexType));
  s.read((char*)&(cols), sizeof(IndexType));

  FileSize headerSize = sizeof(MatrixType) + sizeof(IndexType) * 2;
  if ((matrixType == TlMatrixObject::COOF) || (matrixType == TlMatrixObject::COOS)) {
    s.read(reinterpret_cast<char*>(&numOfItems), sizeof(std::size_t));
    headerSize += sizeof(std::size_t);
  }

  const FileSize estimatedFileSize =
      headerSize +
      TlMatrixUtils::estimateFileSize(
          static_cast<TlMatrixObject::MatrixType>(matrixType), rows, cols, numOfItems);
  FileSize answer = 0;
  if (estimatedFileSize == fileSize) {
    answer = headerSize;
  }
  if (pMatrixType != NULL) {
    *pMatrixType = matrixType;
  }
  if (pNumOfRows != NULL) {
    *pNumOfRows = rows;
  }
  if (pNumOfCols != NULL) {
    *pNumOfCols = cols;
  }
  if (pNumOfItems != NULL) {
    *pNumOfItems = numOfItems;
  }
  return answer;
}

std::size_t TlMatrixUtils::estimateFileSize(
    TlMatrixObject::MatrixType matrixType, const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col, const std::size_t numOfItems) {
  std::size_t answer = 0;

  const std::size_t size_index = sizeof(TlMatrixObject::index_type);
  const std::size_t size_double = sizeof(double);

  switch (matrixType) {
    case TlMatrixObject::RSFD:
    case TlMatrixObject::CSFD:
      answer = size_double * row * col;
      break;

    case TlMatrixObject::RLHD:
    case TlMatrixObject::RUHD:
    case TlMatrixObject::CLHD:
    case TlMatrixObject::CUHD:
      answer = size_double * row * (col + 1) / 2;
      break;

    case TlMatrixObject::COOF:
    case TlMatrixObject::COOS:
      answer = (size_index *2 + size_double) * numOfItems;
      break;

    default:
      std::cerr << "program error. " << std::endl;
      break;
  }

  return answer;
}
