#include "tl_dense_general_matrix_abstract.h"
#include <iostream>
#include "tl_matrix_utils.h"

bool TlDenseGeneralMatrixAbstract::load(const std::string& filePath) {
  bool answer = false;
  MatrixType matrixType;
  TlMatrixObject::index_type row, col;

  const TlMatrixUtils::FileSize headerSize =
      TlMatrixUtils::getHeaderInfo(filePath, &matrixType, &row, &col);
  if (headerSize > 0) {
    this->resize(row, col);

    std::fstream fs;
    fs.open(filePath.c_str(), std::ios::in | std::ios::binary);
    if (!fs.fail()) {
      fs.seekg(headerSize);

      switch (matrixType) {
        case TlMatrixObject::RSFD: {
          double v;
          for (TlMatrixObject::index_type r = 0; r < row; ++r) {
            for (TlMatrixObject::index_type c = 0; c < col; ++c) {
              fs.read(reinterpret_cast<char*>(&v), sizeof(double));
              this->set(r, c, v);
            }
          }
          answer = true;
        } break;

        case TlMatrixObject::CSFD: {
          double v;
          for (TlMatrixObject::index_type c = 0; c < col; ++c) {
            for (TlMatrixObject::index_type r = 0; r < row; ++r) {
              fs.read(reinterpret_cast<char*>(&v), sizeof(double));
              this->set(r, c, v);
            }
          }
          answer = true;
        } break;

        default:
          this->log_.critical(TlUtils::format("not supported format: @%s.%d",
                                              __FILE__, __LINE__));
          break;
      }
    } else {
      this->log_.critical(
          TlUtils::format("cannnot open matrix file: %s @(%s:%d)",
                          filePath.c_str(), __FILE__, __LINE__));
    }

    fs.close();
  }

  return answer;
}

bool TlDenseGeneralMatrixAbstract::save(const std::string& filePath) const {
  bool answer = false;
  std::ofstream fs;
  fs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);
  if (!fs.fail()) {
    const char nType = static_cast<char>(TlMatrixObject::CSFD);
    const TlMatrixObject::index_type row = this->getNumOfRows();
    const TlMatrixObject::index_type col = this->getNumOfCols();

    fs.write(&nType, sizeof(char));
    fs.write(reinterpret_cast<const char*>(&row),
             sizeof(TlMatrixObject::index_type));
    fs.write(reinterpret_cast<const char*>(&col),
             sizeof(TlMatrixObject::index_type));

    for (TlMatrixObject::index_type c = 0; c < col; ++c) {
        for (TlMatrixObject::index_type r = 0; r < row; ++r) {
        const double v = this->get(r, c);
        fs.write(reinterpret_cast<const char*>(&v), sizeof(double));
      }
    }
    fs.flush();

    answer = true;
  }
  fs.close();
}

bool TlDenseGeneralMatrixAbstract::load(const std::string& filePath, double* pBuf,
                                      const TlMatrixObject::size_type length) {
  bool answer = false;
  MatrixType matrixType;
  TlMatrixObject::index_type row, col;

  const TlMatrixUtils::FileSize headerSize =
      TlMatrixUtils::getHeaderInfo(filePath, &matrixType, &row, &col);

  if (headerSize > 0) {
    std::fstream fs;
    fs.open(filePath.c_str(), std::ios::in | std::ios::binary);

    if (!fs.fail()) {
      fs.seekg(headerSize);

      switch (matrixType) {
        case TlMatrixObject::CSFD: {
          fs.read(reinterpret_cast<char*>(pBuf), sizeof(double) * length);
          answer = true;
        } break;

        case TlMatrixObject::RSFD: {
          std::vector<double> tmpBuf(length);
          fs.read(reinterpret_cast<char*>(tmpBuf.data()),
                  sizeof(double) * length);
          TlMatrixUtils::RSFD2CSFD(row, col, tmpBuf.data(), pBuf);
          answer = true;
        } break;

        default:
          this->log_.critical(TlUtils::format("not supported format: @%s.%d",
                                              __FILE__, __LINE__));
          break;
      }
    } else {
      this->log_.critical(
          TlUtils::format("cannnot open matrix file: %s @(%s:%d)",
                          filePath.c_str(), __FILE__, __LINE__));
    }
    fs.close();
  }

  return answer;
}

bool TlDenseGeneralMatrixAbstract::save(
    const std::string& filePath, const double* pBuf,
    const TlMatrixObject::size_type length) const {
  bool answer = false;

  std::ofstream fs;
  fs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);
  if (!fs.fail()) {
    const char nType = static_cast<char>(TlMatrixObject::CSFD);
    const TlMatrixObject::index_type row = this->getNumOfRows();
    const TlMatrixObject::index_type col = this->getNumOfCols();

    fs.write(&nType, sizeof(char));
    fs.write(reinterpret_cast<const char*>(&row),
             sizeof(TlMatrixObject::index_type));
    fs.write(reinterpret_cast<const char*>(&col),
             sizeof(TlMatrixObject::index_type));
    fs.write(reinterpret_cast<const char*>(pBuf), length * sizeof(double));
    fs.flush();

    answer = true;
  }
  fs.close();

  return answer;
}
