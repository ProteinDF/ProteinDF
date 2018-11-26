#include "tl_dense_symmetric_matrix_abstract.h"
#include "tl_matrix_utils.h"
#include <cassert>
#include "TlUtils.h"

bool TlDenseSymmetricMatrixAbstract::load(
    const std::string& filePath, double* pBuf,
    const TlMatrixObject::size_type length) {
  bool answer = false;
  TlMatrixObject::MatrixType matrixType;
  TlMatrixObject::index_type row, col;

  const TlMatrixUtils::FileSize headerSize =
      TlMatrixUtils::getHeaderInfo(filePath, &matrixType, &row, &col);

  if (headerSize > 0) {
    std::ios_base::sync_with_stdio(false);
    std::fstream fs;
    fs.open(filePath.c_str(), std::fstream::binary | std::fstream::in);

    if (!fs.fail()) {
      fs.seekg(headerSize);

      switch (matrixType) {
        case TlMatrixObject::RLHD: {
          const TlMatrixObject::index_type dim = row;
          assert(dim == col);
          // this->resize(dim);

          double buf;
          for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
            for (TlMatrixObject::index_type c = 0; c <= r; ++c) {
              fs.read(reinterpret_cast<char*>(&buf), sizeof(double));
              this->set(r, c, buf);
            }
          }
        } break;

        default: {
          this->log_.critical(TlUtils::format("not supported format: @%s.%d",
                                              __FILE__, __LINE__));
        } break;
      }
      answer = true;
    } else {
      this->log_.critical(
          TlUtils::format("cannnot open matrix file: %s @(%s:%d)",
                          filePath.c_str(), __FILE__, __LINE__));
    }
    fs.close();
  }

  return answer;
}

bool TlDenseSymmetricMatrixAbstract::save(
    const std::string& filePath, const double* pBuf,
    const TlMatrixObject::size_type length) const {
  bool answer = false;

  std::ios_base::sync_with_stdio(false);
  std::ofstream fs;
  fs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);

  if (!fs.fail()) {
    const char nType = static_cast<char>(TlMatrixObject::RLHD);
    const TlMatrixObject::index_type row = this->getNumOfRows();
    const TlMatrixObject::index_type col = this->getNumOfCols();

    fs.write(&nType, sizeof(char));
    fs.write(reinterpret_cast<const char*>(&row),
             sizeof(TlMatrixObject::index_type));
    fs.write(reinterpret_cast<const char*>(&col),
             sizeof(TlMatrixObject::index_type));
    fs.write(reinterpret_cast<const char*>(pBuf), length * sizeof(double));

    fs.close();
    answer = true;
  }

  return answer;
}
