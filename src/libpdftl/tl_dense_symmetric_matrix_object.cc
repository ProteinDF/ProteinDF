#include "tl_dense_symmetric_matrix_object.h"
#include <fstream>
#include <iostream>
#include "tl_matrix_utils.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

TlDenseSymmetricMatrixObject::TlDenseSymmetricMatrixObject(
    TlDenseMatrix_ImplObject* pImpl)
    : pImpl_(pImpl) {}

TlDenseSymmetricMatrixObject::~TlDenseSymmetricMatrixObject() {}

TlMatrixObject::index_type TlDenseSymmetricMatrixObject::getNumOfRows() const {
  return this->pImpl_->getNumOfRows();
}

TlMatrixObject::index_type TlDenseSymmetricMatrixObject::getNumOfCols() const {
  return this->pImpl_->getNumOfCols();
}

void TlDenseSymmetricMatrixObject::resize(
    const TlMatrixObject::index_type dim) {
  // #ifndef NDEBUG
  //   if (newRow != newCol) {
  //     this->log_.warn(
  //         TlUtils::format("The size of the matrix is forced to be limited to
  //         the "
  //                         "number of rows.: (%d != %d)",
  //                         newRow, newCol));
  //   }
  // #endif
  this->pImpl_->resize(dim, dim);
}

double TlDenseSymmetricMatrixObject::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
  return this->pImpl_->get(row, col);
}

void TlDenseSymmetricMatrixObject::set(const TlMatrixObject::index_type row,
                                       const TlMatrixObject::index_type col,
                                       const double value) {
  this->pImpl_->set(row, col, value);
}

void TlDenseSymmetricMatrixObject::add(const TlMatrixObject::index_type row,
                                       const TlMatrixObject::index_type col,
                                       const double value) {
  this->pImpl_->add(row, col, value);
}

// ---------------------------------------------------------------------------
// Operations
// ---------------------------------------------------------------------------
double TlDenseSymmetricMatrixObject::sum() const { return this->pImpl_->sum(); }

double TlDenseSymmetricMatrixObject::trace() const {
  return this->pImpl_->trace();
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseSymmetricMatrixObject::load(const std::string& filePath) {
  bool answer = false;
  MatrixType matrixType;
  TlMatrixObject::index_type row;
  TlMatrixObject::index_type col;

  const TlMatrixUtils::FileSize headerSize =
      TlMatrixUtils::getHeaderInfo(filePath, &matrixType, &row, &col);
  if (headerSize > 0) {
    if (row != col) {
      this->log_.critical(
          TlUtils::format("illegal format: @%s.%d", __FILE__, __LINE__));
    }
    this->resize(row);

    std::fstream fs;
    fs.open(filePath.c_str(), std::ios::in | std::ios::binary);
    if (!fs.fail()) {
      fs.seekg(headerSize);

      switch (matrixType) {
        case TlMatrixObject::RLHD: {
          double v;
          for (TlMatrixObject::index_type r = 0; r < row; ++r) {
            for (TlMatrixObject::index_type c = 0; c <= r; ++c) {
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
  } else {
    this->log_.critical(TlUtils::format("illegal matrix format: %s @(%s:%d)",
                                        filePath.c_str(), __FILE__, __LINE__));
  }

  return answer;
}

bool TlDenseSymmetricMatrixObject::save(const std::string& filePath) const {
  bool answer = false;
  std::ofstream fs;
  fs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);
  if (!fs.fail()) {
    const char nType = static_cast<char>(TlMatrixObject::RLHD);
    const TlMatrixObject::index_type dim = this->getNumOfRows();

    fs.write(&nType, sizeof(char));
    fs.write(reinterpret_cast<const char*>(&dim),
             sizeof(TlMatrixObject::index_type));
    fs.write(reinterpret_cast<const char*>(&dim),
             sizeof(TlMatrixObject::index_type));

    for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
      for (TlMatrixObject::index_type c = 0; c <= r; ++c) {
        const double v = this->get(r, c);
        fs.write(reinterpret_cast<const char*>(&v), sizeof(double));
      }
    }
    fs.flush();
    answer = true;
  } else {
    this->log_.critical(TlUtils::format("cannot write matrix: %s @(%s:%d)",
                                        filePath.c_str(), __FILE__, __LINE__));
  }
  fs.close();

  return answer;
}

#ifdef HAVE_HDF5
bool TlDenseSymmetricMatrixObject::loadHdf5(const std::string& filepath,
                                            const std::string& h5path) {
  TlHdf5Utils h5(filepath);

  int mat_type;
  h5.getAttr(h5path, "type", &mat_type);

  index_type row = 0;
  index_type col = 0;
  h5.getAttr(h5path, "row", &row);
  h5.getAttr(h5path, "col", &col);
  if (row != col) {
    this->log_.critical(
        TlUtils::format("illegal parameter: row(%d) != col(%d) @%s:%d", row,
                        col, __FILE__, __LINE__));
  }
  this->resize(row);
  const TlMatrixObject::index_type dim = row;

  // operate matrix elements
  const std::size_t numOfElements = this->getNumOfElements_RLHD();
  std::vector<double> buf(numOfElements);
  h5.get(h5path, &(buf[0]), numOfElements);

  switch (mat_type) {
    case RLHD: {
      std::size_t i = 0;
      for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
        for (TlMatrixObject::index_type c = 0; c <= r; ++c) {
          this->set(r, c, buf[i]);
          i++;
        }
      }
    } break;
    case CUHD: {
      std::size_t i = 0;
      for (TlMatrixObject::index_type c = 0; c < dim; ++c) {
        for (TlMatrixObject::index_type r = 0; r <= c; ++r) {
          this->set(r, c, buf[i]);
          i++;
        }
      }
    } break;

    default:
      this->log_.critical(TlUtils::format(
          "illegal matrix type for TlDenseSymmetricMatrix_BLAS_Old: %d",
          mat_type));
      break;
  }

  return true;
}

bool TlDenseSymmetricMatrixObject::saveHdf5(const std::string& filepath,
                                            const std::string& h5path) const {
  TlHdf5Utils h5(filepath);

  const index_type row = this->getNumOfRows();
  const index_type col = this->getNumOfCols();
  const TlMatrixObject::index_type dim = row;
  assert(dim == col);

  // matrix elements
  {
    const std::size_t numOfElements = this->getNumOfElements_RLHD();
    std::vector<double> buf(numOfElements);
    std::size_t i = 0;
    for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
      for (TlMatrixObject::index_type c = 0; c <= r; ++c) {
        buf[i] = this->get(r, c);
        i++;
      }
    }

    h5.write(h5path, &(buf[0]), numOfElements);
  }

  h5.setAttr(h5path, "type", static_cast<int>(TlMatrixObject::RLHD));
  h5.setAttr(h5path, "row", row);
  h5.setAttr(h5path, "col", col);

  return true;
}

#endif  // HAVE_HDF5

std::ostream& operator<<(std::ostream& stream,
                         const TlDenseSymmetricMatrixObject& mat) {
  const TlMatrixObject::index_type nNumOfDim =
      mat.getNumOfRows();  // == this->getNumOfCols()

  stream << "\n\n";
  for (TlMatrixObject::index_type ord = 0; ord < nNumOfDim; ord += 10) {
    stream << "       ";
    for (TlMatrixObject::index_type j = ord;
         ((j < ord + 10) && (j < nNumOfDim)); ++j) {
      stream << TlUtils::format("   %5d th", j + 1);
    }
    stream << "\n"
           << " ----";

    for (TlMatrixObject::index_type j = ord;
         ((j < ord + 10) && (j < nNumOfDim)); ++j) {
      stream << "-----------";
    }
    stream << "----\n";

    for (TlMatrixObject::index_type i = 0; i < nNumOfDim; ++i) {
      stream << TlUtils::format(" %5d  ", i + 1);

      for (TlMatrixObject::index_type j = ord;
           ((j < ord + 10) && (j < nNumOfDim)); ++j) {
        if (j > i) {
          stream << "    ----   ";
        } else {
          stream << TlUtils::format(" %10.6lf", mat.get(i, j));
        }
      }
      stream << "\n";
    }
    stream << "\n\n";
  }

  return stream;
}
