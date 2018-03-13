#include "tl_dense_general_matrix_object.h"
#include <fstream>
#include <iostream>
#include "tl_matrix_utils.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseGeneralMatrixObject::TlDenseGeneralMatrixObject(
    TlDenseMatrix_ImplObject* pImpl)
    : pImpl_(pImpl) {}

TlDenseGeneralMatrixObject::~TlDenseGeneralMatrixObject() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::index_type TlDenseGeneralMatrixObject::getNumOfRows() const {
  return this->pImpl_->getNumOfRows();
}

TlMatrixObject::index_type TlDenseGeneralMatrixObject::getNumOfCols() const {
  return this->pImpl_->getNumOfCols();
}

void TlDenseGeneralMatrixObject::resize(const TlMatrixObject::index_type row,
                                        const TlMatrixObject::index_type col) {
  this->pImpl_->resize(row, col);
}

double TlDenseGeneralMatrixObject::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
  return this->pImpl_->get(row, col);
}

void TlDenseGeneralMatrixObject::set(const TlMatrixObject::index_type row,
                                     const TlMatrixObject::index_type col,
                                     const double value) {
  this->pImpl_->set(row, col, value);
}

void TlDenseGeneralMatrixObject::add(const TlMatrixObject::index_type row,
                                     const TlMatrixObject::index_type col,
                                     const double value) {
  this->pImpl_->add(row, col, value);
}

// ---------------------------------------------------------------------------
// Operations
// ---------------------------------------------------------------------------
double TlDenseGeneralMatrixObject::sum() const { return this->pImpl_->sum(); }

double TlDenseGeneralMatrixObject::trace() const {
  return this->pImpl_->trace();
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseGeneralMatrixObject::load(const std::string& filePath) {
  bool answer = false;
  MatrixType matrixType;
  TlMatrixObject::index_type row;
  TlMatrixObject::index_type col;

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
  } else {
    this->log_.critical(TlUtils::format("illegal matrix format: %s @(%s:%d)",
                                        filePath.c_str(), __FILE__, __LINE__));
  }

  return answer;
}

bool TlDenseGeneralMatrixObject::save(const std::string& filePath) const {
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
  } else {
    this->log_.critical(TlUtils::format("cannot write matrix: %s @(%s:%d)",
                                        filePath.c_str(), __FILE__, __LINE__));
  }
  fs.close();

  return answer;
}

#ifdef HAVE_HDF5
bool TlDenseGeneralMatrixObject::loadHdf5(const std::string& filepath,
                                          const std::string& h5path) {
  TlHdf5Utils h5(filepath);

  int matrixType;
  h5.getAttr(h5path, "type", &matrixType);

  TlMatrixObject::index_type row = 0;
  TlMatrixObject::index_type col = 0;
  h5.getAttr(h5path, "row", &row);
  h5.getAttr(h5path, "col", &col);
  this->resize(row, col);

  const TlMatrixObject::size_type numOfElements =
      this->getNumOfRows() * this->getNumOfCols();
  switch (matrixType) {
    case TlMatrixObject::RSFD: {
      std::vector<double> buf(numOfElements);
      h5.get(h5path, &(buf[0]), numOfElements);
      std::size_t i = 0;
      for (TlMatrixObject::index_type r = 0; r < row; ++r) {
        for (TlMatrixObject::index_type c = 0; c < col; ++c) {
          this->set(r, c, buf[i]);
          ++i;
        }
      }
    } break;

    case TlMatrixObject::CSFD: {
      std::vector<double> buf(numOfElements);
      h5.get(h5path, &(buf[0]), numOfElements);
      std::size_t i = 0;
      for (TlMatrixObject::index_type c = 0; c < col; ++c) {
        for (TlMatrixObject::index_type r = 0; r < row; ++r) {
          this->set(r, c, buf[i]);
          ++i;
        }
      }
    } break;

    default:
      this->log_.critical(TlUtils::format("illegal matrix type: %d (%d@%s)",
                                          matrixType, __LINE__, __FILE__));
      break;
  }

  return true;
}

bool TlDenseGeneralMatrixObject::saveHdf5(const std::string& filepath,
                                          const std::string& h5path) const {
  TlHdf5Utils h5(filepath);

  const TlMatrixObject::index_type row = this->getNumOfRows();
  const TlMatrixObject::index_type col = this->getNumOfCols();
  const TlMatrixObject::size_type numOfElements = row * col;
  std::vector<double> buf(numOfElements);
  std::size_t i = 0;
  for (TlMatrixObject::index_type c = 0; c < col; ++c) {
    for (TlMatrixObject::index_type r = 0; r < row; ++r) {
      buf[i] = this->get(r, c);
      ++i;
    }
  }

  h5.write(h5path, &(buf[0]), numOfElements);
  h5.setAttr(h5path, "type", static_cast<int>(TlMatrixObject::CSFD));
  h5.setAttr(h5path, "row", this->getNumOfRows());
  h5.setAttr(h5path, "col", this->getNumOfCols());

  return true;
}

#endif  // HAVE_HDF5

// -----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream,
                         const TlDenseGeneralMatrixObject& mat) {
  const TlMatrixObject::index_type numOfRows = mat.getNumOfRows();
  const TlMatrixObject::index_type numOfCols = mat.getNumOfCols();

  for (TlMatrixObject::index_type ord = 0; ord < numOfCols; ord += 10) {
    stream << "       ";
    for (TlMatrixObject::index_type j = ord;
         ((j < ord + 10) && (j < numOfCols)); ++j) {
      stream << TlUtils::format("   %5d th", j + 1);
    }
    stream << "\n ----";

    for (TlMatrixObject::index_type j = ord;
         ((j < ord + 10) && (j < numOfCols)); ++j) {
      stream << "-----------";
    }
    stream << "----\n";

    for (TlMatrixObject::index_type i = 0; i < numOfRows; ++i) {
      stream << TlUtils::format(" %5d  ", i + 1);

      for (TlMatrixObject::index_type j = ord;
           ((j < ord + 10) && (j < numOfCols)); ++j) {
        stream << TlUtils::format(" %10.6lf", mat.get(i, j));
      }
      stream << "\n";
    }
    stream << "\n\n";
  }

  return stream;
}
