#include "TlUtils.h"
#include "tl_sparse_general_matrix_object.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseGeneralMatrixObject::TlSparseGeneralMatrixObject(
    TlSparseMatrix_ImplObject* pImpl)
    : pImpl_(pImpl) {}

TlSparseGeneralMatrixObject::~TlSparseGeneralMatrixObject() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlMatrixObject::index_type TlSparseGeneralMatrixObject::getNumOfRows() const {
  return this->pImpl_->getNumOfRows();
}

TlMatrixObject::index_type TlSparseGeneralMatrixObject::getNumOfCols() const {
  return this->pImpl_->getNumOfCols();
}

void TlSparseGeneralMatrixObject::resize(const TlMatrixObject::index_type row,
                                        const TlMatrixObject::index_type col) {
  this->pImpl_->resize(row, col);
}

double TlSparseGeneralMatrixObject::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
  return this->pImpl_->get(row, col);
}

void TlSparseGeneralMatrixObject::set(const TlMatrixObject::index_type row,
                                     const TlMatrixObject::index_type col,
                                     const double value) {
  this->pImpl_->set(row, col, value);
}

void TlSparseGeneralMatrixObject::add(const TlMatrixObject::index_type row,
                                     const TlMatrixObject::index_type col,
                                     const double value) {
  this->pImpl_->add(row, col, value);
}

void TlSparseGeneralMatrixObject::mul(const TlMatrixObject::index_type row,
                                     const TlMatrixObject::index_type col,
                                     const double value) {
  this->pImpl_->mul(row, col, value);
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlSparseGeneralMatrixObject::load(const std::string& filePath) {
  // TODO: implement
  return true;
}

bool TlSparseGeneralMatrixObject::save(const std::string& filePath) const {
  // TODO: implement
  return true;
}


// -----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream,
                         const TlSparseGeneralMatrixObject& mat) {
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
