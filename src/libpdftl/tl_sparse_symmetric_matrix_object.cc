#include "tl_sparse_symmetric_matrix_object.h"
#include "TlUtils.h"
#include "tl_sparse_matrix_impl_object.h"

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlSparseSymmetricMatrixObject::TlSparseSymmetricMatrixObject(
    TlSparseMatrix_ImplObject* pImpl)
    : pImpl_(pImpl) {}

TlSparseSymmetricMatrixObject::~TlSparseSymmetricMatrixObject() {}

TlMatrixObject::index_type TlSparseSymmetricMatrixObject::getNumOfRows() const {
    return this->pImpl_->getNumOfRows();
}

TlMatrixObject::index_type TlSparseSymmetricMatrixObject::getNumOfCols() const {
    return this->pImpl_->getNumOfCols();
}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
void TlSparseSymmetricMatrixObject::resize(
    const TlMatrixObject::index_type dim) {
    this->pImpl_->resize(dim, dim);
}

double TlSparseSymmetricMatrixObject::get(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) const {
    return this->pImpl_->get(row, col);
}

void TlSparseSymmetricMatrixObject::set(const TlMatrixObject::index_type row,
                                        const TlMatrixObject::index_type col,
                                        const double value) {
    this->pImpl_->set(row, col, value);
}

void TlSparseSymmetricMatrixObject::add(const TlMatrixObject::index_type row,
                                        const TlMatrixObject::index_type col,
                                        const double value) {
    this->pImpl_->add(row, col, value);
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlSparseSymmetricMatrixObject::load(const std::string& filePath) {
    // TODO: implement
    return true;
}

bool TlSparseSymmetricMatrixObject::save(const std::string& filePath) const {
    // TODO: implement
    return true;
}

// ---------------------------------------------------------------------------
// others
// ---------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream,
                         const TlSparseSymmetricMatrixObject& mat) {
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
