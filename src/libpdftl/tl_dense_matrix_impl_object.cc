#include "tl_dense_matrix_impl_object.h"

double TlDenseMatrix_ImplObject::sum() const {
  const TlMatrixObject::index_type rows = this->getNumOfRows();
  const TlMatrixObject::index_type cols = this->getNumOfCols();

  double answer = 0.0;
  for (TlMatrixObject::index_type r = 0; r < rows; ++r) {
    for (TlMatrixObject::index_type c = 0; c < cols; ++c) {
      answer += this->get(r, c);
    }
  }

  return answer;
}

double TlDenseMatrix_ImplObject::trace() const {
  const TlMatrixObject::index_type dim =
      std::min(this->getNumOfRows(), this->getNumOfCols());

  double answer = 0.0;
  for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
    for (TlMatrixObject::index_type c = 0; c < dim; ++c) {
      answer += this->get(r, c);
    }
  }

  return answer;
}

double TlDenseMatrix_ImplObject::getRMS() const {
  const TlMatrixObject::index_type rows = this->getNumOfRows();
  const TlMatrixObject::index_type cols = this->getNumOfCols();

  double sum2 = 0.0;
  for (TlMatrixObject::index_type r = 0; r < rows; ++r) {
    for (TlMatrixObject::index_type c = 0; c < cols; ++c) {
      const double v = this->get(r, c);
      sum2 += v * v;
    }
  }

  const double elements = this->getNumOfRows() * this->getNumOfCols();
  const double rms = std::sqrt(sum2 / elements);

  return rms;
}

double TlDenseMatrix_ImplObject::getMaxAbsoluteElement(
    TlMatrixObject::index_type* outRow,
    TlMatrixObject::index_type* outCol) const {
  const TlMatrixObject::index_type rows = this->getNumOfRows();
  const TlMatrixObject::index_type cols = this->getNumOfCols();

  double answer = 0.0;
  TlMatrixObject::index_type maxRow = 0;
  TlMatrixObject::index_type maxCol = 0;
  for (TlMatrixObject::index_type r = 0; r < rows; ++r) {
    for (TlMatrixObject::index_type c = 0; c < cols; ++c) {
      const double value = std::fabs(this->get(r, c));
      if (value > answer) {
        maxRow = r;
        maxCol = c;
        answer = value;
      }
    }
  }

  if (outRow != NULL) {
    *outRow = maxRow;
  }
  if (outCol != NULL) {
    *outCol = maxCol;
  }

  return answer;
}
