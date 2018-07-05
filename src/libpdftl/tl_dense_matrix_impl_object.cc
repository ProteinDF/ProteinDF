#include <cmath>

#include "tl_dense_matrix_impl_object.h"
#include "tl_dense_vector_object.h"
#include "TlUtils.h"

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
std::vector<double> TlDenseMatrix_ImplObject::diagonals() const {
  const TlMatrixObject::index_type dim =
      std::min(this->getNumOfRows(), this->getNumOfCols());

  std::vector<double> answer(dim, 0.0);
  for (TlMatrixObject::index_type i = 0; i < dim; ++i) {
    answer[i] = this->get(i, i);
  }

  return answer;
}

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
  for (TlMatrixObject::index_type i = 0; i < dim; ++i) {
    answer += this->get(i, i);
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

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
void TlDenseMatrix_ImplObject::dump(double* buf, const std::size_t size) const {
  const TlMatrixObject::index_type rows = this->getNumOfRows();
  const TlMatrixObject::index_type cols = this->getNumOfCols();

  std::size_t count = 0;
  for (TlMatrixObject::index_type r = 0; r < rows; ++r) {
    for (TlMatrixObject::index_type c = 0; c < cols; ++c) {
      if (count > size) {
        return;
      }

      buf[count] = this->get(r, c);
      ++count;
    }
  }
}

void TlDenseMatrix_ImplObject::restore(const double* buf,
                                       const std::size_t size) {
  const TlMatrixObject::index_type rows = this->getNumOfRows();
  const TlMatrixObject::index_type cols = this->getNumOfCols();

  std::size_t count = 0;
  for (TlMatrixObject::index_type r = 0; r < rows; ++r) {
    for (TlMatrixObject::index_type c = 0; c < cols; ++c) {
      if (count > size) {
        return;
      }

      this->set(r, c, buf[count]);
      ++count;
    }
  }
}

// -----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream,
                         const TlDenseMatrix_ImplObject& mat) {
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
