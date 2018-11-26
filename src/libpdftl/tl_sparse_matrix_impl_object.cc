#include <cmath>

#include "tl_sparse_matrix_impl_object.h"
#include "TlUtils.h"

// -----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream,
                         const TlSparseMatrix_ImplObject& mat) {
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
