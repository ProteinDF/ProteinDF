#ifndef TL_SYMMETRIC_MATRIX_OBJECT_H
#define TL_SYMMETRIC_MATRIX_OBJECT_H

template <class T>
class TlSymmetricMatrixObject : public T {
    template <typename StreamType>
    void print(StreamType& out) const {
        const TlMatrixObject::index_type numOfDim = this->getNumOfRows();
        assert(numOfDim == this->getNumOfCols());

        out << "\n\n";
        for (TlMatrixObject::index_type ord = 0; ord < numOfDim; ord += 10) {
            out << "       ";
            for (TlMatrixObject::index_type j = ord;
                 ((j < ord + 10) && (j < numOfDim)); ++j) {
                out << TlUtils::format("   %5d th", j + 1);
            }
            out << "\n"
                << " ----";

            for (TlMatrixObject::index_type j = ord;
                 ((j < ord + 10) && (j < numOfDim)); ++j) {
                out << "-----------";
            }
            out << "----\n";

            for (TlMatrixObject::index_type i = 0; i < numOfDim; ++i) {
                out << TlUtils::format(" %5d  ", i + 1);

                for (TlMatrixObject::index_type j = ord;
                     ((j < ord + 10) && (j < numOfDim)); ++j) {
                    if (j > i) {
                        out << "    ----   ";
                    } else {
                        out << TlUtils::format(" %10.6lf", this->get(i, j));
                    }
                }
                out << "\n";
            }
            out << "\n\n";
        }
        out.flush();
    }
};

#endif  // TL_SYMMETRIC_MATRIX_OBJECT_H
