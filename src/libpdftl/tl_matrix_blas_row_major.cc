#include "tl_matrix_blas_row_major.h"
#include "tl_matrix_utils.h"

TlMatrixBlasRowMajor::TlMatrixBlasRowMajor(const TlMatrixObject::index_type row,
                                           const TlMatrixObject::index_type col)
    : TlDenseGeneralMatrix_BLAS_old(row, col) {}

TlMatrixBlasRowMajor::~TlMatrixBlasRowMajor() {}

bool TlMatrixBlasRowMajor::save(const std::string& filePath) const {
    bool answer = true;

    std::ofstream ofs;
    ofs.open(filePath.c_str(), std::ofstream::out | std::ofstream::binary);

    const char nType = static_cast<char>(TlMatrixObject::RSFD);
    ofs.write(&nType, sizeof(char));
    ofs.write(reinterpret_cast<const char*>(&this->row_), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&this->col_), sizeof(index_type));

    std::vector<double> buf(this->getNumOfElements());
    TlMatrixUtils::CSFD2RSFD(this->getNumOfRows(), this->getNumOfCols(),
                             this->data_, &(buf[0]));
    ofs.write(reinterpret_cast<char*>(&(buf[0])),
              this->getNumOfElements() * sizeof(double));

    ofs.close();

    return answer;
}
