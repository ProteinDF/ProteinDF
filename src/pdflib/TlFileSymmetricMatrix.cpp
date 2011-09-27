#include <iostream>
#include <iomanip>
#include "TlFileSymmetricMatrix.h"
#include "TlFile.h"

TlFileSymmetricMatrix::TlFileSymmetricMatrix(const std::string& filePath, const int dim, const size_t cacheSize)
        : TlFileMatrix(filePath, dim, dim, false, cacheSize)
{
    // TlFileMatri::xopen() is not called, because base constructer was called using (initialize = false)!
    this->open();
}

TlFileSymmetricMatrix::~TlFileSymmetricMatrix()
{
    // fs_ is closed by parent class.
    //std::cerr << "TlFileSymmetricMatrix::~TlFileSymmetricMatrix() called." << std::endl;
}

void TlFileSymmetricMatrix::open()
{
    if (TlFile::isExist(this->filePath_) == false) {
        // create new file
        this->fs_.open(this->filePath_.c_str(), std::ios::binary | std::ios::trunc | std::ios::out);

        // バッファ機能を停止する
        this->fs_ << std::setiosflags(std::ios::unitbuf);

        const int nType = 2; // means RLHD
        this->fs_.write(reinterpret_cast<const char*>(&nType), sizeof(int));
        this->fs_.write(reinterpret_cast<const char*>(&this->numOfRows_), sizeof(int));
        this->fs_.write(reinterpret_cast<const char*>(&this->numOfCols_), sizeof(int));

        const std::size_t size = std::size_t(this->numOfRows_) * std::size_t(this->numOfCols_ +1) / std::size_t(2);
        // size分ファイルにzeroを埋める
        const long unitSize = this->cacheSize_ / sizeof(double);
        ldiv_t ldivt = ldiv(size, unitSize);
        if (ldivt.quot > 0) {
            double* pBuf = new double[unitSize];
            for (int i = 0; i < unitSize; ++i) {
                pBuf[i] = 0.0;
            }

            for (long i = 0; i < ldivt.quot; ++i) {
                this->fs_.write(reinterpret_cast<const char*>(pBuf), sizeof(double) * unitSize);
            }
            delete[] pBuf;
            pBuf = NULL;
        }
        if (ldivt.rem > 0) {
            double* pBuf = new double[ldivt.rem];
            for (int i = 0; i < ldivt.rem; ++i) {
                pBuf[i] = 0.0;
            }

            this->fs_.write(reinterpret_cast<const char*>(pBuf), sizeof(double) * ldivt.rem);
            delete[] pBuf;
            pBuf = NULL;
        }

        // バッファ機能を再開する
        this->fs_ << std::resetiosflags(std::ios::unitbuf);
        this->fs_.flush();

        this->fs_.close();
    }

    this->fs_.open(this->filePath_.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    this->readHeader();
}

bool TlFileSymmetricMatrix::readHeader()
{
    bool bAnswer = true;

    // get file size
    std::fstream::pos_type nFileSize = 0;
    {
        this->fs_.seekg(0, std::ios_base::beg);
        std::fstream::pos_type begin = this->fs_.tellg();
        this->fs_.seekg(0, std::ios_base::end);
        //std::fstream::pos_type end = this->fs_.tellg();
        this->endPos_ = this->fs_.tellg();

        nFileSize = this->endPos_ - begin;
//     std::cerr << "TlFileSymmetricMatrix::readHeader() b=" << begin << ", e="  << this->endPos_ << std::endl;
    }

//   size_t sz = 0;
//   {
//     FILE* fp = fopen(this->filePath_.c_str(), "rb");
//     fseek(fp, 0, SEEK_END);
//     sz = ftell(fp);
//   }


    bool bCheckVariableType = false;
    this->startPos_ = 0;

    // int case:
    std::fstream::pos_type nEstimatedFileSize_int = 0;
    {
        int nType = 0;
        int nRows = 0;
        int nCols = 0;
        this->fs_.seekg(0, std::ios_base::beg);
        this->fs_.read((char*)&nType, sizeof(int));
        this->fs_.read((char*)&nRows, sizeof(int));
        this->fs_.read((char*)&nCols, sizeof(int));

        nEstimatedFileSize_int = std::ifstream::pos_type(sizeof(int) *3)
            + (std::ifstream::pos_type(nRows) * std::ifstream::pos_type(nCols +1)
               * std::ifstream::pos_type(sizeof(double) / 2));
        if ((nRows == nCols) && (nEstimatedFileSize_int == nFileSize)) {
            this->numOfRows_ = nRows;
            this->numOfCols_ = nCols;
            bCheckVariableType = true;
            this->startPos_ = sizeof(int) * 3;
        }
    }

    // long case:
    std::ifstream::pos_type nEstimatedFileSize_long = 0;
    {
        int nType = 0;
        long nRows = 0;
        long nCols = 0;
        this->fs_.seekg(0, std::ios_base::beg);
        this->fs_.read((char*)&nType, sizeof(int));
        this->fs_.read((char*)&nRows, sizeof(long));
        this->fs_.read((char*)&nCols, sizeof(long));

        nEstimatedFileSize_long = std::ifstream::pos_type(sizeof(int) + sizeof(long) * 2)
            + (std::ifstream::pos_type(nRows) * std::ifstream::pos_type(nCols +1)
               * std::ifstream::pos_type(sizeof(double) / 2));
        if ((nRows == nCols) && (nEstimatedFileSize_long == nFileSize)) {
            this->numOfRows_ = nRows;
            this->numOfCols_ = nCols;
            bCheckVariableType = true;
            this->startPos_ = sizeof(int) + sizeof(long) * 2;
        }
    }

    if (bCheckVariableType == true) {
        this->fs_.seekg(this->startPos_, std::ios_base::beg);
    } else {
        std::cerr << "file size mismatch: " << this->filePath_ << std::endl;
//        << nFileSize << ", i="
//        << nEstimatedFileSize_int << ", l="
//        << nEstimatedFileSize_long
//        << ", file=" << this->filePath_
//        << std::endl;
        bAnswer = false;
//     std::abort();
    }

    return bAnswer;
}

void TlFileSymmetricMatrix::add(const int row, const int col, const double value)
{
    TlFileMatrix::add(row, col, value);
}

void TlFileSymmetricMatrix::add(const TlPartialSymmetricMatrix& psm)
{
    const int startRow = psm.getStartRow();
    const int startCol = psm.getStartCol();
    const int maxRow = psm.getRowRange();
    const int maxCol = psm.getColRange();

    for (int r = 0; r < maxRow; ++r) {
        const int globalRow = startRow + r;
        for (int c = 0; c < maxCol; ++c) {
            const int globalCol = startCol + c;

            // for half matrix
            if (globalRow < globalCol) {
                break;
            }

            TlFileMatrix::add(globalRow, globalCol, psm.getLocal(r, c));
        }
    }
}


TlPartialSymmetricMatrix TlFileSymmetricMatrix::getPartialMatrix(const int startRow, const int startCol,
                                                                 const int range) const
{
    TlPartialSymmetricMatrix matrix(this->getNumOfRows(), startRow, startCol, range);

    const int localRows = matrix.getRowRange();
    const int localCols = matrix.getColRange();
    for (int r = 0; r < localRows; ++r) {
        const int globalRow = startRow + r;
        for (int c = 0; c < localCols; ++c) {
            const int globalCol = startCol + c;

            matrix.set(globalRow, globalCol, this->get(globalRow, globalCol));
        }
    }

    return matrix;
}


TlFileSymmetricMatrix& TlFileSymmetricMatrix::operator*=(const double coef)
{
    // ToDo: 高速化
    const size_type numOfRows = this->getNumOfRows();
    const size_type numOfCols = this->getNumOfCols();
    for (size_type r = 0; r < numOfRows; ++r) {
        for (size_type c = 0; c < numOfCols; ++c) {
            const double value = this->get(r, c) * coef;
            this->set(r, c, value);
        }
    }

    return *this;
}
