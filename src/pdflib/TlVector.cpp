#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>

#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

TlVector::TlVector(size_type size)
        : size_(size), data_(NULL)
{
    this->initialize();
}


TlVector::TlVector(const double* p, size_type size)
        : size_(size), data_(NULL)
{
    this->initialize(false);
    std::copy(p, p +size, this->data_);
}


TlVector::TlVector(const TlVector& rhs)
        : size_(rhs.size_), data_(NULL)
{
    this->initialize(false);
    std::copy(rhs.data_, rhs.data_ + rhs.getSize(), this->data_);
}


TlVector::TlVector(const std::vector<double>& rhs)
        : size_(rhs.size()), data_(NULL)
{
    this->initialize(false);
    std::copy(rhs.begin(), rhs.end(), this->data_);
}


TlVector::~TlVector()
{
    this->destroy();
}


void TlVector::initialize(bool isZeroClear)
{
    const size_type size = this->getSize();
    this->data_ = new double[size];

    if (isZeroClear == true) {
        this->zeroClear();
    }
}


void TlVector::destroy()
{
    this->size_ = 0;
    delete[] this->data_;
    this->data_ = NULL;
}


void TlVector::resize(const size_type size)
{
    double* pOld = this->data_;
    const size_type oldSize = this->size_;

    this->size_ = size;
    this->initialize(true);

    const size_type fillSize = std::min(size, oldSize);
    std::copy(pOld, pOld + fillSize, this->data_);

    delete[] pOld;
    pOld = NULL;
}


void TlVector::push_back(const double value)
{
    const size_type size = this->getSize();
    this->resize(size +1);
    this->data_[size] = value;
}


double TlVector::getMaxAbsoluteElement() const
{
    double answer = 0.0;
    const size_type size = this->getSize();
    for (size_type i = 0; i < size; ++i) {
        answer = std::max(answer, std::fabs(this->get(i)));
    }

    return answer;
}


TlVector& TlVector::dot(const TlVector& rhs)
{
    assert(this->getSize() == rhs.getSize());
    std::transform(rhs.data_, rhs.data_ + rhs.getSize(), this->data_,
                   this->data_, std::multiplies<double>());

    return *this;
}


double TlVector::sum() const
{
    return std::accumulate(this->data_, this->data_ + this->getSize(), 0.0);
}


void TlVector::sortByGrater()
{
    //std::cout << "call TlVector::sortByGraterEqual()" << std::endl;
    std::sort(this->data_, this->data_ + this->getSize(), std::greater<double>());
}


TlVector& TlVector::operator=(const TlVector& rhs)
{
    if (this != &rhs) {
        //this->clear();
        delete[] this->data_;
        this->data_ = NULL;

        this->size_ = rhs.size_;
        this->initialize(false);

        std::copy(rhs.data_, rhs.data_ + this->getSize(), this->data_);
    }

    return (*this);
}


TlVector& TlVector::operator+=(const TlVector& rhs)
{
    assert(this->getSize() == rhs.getSize());

    const size_type size = this->getSize();
#pragma omp parallel for
    for (size_type i = 0; i < size; ++i) {
        this->data_[i] += rhs.get(i);
    }

    return (*this);
}


TlVector& TlVector::operator-=(const TlVector& rhs)
{
    assert(this->getSize() == rhs.getSize());

    const size_type size = this->getSize();
#pragma omp parallel for
    for (size_type i = 0; i < size; ++i) {
        this->data_[i] -= rhs.get(i);
    }

    return (*this);
}


TlVector& TlVector::operator*=(const double& rhs)
{
    const size_type size = this->getSize();
#pragma omp parallel for
    for (size_type i = 0; i < size; ++i) {
        this->data_[i] *= rhs;
    }

    return (*this);
}


TlVector& TlVector::operator/=(const double& rhs)
{
    return (this->operator*=(rhs));
}


TlVector operator+(const TlVector& X, const TlVector& Y)
{
    assert(X.getSize() == Y.getSize());

    TlVector answer = X;
    answer += Y;
    return answer;
}


TlVector operator-(const TlVector& X, const TlVector& Y)
{
    assert(X.getSize() == Y.getSize());

    TlVector answer = X;
    answer -= Y;
    return answer;
}


TlVector operator*(const TlVector& X, const double& Y)
{
    TlVector answer = X;
    answer *= Y;

    return answer;
}


double operator*(const TlVector&X, const TlVector& Y)
{
    assert(X.getSize() == Y.getSize());

    return std::inner_product(X.data_, X.data_ + X.getSize(),
                              Y.data_, 0.0);
}


bool TlVector::isLoadable(std::ifstream& ifs)
{
    if (ifs.is_open() != true) {
        return false;
    }

    ifs.seekg(0, std::ios::beg);
    int nSize = 0;
    const bool bAnswer = TlVector::getHeaderInfo(ifs, &nSize);
    ifs.seekg(0, std::ios::beg);

    return bAnswer;
}


bool TlVector::isLoadable(const std::string& rFilePath)
{
    std::ifstream ifs;
    ifs.open(rFilePath.c_str());
    if (ifs.fail()) {
#ifdef DEBUG
        std::cerr << "[error] TlVector::load(): could not open file. " << rFilePath << std::endl;
#endif // DEBUG
        return false;
    }

    const bool bAnswer = TlVector::isLoadable(ifs);
    ifs.close();
    return bAnswer;
}


bool TlVector::getHeaderInfo(std::ifstream& ifs, int* pnSize)
{
    assert(pnSize != NULL);
    bool bAnswer = true;
    // binary type =====================================================
    // get file size
    std::ifstream::pos_type nFileSize = 0;
    {
        ifs.seekg(0, std::ios_base::beg);
        std::ifstream::pos_type begin = ifs.tellg();
        ifs.seekg(0, std::ios_base::end);
        std::ifstream::pos_type end = ifs.tellg();

        nFileSize = end - begin;
        //std::cerr << "file size = " << nFileSize << " bytes" << std::endl;
    }

    // read header
    {
        bool bCheckVariableType = false;
        std::ifstream::pos_type nStartContentPos = 0;

        //   int case
        {
            int nTrySize = 0;
            ifs.seekg(0, std::ios_base::beg);
            ifs.read((char*)&nTrySize, sizeof(int));

            const std::ifstream::pos_type nEstimatedFileSize = sizeof(int) + (nTrySize * sizeof(double));
            if (nEstimatedFileSize == nFileSize) {
                *pnSize = nTrySize;
                bCheckVariableType = true;
                nStartContentPos = sizeof(int);
            }
        }

        //   long case
        {
            long nTrySize = 0;
            ifs.seekg(0, std::ios_base::beg);
            ifs.read((char*)&nTrySize, sizeof(long));

            const std::ifstream::pos_type nEstimatedFileSize = sizeof(long) + (nTrySize * sizeof(double));
            if (nEstimatedFileSize == nFileSize) {
                *pnSize = nTrySize;
                bCheckVariableType = true;
                nStartContentPos = sizeof(long);
            }
        }

        if (bCheckVariableType == true) {
            ifs.seekg(nStartContentPos, std::ios_base::beg);
        } else {
            bAnswer = false;
        }
    }

    return bAnswer;
}


bool TlVector::load(const std::string& sFilePath)
{
    std::ifstream ifs;

    ifs.open(sFilePath.c_str());
    if (ifs.fail()) {
#ifdef DEBUG
        std::cerr << "[error] TlVector::load(): could not open file. " << sFilePath << std::endl;
#endif // DEBUG
        return false;
    }

    bool bAnswer = this->load(ifs);
    ifs.close();

    if (bAnswer != true) {
        std::cerr << "TlVector::load() failed.: " << sFilePath << std::endl;
        return false;
    }

    return true;
}


bool TlVector::load(std::ifstream& ifs)
{
    int nSize = 0;
    bool bAnswer = TlVector::getHeaderInfo(ifs, &nSize);
    if (bAnswer == true) {
        this->destroy();
        this->size_ = nSize;
        this->initialize();
        ifs.read(reinterpret_cast<char*>(this->data_), sizeof(double) * nSize);
    }

    return bAnswer;
}


bool TlVector::loadText(const std::string& filePath)
{
    bool answer = false;

    std::ifstream ifs;
    ifs.open(filePath.c_str(), std::ios::in);
    if (ifs.fail()) {
#ifdef DEBUG
        std::cerr << "[error] TlVector::loadText(): could not open file. " << filePath << std::endl;
#endif //DEBUG
        return false;
    }

    // load contents
    std::string line = "";
    std::getline(ifs, line); // read 1st line

    if (line == "TEXT") {
        std::string tmp = "";
        ifs >> tmp;
        const int size = std::atoi(tmp.c_str());
        ifs >> tmp; // equal to 'size'
        ifs >> tmp; // equal to '0'
        this->resize(size);
        for (int i = 0; i < size; ++i) {
            ifs >> tmp;
            const double v = std::atof(tmp.c_str());
            (*this)[i] = v;
        }

        answer = true;
    }

    ifs.close();

    return answer;
}


bool TlVector::save(const std::string& sFilePath) const
{
    bool bAnswer = true;

    std::ofstream ofs;
    ofs.open(sFilePath.c_str(), std::ofstream::out | std::ofstream::binary);

    bAnswer = this->save(ofs);

    ofs.close();

    return bAnswer;
}


bool TlVector::save(std::ofstream& ofs) const
{
    bool bAnswer = true;

    const int nSize = this->getSize();

    ofs.write(reinterpret_cast<const char*>(&nSize), sizeof(int));

    const char* p = reinterpret_cast<const char*>(const_cast<TlVector*>(this)->data_);
    ofs.write(p, sizeof(double) * nSize);

    return bAnswer;
}


void TlVector::outputText(std::ostream& os) const
{
    const size_type nSize = this->getSize();

    os << "TEXT\n";
    os << nSize << "\n";
    os << nSize << "\n";
    os << "0\n";

    for (size_type j = 0; j < nSize; j += 10) {
        for (size_type i = j; ((i < j +10) && (i < nSize)); ++i) {
            os << TlUtils::format("  %10.4lf", this->get(i));
        }
        os << std::endl;
    }
    os << std::endl;
}
