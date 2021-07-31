// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include "tl_dense_matrix_arrays_object.h"

#include <algorithm>
#include <cassert>
#include <iostream>

#include "TlLogging.h"
#include "TlMemManager.h"
#include "TlUtils.h"
#include "tl_assert.h"

#define DEFAULT_CHUNK_SIZE 32

TlDenseMatrix_arrays_Object::TlDenseMatrix_arrays_Object(const index_type numOfVectors, const index_type sizeOfVector,
                                                         const int numOfSubunits, const int subunitID,
                                                         const bool isUsingMemManager)
    : TlMatrixObject(RSFD),
      numOfVectors_(0),
      sizeOfVector_(0),
      sizeOfChunk_(DEFAULT_CHUNK_SIZE),
      numOfSubunits_(numOfSubunits),
      subunitID_(subunitID),
      reservedVectorSize_(0),
      isUsingMemManager_(isUsingMemManager),
      numOfLocalChunks_(0) {
    this->resize(numOfVectors, sizeOfVector);
}

TlDenseMatrix_arrays_Object::TlDenseMatrix_arrays_Object(const TlDenseMatrix_arrays_Object& rhs)
    : TlMatrixObject(RSFD),
      numOfVectors_(0),
      sizeOfVector_(0),
      sizeOfChunk_(rhs.sizeOfChunk_),
      numOfSubunits_(rhs.numOfSubunits_),
      subunitID_(rhs.subunitID_),
      reservedVectorSize_(0),
      isUsingMemManager_(rhs.isUsingMemManager_),
      numOfLocalChunks_(0) {
    this->resize(rhs.numOfVectors_, rhs.sizeOfVector_);

    // chunk
    const int numOfLocalChunks = this->numOfLocalChunks_;
    assert(this->chunks_.size() == std::size_t(numOfLocalChunks));
    const std::size_t sizeOfChunkVectors = this->reservedVectorSize_ * this->sizeOfChunk_;
    for (int i = 0; i < numOfLocalChunks; ++i) {
        std::copy(rhs.chunks_[i], rhs.chunks_[i] + sizeOfChunkVectors, this->chunks_[i]);
    }
}

TlDenseMatrix_arrays_Object::~TlDenseMatrix_arrays_Object() {
    this->destroy();
}

TlDenseMatrix_arrays_Object& TlDenseMatrix_arrays_Object::operator=(const TlDenseMatrix_arrays_Object& rhs) {
    if (this != &rhs) {
        this->destroy();

        this->numOfSubunits_ = rhs.numOfSubunits_;
        this->subunitID_ = rhs.subunitID_;
        this->isUsingMemManager_ = rhs.isUsingMemManager_;

        this->resize(rhs.numOfVectors_, rhs.sizeOfVector_);

        const index_type numOfLocalChunks = this->numOfLocalChunks_;
        assert(this->chunks_.size() == std::size_t(numOfLocalChunks));
        const std::size_t sizeOfChunkVectors = this->reservedVectorSize_ * this->sizeOfChunk_;
        for (index_type i = 0; i < numOfLocalChunks; ++i) {
            std::copy(rhs.chunks_[i], rhs.chunks_[i] + sizeOfChunkVectors, this->chunks_[i]);
        }
    }

    return *this;
}

int TlDenseMatrix_arrays_Object::getSizeOfChunk() const {
    return this->sizeOfChunk_;
}

void TlDenseMatrix_arrays_Object::resize(const index_type newNumOfVectors, const index_type newSizeOfVector)

{
    // evaluate chunks
    const int prevNumOfLocalChunks = this->numOfLocalChunks_;
    const int newNumOfLocalChunks =
        this->getNumOfLocalChunks(newNumOfVectors, this->numOfSubunits_, this->sizeOfChunk_);
    this->numOfLocalChunks_ = newNumOfLocalChunks;

    if (newNumOfLocalChunks > prevNumOfLocalChunks) {
        // ベクトル数を拡大
        const index_type reservedVectorSize = this->reservedVectorSize_;

        this->chunks_.resize(newNumOfLocalChunks);
        const std::size_t sizeOfChunkVectors = reservedVectorSize * this->sizeOfChunk_;
        for (int i = prevNumOfLocalChunks; i < newNumOfLocalChunks; ++i) {
            double* pNewChunk = NULL;
            if (this->isUsingMemManager_ == true) {
                TlMemManager& rMemManager = TlMemManager::getInstance();
                pNewChunk = reinterpret_cast<double*>(rMemManager.allocate(sizeof(double) * sizeOfChunkVectors));
            } else {
                pNewChunk = new double[sizeOfChunkVectors];
                std::fill(pNewChunk, pNewChunk + sizeOfChunkVectors, 0.0);
            }
            this->chunks_[i] = pNewChunk;
        }

    } else if (newNumOfLocalChunks < prevNumOfLocalChunks) {
        // ベクトル数を縮小
        if (this->isUsingMemManager_ == true) {
            TlMemManager& rMemManager = TlMemManager::getInstance();
            for (int i = newNumOfLocalChunks; i < prevNumOfLocalChunks; ++i) {
                rMemManager.deallocate((char*)this->chunks_[i]);
                this->chunks_[i] = NULL;
            }
        } else {
            for (int i = newNumOfLocalChunks; i < prevNumOfLocalChunks; ++i) {
                delete[] this->chunks_[i];
                this->chunks_[i] = NULL;
            }
        }
        this->chunks_.resize(newNumOfLocalChunks);
    }

    this->numOfVectors_ = newNumOfVectors;
    this->reserveVectorSize(newSizeOfVector);
    this->sizeOfVector_ = newSizeOfVector;
}

void TlDenseMatrix_arrays_Object::reserveVectorSize(index_type newReservedVectorSize) {
    newReservedVectorSize = std::max(this->sizeOfVector_, newReservedVectorSize);
    const index_type prevReservedVectorSize = this->reservedVectorSize_;

    if (prevReservedVectorSize < newReservedVectorSize) {
        const int numOfLocalChunks = this->numOfLocalChunks_;
        assert(this->chunks_.size() == static_cast<std::size_t>(numOfLocalChunks));
        const std::size_t sizeOfChunkVectors = newReservedVectorSize * this->sizeOfChunk_;
        for (int i = 0; i < numOfLocalChunks; ++i) {
            // allocate new object
            double* pNew = NULL;
            if (this->isUsingMemManager_ == true) {
                TlMemManager& rMemManager = TlMemManager::getInstance();
                pNew = reinterpret_cast<double*>(rMemManager.allocate(sizeof(double) * sizeOfChunkVectors));
            } else {
                pNew = new double[sizeOfChunkVectors];
            }
            std::fill(pNew, pNew + sizeOfChunkVectors, 0.0);

            if (this->chunks_[i] != NULL) {
                // copy old data to new
                const int sizeOfChunk = this->sizeOfChunk_;
                for (int v = 0; v < sizeOfChunk; ++v) {
                    for (int u = 0; u < prevReservedVectorSize; ++u) {
                        const std::size_t prevIndex = v * prevReservedVectorSize + u;
                        const std::size_t newIndex = v * newReservedVectorSize + u;
                        assert(newIndex < sizeOfChunkVectors);
                        pNew[newIndex] = this->chunks_[i][prevIndex];
                    }
                }

                // destroy old object
                if (this->isUsingMemManager_ == true) {
                    TlMemManager& rMemManager = TlMemManager::getInstance();
                    rMemManager.deallocate((char*)this->chunks_[i]);
                } else {
                    delete[] this->chunks_[i];
                }
                this->chunks_[i] = NULL;
            }
            this->chunks_[i] = pNew;
        }

        this->reservedVectorSize_ = newReservedVectorSize;
        this->log_.debug(
            TlUtils::format("reserveVectorSize() this->reservedVectorSize_=%d", this->reservedVectorSize_));
    }
}

int TlDenseMatrix_arrays_Object::getSubunitID(const index_type vectorIndex) const {
    assert((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_));
    int subunitId = -1;
    this->getLocalVectorIndex(vectorIndex, &subunitId);

    return subunitId;
}

TlMatrixObject::index_type TlDenseMatrix_arrays_Object::getNumOfLocalChunks(const index_type numOfVectors,
                                                                            const int numOfSubunits,
                                                                            const int sizeOfChunk) {
    const int numOfVectorsPerUnit = sizeOfChunk * numOfSubunits;
    const int numOfLocalChunks = numOfVectors / numOfVectorsPerUnit + 1;

    return numOfLocalChunks;
}

// defunct
TlMatrixObject::index_type TlDenseMatrix_arrays_Object::getNumOfLocalVectors(const index_type numOfVectors,
                                                                             const int numOfSubunits,
                                                                             const int sizeOfChunk) {
    const index_type numOfLocalChunks =
        TlDenseMatrix_arrays_Object::getNumOfLocalChunks(numOfVectors, numOfSubunits, sizeOfChunk);
    const index_type numOfLocalVectors = sizeOfChunk * numOfLocalChunks;

    return numOfLocalVectors;
}

TlMatrixObject::index_type TlDenseMatrix_arrays_Object::getLocalVectorIndex(const index_type vectorIndex,
                                                                            int* pSubunitId, int* pLocalChunkId,
                                                                            int* pLocalChunkVectorIndex) const {
    // assert((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_));
    TL_ASSERT((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_),
              TlUtils::format("vectorIndex=%d (max:%d)", vectorIndex, this->numOfVectors_));

    const index_type sizeOfBlock = this->sizeOfChunk_ * this->numOfSubunits_;
    const div_t blocks = std::div(vectorIndex, sizeOfBlock);
    const index_type numOfBlocks = blocks.quot;

    const div_t chunks = std::div(blocks.rem, this->sizeOfChunk_);
    const index_type chunkId = chunks.quot;
    const index_type chunkIndex = chunks.rem;

    const index_type localIndex = numOfBlocks * this->sizeOfChunk_ + chunkIndex;

    if (pSubunitId != NULL) {
        *pSubunitId = chunkId;
    }

    if (pLocalChunkId != NULL) {
        *pLocalChunkId = numOfBlocks;
    }

    if (pLocalChunkVectorIndex != NULL) {
        *pLocalChunkVectorIndex = chunkIndex;
    }

    return localIndex;
}

void TlDenseMatrix_arrays_Object::set_to_vm(const index_type vectorIndex, const index_type index, const double value) {
    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    (void)this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkId < static_cast<int>(this->chunks_.size()));
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        const int localIndex = this->reservedVectorSize_ * localChunkVectorIndex + index;
        this->chunks_[localChunkId][localIndex] = value;
    }
}

void TlDenseMatrix_arrays_Object::add_to_vm(const index_type vectorIndex, const index_type index, const double value) {
    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    (void)this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkId < static_cast<int>(this->chunks_.size()));
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        const int localIndex = this->reservedVectorSize_ * localChunkVectorIndex + index;
        this->chunks_[localChunkId][localIndex] += value;
    }
}

double TlDenseMatrix_arrays_Object::get_from_vm(const index_type vectorIndex, const index_type index) const {
    assert((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_));
    assert((0 <= index) && (index < this->sizeOfVector_));

    double answer = 0.0;
    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    (void)this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkId < static_cast<int>(this->chunks_.size()));
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        const int localIndex = this->reservedVectorSize_ * localChunkVectorIndex + index;
        answer = this->chunks_[localChunkId][localIndex];
    }

    return answer;
}

std::vector<double> TlDenseMatrix_arrays_Object::getVector(const index_type vectorIndex) const {
    const index_type vectorSize = this->sizeOfVector_;
    std::vector<double> answer(vectorSize);

    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    (void)this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkId < static_cast<int>(this->chunks_.size()));
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        const int localIndex = this->reservedVectorSize_ * localChunkVectorIndex;

        std::copy(&(this->chunks_[localChunkId][localIndex]), &(this->chunks_[localChunkId][localIndex]) + vectorSize,
                  answer.begin());
    }

    return answer;
}

void TlDenseMatrix_arrays_Object::getVector(const index_type vectorIndex, double* pBuf, const index_type length) const {
    const index_type vectorSize = this->sizeOfVector_;
    const index_type copySize = std::min(length, vectorSize);

    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    (void)this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkId < static_cast<int>(this->chunks_.size()));
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        const int localIndex = this->reservedVectorSize_ * localChunkVectorIndex;
        std::copy(this->chunks_[localChunkId] + localIndex, this->chunks_[localChunkId] + (localIndex + copySize),
                  pBuf);
    }
}

void TlDenseMatrix_arrays_Object::setVector(const index_type vectorIndex, const std::vector<double>& v) {
    assert(v.size() == static_cast<std::size_t>(this->sizeOfVector_));

    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    (void)this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkId < static_cast<int>(this->chunks_.size()));
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        const int localIndex = this->reservedVectorSize_ * localChunkVectorIndex;
        std::copy(v.begin(), v.end(), &(this->chunks_[localChunkId][localIndex]));
    }
}

std::size_t TlDenseMatrix_arrays_Object::getChunk(const index_type vectorIndex, double* pBuf,
                                                  const std::size_t length) const {
    // TODO: length support
    const int sizeOfChunk = this->sizeOfChunk_;
    const int sizeOfVector = this->sizeOfVector_;

    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    assert(subunitId == this->subunitID_);

    std::size_t copiedSize = 0;
    {
        const int sizeOfVectorReserved = this->reservedVectorSize_;
        for (int v = 0; v < sizeOfChunk; ++v) {
            std::copy(&(this->chunks_[localChunkId][v * sizeOfVectorReserved]),
                      &(this->chunks_[localChunkId][v * sizeOfVectorReserved]) + sizeOfVector,
                      &(pBuf[v * sizeOfVector]));
            copiedSize += sizeOfVector;
        }
    }

    return copiedSize;
}

std::string TlDenseMatrix_arrays_Object::getFileName(const std::string& basename, const int subunitID) {
    return TlUtils::format("%s.part%d.mat", basename.c_str(), subunitID);
}

bool TlDenseMatrix_arrays_Object::save(const std::string& basename) const {
    const index_type numOfVectors = this->numOfVectors_;
    const index_type sizeOfVector = this->sizeOfVector_;

    std::ofstream ofs;
    const std::string path = TlDenseMatrix_arrays_Object::getFileName(basename, this->subunitID_);
    ofs.open(path.c_str(), std::ofstream::out | std::ofstream::binary);

    // header
    const int matrixType = TlMatrixObject::ABGD;
    const index_type reservedSizeOfVector = sizeOfVector;
    ofs.write(reinterpret_cast<const char*>(&matrixType), sizeof(char));
    ofs.write(reinterpret_cast<const char*>(&numOfVectors), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&sizeOfVector), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&reservedSizeOfVector), sizeof(index_type));  // common use for array_mmap

    ofs.write(reinterpret_cast<const char*>(&(this->numOfSubunits_)), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&(this->subunitID_)), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&(this->sizeOfChunk_)), sizeof(int));

    {
        const int numOfLocalChunks = this->numOfLocalChunks_;
        assert(this->chunks_.size() == std::size_t(numOfLocalChunks));
        const std::size_t sizeOfVectorReserved = this->reservedVectorSize_;
        const int sizeOfChunk = this->sizeOfChunk_;
        for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
            for (int v = 0; v < sizeOfChunk; ++v) {
                const std::size_t headOfVector = sizeOfVectorReserved * v;
                ofs.write(reinterpret_cast<const char*>(this->chunks_[chunk] + headOfVector),
                          sizeof(double) * sizeOfVector);
            }
        }
    }

    ofs.close();

    return true;
}

bool TlDenseMatrix_arrays_Object::saveSubunitFileWithReservedSizeOfVector(const std::string& path) const {
    const index_type numOfVectors = this->numOfVectors_;
    const index_type sizeOfVector = this->sizeOfVector_;

    std::ofstream ofs;
    ofs.open(path.c_str(), std::ofstream::out | std::ofstream::binary);

    // header
    const int matrixType = TlMatrixObject::ABGD;
    ofs.write(reinterpret_cast<const char*>(&matrixType), sizeof(char));
    ofs.write(reinterpret_cast<const char*>(&numOfVectors), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&sizeOfVector), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&sizeOfVector), sizeof(index_type));  // for reservedSizeOfVector
    ofs.write(reinterpret_cast<const char*>(&(this->numOfSubunits_)), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&(this->subunitID_)), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&(this->sizeOfChunk_)), sizeof(int));

    {
        const int numOfLocalChunks = this->numOfLocalChunks_;
        assert(this->chunks_.size() == std::size_t(numOfLocalChunks));
        const std::size_t sizeOfVectorReserved = this->reservedVectorSize_;
        const int sizeOfChunk = this->sizeOfChunk_;
        for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
            for (int v = 0; v < sizeOfChunk; ++v) {
                const std::size_t headOfVector = sizeOfVectorReserved * v;
                ofs.write(reinterpret_cast<const char*>(this->chunks_[chunk] + headOfVector),
                          sizeof(double) * sizeOfVector);
            }
        }
    }

    ofs.close();

    return true;
}

bool TlDenseMatrix_arrays_Object::saveByTheOtherType(const std::string& basename) const {
    assert(this->numOfSubunits_ == 1);

    const index_type numOfVectors = this->numOfVectors_;
    const index_type sizeOfVector = this->sizeOfVector_;
    const int sizeOfChunk = this->sizeOfChunk_;

    const index_type numOfVectorsB = sizeOfVector;  // swap!
    const index_type sizeOfVectorB = numOfVectors;  // swap!

    std::ofstream ofs;
    const std::string path = TlDenseMatrix_arrays_Object::getFileName(basename, this->subunitID_);
    ofs.open(path.c_str(), std::ofstream::out | std::ofstream::binary);

    // header
    const int matrixType = TlMatrixObject::ABGD;
    ofs.write(reinterpret_cast<const char*>(&matrixType), sizeof(char));
    ofs.write(reinterpret_cast<const char*>(&numOfVectorsB), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&sizeOfVectorB), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&(this->numOfSubunits_)), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&(this->subunitID_)), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&sizeOfChunk), sizeof(int));

    // data
    {
        const int numOfLocalChunks = this->numOfLocalChunks_;
        const int reservedVectorSize = this->reservedVectorSize_;
        double* pBuf = new double[sizeOfVectorB];
        for (index_type j = 0; j < sizeOfVector; ++j) {
            for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
                for (int c = 0; c < sizeOfChunk; ++c) {
                    const std::size_t headOfVector = reservedVectorSize * c;
                    const int i = chunk * sizeOfChunk + c;
                    if (i >= numOfVectors) {
                        break;
                    }
                    pBuf[i] = this->chunks_[chunk][headOfVector + j];
                }
            }
            ofs.write(reinterpret_cast<const char*>(pBuf), sizeof(double) * sizeOfVectorB);
        }

        delete[] pBuf;
        pBuf = NULL;
    }

    ofs.close();

    return true;
}

bool TlDenseMatrix_arrays_Object::isLoadable(const std::string& filepath, index_type* pNumOfVectors,
                                             index_type* pSizeOfVector, int* pNumOfSubunits, int* pSubunitID,
                                             int* pSizeOfChunk) {
    bool answer = false;

    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ifstream::in);
    if (ifs.good()) {
        // read header
        int matrixType = 0;
        index_type numOfVectors = 0;
        index_type sizeOfVector = 0;
        index_type sizeOfReservedVector = 0;
        int numOfSubunits = 0;
        int subunitID = 0;
        int sizeOfChunk = 0;
        ifs.read((char*)&matrixType, sizeof(char));
        ifs.read((char*)&numOfVectors, sizeof(index_type));
        ifs.read((char*)&sizeOfVector, sizeof(index_type));
        ifs.read((char*)&sizeOfReservedVector, sizeof(index_type));  // common use for array_mmap
        ifs.read((char*)&numOfSubunits, sizeof(int));
        ifs.read((char*)&subunitID, sizeof(int));
        ifs.read((char*)&sizeOfChunk, sizeof(int));

        if ((numOfVectors > 0) && (sizeOfVector > 0) && (numOfSubunits > 0) && (subunitID >= 0) &&
            (subunitID < numOfSubunits)) {
            if (pNumOfVectors != NULL) {
                *pNumOfVectors = numOfVectors;
            }
            if (pSizeOfVector != NULL) {
                *pSizeOfVector = sizeOfVector;
            }
            if (pNumOfSubunits != NULL) {
                *pNumOfSubunits = numOfSubunits;
            }
            if (pSubunitID != NULL) {
                *pSubunitID = subunitID;
            }
            if (pSizeOfChunk != NULL) {
                *pSizeOfChunk = sizeOfChunk;
            }
            answer = true;
        }
    }

    return answer;
}

bool TlDenseMatrix_arrays_Object::load(const std::string& basename) {
    return this->load(basename, this->subunitID_);
}

bool TlDenseMatrix_arrays_Object::load(const std::string& basename, const int subunitID) {
    bool answer = false;

    std::string path = basename;
    if (subunitID >= 0) {
        path = TlDenseMatrix_arrays_Object::getFileName(basename, subunitID);
    }

    // std::ifstream ifs;
    // ifs.open(path.c_str(), std::ifstream::in);
    // if (ifs.good()) {
    //     // header
    //     int matrixType = 0;
    //     index_type numOfVectors = 0;
    //     index_type sizeOfVector = 0;
    //     int read_numOfSubunits = 0;
    //     int read_subunitID = 0;
    //     int read_sizeOfChunk = 0;
    //     ifs.read((char*)&matrixType, sizeof(char));
    //     ifs.read((char*)&numOfVectors, sizeof(index_type));
    //     ifs.read((char*)&sizeOfVector, sizeof(index_type));
    //     ifs.read((char*)&read_numOfSubunits, sizeof(int));
    //     ifs.read((char*)&read_subunitID, sizeof(int));
    //     ifs.read((char*)&read_sizeOfChunk, sizeof(int));
    //     this->numOfSubunits_ = read_numOfSubunits;
    //     this->subunitID_ = read_subunitID;
    //     this->resize(numOfVectors, sizeOfVector);

    //     this->numOfSubunits_ = read_numOfSubunits;
    //     this->subunitID_ = read_subunitID;
    //     this->sizeOfChunk_ = read_sizeOfChunk;
    //     this->resize(numOfVectors, sizeOfVector);

    //     const int headerSize = sizeof(char) + sizeof(index_type) * 2 + sizeof(int) * 3;
    //     assert(headerSize == ifs.tellg());

    //     // data
    //     ifs.clear();
    //     ifs.seekg(headerSize, std::ios_base::beg);
    //     assert(headerSize == ifs.tellg());
    //     // std::size_t copied_chunk = 0;
    //     {
    //         const int numOfLocalChunks = this->numOfLocalChunks_;
    //         assert(static_cast<std::size_t>(numOfLocalChunks) == this->chunks_.size());
    //         const int sizeOfChunk = this->sizeOfChunk_;
    //         const int sizeOfVector = this->sizeOfVector_;
    //         const int sizeOfVectorReserved = this->reservedVectorSize_;

    //         std::vector<double> buf(sizeOfVectorReserved);
    //         for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
    //             for (int v = 0; v < sizeOfChunk; ++v) {
    //                 ifs.read((char*)&(buf[0]), sizeof(double) * sizeOfVector);
    //                 // copied_chunk += sizeOfVector;

    //                 std::copy(buf.begin(), buf.begin() + sizeOfVector,
    //                           this->chunks_[chunk] + (v * sizeOfVectorReserved));
    //             }
    //         }
    //     }

    //     answer = true;
    // } else {
    //     log.error(TlUtils::format("cannot open file: %s", path.c_str()));
    // }

    // ifs.close();

    answer = this->loadSubunitFile(path);
    return answer;
}

bool TlDenseMatrix_arrays_Object::loadSubunitFile(const std::string& path) {
    TlLogging& log = TlLogging::getInstance();
    bool answer = false;

    std::ifstream ifs;
    ifs.open(path.c_str(), std::ifstream::in);
    if (ifs.good()) {
        // header
        int matrixType = 0;
        index_type numOfVectors = 0;
        index_type sizeOfVector = 0;
        int read_numOfSubunits = 0;
        int read_subunitID = 0;
        int read_sizeOfChunk = 0;
        ifs.read((char*)&matrixType, sizeof(char));
        ifs.read((char*)&numOfVectors, sizeof(index_type));
        ifs.read((char*)&sizeOfVector, sizeof(index_type));
        ifs.read((char*)&read_numOfSubunits, sizeof(int));
        ifs.read((char*)&read_subunitID, sizeof(int));
        ifs.read((char*)&read_sizeOfChunk, sizeof(int));
        this->numOfSubunits_ = read_numOfSubunits;
        this->subunitID_ = read_subunitID;
        this->resize(numOfVectors, sizeOfVector);

        this->numOfSubunits_ = read_numOfSubunits;
        this->subunitID_ = read_subunitID;
        this->sizeOfChunk_ = read_sizeOfChunk;
        this->resize(numOfVectors, sizeOfVector);

        const int headerSize = sizeof(char) + sizeof(index_type) * 2 + sizeof(int) * 3;
        assert(headerSize == ifs.tellg());

        // data
        ifs.clear();
        ifs.seekg(headerSize, std::ios_base::beg);
        assert(headerSize == ifs.tellg());
        // std::size_t copied_chunk = 0;
        {
            const int numOfLocalChunks = this->numOfLocalChunks_;
            assert(static_cast<std::size_t>(numOfLocalChunks) == this->chunks_.size());
            const int sizeOfChunk = this->sizeOfChunk_;
            const int sizeOfVector = this->sizeOfVector_;
            const int sizeOfVectorReserved = this->reservedVectorSize_;

            std::vector<double> buf(sizeOfVectorReserved);
            for (int chunk = 0; chunk < numOfLocalChunks; ++chunk) {
                for (int v = 0; v < sizeOfChunk; ++v) {
                    ifs.read((char*)&(buf[0]), sizeof(double) * sizeOfVector);
                    // copied_chunk += sizeOfVector;

                    std::copy(buf.begin(), buf.begin() + sizeOfVector,
                              this->chunks_[chunk] + (v * sizeOfVectorReserved));
                }
            }
        }

        answer = true;
    } else {
        log.error(TlUtils::format("cannot open file: %s", path.c_str()));
    }

    ifs.close();

    return answer;
}

std::size_t TlDenseMatrix_arrays_Object::getMemSize() const {
    return this->numOfVectors_ * this->sizeOfVector_;
}

void TlDenseMatrix_arrays_Object::destroy() {
    const index_type numOfLocalChunks = this->numOfLocalChunks_;
    if (this->isUsingMemManager_ == true) {
        TlMemManager& rMemManager = TlMemManager::getInstance();
        for (int i = 0; i < numOfLocalChunks; ++i) {
            rMemManager.deallocate((char*)this->chunks_[i]);
            this->chunks_[i] = NULL;
        }
    } else {
        for (int i = 0; i < numOfLocalChunks; ++i) {
            delete[] this->chunks_[i];
            this->chunks_[i] = NULL;
        }
    }
    this->chunks_.clear();
}
