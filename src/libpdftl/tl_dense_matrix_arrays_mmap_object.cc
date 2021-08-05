#include "tl_dense_matrix_arrays_mmap_object.h"

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>

#include "TlFile.h"
#include "TlLogging.h"
#include "TlMemManager.h"
#include "TlSystem.h"
#include "TlUtils.h"
#include "tl_assert.h"
#include "tl_matrix_utils.h"

#define DEFAULT_CHUNK_SIZE 32
#define CLEAR_BUFSIZE (4096)

TlDenseMatrix_arrays_mmap_Object::TlDenseMatrix_arrays_mmap_Object(
    const std::string& baseFilePath, const index_type numOfVectors, const index_type sizeOfVector,
    const int numOfSubunits, const int subunitID, const index_type reservedSizeOfVector, bool forceCreateNewFile)
    : TlMatrixObject(ABGD),
      numOfVectors_(numOfVectors),
      sizeOfVector_(sizeOfVector),
      reservedSizeOfVector_(reservedSizeOfVector),
      sizeOfChunk_(DEFAULT_CHUNK_SIZE),
      numOfSubunits_(numOfSubunits),
      subunitID_(subunitID) {
    this->filePath_ = TlDenseMatrix_arrays_mmap_Object::getFileName(baseFilePath, this->subunitID_);
    // this->numOfLocalChunks_ = TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(
    //     this->getNumOfVectors(), this->getNumOfSubunits(), this->getSizeOfChunk());

    if (reservedSizeOfVector == 0) {
        this->reservedSizeOfVector_ = this->sizeOfVector_;
    }
    // std::cerr << "reserved size of vector = " << this->reservedSizeOfVector_ << std::endl;

    this->createNewFile(forceCreateNewFile);
    this->openFile();

    // std::cerr << "# subunits: " << this->getNumOfSubunits() << std::endl;
    // std::cerr << "# subunitID: " << this->getSubunitID() << std::endl;
}

TlDenseMatrix_arrays_mmap_Object::TlDenseMatrix_arrays_mmap_Object(const std::string& filePath)
    : TlMatrixObject(ABGD), filePath_(filePath), sizeOfChunk_(DEFAULT_CHUNK_SIZE) {
    // std::cerr << "TlDenseMatrix_arrays_mmap_Object::TlDenseMatrix_arrays_mmap_Object() constructor2 called"
    //           << std::endl;
    this->openFile();

    // std::cerr << "# subunits: " << this->getNumOfSubunits() << std::endl;
    // std::cerr << "# subunitID: " << this->getSubunitID() << std::endl;
}

TlDenseMatrix_arrays_mmap_Object::~TlDenseMatrix_arrays_mmap_Object() {
    this->deleteMmap();
    // this->destroy();
}

int TlDenseMatrix_arrays_mmap_Object::getSizeOfChunk() const {
    return this->sizeOfChunk_;
}

void TlDenseMatrix_arrays_mmap_Object::resize(const index_type newNumOfVectors, const index_type newSizeOfVector) {
    if (newNumOfVectors != this->getNumOfVectors()) {
        // I end up recreating the file, so both the sizeOfVector and the numOfVectors are changed together.

        // move file
        this->deleteMmap();
        const std::string filePath = this->filePath_;
        const std::string backupFilePath = this->filePath_ + ".bak";
        const int ret = TlFile::rename(filePath, backupFilePath);
        assert(ret == 0);

        // create new matrix
        this->numOfVectors_ = newNumOfVectors;
        this->sizeOfVector_ = newSizeOfVector;
        this->reservedSizeOfVector_ = std::max(newSizeOfVector, this->reservedSizeOfVector_);

        this->createNewFile();
        this->openFile();
        this->copyFromBackup();

        TlFile::remove(backupFilePath);
    } else {
        // change sizeOfVector only

        if (newSizeOfVector > this->reservedSizeOfVector_) {
            this->reserveVectorSize(newSizeOfVector);
        }

        if (newSizeOfVector != this->getSizeOfVector()) {
            assert(newSizeOfVector <= this->reservedSizeOfVector_);

            this->sizeOfVector_ = newSizeOfVector;
            this->updateFileHeader();
        }
    }
}

void TlDenseMatrix_arrays_mmap_Object::copyFromBackup() {
    TlDenseMatrix_arrays_mmap_Object backup(this->filePath_ + ".bak");
    const index_type oldNumOfVectors = backup.getNumOfVectors();
    const index_type oldSizeOfVector = backup.getSizeOfVector();
    // const index_type oldReservedSizeOfVector = backup.reservedSizeOfVector_;
    const index_type newNumOfVectors = this->getNumOfVectors();
    const index_type newSizeOfVector = this->getSizeOfVector();
    const index_type newReservedSizeOfVector = this->reservedSizeOfVector_;

    // std::cerr << ">>>> copy From backup: " << std::endl;
    // std::cerr << TlUtils::format("# vec: %d / %d", oldNumOfVectors, newNumOfVectors) << std::endl;
    // std::cerr << TlUtils::format("size of vec: %d / %d", oldSizeOfVector, newSizeOfVector) << std::endl;
    // std::cerr << TlUtils::format("reserved: %d / %d", oldReservedSizeOfVector, newReservedSizeOfVector) << std::endl;

    const int sizeOfChunk = this->getSizeOfChunk();
    const int numOfSubunits = this->getNumOfSubunits();
    const int subunitId = this->getSubunitID();

    const TlMatrixObject::size_type oldNumOfLocalChunks =
        TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(oldNumOfVectors, numOfSubunits, sizeOfChunk);
    const TlMatrixObject::size_type newNumOfLocalChunks =
        TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(newNumOfVectors, numOfSubunits, sizeOfChunk);
    const int copyNumOfLocalChunks = std::min(oldNumOfLocalChunks, newNumOfLocalChunks);

    const TlMatrixObject::size_type copySizeOfVector = std::min(oldSizeOfVector, newSizeOfVector);

    std::size_t count = 0;
    std::vector<double> oldChunk(sizeOfChunk * oldSizeOfVector);
    for (int chunkIndex = 0; chunkIndex < copyNumOfLocalChunks; ++chunkIndex) {
        const index_type globalVectorIndex = sizeOfChunk * (chunkIndex * numOfSubunits + subunitId);
        if (globalVectorIndex > oldNumOfVectors) {
            break;
        }

        backup.getChunk(globalVectorIndex, &(oldChunk[0]), sizeOfChunk * oldSizeOfVector);
        // for (int i = 0; i < sizeOfChunk * oldSizeOfVector; ++i) {
        //     std::cerr << TlUtils::format("[%d] %d: %f", i, count + i, oldChunk[i]) << std::endl;
        // }
        count += sizeOfChunk * oldSizeOfVector;

        const std::size_t copyOutOffset = sizeOfChunk * chunkIndex * newReservedSizeOfVector;
        for (int i = 0; i < sizeOfChunk; ++i) {
            // const std::size_t copyFrom = oldReservedSizeOfVector * i;
            const std::size_t copyFrom = oldSizeOfVector * i;
            const std::size_t copyOut = copyOutOffset + newReservedSizeOfVector * i;
            // std::cerr << "copy out = " << copyOut << std::endl;
            // for (int t = 0; t < copySizeOfVector; ++t) {
            //     std::cerr << TlUtils::format("chunk[%d] = %f", copyFrom + t, oldChunk[copyFrom + t]) << std::endl;
            // }
            std::copy(&(oldChunk[0]) + copyFrom, &(oldChunk[0]) + (copyFrom + copySizeOfVector),
                      this->dataBegin_ + copyOut);
        }
    }
}

void TlDenseMatrix_arrays_mmap_Object::reserveVectorSize(index_type newReservedSizeOfVector) {
    newReservedSizeOfVector = std::max(this->sizeOfVector_, newReservedSizeOfVector);
    const index_type prevReservedSizeOfVector = this->reservedSizeOfVector_;

    // std::cerr << ">>>> reserved: " << std::endl;
    // std::cerr << "reserved1: " << std::endl;
    // std::cerr << TlUtils::format("# vec: %d", this->numOfVectors_) << std::endl;
    // std::cerr << TlUtils::format("size of vec: %d", this->sizeOfVector_) << std::endl;
    // std::cerr << TlUtils::format("reserved: %d", this->reservedSizeOfVector_) << std::endl;

    if (prevReservedSizeOfVector < newReservedSizeOfVector) {
        this->deleteMmap();

        // move file
        const std::string filePath = this->filePath_;
        const std::string backupFilePath = this->filePath_ + ".bak";
        const int ret = TlFile::rename(filePath, backupFilePath);
        assert(ret == 0);

        // create new matrix
        this->reservedSizeOfVector_ = newReservedSizeOfVector;
        this->createNewFile();
        this->openFile();

        this->copyFromBackup();

        // std::cerr << "reserved2: " << std::endl;
        // std::cerr << TlUtils::format("# vec: %d", this->numOfVectors_) << std::endl;
        // std::cerr << TlUtils::format("size of vec: %d", this->sizeOfVector_) << std::endl;
        // std::cerr << TlUtils::format("reserved: %d", this->reservedSizeOfVector_) << std::endl;

        TlFile::remove(backupFilePath);
    }
}

int TlDenseMatrix_arrays_mmap_Object::getSubunitID(const index_type vectorIndex) const {
    assert((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_));
    int subunitId = -1;
    this->getLocalVectorIndex(vectorIndex, &subunitId);

    return subunitId;
}

// static function
TlMatrixObject::index_type TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(const index_type numOfVectors,
                                                                                 const int numOfSubunits,
                                                                                 const int sizeOfChunk) {
    const TlMatrixObject::index_type numOfVectorsPerUnit = sizeOfChunk * numOfSubunits;
    const TlMatrixObject::index_type numOfLocalChunks =
        (numOfVectors + numOfVectorsPerUnit - 1) / numOfVectorsPerUnit;  // round up

    return numOfLocalChunks;
}

// defunct
// TlMatrixObject::index_type TlDenseMatrix_arrays_mmap_Object::getNumOfLocalVectors(const index_type numOfVectors,
//                                                                                   const int numOfSubunits,
//                                                                                   const int sizeOfChunk) {
//     const index_type numOfLocalChunks =
//         TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(numOfVectors, numOfSubunits, sizeOfChunk);
//     const index_type numOfLocalVectors = sizeOfChunk * numOfLocalChunks;

//     return numOfLocalVectors;
// }

std::size_t TlDenseMatrix_arrays_mmap_Object::getLocalVectorIndex(const index_type vectorIndex, int* pSubunitId,
                                                                  int* pLocalChunkId,
                                                                  int* pLocalChunkVectorIndex) const {
    TL_ASSERT((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_),
              TlUtils::format("vectorIndex=%d (max:%d)", vectorIndex, this->numOfVectors_));

    const index_type sizeOfBlock = this->sizeOfChunk_ * this->numOfSubunits_;
    const div_t blocks = std::div(vectorIndex, sizeOfBlock);
    const index_type numOfBlocks = blocks.quot;

    const div_t chunks = std::div(blocks.rem, this->sizeOfChunk_);
    const index_type chunkId = chunks.quot;
    const index_type chunkIndex = chunks.rem;

    // const index_type localIndex = numOfBlocks * this->sizeOfChunk_ + chunkIndex;
    const std::size_t localIndex = (this->sizeOfChunk_ * numOfBlocks + chunkIndex) * this->reservedSizeOfVector_;

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

void TlDenseMatrix_arrays_mmap_Object::set_to_vm(const index_type vectorIndex, const index_type index,
                                                 const double value) {
    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    const std::size_t head = this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        this->dataBegin_[head + index] = value;
    }
}

void TlDenseMatrix_arrays_mmap_Object::add_to_vm(const index_type vectorIndex, const index_type index,
                                                 const double value) {
    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    const std::size_t head = this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        this->dataBegin_[head + index] += value;
    }
}

double TlDenseMatrix_arrays_mmap_Object::get_from_vm(const index_type vectorIndex, const index_type index) const {
    assert((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_));
    assert((0 <= index) && (index < this->sizeOfVector_));

    double answer = 0.0;
    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    const std::size_t head = this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        answer = this->dataBegin_[head + index];

        // std::cerr << TlUtils::format("get (%d, %d) [%d, %d] = %f", vectorIndex, index, head, index, answer)
        //           << std::endl;
    }

    return answer;
}

std::vector<double> TlDenseMatrix_arrays_mmap_Object::getVector(const index_type vectorIndex) const {
    const index_type vectorSize = this->sizeOfVector_;
    std::vector<double> answer(vectorSize);

    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    const std::size_t head = this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        // const int localIndex = this->reservedSizeOfVector_ * localChunkVectorIndex;
        // std::copy(this->dataBegin_ + (this->sizeOfChunk_ * localChunkId + localIndex),
        //           this->dataBegin_ + (this->sizeOfChunk_ * localChunkId + localIndex) + vectorSize,
        //           answer.begin());
        std::copy(this->dataBegin_ + head, this->dataBegin_ + head + vectorSize, answer.begin());
    }

    return answer;
}

std::size_t TlDenseMatrix_arrays_mmap_Object::getVector(const index_type vectorIndex, double* pBuf,
                                                        const std::size_t maxCount) const {
    const std::size_t vectorSize = this->sizeOfVector_;
    std::size_t copySize = 0;

    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    const std::size_t head = this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        // const int localIndex = this->reservedSizeOfVector_ * localChunkVectorIndex;
        // std::copy(this->dataBegin_ + (this->sizeOfChunk_ * localChunkId + localIndex),
        //           this->dataBegin_ + (this->sizeOfChunk_ * localChunkId + localIndex) + copySize, pBuf);

        copySize = std::min(maxCount, vectorSize);
        std::copy(this->dataBegin_ + head, this->dataBegin_ + head + copySize, pBuf);
    }

    return copySize;
}

void TlDenseMatrix_arrays_mmap_Object::setVector(const index_type vectorIndex, const std::vector<double>& v) {
    assert(v.size() == static_cast<std::size_t>(this->sizeOfVector_));

    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    const std::size_t head = this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    if (subunitId == this->subunitID_) {
        assert(localChunkVectorIndex < this->sizeOfChunk_);
        // const int localIndex = this->reservedSizeOfVector_ * localChunkVectorIndex;
        // std::copy(v.begin(), v.end(), this->dataBegin_ + (this->sizeOfChunk_ * localChunkId + localIndex));
        std::copy(v.begin(), v.end(), this->dataBegin_ + head);
    }
}

void TlDenseMatrix_arrays_mmap_Object::setAcrossMultipleVectors(index_type index, const std::valarray<double>& values) {
    const index_type numOfVectors = this->getNumOfVectors();
    TL_ASSERT((values.size() == static_cast<std::size_t>(numOfVectors)),
              TlUtils::format("%ld != %ld @%s,%d", values.size(), numOfVectors, __FILE__, __LINE__));

    for (index_type vectorIndex = 0; vectorIndex < numOfVectors; ++vectorIndex) {
        this->set_to_vm(vectorIndex, index, values[vectorIndex]);
    }
}

std::size_t TlDenseMatrix_arrays_mmap_Object::getChunk(const index_type vectorIndex, double* pBuf,
                                                       const std::size_t length) const {
    TL_ASSERT((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_),
              TlUtils::format("vectorIndex=%d (max:%d)", vectorIndex, this->numOfVectors_));

    const int sizeOfChunk = this->sizeOfChunk_;
    const int sizeOfVector = this->sizeOfVector_;

    int subunitId = 0;
    int localChunkId = 0;
    int localChunkVectorIndex = 0;
    const std::size_t head = this->getLocalVectorIndex(vectorIndex, &subunitId, &localChunkId, &localChunkVectorIndex);
    assert(subunitId == this->subunitID_);

    std::size_t copiedSize = 0;
    {
        const int reservedSizeOfVector = this->reservedSizeOfVector_;
        for (int v = 0; v < sizeOfChunk; ++v) {
            // std::copy(this->dataBegin_ + this->sizeOfChunk_ * localChunkId + (v * reservedSizeOfVector),
            //           this->dataBegin_ + this->sizeOfChunk_ * localChunkId + (v * reservedSizeOfVector) +
            //           sizeOfVector,
            //           &(pBuf[v * sizeOfVector]));
            const std::size_t size = std::min<std::size_t>(sizeOfVector, length - copiedSize);
            std::copy(this->dataBegin_ + head + reservedSizeOfVector * v,
                      this->dataBegin_ + head + reservedSizeOfVector * v + size, &(pBuf[v * sizeOfVector]));
            copiedSize += sizeOfVector;
        }
    }

    return copiedSize;
}

std::string TlDenseMatrix_arrays_mmap_Object::getFileName(const std::string& baseFilePath, const int subunitID) {
    return TlUtils::format("%s.part%d.mat", baseFilePath.c_str(), subunitID);
}

bool TlDenseMatrix_arrays_mmap_Object::isLoadable(const std::string& filepath, index_type* pNumOfVectors,
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
        index_type reservedSizeOfVector = 0;
        int numOfSubunits = 0;
        int subunitID = 0;
        int sizeOfChunk = 0;
        ifs.read((char*)&matrixType, sizeof(char));
        ifs.read((char*)&numOfVectors, sizeof(index_type));
        ifs.read((char*)&sizeOfVector, sizeof(index_type));
        ifs.read((char*)&reservedSizeOfVector, sizeof(index_type));

        ifs.read((char*)&numOfSubunits, sizeof(int));
        ifs.read((char*)&subunitID, sizeof(int));
        ifs.read((char*)&sizeOfChunk, sizeof(int));

        // std::cerr << TlUtils::format("type: %d (%d, %d/%d) [%d/%d] c=%d",
        //                              matrixType, numOfVectors, sizeOfVector, reservedSizeOfVector, subunitID, numOfSubunits, sizeOfChunk)
        //           << std::endl;

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

std::size_t TlDenseMatrix_arrays_mmap_Object::getMemSize() const {
    return this->numOfVectors_ * this->sizeOfVector_;
}

TlMatrixObject::size_type TlDenseMatrix_arrays_mmap_Object::getNumOfReservedElements() const {
    const TlMatrixObject::size_type numOfLocalChunks = TlDenseMatrix_arrays_mmap_Object::getNumOfLocalChunks(
        this->getNumOfVectors(), this->getNumOfSubunits(), this->getSizeOfChunk());
    const TlMatrixObject::size_type elements = (numOfLocalChunks * this->sizeOfChunk_) * this->reservedSizeOfVector_;

    return elements;
}

void TlDenseMatrix_arrays_mmap_Object::createNewFile(const bool isForceCreate) {
    if (TlFile::isExistFile(this->filePath_) && (isForceCreate == true)) {
        TlFile::remove(this->filePath_);
    }

    if (!TlFile::isExistFile(this->filePath_)) {
        this->createNewFile(this->filePath_, this->getNumOfReservedElements());
    }
}

void TlDenseMatrix_arrays_mmap_Object::createNewFile(const std::string& filePath, const std::size_t size) {
    // std::cerr << "TlDenseMatrix_arrays_mmap_Object::createNewFile(): " << filePath << ", " << size << std::endl;

    // create new file
    std::ofstream ofs;
    ofs.open(filePath.c_str(), std::ios::binary | std::ios::trunc | std::ios::out);

    // バッファ機能を停止する
    ofs << std::setiosflags(std::ios::unitbuf);

    (void)this->writeMatrixHeader(&ofs);

    // size分ファイルにzeroを埋める
    std::vector<double> buf(CLEAR_BUFSIZE);
    ldiv_t ldivt = ldiv(size, CLEAR_BUFSIZE);
    for (long i = 0; i < ldivt.quot; ++i) {
        ofs.write(reinterpret_cast<const char*>(&(buf[0])), sizeof(double) * CLEAR_BUFSIZE);
    }
    if (ldivt.rem > 0) {
        ofs.write(reinterpret_cast<const char*>(&(buf[0])), sizeof(double) * ldivt.rem);
    }

    // バッファ機能を再開する
    ofs << std::resetiosflags(std::ios::unitbuf);
    ofs.flush();

    ofs.close();
}

void TlDenseMatrix_arrays_mmap_Object::openFile() {
    // std::cerr << "TlDenseMatrix_arrays_mmap_Object::openFile(): getHeaderInfo" << std::endl;
    this->getHeaderInfo(this->filePath_);
    // std::cerr << "TlDenseMatrix_arrays_mmap_Object::openFile(): new mmap" << std::endl;
    this->newMmap();
    // std::cerr << "TlDenseMatrix_arrays_mmap_Object::openFile(): end" << std::endl;
}

void TlDenseMatrix_arrays_mmap_Object::getHeaderInfo(const std::string& filePath) {
    TlMatrixObject::HeaderInfo headerInfo;
    const bool isLoadable = TlMatrixUtils::getHeaderInfo(filePath, &headerInfo);
    assert(isLoadable == true);
    const std::size_t headerSize = headerInfo.headerSize;

    if ((isLoadable == true) && (headerInfo.matrixType == this->getType())) {
        this->numOfVectors_ = headerInfo.numOfVectors;
        this->sizeOfVector_ = headerInfo.sizeOfVector;
        this->reservedSizeOfVector_ = headerInfo.reservedSizeOfVector;
        TL_ASSERT((this->sizeOfVector_ <= this->reservedSizeOfVector_),
                  TlUtils::format("sizeOfVector, reservedSizeOfVector = %d, %d @%s,%d", this->sizeOfVector_,
                                  this->reservedSizeOfVector_, __FILE__, __LINE__));

        this->numOfSubunits_ = headerInfo.numOfSubunits;
        this->subunitID_ = headerInfo.subunitId;
        this->sizeOfChunk_ = headerInfo.sizeOfChunk;
    } else {
        std::cerr << TlUtils::format("cannot read matrix header: %s@%d", __FILE__, __LINE__) << std::endl;
        std::cerr << TlUtils::format("matrix: type=%d, numOfVectors=%d, sizeOfVector=%d; hs=%ld",
                                     static_cast<int>(headerInfo.matrixType), headerInfo.numOfVectors,
                                     headerInfo.sizeOfVector, headerSize)
                  << std::endl;
    }

    this->headerSize_ = headerSize;
    this->fileSize_ = headerSize + sizeof(double) * this->getNumOfReservedElements();
}

void TlDenseMatrix_arrays_mmap_Object::newMmap() {
    const int fd = open(this->filePath_.c_str(), O_RDWR, S_IRUSR | S_IWUSR);
    if (fd == -1) {
        int errorNo = errno;
        std::string errStr(strerror(errorNo));
        this->log_.critical(TlUtils::format("open error: %s (%s)", this->filePath_.c_str(), errStr.c_str()));
        throw errStr;
    }

    struct stat sb;
    fstat(fd, &sb);
    TL_ASSERT((std::size_t(sb.st_size) == this->fileSize_),
              TlUtils::format("%ld, %ld", std::size_t(sb.st_size), this->fileSize_));

    this->mmapBegin_ = (char*)TlSystem::newMmap(this->fileSize_, fd);
    msync(this->mmapBegin_, this->getNumOfReservedElements(), MS_ASYNC);
    madvise(this->mmapBegin_, this->getNumOfReservedElements(), MADV_WILLNEED);
    this->dataBegin_ = (double*)(this->mmapBegin_ + this->headerSize_);

    // we can close the file after mmap() was called.
    close(fd);
}

std::size_t TlDenseMatrix_arrays_mmap_Object::writeMatrixHeader(std::ofstream* pFs) const {
    std::size_t sizeOfWriteBuffer = 0;

    const char matrixType = static_cast<char>(this->getType());
    const index_type numOfVectors = this->numOfVectors_;
    const index_type sizeOfVector = this->sizeOfVector_;
    const index_type reservedSizeOfVector = this->reservedSizeOfVector_;

    pFs->write(reinterpret_cast<const char*>(&matrixType), sizeof(char));
    pFs->write(reinterpret_cast<const char*>(&numOfVectors), sizeof(index_type));
    pFs->write(reinterpret_cast<const char*>(&sizeOfVector), sizeof(index_type));
    pFs->write(reinterpret_cast<const char*>(&reservedSizeOfVector), sizeof(index_type));

    pFs->write(reinterpret_cast<const char*>(&(this->numOfSubunits_)), sizeof(int));
    pFs->write(reinterpret_cast<const char*>(&(this->subunitID_)), sizeof(int));
    pFs->write(reinterpret_cast<const char*>(&(this->sizeOfChunk_)), sizeof(int));

    sizeOfWriteBuffer = sizeof(char) + sizeof(index_type) * 3 + sizeof(int) * 3;

    return sizeOfWriteBuffer;
}

void TlDenseMatrix_arrays_mmap_Object::updateFileHeader() {
    const char matrixType = static_cast<char>(this->getType());
    const index_type numOfVectors = this->numOfVectors_;
    const index_type sizeOfVector = this->sizeOfVector_;
    const index_type reservedSizeOfVector = this->reservedSizeOfVector_;
    const int numOfSubunits = this->numOfSubunits_;
    const int subunitID = this->subunitID_;
    const int sizeOfChunk = this->sizeOfChunk_;

    std::size_t offset = 0;
    char const* p = NULL;

    p = reinterpret_cast<const char*>(&matrixType);
    std::copy(p, p + sizeof(char), this->mmapBegin_ + offset);
    offset += sizeof(char);

    p = reinterpret_cast<const char*>(&numOfVectors);
    std::copy(p, p + sizeof(index_type), this->mmapBegin_ + offset);
    offset += sizeof(index_type);

    p = reinterpret_cast<const char*>(&sizeOfVector);
    std::copy(p, p + sizeof(index_type), this->mmapBegin_ + offset);
    offset += sizeof(index_type);

    p = reinterpret_cast<const char*>(&reservedSizeOfVector);
    std::copy(p, p + sizeof(index_type), this->mmapBegin_ + offset);
    offset += sizeof(index_type);

    p = reinterpret_cast<const char*>(&numOfSubunits);
    std::copy(p, p + sizeof(int), this->mmapBegin_ + offset);
    offset += sizeof(int);

    p = reinterpret_cast<const char*>(&subunitID);
    std::copy(p, p + sizeof(int), this->mmapBegin_ + offset);
    offset += sizeof(int);

    p = reinterpret_cast<const char*>(&sizeOfChunk);
    std::copy(p, p + sizeof(int), this->mmapBegin_ + offset);
    offset += sizeof(int);
}

void TlDenseMatrix_arrays_mmap_Object::syncMmap() {
    int ret = msync(this->mmapBegin_, this->fileSize_, MS_SYNC);
    if (ret != 0) {
        std::cerr << "msync errored. " << ret << std::endl;
    }
}

void TlDenseMatrix_arrays_mmap_Object::deleteMmap() {
    this->syncMmap();
    munmap(this->mmapBegin_, this->fileSize_);

    this->mmapBegin_ = NULL;
}
