#include <iostream>
#include <cassert>
#include <sstream>

#include "GridDataManager.h"
#include "TlFile.h"
#include "TlVector.h"
#include "TlUtils.h"

const size_t GridDataManager::headerSize =
    sizeof(int) + sizeof(GridDataManager::ChunkType) + sizeof(size_t);

const std::string GridDataManager::chunkTypeStr[] = {
    "UNDEFINED",
    "COORD_X",
    "COORD_Y",
    "COORD_Z",
    "GRID_WEIGHT",
    "DENSITY",
    "DENSITY_ALPHA",
    "DENSITY_BETA",
    "GRAD_DENSITY_X",
    "GRAD_DENSITY_Y",
    "GRAD_DENSITY_Z",
    "GRAD_DENSITY_X_ALPHA",
    "GRAD_DENSITY_Y_ALPHA",
    "GRAD_DENSITY_Z_ALPHA",
    "GRAD_DENSITY_X_BETA",
    "GRAD_DENSITY_Y_BETA",
    "GRAD_DENSITY_Z_BETA"
};


const std::size_t GridDataManager::bufferSize_ = 128UL * 1024UL * 1024UL; // 128 MB

GridDataManager::GridDataManager(const std::string& filePath)
        : filePath_(filePath)
{
    this->open();
}


GridDataManager::~GridDataManager()
{
    AtomChunkDataType::iterator atomItrEnd = this->atomChunkData_.end();
    for (AtomChunkDataType::iterator atomItr = this->atomChunkData_.begin();
            atomItr != atomItrEnd; ++atomItr) {
        const int atomIndex = atomItr->first;
        ChunkDataType::iterator chunkItrEnd = atomItr->second.end();
        for (ChunkDataType::iterator chunkItr = atomItr->second.begin();
                chunkItr != chunkItrEnd; ++chunkItr) {
            const ChunkType chunkType = chunkItr->first;
            this->writeChunk(atomIndex, chunkType);
        }
    }

    this->fs_.flush();
    this->fs_.close();
}


void GridDataManager::open()
{
    if (TlFile::isExist(this->filePath_) == false) {
        // create new file
        this->fs_.open(this->filePath_.c_str(),
                       std::ios::binary | std::ios::trunc | std::ios::out);
        this->fs_.close();
    }

    this->fs_.open(this->filePath_.c_str(),
                   std::ios::binary | std::ios::in | std::ios::out);
    std::ios_base::iostate state = this->fs_.rdstate();
    if ((state & std::ios_base::eofbit) || (state & std::ios_base::failbit)) {
        this->fs_.clear();
    }
    this->fs_.setf(std::ios::unitbuf);
    assert(!this->fs_.fail());

    this->readHeaders();
}


void GridDataManager::readHeaders()
{
    // file size
    this->fs_.seekg(0, std::ios::end);
    std::fstream::pos_type maxPos = this->fs_.tellg();
    this->fs_.seekg(0, std::ios::beg);
    maxPos -= this->fs_.tellg();

    //
    int atomIndex = 0;
    ChunkType chunkType = UNDEFINED;
    size_t bufferSize = 0;
    while (this->fs_.eof() != true) {
        const std::fstream::pos_type filePosition = this->fs_.tellg();
        if (filePosition == (std::fstream::pos_type) -1) {
            break;
        }

        this->fs_.read(reinterpret_cast<char*>(&atomIndex), sizeof(int));
        this->fs_.read(reinterpret_cast<char*>(&chunkType), sizeof(ChunkType));
        this->fs_.read(reinterpret_cast<char*>(&bufferSize), sizeof(size_t));

        if (int(chunkType) > 0) {
            ChunkHeader ch;
            ch.filePosition = filePosition;
            ch.bufferSize = bufferSize;
            this->atomChunkHeader_[atomIndex][chunkType] = ch;
//       std::cout << "readHeaders: atom=" << atomIndex
//      << ", chunk=" << chunkType
//      << ", pos=" << filePosition
//      << ", buf=" << bufferSize
//      << std::endl;
        }

        std::fstream::pos_type newPosition = filePosition;
        newPosition += GridDataManager::headerSize;
        newPosition += bufferSize;
        if (newPosition >= maxPos) {
            break;
        }
        this->fs_.seekg(newPosition, std::ios::beg);
    }

    this->fs_.clear();
}

void GridDataManager::decomposeGridInfo(const std::vector<GridDataManager::GridInfo>& grids,
                                        std::vector<double>* pCoordX,
                                        std::vector<double>* pCoordY,
                                        std::vector<double>* pCoordZ,
                                        std::vector<double>* pWeight)
{
    assert(pCoordX != NULL);
    assert(pCoordY != NULL);
    assert(pCoordZ != NULL);
    assert(pWeight != NULL);

    const size_t counts = grids.size();
    pCoordX->resize(counts);
    pCoordY->resize(counts);
    pCoordZ->resize(counts);
    pWeight->resize(counts);

    for (size_t i = 0; i < counts; ++i) {
        const TlPosition pos = grids[i].position;
        (*pCoordX)[i] = pos.x();
        (*pCoordY)[i] = pos.y();
        (*pCoordZ)[i] = pos.z();
        (*pWeight)[i] = grids[i].weight;
    }
}

std::vector<GridDataManager::GridInfo>
GridDataManager::composeGridInfo(const std::vector<double>& coordX,
                                 const std::vector<double>& coordY,
                                 const std::vector<double>& coordZ,
                                 const std::vector<double>& weight)
{
    const size_t counts = coordX.size();
    assert(counts == coordY.size());
    assert(counts == coordZ.size());
    assert(counts == weight.size());

    std::vector<GridInfo> grids(counts);
    for (size_t i = 0; i < counts; ++i) {
        grids[i].position = TlPosition(coordX[i], coordY[i], coordZ[i]);
        grids[i].weight = weight[i];
    }

    return grids;
}

void GridDataManager::setData(const int atomIndex, const ChunkType chunkType,
                              const std::vector<double>& data)
{
    ChunkData chunkData;
    chunkData.data = data;
    chunkData.isStored = false;

//   std::cerr << TlUtils::format("GridDataManager::setData() atom=%d, chunktype=%d, datasize=%d",
//                 atomIndex, (int)chunkType, data.size())
//      << std::endl;

    this->atomChunkData_[atomIndex][chunkType] = chunkData;

    // flush
    this->unflushedDataSize_ += data.size() * sizeof(double);
    if (this->unflushedDataSize_ > this->bufferSize_) {
        //std::cerr << "GridDataManager::setData() data flush." << std::endl;
        AtomChunkDataType::iterator atomItrEnd = this->atomChunkData_.end();
        for (AtomChunkDataType::iterator atomItr = this->atomChunkData_.begin();
                atomItr != atomItrEnd; ++atomItr) {
            const int atomIndex = atomItr->first;
            ChunkDataType::iterator chunkItrEnd = atomItr->second.end();
            for (ChunkDataType::iterator chunkItr = atomItr->second.begin();
                    chunkItr != chunkItrEnd; ++chunkItr) {
                const ChunkType chunkType = chunkItr->first;
                this->writeChunk(atomIndex, chunkType);
            }
        }

        this->atomChunkData_.clear();
        this->unflushedDataSize_ = 0;
    }
}


std::vector<double> GridDataManager::getData(const int atomIndex,
                                             const ChunkType chunkType) const
{
    std::vector<double> answer;

    // まだディスクから読み取れていない場合は読み取る。
    {
        AtomChunkDataType::const_iterator atomChunk = this->atomChunkData_.find(atomIndex);
        if ((atomChunk == this->atomChunkData_.end()) ||
            (atomChunk->second.find(chunkType) == atomChunk->second.end())) {
            this->readChunk(atomIndex, chunkType);
        }
    }

    bool readFlag = false;
    AtomChunkDataType::const_iterator atomChunk = this->atomChunkData_.find(atomIndex);
    if (atomChunk != this->atomChunkData_.end()) {
        ChunkDataType::const_iterator chunk = atomChunk->second.find(chunkType);
        if (chunk != atomChunk->second.end()) {
            answer = chunk->second.data;
            readFlag = true;
        }
    }
    if (readFlag == false) {
        std::cerr << "could not read atom=" << atomIndex << ", type=" << chunkType << std::endl;
    }

    return answer;
}


void GridDataManager::readChunk(const int atomIndex, const ChunkType chunkType) const
{
    std::vector<double> answer;

    std::fstream::pos_type pos = -1;
    {
        AtomChunkHeaderType::const_iterator atomItr = this->atomChunkHeader_.find(atomIndex);
        if (atomItr != this->atomChunkHeader_.end()) {
            ChunkHeaderType::const_iterator chunkItr = atomItr->second.find(chunkType);
            if (chunkItr != atomItr->second.end()) {
                pos = chunkItr->second.filePosition;
            }
        }
    }
    if (pos == (std::fstream::pos_type)-1) {
//     std::cerr << TlUtils::format("GridDataManager::readChunk() unreadable header. atomIndex=%d, chunkType=%d",
//               atomIndex, chunkType)
//        << std::endl;
        return;
    }

    this->fs_.seekg(pos, std::ios::beg);

    // header
    int diskAtomIndex = 0;
    ChunkType diskChunkType = UNDEFINED;
    size_t bufferSize = 0;
    this->fs_.read(reinterpret_cast<char*>(&diskAtomIndex), sizeof(int));
    this->fs_.read(reinterpret_cast<char*>(&diskChunkType), sizeof(ChunkType));
    this->fs_.read(reinterpret_cast<char*>(&bufferSize), sizeof(size_t));
    //   std::cout << "readChunk: atom=" << diskAtomIndex
    //        << ", chunk=" << diskChunkType
    //        << ", bufferSize=" << bufferSize
    //        << std::endl;
    assert(atomIndex == diskAtomIndex);
    assert(chunkType == diskChunkType);
    assert(bufferSize == const_cast<GridDataManager*>(this)->atomChunkHeader_[atomIndex][chunkType].bufferSize);

    // contents
    size_t count = bufferSize / sizeof(double);
    std::vector<double> data(count, 0.0);
    this->fs_.read(reinterpret_cast<char*>(&(data[0])), sizeof(double) * count);

    //this->setData(atomIndex, chunkType, data);
    {
        ChunkData chunkData;
        chunkData.data = data;
        chunkData.isStored = true;

        this->atomChunkData_[atomIndex][chunkType] = chunkData;
    }
}


void GridDataManager::writeChunk(int atomIndex, ChunkType chunkType)
{
    if ((this->atomChunkData_.find(atomIndex) == this->atomChunkData_.end()) ||
            (this->atomChunkData_[atomIndex].find(chunkType) == this->atomChunkData_[atomIndex].end()) ||
            (this->atomChunkData_[atomIndex][chunkType].isStored == true)) {
        // need not write!
        return;
    }
    size_t needBufferSize = sizeof(double) * this->atomChunkData_[atomIndex][chunkType].data.size();

    // start position
    std::fstream::pos_type startPos;
    {
        this->fs_.seekg(0, std::ios::end);
        startPos = this->fs_.tellg();  // 書き込み位置を設定
        if (startPos == (std::fstream::pos_type)-1) {
            startPos = 0;
        }
    }

    AtomChunkHeaderType::const_iterator pAtom = this->atomChunkHeader_.find(atomIndex);
    if (pAtom != this->atomChunkHeader_.end()) {
        ChunkHeaderType::const_iterator pChunkType = this->atomChunkHeader_[atomIndex].find(chunkType);
        if (pChunkType != this->atomChunkHeader_[atomIndex].end()) {
            size_t bufferSize = this->atomChunkHeader_[atomIndex][chunkType].bufferSize;
            if (bufferSize < needBufferSize) {
                // 現ヘッダの無効化
                this->fs_.seekp(this->atomChunkHeader_[atomIndex][chunkType].filePosition);
                this->fs_.write(reinterpret_cast<char*>(&atomIndex), sizeof(int));
                ChunkType invalid = INVALID;
                this->fs_.write(reinterpret_cast<char*>(&invalid), sizeof(ChunkType));
                this->fs_.write(reinterpret_cast<char*>(&bufferSize), sizeof(size_t));
            } else {
                // 上書き用に書き込み位置を変更
                startPos = this->atomChunkHeader_[atomIndex][chunkType].filePosition;
            }
        }
    }

    // header
//   std::cout << "writeChunk: atom=" << atomIndex
//      << ", chunk=" << chunkType
//      << ", pos=" << startPos
//      << ", buf=" << needBufferSize
//      << std::endl;
    std::ios_base::iostate state = this->fs_.rdstate();
    if ((state & std::ios_base::eofbit) || (state & std::ios_base::failbit)) {
        this->fs_.clear();
    }
    assert(!this->fs_.fail());

    this->fs_.seekp(startPos, std::ios::beg);
    this->fs_.write(reinterpret_cast<char*>(&atomIndex), sizeof(int));
    this->fs_.write(reinterpret_cast<char*>(&chunkType), sizeof(ChunkType));
    this->fs_.write(reinterpret_cast<char*>(&needBufferSize), sizeof(size_t));

    // contents
    this->fs_.write(reinterpret_cast<char*>(&(this->atomChunkData_[atomIndex][chunkType].data[0])),
                    sizeof(double) * this->atomChunkData_[atomIndex][chunkType].data.size());
    this->fs_.flush();

    // update
    ChunkHeader chunkHeader;
    chunkHeader.filePosition = startPos;
    chunkHeader.bufferSize = needBufferSize;
    this->atomChunkHeader_[atomIndex][chunkType] = chunkHeader;
    this->atomChunkData_[atomIndex][chunkType].isStored = true;
}

std::string GridDataManager::str() const
{
    std::string str = "";

    AtomChunkHeaderType::const_iterator pAtomEnd = this->atomChunkHeader_.end();
    for (AtomChunkHeaderType::const_iterator pAtom = this->atomChunkHeader_.begin();
            pAtom != pAtomEnd; ++pAtom) {
        const int atom = pAtom->first;
        ChunkHeaderType::const_iterator pChunkEnd = pAtom->second.end();
        for (ChunkHeaderType::const_iterator pChunk = pAtom->second.begin();
                pChunk != pChunkEnd; ++pChunk) {
            const ChunkType chunkType = pChunk->first;
            //const ChunkHeader chunkHeader = pChunk->second;
            if ((0 <= chunkType) && (chunkType <= 16)) {
                str += TlUtils::format(">>>> atom=%d, chunkType=%s\n",
                                       atom, chunkTypeStr[chunkType].c_str());
                const TlVector v(this->getData(atom, chunkType));
                std::stringstream ss;
                v.print(ss);
                str += ss.str();
                str += "\n";
            }

        }
    }

    return str;
}



