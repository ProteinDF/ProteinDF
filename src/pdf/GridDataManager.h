#ifndef GRIDDATAMANAGER_H
#define GRIDDATAMANAGER_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include "TlPosition.h"

class GridDataManager {
public:
    struct GridInfo {
public:
        GridInfo(const TlPosition& p = TlPosition(), double w = 0.0)
                : position(p), weight(w) {
        }

public:
        TlPosition position;
        double weight;
    };

    enum ChunkType {
        INVALID = -1,
        UNDEFINED = 0,
        COORD_X = 1,
        COORD_Y,
        COORD_Z,
        GRID_WEIGHT,
        DENSITY,
        DENSITY_ALPHA,
        DENSITY_BETA,
        GRAD_DENSITY_X,
        GRAD_DENSITY_Y,
        GRAD_DENSITY_Z,
        GRAD_DENSITY_X_ALPHA,
        GRAD_DENSITY_Y_ALPHA,
        GRAD_DENSITY_Z_ALPHA,
        GRAD_DENSITY_X_BETA,
        GRAD_DENSITY_Y_BETA,
        GRAD_DENSITY_Z_BETA
    };

    struct ChunkHeader {
public:
        std::fstream::pos_type filePosition;
        std::size_t bufferSize; // = sizeof(オブジェクト) * 数
    };

    /// chunkのデータ構造
    struct ChunkData {
public:
        ChunkData() : isStored(false) {
        }
public:
        std::vector<double> data;
        bool isStored;
    };

    typedef std::map<ChunkType, ChunkHeader> ChunkHeaderType;
    typedef std::map<int, ChunkHeaderType> AtomChunkHeaderType;

    typedef std::map<ChunkType, ChunkData> ChunkDataType;
    typedef std::map<int, ChunkDataType> AtomChunkDataType;

public:
    GridDataManager(const std::string& fileName);
    ~GridDataManager();

public:
    void setData(int atomIndex, ChunkType chunkType,
                 const std::vector<double>& data);

    std::vector<double> getData(int atomIndex,
                                ChunkType chunkType) const;

    // for DfGenerateGrid
    static void decomposeGridInfo(const std::vector<GridInfo>& grids,
                                  std::vector<double>* pCoordX,
                                  std::vector<double>* pCoordY,
                                  std::vector<double>* pCoordZ,
                                  std::vector<double>* pWeight);
    static std::vector<GridInfo> composeGridInfo(const std::vector<double>& coordX,
                                                 const std::vector<double>& coordY,
                                                 const std::vector<double>& coordZ,
                                                 const std::vector<double>& weight);
    std::string str() const;

private:
    void open();
    void readHeaders();

    void readChunk(int atomIndex,  ChunkType chunkType) const;

    void writeChunk(int atomIndex, ChunkType chunkType);

private:
    static const size_t headerSize;

    /// 保存先ファイル名
    std::string filePath_;

    mutable std::fstream fs_;

    // atomChunkHeader_[atomIndex][chunkType] = chunkHeader
    AtomChunkHeaderType atomChunkHeader_;

    // contents_[atomIndex][chunkType] = chunkData
    mutable AtomChunkDataType atomChunkData_;

    static const std::string chunkTypeStr[];

    static const std::size_t bufferSize_;
    std::size_t unflushedDataSize_;
};

#endif // GRIDDATAMANAGER_H
