#include <iostream>
#include <fstream>
#include "makeField_common.h"
#include "TlSerializeData.h"
#include "Fl_Geometry.h"

void getDefaultSize(const TlSerializeData& param,
                    TlPosition* pStartPos,
                    TlPosition* pEndPos)
{
    // 計算サイズの決定
    Fl_Geometry flGeom(param["model"]["coordinates"]); // 単位はa.u.
    double min_x = 0.0;
    double min_y = 0.0;
    double min_z = 0.0;
    double max_x = 0.0;
    double max_y = 0.0;
    double max_z = 0.0;
    const std::size_t numOfAtoms = flGeom.getNumOfAtoms();
    // define start value
    if (numOfAtoms > 0) {
        const TlPosition p = flGeom.getCoordinate(0);
        min_x = max_x = p.x();
        min_y = max_y = p.y();
        min_z = max_z = p.z();
    }
    for (std::size_t i = 0; i < numOfAtoms; ++i) {
        const TlPosition pos = flGeom.getCoordinate(i);
        const double x = pos.x();
        const double y = pos.y();
        const double z = pos.z();
        min_x = std::min(min_x, x);
        min_y = std::min(min_y, y);
        min_z = std::min(min_z, z);
        max_x = std::max(max_x, x);
        max_y = std::max(max_y, y);
        max_z = std::max(max_z, z);
    }

    const double buffer = 10.0 / AU_PER_ANG; // 10A
    min_x -= buffer;
    max_x += buffer;
    min_y -= buffer;
    max_y += buffer;
    min_z -= buffer;
    max_z += buffer;

    if (pStartPos != NULL) {
        *pStartPos = TlPosition(min_x, min_y, min_z);
    }
    if (pEndPos != NULL) {
        *pEndPos = TlPosition(max_x, max_y, max_z);
    }
}


void compensateRange(TlPosition* pStartPos,
                     TlPosition* pEndPos)
{
    assert(pStartPos != NULL);
    assert(pEndPos != NULL);

    double sx = pStartPos->x();
    double sy = pStartPos->y();
    double sz = pStartPos->z();
    double ex = pEndPos->x();
    double ey = pEndPos->y();
    double ez = pEndPos->z();

    compensateRange(&sx, &ex);
    compensateRange(&sy, &ey);
    compensateRange(&sz, &ez);

    pStartPos->moveTo(sx, sy, sz);
    pEndPos->moveTo(ex, ey, ez);
}


void compensateRange(double* pMinValue,
                     double* pMaxValue)
{
    assert(pMinValue != NULL);
    assert(pMaxValue != NULL);

    const double distance = *pMaxValue - *pMinValue;
    if (distance < MIN_LENGTH) {
        const double center = (*pMinValue + *pMaxValue) / 2.0;
        const double half_length = MIN_LENGTH / 2.0;
        *pMinValue = center - half_length;
        *pMaxValue = center + half_length;
    }
}



std::vector<TlPosition> makeGrids(const TlPosition& startPos,
                                  const TlPosition& endPos,
                                  const TlPosition& gridPitch,
                                  int* pNumOfGridX,
                                  int* pNumOfGridY,
                                  int* pNumOfGridZ)
{
    const TlPosition range = endPos - startPos;
    const int numOfGridX = range.x() / gridPitch.x();
    const int numOfGridY = range.y() / gridPitch.y();
    const int numOfGridZ = range.z() / gridPitch.z();

    const double sx = startPos.x();
    const double sy = startPos.y();
    const double sz = startPos.z();
    const double dx = gridPitch.x();
    const double dy = gridPitch.y();
    const double dz = gridPitch.z();

    // i=x, j=y, k=zなので、AVS FieldDataの扱いとしては
    // はじめにi(=x)が増加するように計算しなければならない    
    std::size_t index = 0;
    std::vector<TlPosition> answer(numOfGridX * numOfGridY * numOfGridZ);
    for (int z = 0; z < numOfGridZ; ++z) {
        const double pz = sz + dz * z;

        for (int y = 0; y < numOfGridY; ++y) {
            const double py = sy + dy * y;
            
            for (int x = 0; x < numOfGridX; ++x) {
                const double px = sx + dx * x;
                
                answer[index] = TlPosition(px, py, pz);
                ++index;
            }
        }
    }

    if (pNumOfGridX != NULL) {
        *pNumOfGridX = numOfGridX;
    }
    if (pNumOfGridY != NULL) {
        *pNumOfGridY = numOfGridY;
    }
    if (pNumOfGridZ != NULL) {
        *pNumOfGridZ = numOfGridZ;
    }
    
    return answer;
}


void saveFieldData(const std::size_t numOfGridX,
                   const std::size_t numOfGridY,
                   const std::size_t numOfGridZ,
                   const std::vector<TlPosition>& grids,
                   const std::vector<std::vector<double> >& data,
                   const std::string& label,
                   const std::string& filePath)
{
    static const int nspace = 3;
    const int veclen = data.size();
    const std::size_t numOfGrids = grids.size();
    
    // data
    float* pData = new float[numOfGrids * veclen];
    for (int v = 0; v < veclen; ++v) {
        float max_value = data[v][0];
        float min_value = max_value;
        for (std::size_t i = 0; i < numOfGrids; ++i) {
            const float value = data[v][i];
            pData[i * veclen + v] = value;
            max_value = std::max(max_value, value);
            min_value = std::min(min_value, value);
        }
        std::cerr << "veclen: " << v << std::endl;
        std::cerr << " max value = " << max_value << std::endl;
        std::cerr << " min value = " << min_value << std::endl;
    }

    // coordination
    float* pCoord = new float[numOfGrids * 3];
    for (std::size_t i = 0; i < numOfGrids; ++i) {
        const TlPosition& p = grids[i];
        pCoord[               + i] = p.x();
        pCoord[numOfGrids     + i] = p.y();
        pCoord[numOfGrids * 2 + i] = p.z();
    }


    std::ofstream ofs;
    // header
    ofs.open(filePath.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
    ofs << "# AVS field file\n";
    ofs << "ndim = 3\n";
    ofs << TlUtils::format("dim1 = %ld\n", numOfGridX);
    ofs << TlUtils::format("dim2 = %ld\n", numOfGridY);
    ofs << TlUtils::format("dim3 = %ld\n", numOfGridZ);
    ofs << TlUtils::format("nspace = %d\n", nspace);
    ofs << TlUtils::format("veclen = %d\n", veclen);
    ofs << "data = float\n";
    ofs << "field = irregular\n";
    ofs << TlUtils::format("label = %s\n", label.c_str());
    ofs << "\f\f";

    // 本体
    ofs.write((char*)pData, sizeof(float) * numOfGrids * veclen);
    ofs.write((char*)pCoord, sizeof(float) * numOfGrids * 3);

    ofs.close();

    delete[] pData;
    pData = NULL;
    delete[] pCoord;
    pCoord = NULL;
}


void saveCubeData(const std::vector<TlAtom>& atoms,
                  const std::size_t numOfGridX,
                  const std::size_t numOfGridY,
                  const std::size_t numOfGridZ,
                  const TlPosition& startPos,
                  const TlPosition& gridPitch,
                  const std::vector<double>& data,
                  const std::string& label,
                  const std::string& filePath)
{
    const std::size_t numOfGrids = data.size();
    const int numOfAtoms = atoms.size();
    
    std::ofstream ofs;
    ofs.open(filePath.c_str(), std::ios::out | std::ios::trunc);
    // header
    ofs << "comemnt line:\n";
    ofs << TlUtils::format("comment line:%s\n", label.c_str());
    ofs << TlUtils::format("%5d%12.6f%12.6f%12.6f\n", numOfAtoms, startPos.x(), startPos.y(), startPos.z()); // 開始点
    ofs << TlUtils::format("%5d%12.6f%12.6f%12.6f\n", numOfGridZ, 0.0, 0.0, gridPitch.z());
    ofs << TlUtils::format("%5d%12.6f%12.6f%12.6f\n", numOfGridY, 0.0, gridPitch.y(), 0.0);
    ofs << TlUtils::format("%5d%12.6f%12.6f%12.6f\n", numOfGridX, gridPitch.x(), 0.0, 0.0);

    // 原子
    for (int i = 0; i < numOfAtoms; ++i) {
        ofs << TlUtils::format("%5d%12.6f%12.6f%12.6f%12.6f\n",
                               TlAtom::getElementNumber(atoms[i].getSymbol()),
                               atoms[i].getCharge(),
                               atoms[i].getPosition().x(), atoms[i].getPosition().y(), atoms[i].getPosition().z());
    }

    // 物理量
    const int numOfOutputGridGroups = numOfGrids / 6;
    for (int i = 0; i < numOfOutputGridGroups; ++i) {
        const int base = 6 * i;
        ofs << TlUtils::format("% 12e % 12e % 12e % 12e % 12e % 12e\n",
                               data[base], data[base +1], data[base +2],
                               data[base +3], data[base +4], data[base +5]);
    }
    const int numOfOutputRests = numOfGrids - numOfOutputGridGroups * 6;
    for (int i = 0; i < numOfOutputRests; ++i) {
        ofs << TlUtils::format("% 12e ", data[numOfOutputGridGroups * 6 + i]);
    }
    ofs << "\n";
    
    ofs.close();
}


