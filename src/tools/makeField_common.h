#ifndef MAKEFIELD_COMMON_H
#define MAKEFIELD_COMMON_H

#include <vector>
#include <string>

#include "TlSerializeData.h"
#include "TlPosition.h"
#include "TlAtom.h"

static const double ANG_PER_AU = 0.529177249;
static const double AU_PER_ANG = 1.0 / ANG_PER_AU;
static const double MIN_LENGTH = 4.0;
static const double GRID_PITCH = 0.5 * AU_PER_ANG; // 0.5A

void getDefaultSize(const TlSerializeData& param,
                    TlPosition* pStartPos,
                    TlPosition* pEndPos);

void compensateRange(TlPosition* pStartPos,
                     TlPosition* pEndPos);
void compensateRange(double* pMinValue,
                     double* pMaxValue);

std::vector<TlPosition> makeGrids(const TlPosition& startPos,
                                  const TlPosition& endPos,
                                  const TlPosition& gridPitch,
                                  int* pNumOfGridX = NULL,
                                  int* pNumOfGridY = NULL,
                                  int* pNumOfGridZ = NULL);

void saveFieldData(const std::size_t numOfGridX,
                   const std::size_t numOfGridY,
                   const std::size_t numOfGridZ,
                   const std::vector<TlPosition>& grids,
                   const std::vector<std::vector<double> >& data,
                   const std::string& label,
                   const std::string& filePath);

void saveCubeData(const std::vector<TlAtom>& atoms,
                  const std::size_t numOfGridX,
                  const std::size_t numOfGridY,
                  const std::size_t numOfGridZ,
                  const TlPosition& startPos,
                  const TlPosition& gridPitch,
                  const std::vector<double>& data,
                  const std::string& label,
                  const std::string& filePath);

#endif // MAKEFIELD_COMMON_H
