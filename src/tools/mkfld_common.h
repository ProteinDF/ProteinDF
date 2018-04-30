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

#ifndef MAKEFIELD_COMMON_H
#define MAKEFIELD_COMMON_H

#include <string>
#include <vector>

#include "TlAtom.h"
#include "TlPosition.h"
#include "TlSerializeData.h"

static const double ANG_PER_AU = 0.529177249;
static const double AU_PER_ANG = 1.0 / ANG_PER_AU;
static const double MIN_LENGTH = 4.0;
static const double GRID_PITCH = 0.5 * AU_PER_ANG;  // 0.5A

void getDefaultSize(const TlSerializeData& param, TlPosition* pStartPos,
                    TlPosition* pEndPos);

void compensateRange(TlPosition* pStartPos, TlPosition* pEndPos);
void compensateRange(double* pMinValue, double* pMaxValue);

std::vector<TlPosition> makeGrids(const TlPosition& startPos,
                                  const TlPosition& endPos,
                                  const TlPosition& gridPitch,
                                  int* pNumOfGridX = NULL,
                                  int* pNumOfGridY = NULL,
                                  int* pNumOfGridZ = NULL);

void saveFieldData(const std::size_t numOfGridX, const std::size_t numOfGridY,
                   const std::size_t numOfGridZ,
                   const std::vector<TlPosition>& grids,
                   const std::vector<std::vector<double> >& data,
                   const std::string& label, const std::string& filePath);

/// output Gaussian cube file format
///
/// startPos, gridPitchはBohr単位で入力すること。
void saveCubeData(const std::vector<TlAtom>& atoms,
                  const std::size_t numOfGridX, const std::size_t numOfGridY,
                  const std::size_t numOfGridZ, const TlPosition& startPos,
                  const TlPosition& gridPitch, const std::vector<double>& data,
                  const std::string& label, const std::string& filePath);

#endif  // MAKEFIELD_COMMON_H
