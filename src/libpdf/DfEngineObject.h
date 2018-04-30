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

#ifndef DFENGINEOBJECT_H
#define DFENGINEOBJECT_H

#include <vector>
#include "TlLogging.h"
#include "TlOrbitalInfoObject.h"

// s=0, p=1, d=2
#define AM_LEVEL (2)
// 2次微分まで対応
#define AM_GRAD_LEVEL (2)

class DfEngineObject {
 public:
  typedef int index_type;

 public:
  DfEngineObject()
      : primitiveLevelThreshold_(0.0), log_(TlLogging::getInstance()) {}

  void setPrimitiveLevelThreshold(const double threshold) {
    this->primitiveLevelThreshold_ = threshold;
  }
  double getPrimitiveLevelThreshold() const {
    return this->primitiveLevelThreshold_;
  }

  virtual void calc(const int diff1, const TlOrbitalInfoObject& orbInfo1,
                    const index_type shell1, const int diff2,
                    const TlOrbitalInfoObject& orbInfo2,
                    const index_type shell2, const int diff3,
                    const TlOrbitalInfoObject& orbInfo3,
                    const index_type shell3, const int diff4,
                    const TlOrbitalInfoObject& orbInfo4,
                    const index_type shell4) = 0;

 public:
  virtual ~DfEngineObject(){};

 public:
  virtual double value(const index_type index) const = 0;

 protected:
  // see J. Chem. Phys., 105, 2726 (1996), eq33
  double primitiveLevelThreshold_;

  TlLogging& log_;
};

#endif  // DFENGINEOBJECT_H
