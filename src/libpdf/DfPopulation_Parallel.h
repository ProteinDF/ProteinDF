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

#ifndef DFPOPULATION_PARALLEL_H
#define DFPOPULATION_PARALLEL_H

#include "DfPopulation.h"
#include "TlCommunicate.h"

class TlDenseSymmetricMatrix_Lapack;
class TlDenseSymmetricMatrix_Scalapack;

class DfPopulation_Parallel : public DfPopulation {
   public:
    DfPopulation_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfPopulation_Parallel();

   public:
    template <class SymmetricMatrixType>
    double getSumOfElectrons(const SymmetricMatrixType& P);

   protected:
    virtual void calcPop(const int iteration);
};

template <class SymmetricMatrixType>
double DfPopulation_Parallel::getSumOfElectrons(const SymmetricMatrixType& P) {
    double answer = 0.0;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        answer = DfPopulation::getSumOfElectrons(P);
    }
    rComm.broadcast(answer);

    return answer;
}

#endif  // DFPOPULATION_PARALLEL_H
