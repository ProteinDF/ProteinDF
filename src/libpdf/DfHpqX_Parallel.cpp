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

#include "DfHpqX_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"
#include "tl_dense_symmetric_matrix_blacs.h"

DfHpqX_Parallel::DfHpqX_Parallel(TlSerializeData* pPdfParam)
    : DfHpqX(pPdfParam) {}

DfHpqX_Parallel::~DfHpqX_Parallel() {}

void DfHpqX_Parallel::logger(const std::string& str) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfHpqX::logger(str);
  }
}

void DfHpqX_Parallel::getHpqD(TlDenseSymmetricMatrix_blacs* pHpq,
                              TlDenseSymmetricMatrix_blacs* pHpq2) {
  const int numOfAOs = this->m_nNumOfAOs;

  // make coordinates
  const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);
  const int numOfAtoms = flGeom.getNumOfAtoms();
  const int numOfDummyAtoms = flGeom.getNumOfDummyAtoms();
  const int numOfRealAtoms = numOfAtoms - numOfDummyAtoms;
  std::vector<TlAtom> Cs(numOfRealAtoms);
  std::vector<TlAtom> Xs(numOfDummyAtoms);
  std::size_t realAtomIndex = 0;
  std::size_t dummyAtomIndex = 0;
  for (int i = 0; i < numOfAtoms; ++i) {
    const std::string atomName = flGeom.getAtomSymbol(i);
    const TlPosition p = flGeom.getCoordinate(i);
    const double charge = flGeom.getCharge(i);
    const TlAtom atom(atomName, p, charge);
    if (atomName == "X") {
      Xs[dummyAtomIndex] = atom;
      ++dummyAtomIndex;
    } else {
      Cs[realAtomIndex] = atom;
      ++realAtomIndex;
    }
  }
  assert(realAtomIndex == Cs.size());
  assert(dummyAtomIndex == Xs.size());

  DfHpqEngine engine;

  pHpq->resize(numOfAOs);
  pHpq2->resize(numOfAOs);
  TlSparseSymmetricMatrix tmpHpq(numOfAOs);
  TlSparseSymmetricMatrix tmpHpq2(numOfAOs);

  this->createEngines();
  DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
  std::vector<DfTaskCtrl::Task2> taskList;
  bool hasTask = pTaskCtrl->getQueue2(this->orbitalInfo_, true,
                                      this->grainSize_, &taskList, true);
  while (hasTask == true) {
    this->getHpq_part(this->orbitalInfo_, taskList, Cs, Xs, &tmpHpq, &tmpHpq2);

    hasTask = pTaskCtrl->getQueue2(this->orbitalInfo_, true, this->grainSize_,
                                   &taskList);
  }

  if (this->chargeExtrapolateNumber_ > 0) {
    tmpHpq2 /= static_cast<double>(this->chargeExtrapolateNumber_);
  }

  this->loggerTime(" finalize");
  pHpq->mergeSparseMatrix(tmpHpq);
  pHpq2->mergeSparseMatrix(tmpHpq2);

  delete pTaskCtrl;
  pTaskCtrl = NULL;
  this->destroyEngines();
}

DfTaskCtrl* DfHpqX_Parallel::getDfTaskCtrlObject() const {
  DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl_Parallel(this->pPdfParam_);
  return pDfTaskCtrl;
}

void DfHpqX_Parallel::finalize(TlDenseSymmetricMatrix_BLAS_Old* pHpq,
                               TlDenseSymmetricMatrix_BLAS_Old* pHpq2) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  rComm.allReduce_SUM(*pHpq);
  rComm.allReduce_SUM(*pHpq2);
}
