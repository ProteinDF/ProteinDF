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

#include "DfIntegrals_Parallel.h"
#include "DfCD_Parallel.h"
#include "DfEriX_Parallel.h"
#include "DfGenerateGrid_Parallel.h"
#include "DfHpqX_Parallel.h"
#include "DfInvMatrix_Parallel.h"
#include "DfOverlapX_Parallel.h"
#include "DfXMatrix_Parallel.h"
#include "Fl_Geometry.h"
#include "TlCommunicate.h"

DfIntegrals_Parallel::DfIntegrals_Parallel(TlSerializeData* pParam)
    : DfIntegrals(pParam) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  DfObject::rank_ = rComm.getRank();
  //     std::cerr << TlUtils::format("[%d] DfIntegrals_Parallel constructor
  //     called.", rComm.getRank())
  //               << std::endl;
}

DfIntegrals_Parallel::~DfIntegrals_Parallel() {
  //     TlCommunicate& rComm = TlCommunicate::getInstance();
  //     std::cerr << TlUtils::format("[%d] DfIntegrals_Parallel destructor
  //     called.", rComm.getRank())
  //               << std::endl;
}

void DfIntegrals_Parallel::logger(const std::string& str) const {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  if (rComm.isMaster() == true) {
    DfIntegrals::logger(str);
  }
}

DfCD* DfIntegrals_Parallel::getDfCDObject() {
  DfCD* pDfCD = new DfCD_Parallel(this->pPdfParam_);
  return pDfCD;
}

DfXMatrix* DfIntegrals_Parallel::getDfXMatrixObject() {
  DfXMatrix* pDfXMatrix = new DfXMatrix_Parallel(this->pPdfParam_);

  return pDfXMatrix;
}

DfInvMatrix* DfIntegrals_Parallel::getDfInvMatrixObject() {
  DfInvMatrix* pDfInvMatrix = new DfInvMatrix_Parallel(this->pPdfParam_);

  return pDfInvMatrix;
}

DfGenerateGrid* DfIntegrals_Parallel::getDfGenerateGridObject() {
  DfGenerateGrid* pDfGenerateGrid =
      new DfGenerateGrid_Parallel(this->pPdfParam_);

  return pDfGenerateGrid;
}

void saveInvSquareVMatrix(const TlDenseSymmetricMatrix_Lapack& v) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster()) {
    v.save("fl_Work/fl_Mtr_invSquareV.matrix");
  }
  rComm.barrier();
}

void DfIntegrals_Parallel::outputStartTitle(const std::string& stepName,
                                            const char lineChar) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfIntegrals::outputStartTitle(stepName, lineChar);
  }
}

void DfIntegrals_Parallel::outputEndTitle(const std::string& stepName,
                                          const char lineChar) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    DfIntegrals::outputEndTitle(stepName, lineChar);
  }
}

void DfIntegrals_Parallel::createHpqMatrix() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->createHpqMatrix_ScaLAPACK();
    return;
  }
#endif  // HAVE_SCALAPACK

  this->createHpqMatrix_LAPACK();
}

void DfIntegrals_Parallel::createHpqMatrix_LAPACK() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  // const std::size_t needMem =
  //     this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) * sizeof(double);
  // if (this->procMaxMemSize_ < needMem) {
  //   this->logger(" H_pq is build up on disk.\n");
  //   TlMatrix::useMemManager(true);
  // } else {
  //   this->logger(" H_pq is build up on memory.\n");
  //   TlMatrix::useMemManager(false);
  // }
  TlDenseSymmetricMatrix_Lapack Hpq(this->m_nNumOfAOs);
  TlDenseSymmetricMatrix_Lapack Hpq2(this->m_nNumOfAOs);

  DfHpqX_Parallel dfHpq(this->pPdfParam_);
  dfHpq.getHpq(&Hpq, &Hpq2);
  // if (this->isUseNewEngine_ == true) {
  //     this->logger(" use new engine\n");
  //     DfHpqX_Parallel dfHpq(this->pPdfParam_);
  //     dfHpq.getHpq(&Hpq, &Hpq2);
  // } else {
  //     DfHpq_Parallel dfHpq(this->pPdfParam_);
  //     dfHpq.getHpq(&Hpq, &Hpq2);
  // }

  if (rComm.isMaster() == true) {
    this->saveHpqMatrix(Hpq);
    this->saveHpq2Matrix(Hpq2);
  }
}

void DfIntegrals_Parallel::createHpqMatrix_ScaLAPACK() {
  TlDenseSymmetricMatrix_Scalapack Hpq(this->m_nNumOfAOs);
  TlDenseSymmetricMatrix_Scalapack Hpq2(this->m_nNumOfAOs);

  this->logger(" Hpq build using distribute matrix.\n");
  DfHpqX_Parallel dfHpq(this->pPdfParam_);
  dfHpq.getHpqD(&Hpq, &Hpq2);
  // if (this->isUseNewEngine_ == true) {
  //     this->logger(" use new engine\n");
  //     DfHpqX_Parallel dfHpq(this->pPdfParam_);
  //     dfHpq.getHpqD(&Hpq, &Hpq2);
  // } else {
  //     DfHpq_Parallel dfHpq(this->pPdfParam_);
  //     dfHpq.getHpq(&Hpq, &Hpq2);
  // }

  this->saveHpqMatrix(Hpq);
  this->saveHpq2Matrix(Hpq2);
}

void DfIntegrals_Parallel::createOverlapMatrix() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->createOverlapMatrix_ScaLAPACK();
    return;
  }
#endif  // HAVE_SCALAPACK

  this->createOverlapMatrix_LAPACK();
}

void DfIntegrals_Parallel::createOverlapMatrix_LAPACK() {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  unsigned int calcState =
      (*this->pPdfParam_)["control"]["integrals_state"].getUInt();
  DfOverlapX_Parallel dfOverlapX(this->pPdfParam_);

  // Spq
  if ((calcState & DfIntegrals::Spq) == 0) {
    this->outputStartTitle("Spq");

    std::size_t needMem =
        this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2 * sizeof(double);
    needMem += rComm.getWorkMemSize();
    // if (this->procMaxMemSize_ < needMem) {
    //   this->logger(" S_(p q) is build on disk.\n");
    //   TlMatrix::useMemManager(true);
    // } else {
    //   TlMatrix::useMemManager(false);
    //   this->logger(" S_(p q) is build on memory.\n");
    // }

    TlDenseSymmetricMatrix_Lapack Spq(this->m_nNumOfAOs);
    dfOverlapX.getSpq(&Spq);
    // if (this->isUseNewEngine_ == true) {
    //     this->logger(" use new engine.\n");
    //     dfOverlapX.getSpq(&Spq);
    // } else {
    //     dfOverlap.getSpq(&Spq);
    // }

    if (rComm.isMaster() == true) {
      this->saveSpqMatrix(Spq);
    }
    rComm.barrier();

    this->outputEndTitle();

    calcState |= DfIntegrals::Spq;
    (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
    this->saveParam();
  }

  if (this->K_engine_ == K_ENGINE_RI_K) {
    // Sgd
    if ((calcState & DfIntegrals::Sgd) == 0) {
      if (this->m_bIsXCFitting == true) {
        this->outputStartTitle("Sgd");

        // const std::size_t needMem =
        //     this->numOfAuxXC_ * (this->numOfAuxXC_ + 1) / 2 * sizeof(double);
        // if (this->procMaxMemSize_ < needMem) {
        //   this->logger(" S_(gamma delta) is build on disk.\n");
        //   TlMatrix::useMemManager(true);
        // } else {
        //   this->logger(" S_(gamma delta) is build on memory.\n");
        //   TlMatrix::useMemManager(false);
        // }

        TlDenseSymmetricMatrix_Lapack Sgd(this->numOfAuxXC_);
        dfOverlapX.getSgd(&Sgd);

        if (rComm.isMaster() == true) {
          this->saveSgdMatrix(Sgd);
        }
        rComm.barrier();

        this->outputEndTitle();
      }

      calcState |= DfIntegrals::Sgd;
      (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
      this->saveParam();
    }
  }

  if (this->J_engine_ == J_ENGINE_RI_J) {
    // Sab2
    if ((calcState & DfIntegrals::Sab2) == 0) {
      this->outputStartTitle("Sab2");

      // const std::size_t needMem =
      //     this->m_nNumOfAux * (this->m_nNumOfAux + 1) / 2 * sizeof(double);
      // if (this->procMaxMemSize_ < needMem) {
      //   this->logger(" S_(alpha beta) is build on disk.\n");
      //   TlMatrix::useMemManager(true);
      // } else {
      //   this->logger(" S_(alpha beta) is build on memory.\n");
      //   TlMatrix::useMemManager(false);
      // }

      TlDenseSymmetricMatrix_Lapack Sab2(this->m_nNumOfAux);
      dfOverlapX.getSab(&Sab2);
      // if (this->isUseNewEngine_ == true) {
      //     this->logger(" use new engine.\n");
      //     dfOverlapX.getSab(&Sab2);
      // } else {
      //     dfOverlap.getSab2(&Sab2);
      // }

      if (rComm.isMaster() == true) {
        this->saveSab2Matrix(Sab2);
      }
      rComm.barrier();

      this->outputEndTitle();

      calcState |= DfIntegrals::Sab2;
      (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
      this->saveParam();
    }

    // Na
    if ((calcState & DfIntegrals::Na) == 0) {
      this->outputStartTitle("N_alpha");

      // const std::size_t needMem = this->m_nNumOfAux * sizeof(double);
      // if (this->procMaxMemSize_ < needMem) {
      //   this->logger(" [alpha] is build on disk.\n");
      //   TlMatrix::useMemManager(true);
      // } else {
      //   this->logger(" [alpha] is build on memory.\n");
      //   TlMatrix::useMemManager(false);
      // }

      TlDenseVector_Lapack Na(this->m_nNumOfAux);
      dfOverlapX.getNalpha(&Na);
      this->saveNalpha(Na);

      this->outputEndTitle();

      calcState |= DfIntegrals::Na;
      (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
      this->saveParam();
    }
  }
}

void DfIntegrals_Parallel::createOverlapMatrix_ScaLAPACK() {
  // TlCommunicate& rComm = TlCommunicate::getInstance();
  unsigned int calcState =
      (*this->pPdfParam_)["control"]["integrals_state"].getUInt();
  DfOverlapX_Parallel dfOverlapX(this->pPdfParam_);

  // Spq
  if ((calcState & DfIntegrals::Spq) == 0) {
    this->outputStartTitle("Spq");

    this->logger(" S_(p q) is build using distribute matrix.\n");
    TlDenseSymmetricMatrix_Scalapack Spq(this->m_nNumOfAOs);
    dfOverlapX.getSpqD(&Spq);
    // if (this->isUseNewEngine_ == true) {
    //     this->logger(" use new engine.\n");
    //     dfOverlapX.getSpqD(&Spq);
    // } else {
    //     dfOverlap.getSpq(&Spq);
    // }
    this->saveSpqMatrix(Spq);

    this->outputEndTitle();

    calcState |= DfIntegrals::Spq;
    (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
    this->saveParam();
  }

  if (this->K_engine_ == K_ENGINE_RI_K) {
    // Sgd
    if ((calcState & DfIntegrals::Sgd) == 0) {
      if (this->m_bIsXCFitting == true) {
        this->outputStartTitle("Sgd");

        TlDenseSymmetricMatrix_Scalapack Sgd(this->numOfAuxXC_);
        dfOverlapX.getSgd(&Sgd);
        this->saveSgdMatrix(Sgd);

        this->outputEndTitle();
      }

      calcState |= DfIntegrals::Sgd;
      (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
      this->saveParam();
    }
  }

  if (this->J_engine_ == J_ENGINE_RI_J) {
    // Sab2
    if ((calcState & DfIntegrals::Sab2) == 0) {
      this->outputStartTitle("Sab2");

      this->logger(" S_(alpha beta) is build using on distribute matrix.\n");
      TlDenseSymmetricMatrix_Scalapack Sab2(this->m_nNumOfAux);
      dfOverlapX.getSabD(&Sab2);
      // if (this->isUseNewEngine_ == true) {
      //     this->logger(" use new engine.\n");
      //     dfOverlapX.getSabD(&Sab2);
      // } else {
      //     dfOverlap.getSab2(&Sab2);
      // }
      this->saveSab2Matrix(Sab2);

      this->outputEndTitle();

      calcState |= DfIntegrals::Sab2;
      (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
      this->saveParam();
    }

    // Na
    if ((calcState & DfIntegrals::Na) == 0) {
      this->outputStartTitle("N_alpha");

      // const std::size_t needMem = this->m_nNumOfAux * sizeof(double);
      // if (this->procMaxMemSize_ < needMem) {
      //   this->logger(" [alpha] is build on disk.\n");
      //   TlMatrix::useMemManager(true);
      // } else {
      //   this->logger(" [alpha] is build on memory.\n");
      //   TlMatrix::useMemManager(false);
      // }

      TlDenseVector_Lapack Na(this->m_nNumOfAux);
      dfOverlapX.getNalpha(&Na);
      this->saveNalpha(Na);

      this->outputEndTitle();

      calcState |= DfIntegrals::Na;
      (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
      this->saveParam();
    }
  }
}

void DfIntegrals_Parallel::createERIMatrix() {
#ifdef HAVE_SCALAPACK
  if (this->m_bUsingSCALAPACK == true) {
    this->createERIMatrix_ScaLAPACK();
    return;
  }
#endif  // HAVE_SCALAPACK

  this->createERIMatrix_LAPACK();
}

void DfIntegrals_Parallel::createERIMatrix_LAPACK() {
  unsigned int calcState =
      (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

  if (this->J_engine_ == J_ENGINE_RI_J) {
    if ((calcState & DfIntegrals::Sab) == 0) {
      this->outputStartTitle("Sab");

      // const std::size_t needMem =
      //     this->m_nNumOfAux * (this->m_nNumOfAux + 1) / 2 * sizeof(double);
      // if (this->procMaxMemSize_ < needMem) {
      //   this->logger(" <alpha|beta> is build on disk.\n");
      //   TlMatrix::useMemManager(true);
      // } else {
      //   this->logger(" <alpha|beta> is build on memory.\n");
      //   TlMatrix::useMemManager(false);
      // }

      TlDenseSymmetricMatrix_Lapack Sab(this->m_nNumOfAux);

      DfEriX_Parallel dfEri(this->pPdfParam_);
      dfEri.getJab(&Sab);
      // if (this->isUseNewEngine_ == true) {
      //     this->logger(" use new engine.\n");
      //     DfEriX_Parallel dfEri(this->pPdfParam_);
      //     dfEri.getJab(&Sab);
      // } else {
      //     DfEri_Parallel dfEri(this->pPdfParam_);
      //     dfEri.getSab(&Sab);
      // }

      TlCommunicate& rComm = TlCommunicate::getInstance();
      if (rComm.isMaster() == true) {
        this->saveSabMatrix(Sab);
      }

      this->outputEndTitle();

      calcState |= DfIntegrals::Sab;
      (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
      this->saveParam();
    }
  }
}

void DfIntegrals_Parallel::createERIMatrix_ScaLAPACK() {
  unsigned int calcState =
      (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

  if (this->J_engine_ == J_ENGINE_RI_J) {
    if ((calcState & DfIntegrals::Sab) == 0) {
      this->outputStartTitle("Sab");

      DfEriX_Parallel dfEri(this->pPdfParam_);
      TlDenseSymmetricMatrix_Scalapack Jab(this->m_nNumOfAux);
      dfEri.getJab(&Jab);
      this->saveSabMatrix(Jab);

      this->outputEndTitle();

      calcState |= DfIntegrals::Sab;
      (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
      this->saveParam();
    }
  }
}
