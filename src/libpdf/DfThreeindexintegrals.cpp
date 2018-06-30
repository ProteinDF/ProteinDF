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

#include "DfThreeindexintegrals.h"
#include "CnError.h"
#include "DfEriX.h"
#include "DfOverlapX.h"
#include "Fl_Geometry.h"
#include "TlUtils.h"
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_dense_vector_blas.h"

DfThreeindexintegrals::DfThreeindexintegrals(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {}

DfThreeindexintegrals::~DfThreeindexintegrals() {}

void DfThreeindexintegrals::DfThreeindexintegralsMain() {
  // DIRECT SCHEME
  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      if (this->m_bIsXCFitting == true) {
        this->mainDIRECT_RKS(this->m_nIteration);  // RKS
      } else {
        this->mainDIRECT_RKS2(this->m_nIteration);  // RKS
      }
    } break;

    case METHOD_UKS:
      this->mainDIRECT_UKS(this->m_nIteration);  // UKS alpha spin
      break;

    case METHOD_ROKS:
      this->mainDIRECT_ROKS(this->m_nIteration);  // ROKS
      break;

    default:
      CnErr.abort("DfThreeindexintegrals", "", "main", "runtype is illegal");
      break;
  }
}

TlDenseSymmetricMatrix_BLAS_Old DfThreeindexintegrals::getPMatrix(
    const RUN_TYPE runType, int iteration) {
  assert(iteration >= 0);

  TlDenseSymmetricMatrix_BLAS_Old P;
  P.load(this->getPpqMatrixPath(runType, iteration));

  return P;
}

void DfThreeindexintegrals::mainDIRECT_RKS(int iteration) {
  this->log_.info("Direct scheme method is employed");

  // read Rou
  TlVector_BLAS currRho;
  currRho.load(this->getRhoPath(RUN_RKS, iteration));
  assert(this->m_nNumOfAux == currRho.getSize());

  // read Myu
  TlVector_BLAS currMyu;
  currMyu.load(this->getMyuPath(RUN_RKS, iteration));
  assert(this->numOfAuxXC_ == currMyu.getSize());

  // read Epsilon
  // energy for xc energy term (compair myu*Ppq*Pqa with 4/3*Ex1)
  TlVector_BLAS currEps;
  currEps.load("fl_Work/fl_Vct_Epsilon" + TlUtils::xtos(iteration));
  assert(this->numOfAuxXC_ == currEps.getSize());

  // read previous Rou, Myu, Eps & calculate delta
  if (iteration > 1) {  // iteration != 1
    {
      TlVector_BLAS prevRho(this->m_nNumOfAux);
      prevRho.load(this->getRhoPath(RUN_RKS, iteration - 1));
      assert(this->m_nNumOfAux == prevRho.getSize());

      currRho -= prevRho;
    }

    {
      TlVector_BLAS prevMyu;
      prevMyu.load(this->getMyuPath(RUN_RKS, iteration - 1));
      assert(this->numOfAuxXC_ == prevMyu.getSize());

      currMyu -= prevMyu;
    }

    {
      TlVector_BLAS prevEps;
      prevEps.load("fl_Work/fl_Vct_Epsilon" + TlUtils::xtos(iteration - 1));
      assert(this->numOfAuxXC_ == prevEps.getSize());

      currEps -= prevEps;
    }
  }

  // prepare Fock
  TlDenseSymmetricMatrix_BLAS_Old F(this->m_nNumOfAOs);

  // クーロン項
  // [pq|alpha]
  {
    DfEriX dferi(this->pPdfParam_);
    // dferi.getdeltaHpqA(currRho, F);
    dferi.getJ(currRho, &F);
  }

  // prepare E
  TlDenseSymmetricMatrix_BLAS_Old E = F;

  // 交換相関項
  if (this->m_sXCFunctional != "hf") {
    // myu * [pq gamma] for Fock
    // eps * [pq gamma] for E
    {
      DfOverlapX dfovr(this->pPdfParam_);
      // dfovr.getdeltaHpqG(currMyu, currEps, F, E);
      dfovr.get_pqg(currMyu, currEps, &F, &E);
    }
  } else {
    // Ex using HF
    TlDenseSymmetricMatrix_BLAS_Old P(this->m_nNumOfAOs);
    if (iteration == 1) {
      P = this->getPMatrix(RUN_RKS, iteration - 1);
    } else {
      P = this->getPMatrix(RUN_RKS, iteration - 1);
      {
        const TlDenseSymmetricMatrix_BLAS_Old prevP =
            this->getPMatrix(RUN_RKS, iteration - 2);
        P -= prevP;
      }
    }

    //     std::cerr << "DfThreeindexintegral using RIHF." << std::endl;
    //     DfEri2 dfEri2(this->m_flGbi, iteration);
    //     TlDenseSymmetricMatrix_BLAS_Old K = dfEri2.getKMatrix(P);
    //     K *= -0.5;
    //     std::cerr << "K(RI) >>>>" << std::endl;
    //     K.print(std::cerr);

    //     std::cerr << "DfThreeindexintegral using HF." << std::endl;
    // DfTwoElectronIntegral dfTEI(this->pPdfParam_);
    DfEriX dfEri(this->pPdfParam_);
    TlDenseSymmetricMatrix_BLAS_Old K(this->m_nNumOfAOs);
    // dfTEI.getContractKMatrixByIntegralDriven(P, &K);
    dfEri.getK(P, &K);

    // TlDenseSymmetricMatrix_BLAS_Old K = dfTEI.getContractKMatrixByRTmethod(P);
    //     std::cerr << "K(HF) >>>>" << std::endl;
    //     K.print(std::cerr);

    F += K;
    E += (0.5 * K);
  }

  DfObject::saveFpqMatrix(RUN_RKS, iteration, F);
  E.save("fl_Work/fl_Mtr_Epqtmp" + TlUtils::xtos(iteration));
}

void DfThreeindexintegrals::mainDIRECT_RKS2(int iteration) {
  this->log_.info("Direct scheme method is employed");

  // read Rou
  TlVector_BLAS currRho;
  currRho.load("fl_Work/fl_Vct_Rou" + TlUtils::xtos(iteration));
  assert(this->m_nNumOfAux == currRho.getSize());

  // read Myu
  //   TlVector_BLAS currMyu;
  //   currMyu.load("fl_Work/fl_Vct_Myu" + TlUtils::xtos(iteration));
  //   assert(this->numOfAuxXC_ == currMyu.getSize());

  // read Epsilon
  // energy for xc energy term (compair myu*Ppq*Pqa with 4/3*Ex1)
  //   TlVector_BLAS currEps;
  //   currEps.load("fl_Work/fl_Vct_Epsilon" + TlUtils::xtos(iteration));
  //   assert(this->numOfAuxXC_ == currEps.getSize());

  // read previous Rou, Myu, Eps & calculate delta
  if (iteration > 1) {  // iteration != 1
    {
      TlVector_BLAS prevRho(this->m_nNumOfAux);
      prevRho.load("fl_Work/fl_Vct_Rou" + TlUtils::xtos(iteration - 1));
      assert(this->m_nNumOfAux == prevRho.getSize());

      currRho -= prevRho;
    }

    //     {
    //       TlVector_BLAS prevMyu;
    //       prevMyu.load("fl_Work/fl_Vct_Myu" + TlUtils::xtos(iteration -1));
    //       assert(this->numOfAuxXC_ == prevMyu.getSize());

    //       currMyu -= prevMyu;
    //     }

    //     {
    //       TlVector_BLAS prevEps;
    //       prevEps.load("fl_Work/fl_Vct_Epsilon" + TlUtils::xtos(iteration
    //       -1)); assert(this->numOfAuxXC_ == prevEps.getSize());

    //       currEps -= prevEps;
    //     }
  }

  // prepare Fock
  TlDenseSymmetricMatrix_BLAS_Old F(this->m_nNumOfAOs);

  // クーロン項
  // [pq|alpha]
  {
    DfEriX dferi(this->pPdfParam_);
    // dferi.getdeltaHpqA(currRho, F);
    dferi.getJ(currRho, &F);
  }

  // prepare E
  TlDenseSymmetricMatrix_BLAS_Old E = F;

  // 交換相関項
  if (this->m_sXCFunctional != "hf") {
    //     TlDenseSymmetricMatrix_BLAS_Old P(this->m_nNumOfAOs);
    //     if (iteration == 1){
    //       P = this->getPMatrix("rks", iteration -1);
    //     } else {
    //       P = this->getPMatrix("rks", iteration -1);
    //       {
    //  const TlDenseSymmetricMatrix_BLAS_Old prevP = this->getPMatrix("rks",
    //  iteration
    //  -2);
    //  P -= prevP;
    //       }
    //     }

    // new method
    //     {
    //       TlDenseSymmetricMatrix_BLAS_Old Fxc; //(this->m_nNumOfAOs);
    //       DfXCFunctional dfXCFunctional(this->m_flGbi);

    //       const TlDenseSymmetricMatrix_BLAS_Old P = this->getPMatrix("rks",
    //       iteration
    //       -1);

    //       Fxc = dfXCFunctional.getFxc(P);
    //       Fxc.save("fl_Work/fl_Mtr_Fxc.matrix." + TlUtils::xtos(iteration));

    //       //       if (iteration != 1){
    //       //     TlDenseSymmetricMatrix_BLAS_Old Fxc_prev;
    //       //     Fxc_prev.load("fl_Work/fl_Mtr_Fxc.matrix." +
    //       TlUtils::xtos(iteration -1));
    //       //     Fxc -= Fxc_prev;
    //       //       }

    //       F += Fxc;
    //     }

    //     DfCalcGridX dfCalcGrid(this->m_flGbi);
    //     dfCalcGrid.calcXCIntegForFock(P, F);
    // CnErr.abort("not implemented yet!! stop.");
  } else {
    // Ex using HF
    TlDenseSymmetricMatrix_BLAS_Old P(this->m_nNumOfAOs);
    if (iteration == 1) {
      P = this->getPMatrix(RUN_RKS, iteration - 1);
    } else {
      P = this->getPMatrix(RUN_RKS, iteration - 1);
      {
        const TlDenseSymmetricMatrix_BLAS_Old prevP =
            this->getPMatrix(RUN_RKS, iteration - 2);
        P -= prevP;
      }
    }

    //     std::cerr << "delta P >>>>" << std::endl;
    //     P.print(std::cerr);

    //     std::cerr << "DfThreeindexintegral using RIHF." << std::endl;
    //     DfEri2 dfEri2(this->m_flGbi, iteration);
    //     TlDenseSymmetricMatrix_BLAS_Old K = dfEri2.getKMatrix(P);
    //     K *= -0.5;
    //     std::cerr << "K(RI) >>>>" << std::endl;
    //     K.print(std::cerr);

    //     std::cerr << "DfThreeindexintegral using HF." << std::endl;
    // DfTwoElectronIntegral dfTEI(this->pPdfParam_);
    DfEriX dfEri(this->pPdfParam_);
    TlDenseSymmetricMatrix_BLAS_Old K(this->m_nNumOfAOs);
    // dfTEI.getContractKMatrixByIntegralDriven(P, &K);
    dfEri.getK(P, &K);

    // TlDenseSymmetricMatrix_BLAS_Old K = dfTEI.getContractKMatrixByRTmethod(P);
    //     std::cerr << "K(HF) >>>>" << std::endl;
    //     K.print(std::cerr);

    F += K;
    E += (0.5 * K);
  }

  DfObject::saveFpqMatrix(RUN_RKS, iteration, F);
  E.save("fl_Work/fl_Mtr_Epqtmp" + TlUtils::xtos(iteration));
}

void DfThreeindexintegrals::mainDIRECT_UKS(int iteration) {
  this->log_.info("Direct scheme method is employed");

  // read Rho
  TlVector_BLAS currRho;
  {
    TlVector_BLAS currRhoA;
    TlVector_BLAS currRhoB;
    currRhoA.load("fl_Work/fl_Vct_Roua" + TlUtils::xtos(iteration));
    currRhoB.load("fl_Work/fl_Vct_Roub" + TlUtils::xtos(iteration));
    assert(currRhoA.getSize() == currRhoB.getSize());
    assert(this->m_nNumOfAux == currRhoA.getSize());

    currRho = currRhoA + currRhoB;
  }

  // read Myu
  TlVector_BLAS currMyuA;
  TlVector_BLAS currMyuB;
  currMyuA.load("fl_Work/fl_Vct_Myua" + TlUtils::xtos(iteration));
  currMyuB.load("fl_Work/fl_Vct_Myub" + TlUtils::xtos(iteration));
  assert(currMyuA.getSize() == currMyuB.getSize());
  assert(this->numOfAuxXC_ == currMyuA.getSize());

  // read Epsilon
  // energy for xc energy term (compair myu*Ppq*Pqa with 4/3*Ex1)
  TlVector_BLAS currEpsA;
  TlVector_BLAS currEpsB;
  if ("xalpha" == m_sXCFunctional || "gxalpha" == m_sXCFunctional) {
    currEpsA.load("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(iteration));
    currEpsB.load("fl_Work/fl_Vct_Epsilonb" + TlUtils::xtos(iteration));
    assert(currEpsA.getSize() == currEpsB.getSize());
  } else {
    currEpsA.load("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(iteration));
  }
  assert(this->numOfAuxXC_ == currEpsA.getSize());

  // read previous
  if (iteration > 1) {
    {
      TlVector_BLAS prevRhoA;
      TlVector_BLAS prevRhoB;
      prevRhoA.load("fl_Work/fl_Vct_Roua" + TlUtils::xtos(iteration - 1));
      prevRhoB.load("fl_Work/fl_Vct_Roub" + TlUtils::xtos(iteration - 1));
      assert(prevRhoA.getSize() == prevRhoB.getSize());
      assert(this->m_nNumOfAux == prevRhoA.getSize());

      currRho -= (prevRhoA + prevRhoB);
    }

    {
      TlVector_BLAS prevMyuA;
      TlVector_BLAS prevMyuB;
      prevMyuA.load("fl_Work/fl_Vct_Myua" + TlUtils::xtos(iteration - 1));
      prevMyuB.load("fl_Work/fl_Vct_Myub" + TlUtils::xtos(iteration - 1));
      assert(prevMyuA.getSize() == prevMyuB.getSize());
      assert(this->numOfAuxXC_ == prevMyuA.getSize());

      currMyuA -= prevMyuA;
      currMyuB -= prevMyuB;
    }

    {
      if ("xalpha" == m_sXCFunctional || "gxalpha" == m_sXCFunctional) {
        TlVector_BLAS prevEpsA;
        TlVector_BLAS prevEpsB;
        prevEpsA.load("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(iteration - 1));
        prevEpsB.load("fl_Work/fl_Vct_Epsilonb" + TlUtils::xtos(iteration - 1));
        assert(prevEpsA.getSize() == prevEpsB.getSize());
        assert(this->numOfAuxXC_ == prevEpsA.getSize());

        currEpsA -= prevEpsA;
        currEpsB -= prevEpsB;
      } else {
        TlVector_BLAS prevEpsA;
        prevEpsA.load("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(iteration - 1));
        assert(this->numOfAuxXC_ == prevEpsA.getSize());

        currEpsA -= prevEpsA;
      }
    }
  }

  TlDenseSymmetricMatrix_BLAS_Old F(this->m_nNumOfAOs);

  {
    DfEriX dferi(this->pPdfParam_);
    // dferi.getdeltaHpqA(currRho, F);
    dferi.getJ(currRho, &F);
  }

  F.save("fl_Work/fl_Mtr_Epqtmp" + TlUtils::xtos(iteration));

  TlDenseSymmetricMatrix_BLAS_Old E = F;

  {
    DfOverlapX dfovr(this->pPdfParam_);
    dfovr.get_pqg(currMyuA, currMyuB, &F, &E);
  }

  F.save("fl_Work/fl_Mtr_Fpq.matrix.uks-alpha" + TlUtils::xtos(iteration));
  E.save("fl_Work/fl_Mtr_Fpq.matrix.uks-beta" + TlUtils::xtos(iteration));

  E.load("fl_Work/fl_Mtr_Epqtmp" + TlUtils::xtos(iteration));

  {
    DfOverlapX dfovr(this->pPdfParam_);
    if (m_sXCFunctional == "xalpha" || m_sXCFunctional == "gxalpha") {
      dfovr.get_pqg(currEpsA, &E);
      dfovr.get_pqg(currEpsB, &E);
    } else {
      dfovr.get_pqg(currEpsA, &E);
    }
  }

  E.save("fl_Work/fl_Mtr_Epqtmp" + TlUtils::xtos(iteration));
}

void DfThreeindexintegrals::mainDIRECT_ROKS(int iteration) {
  this->log_.info("ROKS Direct scheme method is employed");

  // read Rou
  TlVector_BLAS currRho;
  {
    TlVector_BLAS currRhoA;
    currRhoA.load("fl_Work/fl_Vct_Roua" + TlUtils::xtos(iteration));
    TlVector_BLAS currRhoB;
    currRhoB.load("fl_Work/fl_Vct_Roub" + TlUtils::xtos(iteration));

    assert(currRhoA.getSize() == currRhoB.getSize());
    assert(this->m_nNumOfAux == currRhoA.getSize());

    currRho = currRhoA + currRhoB;
  }

  // read Myu
  TlVector_BLAS currMyuA;
  currMyuA.load("fl_Work/fl_Vct_Myua" + TlUtils::xtos(iteration));
  TlVector_BLAS currMyuB;
  currMyuB.load("fl_Work/fl_Vct_Myub" + TlUtils::xtos(iteration));

  assert(currMyuA.getSize() == currMyuB.getSize());
  assert(this->numOfAuxXC_ == currMyuA.getSize());

  // read Epsilon
  // energy for xc energy term (compair myu*Ppq*Pqa with 4/3*Ex1)
  TlVector_BLAS currEpsA;
  TlVector_BLAS currEpsB;

  if (m_sXCFunctional == "xalpha" || m_sXCFunctional == "gxalpha") {
    currEpsA.load("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(iteration));
    currEpsB.load("fl_Work/fl_Vct_Epsilonb" + TlUtils::xtos(iteration));

    assert(currEpsA.getSize() == currEpsB.getSize());
    assert(this->numOfAuxXC_ == currEpsA.getSize());
  } else {
    currEpsA.load("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(iteration));
    assert(this->numOfAuxXC_ == currEpsA.getSize());
  }

  // read previous
  if (iteration > 1) {
    {
      TlVector_BLAS prevRhoA;
      prevRhoA.load("fl_Work/fl_Vct_Roua" + TlUtils::xtos(iteration - 1));
      TlVector_BLAS prevRhoB;
      prevRhoB.load("fl_Work/fl_Vct_Roub" + TlUtils::xtos(iteration - 1));

      assert(prevRhoA.getSize() == prevRhoB.getSize());
      assert(this->m_nNumOfAux == prevRhoA.getSize());

      currRho -= (prevRhoA + prevRhoB);
    }

    {
      TlVector_BLAS prevMyuA;
      prevMyuA.load("fl_Work/fl_Vct_Myua" + TlUtils::xtos(iteration - 1));
      TlVector_BLAS prevMyuB;
      prevMyuB.load("fl_Work/fl_Vct_Myub" + TlUtils::xtos(iteration - 1));

      assert(prevMyuA.getSize() == prevMyuB.getSize());
      assert(this->numOfAuxXC_ == prevMyuA.getSize());

      currMyuA -= prevMyuA;
      currMyuB -= prevMyuB;
    }

    {
      if (m_sXCFunctional == "xalpha" || m_sXCFunctional == "gxalpha") {
        TlVector_BLAS prevEpsA;
        TlVector_BLAS prevEpsB;
        prevEpsA.load("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(iteration - 1));
        prevEpsB.load("fl_Work/fl_Vct_Epsilonb" + TlUtils::xtos(iteration - 1));

        assert(prevEpsA.getSize() == prevEpsB.getSize());
        assert(this->numOfAuxXC_ == prevEpsA.getSize());

        currEpsA -= prevEpsA;
        currEpsB -= prevEpsB;
      } else {
        TlVector_BLAS prevEpsA;
        prevEpsA.load("fl_Work/fl_Vct_Epsilona" + TlUtils::xtos(iteration - 1));
        assert(this->numOfAuxXC_ == prevEpsA.getSize());

        currEpsA -= prevEpsA;
      }
    }
  }

  // construct kohn-sham matrix of closed part ( F1 )
  this->log_.info("construct Kohn-Sham matrix F1");
  {
    TlDenseSymmetricMatrix_BLAS_Old F(this->m_nNumOfAOs);
    TlDenseSymmetricMatrix_BLAS_Old E(this->m_nNumOfAOs);

    {
      // integral object generated
      // then add coulomb contribution of rho to F ?
      {
        DfEriX dferi(this->pPdfParam_);
        dferi.getJ(currRho, &F);
      }

      E = F;

      DfOverlapX dfovr(this->pPdfParam_);
      if (m_sXCFunctional == "xalpha" || m_sXCFunctional == "gxalpha") {
        TlVector_BLAS myu = 0.5 * currMyuA;

        dfovr.get_pqg(myu, currEpsA, &F, &E);

        myu = 0.5 * currMyuB;
        dfovr.get_pqg(myu, currEpsB, &F, &E);
      } else {
        TlVector_BLAS myu = 0.5 * currMyuA;
        dfovr.get_pqg(myu, &F);

        myu = 0.5 * currMyuB;
        dfovr.get_pqg(myu, currEpsA, &F, &E);
      }
    }

    // temporarily write F matrix to a file to calculate E
    E.save("fl_Work/fl_Mtr_Epqtmp" + TlUtils::xtos(iteration));

    // temporarily write F matrix to a file
    F.save("fl_Work/fl_Mtr_F1pqtmp");
  }

  // construct kohn-sham matrix of open part ( F2 )
  this->log_.info("construct Kohn-Sham matrix F2");
  {
    TlDenseSymmetricMatrix_BLAS_Old F(this->m_nNumOfAOs);
    TlDenseSymmetricMatrix_BLAS_Old E(this->m_nNumOfAOs);

    // integral object generated
    // then add coulomb contribution of rho to F ?
    {
      TlVector_BLAS tmpRho = 0.5 * currRho;
      DfEriX dferi(this->pPdfParam_);
      dferi.getJ(tmpRho, &F);
    }

    // integral object generated
    // then add xc potential contribution of myu-alpha to F ?
    {
      TlVector_BLAS myu = 0.5 * currMyuA;
      DfOverlapX dfovr(this->pPdfParam_);
      dfovr.get_pqg(currMyuA, &F);
    }

    // integral object generated
    // then add xc potential contribution of myu-beta to F ?
    {
      // Log << "■不対電子はαだけにあること！■\n";
    }

    F.save("fl_Work/fl_Mtr_F2pqtmp");
  }
}
