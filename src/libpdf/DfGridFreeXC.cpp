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

#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#include "DfCD.h"
#include "DfFunctional_B3LYP.h"
#include "DfFunctional_B88LYP.h"
#include "DfFunctional_Becke88.h"
#include "DfFunctional_HFS.h"
#include "DfFunctional_SVWN.h"
#include "DfGridFreeXC.h"
#include "DfOverlapX.h"
#include "DfXCFunctional.h"
#include "DfXMatrix.h"
#include "Fl_Geometry.h"
#include "TlSystem.h"
#include "TlTime.h"
#include "TlUtils.h"

const int DfGridFreeXC::MAX_SHELL_TYPE = 2 + 1;
const double DfGridFreeXC::ONE_THIRD = 1.0 / 3.0;

DfGridFreeXC::DfGridFreeXC(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      pOvpEngines_(NULL),
      orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]) {
  this->numOfPQs_ = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

  this->tau_ = 1.0E-10;
  if ((*pPdfParam)["gridfree/CDAM_tau"].getStr().empty() != true) {
    this->tau_ = (*pPdfParam)["grid_free/CDAM_tau"].getDouble();
  }

  this->epsilon_ = 1.0E-4;
  if ((*pPdfParam)["gridfree/CD_epsilon"].getStr().empty() != true) {
    this->epsilon_ = (*pPdfParam)["grid_free/CD_epsilon"].getDouble();
  }

  this->isCanonicalOrthogonalize_ = true;
  if ((*pPdfParam)["gridfree/orthogonalize_method"].getStr().empty() != true) {
    const std::string method = TlUtils::toUpper(
        (*pPdfParam)["gridfree/orthogonalize_method"].getStr());
    if (method == "LOWDIN") {
      this->isCanonicalOrthogonalize_ = false;
    }
  }

  this->GfVEigvalVtrPath_ = "";
  if ((*pPdfParam)["gridfree/save_v_eigval"].getBoolean()) {
    this->GfVEigvalVtrPath_ = DfObject::getGfVEigvalVtrPath();
  }

  this->debugSaveM_ = (*pPdfParam)["debug/DfGridFreeXC/saveM"].getBoolean();
  if (this->debugSaveM_) {
    this->log_.info("using GAMESS formula");
  }
}

DfGridFreeXC::~DfGridFreeXC() {}

DfOverlapX* DfGridFreeXC::getDfOverlapObject() {
  DfOverlapX* pDfOverlapX = new DfOverlapX(this->pPdfParam_);
  return pDfOverlapX;
}

DfXMatrix* DfGridFreeXC::getDfXMatrixObject() {
  DfXMatrix* pDfXMatrix = new DfXMatrix(this->pPdfParam_);
  return pDfXMatrix;
}

// before SCF ==================================================================
void DfGridFreeXC::preprocessBeforeSCF() {
  this->preprocessBeforeSCF_templ<DfOverlapX, DfXMatrix,
                                  TlDenseSymmetricMatrix_BLAS_Old,
                                  TlDenseGeneralMatrix_BLAS_old>();
}

// in SCF ======================================================================
void DfGridFreeXC::buildFxc() {
  const DfXCFunctional xcFunc(this->pPdfParam_);
  if (xcFunc.getXcType() == DfXCFunctional::HF) {
    // need not pure-DFT term
    return;
  }

  const DfXCFunctional::FUNCTIONAL_TYPE funcType = xcFunc.getFunctionalType();
  switch (funcType) {
    case DfXCFunctional::LDA:
      this->buildFxc_LDA();
      break;

    case DfXCFunctional::GGA:
      this->buildFxc_GGA();
      break;

    default:
      this->log_.critical("unknown XC functional type. stop.");
      CnErr.abort();
      break;
  }
}

void DfGridFreeXC::buildFxc_LDA() {
  this->log_.info("DfGridFreeXC::buildFxc_LDA()");
  this->buildFxc_LDA_method<DfOverlapX, DfCD, TlDenseSymmetricMatrix_BLAS_Old,
                            TlDenseGeneralMatrix_BLAS_old>();
}

void DfGridFreeXC::createEngines() {
  assert(this->pOvpEngines_ == NULL);

  this->log_.info(
      TlUtils::format("create ERI engine: %d", this->numOfThreads_));
  this->pOvpEngines_ = new DfOverlapEngine[this->numOfThreads_];
}

void DfGridFreeXC::destroyEngines() {
  this->log_.info("delete OpenMP ERI engine");
  if (this->pOvpEngines_ != NULL) {
    delete[] this->pOvpEngines_;
  }
  this->pOvpEngines_ = NULL;
}

DfTaskCtrl* DfGridFreeXC::getDfTaskCtrlObject() const {
  DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
  return pDfTaskCtrl;
}

void DfGridFreeXC::finalize(TlDenseSymmetricMatrix_BLAS_Old* pMtx) {
  // do nothing
}

void DfGridFreeXC::get_F_lamda(const TlVector_BLAS lamda,
                               TlMatrixObject* pF_lamda,
                               TlMatrixObject* pE_lamda) {
  const int dim = lamda.getSize();
  assert(pF_lamda->getNumOfRows() == dim);
  assert(pF_lamda->getNumOfCols() == dim);
  assert(pE_lamda->getNumOfRows() == dim);
  assert(pE_lamda->getNumOfCols() == dim);

  DfFunctional_LDA* pFunc = NULL;
  std::string checkXC = this->m_sXCFunctional;
  if (checkXC == "SVWN") {
    pFunc = new DfFunctional_SVWN();
  } else if (checkXC == "HFS") {
    pFunc = new DfFunctional_HFS();
  } else {
    this->log_.critical(
        TlUtils::format("not support functional: %s", checkXC.c_str()));
    abort();
  }

  double fv_a = 0.0;
  double fv_b = 0.0;
  for (int i = 0; i < dim; ++i) {
    const double v = lamda.get(i);
    if (v > 1.0E-16) {
      pFunc->getDerivativeFunctional(v, v, &fv_a, &fv_b);
      pF_lamda->set(i, i, fv_a);

      const double f = pFunc->getFunctional(v, v) / (2.0 * v);
      pE_lamda->set(i, i, f);
    }
  }

  delete pFunc;
  pFunc = NULL;
}

void DfGridFreeXC::getM_exact(const TlDenseSymmetricMatrix_BLAS_Old& P,
                              TlDenseSymmetricMatrix_BLAS_Old* pM) {
  assert(pM != NULL);
  TlDenseGeneralMatrix_BLAS_old M(this->m_nNumOfAOs, this->m_nNumOfAOs);
  pM->resize(this->m_nNumOfAOs);

  DfOverlapEngine engine;

  const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                  (*(this->pPdfParam_))["basis_set"]);

  const ShellArrayTable shellArrayTable =
      this->makeShellArrayTable(orbitalInfo);
  // const ShellPairArrayTable shellPairArrayTable =
  // this->getShellPairArrayTable(shellArrayTable);

  for (int shellTypeP = MAX_SHELL_TYPE - 1; shellTypeP >= 0; --shellTypeP) {
    const int maxStepsP = 2 * shellTypeP + 1;
    const ShellArray shellArrayP = shellArrayTable[shellTypeP];
    ShellArray::const_iterator pItEnd = shellArrayP.end();

    for (int shellTypeQ = MAX_SHELL_TYPE - 1; shellTypeQ >= 0; --shellTypeQ) {
      const int maxStepsQ = 2 * shellTypeQ + 1;
      const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
      ShellArray::const_iterator qItEnd = shellArrayQ.end();

      for (int shellTypeR = MAX_SHELL_TYPE - 1; shellTypeR >= 0; --shellTypeR) {
        const int maxStepsR = 2 * shellTypeR + 1;
        const ShellArray shellArrayR = shellArrayTable[shellTypeR];
        ShellArray::const_iterator rItEnd = shellArrayR.end();

        for (int shellTypeS = MAX_SHELL_TYPE - 1; shellTypeS >= 0;
             --shellTypeS) {
          const int maxStepsS = 2 * shellTypeS + 1;
          const ShellArray shellArrayS = shellArrayTable[shellTypeS];
          ShellArray::const_iterator sItEnd = shellArrayS.end();

          const DfOverlapEngine::Query query(0, 0, 0, 0, shellTypeP, shellTypeQ,
                                             shellTypeR, shellTypeS);

          for (ShellArray::const_iterator pIt = shellArrayP.begin();
               pIt != pItEnd; ++pIt) {
            const index_type shellIndexP = *pIt;
            // const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            // const DfOverlapEngine::PGTOs pgtosP =
            // DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);

            for (ShellArray::const_iterator qIt = shellArrayQ.begin();
                 qIt != qItEnd; ++qIt) {
              const index_type shellIndexQ = *qIt;
              // const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
              // const DfOverlapEngine::PGTOs pgtosQ =
              // DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);

              for (ShellArray::const_iterator rIt = shellArrayR.begin();
                   rIt != rItEnd; ++rIt) {
                const index_type shellIndexR = *rIt;
                // const TlPosition posR = orbitalInfo.getPosition(shellIndexR);
                // const DfOverlapEngine::PGTOs pgtosR =
                // DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexR);

                for (ShellArray::const_iterator sIt = shellArrayS.begin();
                     sIt != sItEnd; ++sIt) {
                  const index_type shellIndexS = *sIt;
                  // const TlPosition posS =
                  // orbitalInfo.getPosition(shellIndexS); const
                  // DfOverlapEngine::PGTOs pgtosS =
                  // DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexS);

                  // engine.calc0(query,
                  //              posP, posQ, posR, posS,
                  //              pgtosP, pgtosQ, pgtosR, pgtosS);
                  engine.calc(0, orbitalInfo, shellIndexP, 0, orbitalInfo,
                              shellIndexQ, 0, orbitalInfo, shellIndexR, 0,
                              orbitalInfo, shellIndexS);

                  int index = 0;
                  for (int i = 0; i < maxStepsP; ++i) {
                    const int indexP = shellIndexP + i;

                    for (int j = 0; j < maxStepsQ; ++j) {
                      const int indexQ = shellIndexQ + j;

                      for (int k = 0; k < maxStepsR; ++k) {
                        const int indexR = shellIndexR + k;

                        for (int l = 0; l < maxStepsS; ++l) {
                          const int indexS = shellIndexS + l;

                          const double P_rs = P.get(indexR, indexS);
                          const double value = engine.WORK[index];
                          M.add(indexP, indexQ, P_rs * value);

                          ++index;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  *pM = M;
}

DfGridFreeXC::ShellArrayTable DfGridFreeXC::makeShellArrayTable(
    const TlOrbitalInfoObject& orbitalInfo) {
  ShellArrayTable shellArrayTable(MAX_SHELL_TYPE);
  const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();

  index_type shellIndex = 0;
  while (shellIndex < maxShellIndex) {
    // shellType: 0=s, 1=p, 2=d
    const int shellType = orbitalInfo.getShellType(shellIndex);
    const int steps = 2 * shellType + 1;

    shellArrayTable[shellType].push_back(shellIndex);

    shellIndex += steps;
  }

  return shellArrayTable;
}

DfGridFreeXC::ShellPairArrayTable DfGridFreeXC::getShellPairArrayTable(
    const ShellArrayTable& shellArrayTable) {
  ShellPairArrayTable shellPairArrayTable(MAX_SHELL_TYPE * MAX_SHELL_TYPE);

  for (int shellTypeP = MAX_SHELL_TYPE - 1; shellTypeP >= 0; --shellTypeP) {
    const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
    ShellArray::const_iterator pItEnd = shellArrayP.end();

    for (int shellTypeR = MAX_SHELL_TYPE - 1; shellTypeR >= 0; --shellTypeR) {
      const ShellArray& shellArrayR = shellArrayTable[shellTypeR];
      ShellArray::const_iterator rItEnd = shellArrayR.end();

      const int shellPairType_PR = shellTypeP * MAX_SHELL_TYPE + shellTypeR;
      for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd;
           ++pIt) {
        const index_type indexP = *pIt;

        for (ShellArray::const_iterator rIt = shellArrayR.begin();
             rIt != rItEnd; ++rIt) {
          const index_type indexR = *rIt;

          if (indexP >= indexR) {
            ShellPair shellPair(indexP, indexR);
            shellPairArrayTable[shellPairType_PR].push_back(shellPair);
          }
        }
      }
    }
  }

  return shellPairArrayTable;
}

// TlDenseSymmetricMatrix_BLAS_Old DfGridFreeXC::getPMatrix()
// {
//     TlDenseSymmetricMatrix_BLAS_Old P =
//     this->getPpqMatrix<TlDenseSymmetricMatrix_BLAS_Old>(RUN_RKS,
//     this->m_nIteration -1); return P;
// }

TlDenseGeneralMatrix_BLAS_old DfGridFreeXC::getL() {
  // TlDenseGeneralMatrix_BLAS_old L =
  // DfObject::getLMatrix<TlDenseGeneralMatrix_BLAS_old>();
  TlDenseGeneralMatrix_BLAS_old L;
  L.load("GF_L.mat");

  return L;
}

DfGridFreeXC::PQ_PairArray DfGridFreeXC::getI2PQ() {
  std::string filepath = this->getI2pqVtrPath();
  std::ifstream ifs;
  ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
  if (ifs.fail()) {
    abort();
  }

  std::size_t size = 0;
  ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

  PQ_PairArray answer(size);
  index_type shellIndex1 = 0;
  index_type shellIndex2 = 0;
  for (std::size_t i = 0; i < size; ++i) {
    ifs.read(reinterpret_cast<char*>(&shellIndex1), sizeof(index_type));
    ifs.read(reinterpret_cast<char*>(&shellIndex2), sizeof(index_type));
    answer[i] = IndexPair2(shellIndex1, shellIndex2);
  }

  ifs.close();
  return answer;
}

void DfGridFreeXC::divideCholeskyBasis(const index_type numOfCBs,
                                       index_type* pStart, index_type* pEnd) {
  *pStart = 0;
  *pEnd = numOfCBs;
}

TlDenseSymmetricMatrix_BLAS_Old DfGridFreeXC::getCholeskyVector(
    const TlVector_BLAS& L_col, const PQ_PairArray& I2PQ) {
  const index_type numOfItilde = L_col.getSize();
  TlDenseSymmetricMatrix_BLAS_Old answer(this->m_nNumOfAOs);
  for (index_type i = 0; i < numOfItilde; ++i) {
    answer.set(I2PQ[i].index1(), I2PQ[i].index2(), L_col[i]);
  }

  return answer;
}

// -----------------------------------------------------------------------------
void DfGridFreeXC::buildFxc_GGA() {
  this->log_.info("DfGridFreeXC::buildFxc_GGA()");
  this->buildFxc_GGA_method<DfOverlapX, DfCD, TlDenseSymmetricMatrix_BLAS_Old,
                            TlDenseGeneralMatrix_BLAS_old>();
}

// void DfGridFreeXC::buildFxc_GGA()
// {
//     this->log_.info("build Fxc by grid-free method: functional type is
//     GGA.");

//     std::string basisset_param = "basis_set";
//     if (this->isDedicatedBasisForGridFree_) {
//         basisset_param = "basis_set_gridfree";
//     }
//     const TlOrbitalInfo orbitalInfo_GF((*(this->pPdfParam_))["coordinates"],
//                                        (*(this->pPdfParam_))[basisset_param]);
//     const index_type numOfAOs = this->m_nNumOfAOs;
//     const index_type numOfGfOrbs = orbitalInfo_GF.getNumOfOrbitals();
//     this->log_.info(TlUtils::format("AOs = %d", numOfAOs));
//     this->log_.info(TlUtils::format("auxAOs for GF = %d", numOfGfOrbs));

//     // RKS
//     const TlDenseSymmetricMatrix_BLAS_Old PA = 0.5 *
//     DfObject::getPpqMatrix<TlDenseSymmetricMatrix_BLAS_Old>(RUN_RKS,
//                                                                                  this->m_nIteration -1);
//     assert(PA.getNumOfRows() == numOfAOs);

//     TlDenseSymmetricMatrix_BLAS_Old M;
//     if (this->XC_engine_ == XC_ENGINE_GRIDFREE_CD) {
//         this->log_.info("begin to create M matrix based on CD.");
//         {
//             DfCD dfCD(this->pPdfParam_);
//             dfCD.getM(PA, &M);
//             // M *= 2.0;
//         }
//     } else {
//         this->log_.info("begin to create M matrix using 4-center overlap.");
//         DfOverlapX dfOvp(this->pPdfParam_);
//         if (this->isDedicatedBasisForGridFree_) {
//             dfOvp.getM_A(PA, &M);
//         } else {
//             dfOvp.getM(PA, &M);
//         }
//     }
//     if (this->debugSaveM_) {
//         M.save(TlUtils::format("fl_Work/debug_M.%d.mat",
//         this->m_nIteration));
//     }

//     this->log_.info("begin to generate Fxc using grid-free method.");
//     // M~(=V^t * M * V) および SVU(=SVU, but now only SV)の作成
//     TlDenseGeneralMatrix_BLAS_old S;
//     TlDenseGeneralMatrix_BLAS_old V;
//     if (this->isDedicatedBasisForGridFree_) {
//         S = DfObject::getGfStildeMatrix<TlDenseGeneralMatrix_BLAS_old>();
//         V = DfObject::getGfVMatrix<TlDenseGeneralMatrix_BLAS_old>();
//     } else {
//         S = DfObject::getSpqMatrix<TlDenseSymmetricMatrix_BLAS_Old>();
//         V = DfObject::getXMatrix<TlDenseGeneralMatrix_BLAS_old>();
//     }

//     const index_type numOfGFOrthNormBasis = V.getNumOfCols();
//     this->log_.info(TlUtils::format("orthonormal basis = %d",
//     numOfGFOrthNormBasis)); TlDenseGeneralMatrix_BLAS_old Vt = V;
//     Vt.transposeInPlace();

//     TlDenseGeneralMatrix_BLAS_old St = S;
//     St.transposeInPlace();

//     TlDenseSymmetricMatrix_BLAS_Old Mtilde = Vt * M * V;
//     //Mtilde.save("Mtilde.mat");

//     TlVector_BLAS lambda;
//     TlDenseGeneralMatrix_BLAS_old U;
//     Mtilde.diagonal(&lambda, &U);
//     //lambda.save("lambda.vct");
//     //U.save("U.mat");
//     TlDenseGeneralMatrix_BLAS_old Ut = U;
//     Ut.transposeInPlace();
//     assert(lambda.getSize() == numOfGFOrthNormBasis);
//     assert(U.getNumOfRows() == numOfGFOrthNormBasis);
//     assert(U.getNumOfCols() == numOfGFOrthNormBasis);

//     // M[rho^(-1/3)]
//     TlDenseSymmetricMatrix_BLAS_Old Mtilde_13;
//     {
//         TlDenseSymmetricMatrix_BLAS_Old lambda_13(numOfGFOrthNormBasis);
//         for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
//             double v_13 = 0.0;
//             const double v = lambda.get(i);
//             if (v > 1.0E-16) {
//                 v_13 = std::pow(v, - ONE_THIRD);
//             }
//             lambda_13.set(i, i, v_13);
//         }
//         Mtilde_13 = U * lambda_13 * Ut;
//     }
//     //Mtilde_13.save("Mtilde_13.mat");

//     // GGA用gradient
//     TlDenseGeneralMatrix_BLAS_old Gx =
//     this->getDipoleVelocityIntegralsXMatrix<TlDenseGeneralMatrix_BLAS_old>();
//     TlDenseGeneralMatrix_BLAS_old
//     Gy =
//     this->getDipoleVelocityIntegralsYMatrix<TlDenseGeneralMatrix_BLAS_old>();
//     TlDenseGeneralMatrix_BLAS_old Gz =
//     this->getDipoleVelocityIntegralsZMatrix<TlDenseGeneralMatrix_BLAS_old>(); Gx
//     *=
//     -1.0; Gy
//     *= -1.0; Gz *= -1.0; const TlDenseGeneralMatrix_BLAS_old DX = Vt * Gx * V;
//     const
//     TlDenseGeneralMatrix_BLAS_old DY = Vt * Gy * V; const
//     TlDenseGeneralMatrix_BLAS_old DZ = Vt *
//     Gz * V;
//     // DX.save("DX.mat");
//     // DY.save("DY.mat");
//     // DZ.save("DZ.mat");

//     TlDenseGeneralMatrix_BLAS_old DXt = DX;
//     DXt.transposeInPlace();
//     TlDenseGeneralMatrix_BLAS_old DYt = DY;
//     DYt.transposeInPlace();
//     TlDenseGeneralMatrix_BLAS_old DZt = DZ;
//     DZt.transposeInPlace();
//     const TlDenseGeneralMatrix_BLAS_old RTX = 3.0*(DXt * Mtilde_13 + Mtilde_13 *
//     DX);
//     const TlDenseGeneralMatrix_BLAS_old RTY = 3.0*(DYt * Mtilde_13 + Mtilde_13 *
//     DY);
//     const TlDenseGeneralMatrix_BLAS_old RTZ = 3.0*(DZt * Mtilde_13 + Mtilde_13 *
//     DZ);
//     TlDenseGeneralMatrix_BLAS_old RTXt = RTX;
//     RTXt.transposeInPlace();
//     TlDenseGeneralMatrix_BLAS_old RTYt = RTY;
//     RTYt.transposeInPlace();
//     TlDenseGeneralMatrix_BLAS_old RTZt = RTZ;
//     RTZt.transposeInPlace();

//     // RX2 := M[{nabla rho / rho^(-4/3)}^2]
//     const TlDenseSymmetricMatrix_BLAS_Old RX2 = RTXt*RTX + RTYt*RTY + RTZt*RTZ;
//     // RTX.save("RTX.mat");
//     // RTY.save("RTY.mat");
//     // RTZ.save("RTZ.mat");
//     // RX2.save("RX2.mat");

//     // RZ2 := [nabla rho / rho^(4/3)] * (DX + DY + DZ)
//     const TlDenseGeneralMatrix_BLAS_old RZ2 = RTX*DX + RTY*DY + RTZ*DZ;
//     //RZ2.save("RZ2.mat");
//     TlDenseGeneralMatrix_BLAS_old RZ2t = RZ2;
//     RZ2t.transposeInPlace();

//     TlVector_BLAS x2;
//     TlDenseGeneralMatrix_BLAS_old Ux2;
//     RX2.diagonal(&x2, &Ux2);
//     //x2.save("x2.vct");
//     //Ux2.save("Ux2.mat");
//     TlDenseGeneralMatrix_BLAS_old Ux2t = Ux2;
//     Ux2t.transposeInPlace();

//     // ------------------
//     assert(lambda.getSize() == numOfGFOrthNormBasis);
//     // TlVector_BLAS rhoAs(numOfGfOrbs);
//     // TlVector_BLAS xAs(numOfGfOrbs);
//     TlVector_BLAS rhoAs(numOfGFOrthNormBasis);
//     TlVector_BLAS xAs(numOfGFOrthNormBasis);
//     //for (index_type i = 0; i < numOfGfOrbs; ++i) {
//     for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
//         const double rho_value = lambda[i];
//         const double rho = (rho_value > 1.0E-16) ? rho_value : 0.0;
//         rhoAs[i] = rho;

//         const double x2_value = x2.get(i);
//         const double x = (x2_value > 1.0E-16) ? std::sqrt(x2_value) : 0.0;
//         xAs[i] = x;
//     }
//     const TlVector_BLAS rhoBs = rhoAs;
//     const TlVector_BLAS xBs = xAs;

//     DfFunctional_GGA* pFunc = this->getFunctionalGGA();
//     // Fxc -------------------------------------------------------------
//     TlDenseSymmetricMatrix_BLAS_Old FxcA(numOfAOs); // alpha spin
//     TlDenseSymmetricMatrix_BLAS_Old FxcB(numOfAOs); // beta spin
//     {
//         DerivativeFunctionalSets dfs =
//         pFunc->getDerivativeFunctional_GF(rhoAs, rhoBs, xAs, xBs);

//         TlVector_BLAS rhoAA43(numOfGFOrthNormBasis);
//         TlVector_BLAS rhoAB43(numOfGFOrthNormBasis);
//         TlVector_BLAS rhoBB43(numOfGFOrthNormBasis);
//         for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
//             const double rhoA = lambda[i];
//             if (rhoA > 1.0E-16) {
//                 const double rhoA43 = std::pow(rhoA, 4.0/3.0);
//                 const double rhoB43 = rhoA43; // RKS
//                 rhoAA43[i] = rhoA43;
//                 rhoBB43[i] = rhoB43;
//                 rhoAB43[i] = std::sqrt(rhoA43 * rhoB43);
//             }
//         }

//         TlDenseGeneralMatrix_BLAS_old FxcA_tilde(numOfGFOrthNormBasis,
//         numOfGFOrthNormBasis);
//         TlDenseGeneralMatrix_BLAS_old FxcB_tilde(numOfGFOrthNormBasis,
//         numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_RAR(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_RAX(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_RBR(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_RBX(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_GAAR(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_GAAX(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_GABR(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_GABX(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_GBBR(numOfGFOrthNormBasis);
//         TlDenseSymmetricMatrix_BLAS_Old diag_GBBX(numOfGFOrthNormBasis);
//         const int numOfTerms = pFunc->getNumOfDerivativeFunctionalTerms();
//         for (int term = 0; term < numOfTerms; ++term) {
//             for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
//                 diag_RAR(i, i) = dfs.rFrRhoA_R(term, i);
//                 diag_RAX(i, i) = dfs.rFrRhoA_X(term, i);
//                 diag_RBR(i, i) = dfs.rFrRhoB_R(term, i);
//                 diag_RBX(i, i) = dfs.rFrRhoB_X(term, i);

//                 diag_GAAR(i, i) = dfs.rFrGAA_R(term, i) * rhoAA43[i];
//                 diag_GAAX(i, i) = dfs.rFrGAA_X(term, i);
//                 diag_GABR(i, i) = dfs.rFrGAB_R(term, i) * rhoAB43[i];
//                 diag_GABX(i, i) = dfs.rFrGAB_X(term, i);
//                 diag_GBBR(i, i) = dfs.rFrGBB_R(term, i) * rhoBB43[i];
//                 diag_GBBX(i, i) = dfs.rFrGBB_X(term, i);
//             }

//             // alpha spin ------------
//             {
//                 const TlDenseSymmetricMatrix_BLAS_Old Fxc_RR = U * diag_RAR * Ut;
//                 const TlDenseSymmetricMatrix_BLAS_Old Fxc_RX = Ux2 * diag_RAX *
//                 Ux2t;
//                 TlDenseGeneralMatrix_BLAS_old Fxc_tilde_term1 = 0.5 * (Fxc_RR *
//                 Fxc_RX
//                 +
//                 Fxc_RX * Fxc_RR); FxcA_tilde += Fxc_tilde_term1;
//             }

//             TlDenseGeneralMatrix_BLAS_old Fxc_GAA;
//             {
//                 const TlDenseGeneralMatrix_BLAS_old Fxc_GAAR = U * diag_GAAR *
//                 Ut;
//                 const TlDenseGeneralMatrix_BLAS_old Fxc_GAAX = Ux2 * diag_GAAX *
//                 Ux2t;
//                 Fxc_GAA = 2.0 * 0.5 * (Fxc_GAAR * Fxc_GAAX + Fxc_GAAX *
//                 Fxc_GAAR);
//             }

//             TlDenseGeneralMatrix_BLAS_old Fxc_GAB;
//             {
//                 const TlDenseSymmetricMatrix_BLAS_Old Fxc_GABR = U * diag_GABR *
//                 Ut;
//                 const TlDenseSymmetricMatrix_BLAS_Old Fxc_GABX = Ux2 * diag_GABX
//                 *
//                 Ux2t;
//                 Fxc_GAB = 0.5 * (Fxc_GABR * Fxc_GABX + Fxc_GABX * Fxc_GABR);
//             }

//             TlDenseGeneralMatrix_BLAS_old Fxc_GBB;
//             {
//                 const TlDenseSymmetricMatrix_BLAS_Old Fxc_GBBR = U * diag_GBBR *
//                 Ut;
//                 const TlDenseSymmetricMatrix_BLAS_Old Fxc_GBBX = Ux2 * diag_GBBX
//                 *
//                 Ux2t;
//                 Fxc_GBB = 2.0 * 0.5 * (Fxc_GBBR * Fxc_GBBX + Fxc_GBBX *
//                 Fxc_GBBR);
//             }

//             {
//                 TlDenseGeneralMatrix_BLAS_old FxcA_term2 = Fxc_GAA + Fxc_GAB;
//                 TlDenseGeneralMatrix_BLAS_old FxcA_term2t = FxcA_term2;
//                 FxcA_term2t.transposeInPlace();
//                 TlDenseGeneralMatrix_BLAS_old FxcA_tilde_term2 = FxcA_term2t *
//                 RZ2 +
//                 RZ2t *
//                 FxcA_term2; FxcA_tilde += FxcA_tilde_term2;
//             }
//             {
//                 TlDenseGeneralMatrix_BLAS_old FxcB_term2 = Fxc_GBB + Fxc_GAB;
//                 TlDenseGeneralMatrix_BLAS_old FxcB_term2t = FxcB_term2;
//                 FxcB_term2t.transposeInPlace();
//                 TlDenseGeneralMatrix_BLAS_old FxcB_tilde_term2 = FxcB_term2t *
//                 RZ2 +
//                 RZ2t *
//                 FxcB_term2; FxcB_tilde += FxcB_tilde_term2;
//             }
//         }
//         FxcA = S * V * FxcA_tilde * Vt * St;
//         FxcB = S * V * FxcB_tilde * Vt * St;
//     }
//     DfObject::saveFxcMatrix(RUN_RKS, this->m_nIteration,
//     TlDenseSymmetricMatrix_BLAS_Old(FxcA));

//     // Exc -------------------------------------------------------------
//     TlDenseSymmetricMatrix_BLAS_Old ExcA(numOfAOs);
//     //TlDenseSymmetricMatrix_BLAS_Old ExcB(numOfAOs);
//     {
//         const FunctionalSets fs = pFunc->getFunctional_GF(rhoAs, rhoBs, xAs,
//         xBs);

//         TlDenseGeneralMatrix_BLAS_old ExcA_tilde(numOfGFOrthNormBasis,
//         numOfGFOrthNormBasis);
//         //TlDenseGeneralMatrix_BLAS_old ExcB_tilde(numOfGFOrthNormBasis,
//         numOfGFOrthNormBasis);
//         {
//             TlDenseSymmetricMatrix_BLAS_Old diag_AR(numOfGFOrthNormBasis);
//             TlDenseSymmetricMatrix_BLAS_Old diag_AX(numOfGFOrthNormBasis);
//             //TlDenseSymmetricMatrix_BLAS_Old diag_BR(numOfGFOrthNormBasis);
//             //TlDenseSymmetricMatrix_BLAS_Old diag_BX(numOfGFOrthNormBasis);
//             const int numOfTerms = pFunc->getNumOfFunctionalTerms();
//             for (int term = 0; term < numOfTerms; ++term) {
//                 for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
//                     const double rho = rhoAs[i] + rhoBs[i];
//                     if (rho > 1.0E-16) {
//                         const double inv_rho = 1.0 / rho;
//                         diag_AR(i, i) = fs.FA_termR(term, i) * inv_rho;
//                         diag_AX(i, i) = fs.FA_termX(term, i);
//                         // diag_BR(i, i) = fs.FB_termR(term, i) * inv_rho;
//                         // diag_BX(i, i) = fs.FB_termX(term, i);
//                     }
//                 }

//                 const TlDenseSymmetricMatrix_BLAS_Old ExcA_R = U * diag_AR * Ut;
//                 const TlDenseSymmetricMatrix_BLAS_Old ExcA_X = Ux2 * diag_AX *
//                 Ux2t;
//                 const TlDenseGeneralMatrix_BLAS_old ExcA_tilde_term = ExcA_R *
//                 ExcA_X +
//                 ExcA_X * ExcA_R;
//                 // const TlDenseSymmetricMatrix_BLAS_Old ExcB_R = U * diag_BR *
//                 Ut;
//                 // const TlDenseSymmetricMatrix_BLAS_Old ExcB_X = Ux2 * diag_BX *
//                 Ux2t;
//                 // const TlDenseGeneralMatrix_BLAS_old ExcB_tilde_term = ExcB_R *
//                 ExcB_X +
//                 ExcB_X
//                 * ExcB_R;

//                 ExcA_tilde += ExcA_tilde_term;
//                 // ExcB_tilde += ExcB_tilde_term;
//             }
//         }

//         ExcA = S * V * ExcA_tilde * Vt * St;
//         // ExcB = S * V * ExcB_tilde * Vt * St;
//         // Exc *= 2.0; // means RKS
//     }
//     DfObject::saveExcMatrix(RUN_RKS, this->m_nIteration,
//     TlDenseSymmetricMatrix_BLAS_Old(ExcA));
//     // DfObject::saveExcMatrix(RUN_UKS_ALPHA, this->m_nIteration,
//     TlDenseSymmetricMatrix_BLAS_Old(ExcA));
//     // DfObject::saveExcMatrix(RUN_UKS_BETA,  this->m_nIteration,
//     TlDenseSymmetricMatrix_BLAS_Old(ExcB));

//     delete pFunc;
//     pFunc = NULL;
// }

DfFunctional_GGA* DfGridFreeXC::getFunctionalGGA() {
  DfFunctional_GGA* pFunc = NULL;

  DfXCFunctional xcFunc(this->pPdfParam_);
  const DfXCFunctional::XC_TYPE xcType = xcFunc.getXcType();
  switch (xcType) {
    case DfXCFunctional::HFB:
      pFunc = new DfFunctional_Becke88();
      break;

    case DfXCFunctional::BLYP:
      pFunc = new DfFunctional_B88LYP();
      break;

    case DfXCFunctional::B3LYP:
      pFunc = new DfFunctional_B3LYP();
      break;

    default:
      std::cerr << "DfGridFreeXC::getFunctionalGGA() " << xcType << std::endl;
      pFunc = NULL;
      abort();
      break;
  }

  return pFunc;
}

TlDenseGeneralMatrix_BLAS_old DfGridFreeXC::getForce() {
  const RUN_TYPE runType = RUN_RKS;
  const int itr = this->m_nIteration;

  const TlDenseSymmetricMatrix_BLAS_Old S =
      this->getSpqMatrix<TlDenseSymmetricMatrix_BLAS_Old>();
  const TlDenseGeneralMatrix_BLAS_old X =
      this->getXMatrix<TlDenseGeneralMatrix_BLAS_old>();
  TlDenseGeneralMatrix_BLAS_old Xt = X;
  Xt.transposeInPlace();

  // TlDenseSymmetricMatrix_BLAS_Old P = 0.5 *
  // this->getPpqMatrix<TlDenseSymmetricMatrix_BLAS_Old>(RUN_RKS,
  // itr);
  const TlDenseGeneralMatrix_BLAS_old C =
      this->getCMatrix<TlDenseGeneralMatrix_BLAS_old>(runType, itr);
  TlDenseGeneralMatrix_BLAS_old Ct = C;
  Ct.transposeInPlace();

  // TlDenseGeneralMatrix_BLAS_old CCt = C * Ct;
  // CCt.save("CCt.mat");

  // Exc =====================================================================
  TlDenseSymmetricMatrix_BLAS_Old Exc =
      this->getFxcMatrix<TlDenseSymmetricMatrix_BLAS_Old>(RUN_RKS, itr);
  Exc.save("GF_Exc.mat");

  // dS ======================================================================
  TlDenseGeneralMatrix_BLAS_old dx, dy, dz;
  DfOverlapX dfOvp(this->pPdfParam_);
  dfOvp.getGradient(orbitalInfo_, &dx, &dy, &dz);
  // dx.save("GF_dx.mat");
  // dy.save("GF_dy.mat");
  // dz.save("GF_dz.mat");

  // 規格直交化 ==============================================================
  // TlDenseGeneralMatrix_BLAS_old o_Exc = Ct * Exc;
  TlDenseGeneralMatrix_BLAS_old o_Exc = Xt * Exc * X;
  o_Exc.save("GF_o_Exc.mat");

  // TlDenseGeneralMatrix_BLAS_old o_dx = Ct * dx;
  // TlDenseGeneralMatrix_BLAS_old o_dy = Ct * dy;
  // TlDenseGeneralMatrix_BLAS_old o_dz = Ct * dz;
  TlDenseGeneralMatrix_BLAS_old o_dx = Xt * dx * X;
  TlDenseGeneralMatrix_BLAS_old o_dy = Xt * dy * X;
  TlDenseGeneralMatrix_BLAS_old o_dz = Xt * dz * X;
  // TlDenseGeneralMatrix_BLAS_old o_dx = dx * C;
  // TlDenseGeneralMatrix_BLAS_old o_dy = dy * C;
  // TlDenseGeneralMatrix_BLAS_old o_dz = dz * C;
  // o_dx.save("GF_o_dx.mat");
  // o_dy.save("GF_o_dy.mat");
  // o_dz.save("GF_o_dz.mat");

  // 積 ======================================================================
  TlDenseGeneralMatrix_BLAS_old o_dxt = o_dx;
  o_dxt.transposeInPlace();
  TlDenseGeneralMatrix_BLAS_old o_dyt = o_dy;
  o_dyt.transposeInPlace();
  TlDenseGeneralMatrix_BLAS_old o_dzt = o_dz;
  o_dzt.transposeInPlace();

  TlDenseGeneralMatrix_BLAS_old o_dx_Exc = o_dxt * o_Exc;
  TlDenseGeneralMatrix_BLAS_old o_dy_Exc = o_dyt * o_Exc;
  TlDenseGeneralMatrix_BLAS_old o_dz_Exc = o_dzt * o_Exc;
  // TlDenseGeneralMatrix_BLAS_old o_dx_Exc = o_dx * o_Exc;
  // TlDenseGeneralMatrix_BLAS_old o_dy_Exc = o_dy * o_Exc;
  // TlDenseGeneralMatrix_BLAS_old o_dz_Exc = o_dz * o_Exc;
  // o_dx_Exc.save("GF_o_dxt_Exc.mat");
  // o_dy_Exc.save("GF_o_dyt_Exc.mat");
  // o_dz_Exc.save("GF_o_dzt_Exc.mat");

  // 規格直交化から戻す
  TlDenseGeneralMatrix_BLAS_old SX = S * X;
  TlDenseGeneralMatrix_BLAS_old SXt = SX;
  SXt.transposeInPlace();
  o_dx_Exc = SX * o_dx_Exc * SXt;
  o_dy_Exc = SX * o_dy_Exc * SXt;
  o_dz_Exc = SX * o_dz_Exc * SXt;

  const index_type numOfAOs = this->m_nNumOfAOs;
  const index_type numOfMOs = this->m_nNumOfMOs;

  // 右からCをかける
  // o_dx_Exc *= C;
  // o_dy_Exc *= C;
  // o_dz_Exc *= C;
  TlDenseGeneralMatrix_BLAS_old C_mo = C;
  {
    TlDenseSymmetricMatrix_BLAS_Old E(numOfAOs);
    TlVector_BLAS currOcc;
    currOcc.load(this->getOccupationPath(runType));
    for (index_type i = 0; i < numOfMOs; ++i) {
      if (std::fabs(currOcc[i] - 2.0) < 1.0E-5) {
        E.set(i, i, 1.0);
      }
    }
    C_mo *= E;
  }
  o_dx_Exc *= C_mo;
  o_dy_Exc *= C_mo;
  o_dz_Exc *= C_mo;

  // 左からCをかける
  const index_type numOfAtoms = this->m_nNumOfAtoms;
  TlDenseGeneralMatrix_BLAS_old force(numOfAtoms, 3);
  TlVector_BLAS Hx(numOfAOs), Hy(numOfAOs), Hz(numOfAOs);
  {
    TlVector_BLAS currOcc;
    currOcc.load(this->getOccupationPath(runType));
    for (int i = 0; i < numOfMOs; ++i) {
      if (std::fabs(currOcc[i] - 2.0) < 1.0E-5) {
        for (int j = 0; j < numOfAOs; ++j) {
          const double vx = C.get(j, i) * o_dx_Exc.get(j, i);
          Hx.add(j, vx);
          const double vy = C.get(j, i) * o_dy_Exc.get(j, i);
          Hy.add(j, vy);
          const double vz = C.get(j, i) * o_dz_Exc.get(j, i);
          Hz.add(j, vz);
        }
      }
    }
    Hx.save("GF_Hx.vct");
    Hy.save("GF_Hy.vct");
    Hz.save("GF_Hz.vct");
  }
  for (int i = 0; i < numOfAOs; ++i) {
    const index_type atomIndex = this->orbitalInfo_.getAtomIndex(i);
    force.add(atomIndex, 0, Hx.get(i));
    force.add(atomIndex, 1, Hy.get(i));
    force.add(atomIndex, 2, Hz.get(i));
  }

  force *= -4.0;
  force *= 2.0;  // rks

  force.save("GF_force.mat");

  // calc center of atoms
  // const Fl_Geometry flGeom((*(this->pPdfParam_))["coordinates"]);
  // TlPosition wc;
  // double sum_w = 0.0;
  // for (int i = 0; i < numOfAtoms; ++i) {
  //     const TlAtom atom = flGeom.getAtom(i);
  //     const TlPosition p = atom.getPosition();
  //     const double weight = atom.getStdWeight();
  //     std::cerr << TlUtils::format("(% f, %f, %f), %f", p.x(), p.y(), p.z(),
  //     weight) << std::endl; wc += weight * p; sum_w += weight;
  // }
  // wc *= -1.0 / sum_w;
  // std::cerr << TlUtils::format("wc: % f, % f, % f", wc.x(), wc.y(), wc.z())
  // << std::endl;

  // for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
  //     force.add(atomIndex, 0, wc.x());
  //     force.add(atomIndex, 1, wc.y());
  //     force.add(atomIndex, 2, wc.z());
  // }
  // force.save("GF_force2.mat");

  return force;
}

TlDenseGeneralMatrix_BLAS_old DfGridFreeXC::selectGradMat(
    const TlDenseGeneralMatrix_BLAS_old& input, const int atomIndex) {
  const index_type numOfAOs = this->m_nNumOfAOs;
  assert(input.getNumOfRows() == numOfAOs);
  assert(input.getNumOfCols() == numOfAOs);
  TlDenseGeneralMatrix_BLAS_old output(numOfAOs, numOfAOs);
  for (index_type p = 0; p < numOfAOs; ++p) {
    if (this->orbitalInfo_.getAtomIndex(p) == atomIndex) {
      for (index_type q = 0; q < numOfAOs; ++q) {
        output.set(p, q, input.get(p, q));
        // if (this->orbitalInfo_.getAtomIndex(q) == atomIndex) {
        //     output.set(p, q, input.get(p, q));
        // }
      }
    }
  }
  return output;
}
