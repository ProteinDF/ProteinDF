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

#include <cctype>
#include <sstream>

#include "DfInitialGuessHarris.h"
#include "DfInitialGuessHuckel.h"
#include "DfPreScf.h"
#include "tl_dense_general_matrix_lapack.h"

#include "CnError.h"
#include "DfDmatrix.h"
#include "TlStringTokenizer.h"

DfPreScf::DfPreScf(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {}

DfPreScf::~DfPreScf() {}

// ・占有数を持つファイルを作る
// ・必要なら初期ＭＯのファイルと、規格直交表現の初期ＭＯを作り、ファイルに保存する
// ・必要なら密度行列からRou1を作る
// ・Rou1からMyu1, Nyu1を作る
void DfPreScf::prepareGuess() {
  switch (this->initialGuessType_) {
    case GUESS_LCAO: {
      switch (this->m_nMethodType) {
        case METHOD_RKS:
          this->createInitialGuessUsingLCAO(RUN_RKS);
          break;

        case METHOD_UKS:
          this->createInitialGuessUsingLCAO(RUN_UKS_ALPHA);
          this->createInitialGuessUsingLCAO(RUN_UKS_BETA);
          break;

        case METHOD_ROKS:
          this->createInitialGuessUsingLCAO(RUN_ROKS);
          break;

        default:
          CnErr.abort();
          break;
      }
    } break;

    case GUESS_RHO:
    case GUESS_DENSITY: {
      switch (this->m_nMethodType) {
        case METHOD_RKS:
          this->createOccupation(RUN_RKS);
          break;

        case METHOD_UKS:
          this->createOccupation(RUN_UKS_ALPHA);
          this->createOccupation(RUN_UKS_BETA);
          break;

        case METHOD_ROKS:
          this->createOccupation(RUN_ROKS);
          break;

        default:
          CnErr.abort();
          break;
      }
    } break;

    case GUESS_HUCKEL:  // go below
    case GUESS_CORE: {
      DfInitialGuessHuckel huckel(this->pPdfParam_);
    } break;

    case GUESS_HARRIS: {
      switch (this->m_nMethodType) {
        case METHOD_RKS: {
          DfInitialGuessHarris harris(this->pPdfParam_);
          harris.main();
        } break;

        case METHOD_UKS:
          this->logger(
              "Sorry. harris method is not supported except RKS. stop.\n");
          CnErr.abort();
          break;

        case METHOD_ROKS:
          this->logger(
              "Sorry. harris method is not supported except RKS. stop.\n");
          CnErr.abort();
          break;

        default:
          CnErr.abort();
          break;
      }
    } break;

    default: {
      this->logger("**** Error :: DfPreScf \n");
      this->logger(" Inputted scf-start-guess is wrong.\n");
      CnErr.abort();
    } break;
  }
}

// case of "scf-start-guess = lcao"
void DfPreScf::createInitialGuessUsingLCAO(const RUN_TYPE runType) {
  // read guess lcao
  const TlDenseGeneralMatrix_Lapack LCAO =
      this->getLCAO<TlDenseGeneralMatrix_Lapack>(runType);
  this->saveC0(runType, LCAO);

  // read guess occupation
  const TlDenseVector_Lapack aOccupation = this->getOccupation(runType);
  this->saveOccupation(runType, aOccupation);

  // output guess lcao in orthonormal basis to a files in fl_Work directory
  this->buildCprime0(runType, LCAO);

  {
    TlSerializeData tmpParam = *(this->pPdfParam_);
    tmpParam["orbital-overlap-correspondence-method"] = "keep";
    tmpParam["num_of_iterations"] = 0;
    DfDmatrix dfDmatrix(&tmpParam);
    dfDmatrix.DfDmatrixMain();  // RKS only?
  }
}

std::vector<int> DfPreScf::getLevel(std::string sLevel) {
  std::vector<int> answer;
  answer.clear();

  // "-" を " - " に置換
  TlUtils::replace(sLevel, "-", " - ");

  TlStringTokenizer token(sLevel);
  bool bRegionMode = false;
  int nPrevIndex = 0;
  while (token.hasMoreTokens()) {
    std::string tmp = token.nextToken();

    if (tmp == "nil") {
      continue;
    }

    if (tmp == "-") {
      if (nPrevIndex != 0) {
        bRegionMode = true;
        continue;
      } else {
        abort();
        // CnErr.abort("DfPreScf", "", "putdoclevel", "syntax error.");
      }
    }

    int nIndex = atoi(tmp.c_str());
    // std::cerr << "nIndex = " << nIndex << std::endl;
    if (nIndex > 0) {
      if (bRegionMode == true) {
        // 数字が xx - yy の形で入力
        for (int i = nPrevIndex + 1; i <= nIndex; i++) {
          answer.push_back(i);
        }
        bRegionMode = false;
        nPrevIndex = 0;
      } else {
        // 数字が単独で入力
        answer.push_back(nIndex);
        nPrevIndex = nIndex;
        continue;
      }
    } else {
      abort();
      // CnErr.abort("DfPreScf", "", "putdoclevel", "syntax error.");
    }
  }

  return answer;
}

// memo
// vectorクラスを使っているので、書き換える
TlDenseVector_Lapack DfPreScf::getOccupation(const RUN_TYPE runType) {
  TlDenseVector_Lapack occupation;
  const std::string sFile =
      std::string("./guess.occ.") + this->m_sRunTypeSuffix[runType];
  occupation.loadText(sFile.c_str());
  if (this->m_nNumOfMOs < occupation.getSize()) {
    this->logger("Occupation vector is shrinked.");
    occupation.resize(this->m_nNumOfMOs);
  }

  return occupation;
}

void DfPreScf::saveOccupation(const RUN_TYPE runType,
                              const TlDenseVector_Lapack& rOccupation) {
  const std::string sOccFileName = this->getOccupationPath(runType);
  rOccupation.save(sOccFileName);
}

// case of "scf-start-guess = rho"
void DfPreScf::createOccupation(const RUN_TYPE runType) {
  const TlSerializeData& pdfParam = *(this->pPdfParam_);

  // construct guess occupations
  TlDenseVector_Lapack guess_occ(this->m_nNumOfMOs);
  switch (runType) {
    case RUN_RKS: {
      std::vector<int> docLevel =
          this->getLevel(pdfParam["method/nsp/occlevel"].getStr());
      for (std::vector<int>::const_iterator p = docLevel.begin();
           p != docLevel.end(); p++) {
        guess_occ.set(*p - 1, 2.0);
      }
    } break;

    case RUN_UKS_ALPHA: {
      std::vector<int> aoocLevel =
          this->getLevel(pdfParam["method/sp/alpha-spin-occlevel"].getStr());
      for (std::vector<int>::const_iterator p = aoocLevel.begin();
           p != aoocLevel.end(); p++) {
        guess_occ.set(*p - 1, 1.0);
      }
    } break;

    case RUN_UKS_BETA: {
      std::vector<int> boocLevel =
          this->getLevel(pdfParam["method/sp/beta-spin-occlevel"].getStr());
      for (std::vector<int>::const_iterator p = boocLevel.begin();
           p != boocLevel.end(); p++) {
        guess_occ.set(*p - 1, 1.0);
      }
    } break;

    case RUN_ROKS: {
      std::vector<int> docLevel =
          this->getLevel(pdfParam["method/roks/closed-shell"].getStr());
      for (std::vector<int>::const_iterator p = docLevel.begin();
           p != docLevel.end(); p++) {
        guess_occ.set(*p - 1, 2.0);
      }

      std::vector<int> socLevel =
          this->getLevel(pdfParam["method/roks/open-shell"].getStr());
      for (std::vector<int>::const_iterator p = socLevel.begin();
           p != socLevel.end(); p++) {
        // nsoc
        guess_occ.set(*p - 1, 1.0);
      }
    } break;

    default:
      CnErr.abort();
      break;
  }

  // output occupation number to a files in fl_Work directory
  const std::string sOccFileName = this->getOccupationPath(runType);
  guess_occ.save(sOccFileName);
}
