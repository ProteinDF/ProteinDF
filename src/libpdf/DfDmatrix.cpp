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

#include <cmath>
#include <ios>

#include "DfDmatrix.h"
#include "TlFile.h"
#include "TlStringTokenizer.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_blas_old.h"
#include "tl_dense_symmetric_matrix_blas_old.h"
#include "tl_dense_vector_blas.h"

/*********************************************************
MO_OVERLAP_ITER:
軌道の重なりの対応を使用しはじめるiteration回数を指定する。
**********************************************************/

DfDmatrix::DfDmatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {
  const TlSerializeData& pdfParam = *pPdfParam;

  this->orbitalCorrespondenceMethod_ = OCM_NONE;
  const bool isOrbitalCorrespondence =
      pdfParam["orbital-correspondence"].getBoolean();
  const int startItr = pdfParam["orbital-correspondence-start"].getInt();
  if (isOrbitalCorrespondence == true) {
    const std::string method =
        TlUtils::toUpper(pdfParam["orbital-correspondence-method"].getStr());
    if ((method == "MO-OVERLAP") && (this->m_nIteration >= startItr)) {
      this->orbitalCorrespondenceMethod_ = OCM_OVERLAP;
    } else if ((method == "MO-PROJECTION") &&
               (this->m_nIteration >= startItr)) {
      this->orbitalCorrespondenceMethod_ = OCM_PROJECTION;
    }
  }
}

DfDmatrix::~DfDmatrix() {}

void DfDmatrix::DfDmatrixMain() {
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->main(RUN_RKS);
      break;

    case METHOD_UKS:
      this->main(RUN_UKS_ALPHA);
      this->main(RUN_UKS_BETA);
      break;

    case METHOD_ROKS:
      this->main(RUN_ROKS_CLOSED);
      this->main(RUN_ROKS_OPEN);
      break;

    default:
      CnErr.abort();
      break;
  }
}

void DfDmatrix::main(const DfObject::RUN_TYPE runType) {
  // occupation
  TlVector_BLAS currOcc;
  switch (this->orbitalCorrespondenceMethod_) {
    case OCM_OVERLAP:
      this->log_.info(" orbital correspondence method: MO-overlap");
      currOcc =
          this->getOccupationUsingOverlap<TlDenseGeneralMatrix_BLAS_old>(runType);
      currOcc.save(this->getOccupationPath(runType));
      break;

    case OCM_PROJECTION:
      this->log_.info(" orbital correspondence method: MO-projection");
      currOcc = this->getOccupationUsingProjection<TlDenseGeneralMatrix_BLAS_old,
                                                   TlDenseSymmetricMatrix_BLAS_Old>(
          runType);
      currOcc.save(this->getOccupationPath(runType));
      break;

    default:
      this->log_.info(" orbital correspondence method: none");
      currOcc.load(this->getOccupationPath(runType));
      break;
  }

  this->generateDensityMatrix<TlDenseGeneralMatrix_BLAS_old,
                              TlDenseSymmetricMatrix_BLAS_Old>(runType, currOcc);
}

// TlVector_BLAS DfDmatrix::getOccupation(const DfObject::RUN_TYPE runType)
// {
//     const std::string sFileName = this->getOccupationPath(runType);

//     TlVector_BLAS occ;
//     occ.load(sFileName);
//     assert(occ.getSize() == this->m_nNumOfMOs);

//     return occ;
// }

void DfDmatrix::checkOccupation(const TlVector_BLAS& prevOcc,
                                const TlVector_BLAS& currOcc) {
  const double xx = prevOcc.sum();
  const double yy = currOcc.sum();

  if (std::fabs(xx - yy) > 1.0e-10) {
    this->log_.error("SUM pre_occ != SUM crr_occ");
    this->log_.error(TlUtils::format(
        " SUM pre_occ is %10.4lf,  SUM crr_occ is %10.4lf\n", xx, yy));
    this->log_.error("previous occupation");
    {
      std::stringstream ss;
      prevOcc.print(ss);
      this->log_.error(ss.str());
    }
    this->log_.error("current occupation");
    {
      std::stringstream ss;
      currOcc.print(ss);
      this->log_.error(ss.str());
    }

    CnErr.abort("DfDmatrix", "", "", "SUM error for occ !!");
  }
}

void DfDmatrix::printOccupation(const TlVector_BLAS& occ) {
  std::stringstream ss;
  occ.print(ss);
  this->log_.info(ss.str());
}

// print out Two Vectors' elements
void DfDmatrix::printTwoVectors(const TlVector_BLAS& a, const TlVector_BLAS& b,
                                const std::string& title, int pnumcol) {
  assert(a.getSize() == b.getSize());
  this->log_.info(TlUtils::format("\n\n       %s\n\n", title.c_str()));

  this->log_.info("       two vectors");
  const index_type number_of_emt = a.getSize();
  for (int ord = 0; ord < number_of_emt; ord += pnumcol) {
    this->log_.info("       ");
    for (int j = ord; ((j < ord + pnumcol) && (j < number_of_emt)); ++j) {
      this->log_.info(TlUtils::format("   %5d th", j + 1));
    }
    this->log_.info("\n     ");

    for (int j = ord; ((j < ord + pnumcol) && (j < number_of_emt)); ++j) {
      this->log_.info("-----------");
    }
    this->log_.info("----\n       ");

    for (int j = ord; j < ord + pnumcol && j < number_of_emt; ++j) {
      const double aj = a[j];
      this->log_.info(TlUtils::format(" %6.0lf    ", aj + 1));
    }
    this->log_.info("\n\n       ");

    for (int j = ord; ((j < ord + pnumcol) && (j < number_of_emt)); ++j) {
      const double bj = b[j];
      this->log_.info(TlUtils::format(" %10.6lf", bj));
    }
    this->log_.info("\n\n");
  }
}

// TlVector_BLAS DfDmatrix::createOccupation(const DfObject::RUN_TYPE runType)
// {
//     const TlSerializeData& pdfParam = *(this->pPdfParam_);

//     // construct guess occupations
//     TlVector_BLAS occ(this->m_nNumOfMOs);
//     switch (runType) {
//     case RUN_RKS: {
//         std::vector<int> docLevel =
//         this->getLevel(pdfParam["method/rks/occlevel"].getStr()); for
//         (std::vector<int>::const_iterator p = docLevel.begin(); p !=
//         docLevel.end(); p++) {
//             occ[*p -1] = 2.0;
//         }
//     }
//     break;

//     case RUN_UKS_ALPHA: {
//         std::vector<int> aoocLevel =
//         this->getLevel(pdfParam["method/uks/alpha_spin_occlevel"].getStr());
//         for (std::vector<int>::const_iterator p = aoocLevel.begin(); p !=
//         aoocLevel.end(); p++) {
//             occ[*p -1] = 1.0;
//         }
//     }
//     break;

//     case RUN_UKS_BETA: {
//         std::vector<int> boocLevel =
//         this->getLevel(pdfParam["method/uks/beta_spin_occlevel"].getStr());
//         for (std::vector<int>::const_iterator p = boocLevel.begin(); p !=
//         boocLevel.end(); p++) {
//             occ[*p -1] = 1.0;
//         }
//     }
//     break;
//     case RUN_ROKS: {
//         std::vector<int> docLevel =
//         this->getLevel(pdfParam["method/roks/closed-shell"].getStr()); for
//         (std::vector<int>::const_iterator p = docLevel.begin(); p !=
//         docLevel.end(); p++) {
//             occ[*p -1] = 2.0;
//         }

//         std::vector<int> socLevel =
//         this->getLevel(pdfParam["method/roks/open-shell"].getStr()); for
//         (std::vector<int>::const_iterator p = socLevel.begin(); p !=
//         socLevel.end(); p++) {
//             // nsoc
//             occ[*p -1] = 1.0;
//         }
//     }
//     break;

//     default:
//         std::cerr << "program error. (DfDmatrix::createOccupation)" <<
//         std::endl; break;
//     }

//     return occ;
// }

// std::vector<int> DfDmatrix::getLevel(std::string sLevel)
// {
//     std::vector<int> answer;
//     answer.clear();

//     // "-" を " - " に置換
//     TlUtils::replace(sLevel, "-", " - ");

//     TlStringTokenizer token(sLevel);
//     bool bRegionMode = false;
//     int nPrevIndex = 0;
//     while (token.hasMoreTokens()) {
//         std::string tmp = token.nextToken();

//         if (tmp == "nil") {
//             continue;
//         }

//         if (tmp == "-") {
//             if (nPrevIndex != 0) {
//                 bRegionMode = true;
//                 continue;
//             } else {
//                 abort();
//                 //CnErr.abort("DfPreScf", "", "putdoclevel", "syntax
//                 error.");
//             }
//         }

//         int nIndex = atoi(tmp.c_str());
//         //std::cerr << "nIndex = " << nIndex << std::endl;
//         if (nIndex > 0) {
//             if (bRegionMode == true) {
//                 // 数字が xx - yy の形で入力
//                 for (int i = nPrevIndex +1; i <= nIndex; i++) {
//                     answer.push_back(i);
//                 }
//                 bRegionMode = false;
//                 nPrevIndex = 0;
//             } else {
//                 // 数字が単独で入力
//                 answer.push_back(nIndex);
//                 nPrevIndex = nIndex;
//                 continue;
//             }
//         } else {
//             abort();
//             //CnErr.abort("DfPreScf", "", "putdoclevel", "syntax error.");
//         }
//     }

//     return answer;
// }
