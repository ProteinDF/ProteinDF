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
#include <cmath>
#include <ios>

#include "DfDmatrix.h"
#include "TlFile.h"
#include "TlStringTokenizer.h"
#include "TlUtils.h"
#include "common.h"

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
        const std::string method = TlUtils::toUpper(
            pdfParam["orbital-correspondence-method"].getStr());
        if ((method == "MO-OVERLAP") && (this->m_nIteration >= startItr)) {
            this->orbitalCorrespondenceMethod_ = OCM_OVERLAP;
        } else if ((method == "MO-PROJECTION") &&
                   (this->m_nIteration >= startItr)) {
            this->orbitalCorrespondenceMethod_ = OCM_PROJECTION;
        }
    }
}

DfDmatrix::~DfDmatrix() {}

void DfDmatrix::run() {
    switch (this->linearAlgebraPackage_) {
        case LAP_LAPACK: {
            this->log_.info("Linear Algebra Package: LAPACK");
            this->run_impl<TlDenseGeneralMatrix_Lapack,
                           TlDenseSymmetricMatrix_Lapack,
                           TlDenseVector_Lapack>();
        } break;

#ifdef HAVE_EIGEN
        case LAP_EIGEN: {
            this->log_.info("Linear Algebra Package: Eigen");
            this->run_impl<TlDenseGeneralMatrix_Eigen,
                           TlDenseSymmetricMatrix_Eigen, TlDenseVector_Eigen>();
        } break;
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
        case LAP_VIENNACL: {
            this->log_.info("Linear Algebra Package: ViennaCL");
            this->run_impl<TlDenseGeneralMatrix_ViennaCL,
                           TlDenseSymmetricMatrix_ViennaCL,
                           TlDenseVector_ViennaCL>();
        } break;
#endif  // HAVE_VIENNACL

        default:
            CnErr.abort(
                TlUtils::format("program error: @%s,%d", __FILE__, __LINE__));
    }
}

void DfDmatrix::checkOccupation(const TlDenseVectorObject& prevOcc,
                                const TlDenseVectorObject& currOcc) {
    const double xx = prevOcc.sum();
    const double yy = currOcc.sum();

    if (std::fabs(xx - yy) > 1.0e-10) {
        this->log_.error("SUM pre_occ != SUM crr_occ");
        this->log_.error(TlUtils::format(
            " SUM pre_occ is %10.4lf,  SUM crr_occ is %10.4lf\n", xx, yy));
        this->log_.error("previous occupation");
        {
            std::stringstream ss;
            ss << prevOcc << std::endl;
            this->log_.error(ss.str());
        }
        this->log_.error("current occupation");
        {
            std::stringstream ss;
            ss << currOcc << std::endl;
            this->log_.error(ss.str());
        }

        CnErr.abort("DfDmatrix", "", "", "SUM error for occ !!");
    }
}

void DfDmatrix::printOccupation(const TlDenseVectorObject& occ) {
    std::stringstream ss;
    ss << occ << std::endl;
    this->log_.info(ss.str());
}

// print out Two Vectors' elements
void DfDmatrix::printTwoVectors(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::string& title, int pnumcol) {
    assert(a.size() == b.size());
    this->log_.info(TlUtils::format("\n\n       %s\n\n", title.c_str()));

    this->log_.info("       two vectors");
    const index_type number_of_emt = a.size();
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

// ----------------------------------------------------------------------------
