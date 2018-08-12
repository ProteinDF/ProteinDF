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

#include "DfQclo.h"
#include <stdlib.h>
#include "DfDiagonal.h"
#include "DfLevelshift.h"
#include "DfTransFmatrix.h"
#include "DfTransatob.h"
#include "Fl_Fragment.h"
#include "Fl_Tbl_Fragment.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

DfQclo::DfQclo(TlSerializeData* pPdfParam, int num_iter, bool bExecDiis)
    : DfObject(pPdfParam),
      m_bExecDiis(bExecDiis),
      FlFrag(Fl_Geometry((*pPdfParam)["coordinates"])) {
  //     TlLogX& Log = TlLogX::getInstance();
  //     this->number_iteration = num_iter;
  //     this->scftype         = this->m_flGbi["SCF"]["method"];
  //     this->number_ao_basis =
  //     atoi(this->m_flGbi["SCF"]["control-norb"].c_str());

  this->number_fragment = FlFrag.getNumOfFragments();

  // output informations
  //     Log << "number_iteration       = " << this->number_iteration << "\n";
  //     Log << "number_ao_basis        = " << this->number_ao_basis  << "\n";
  //     Log << "number_fragment        = " << this->number_fragment  << "\n";
}

DfQclo::~DfQclo() {}

void DfQclo::DfQcloMain() {
  const int maxNumOfFragments = this->number_fragment;
  for (int i = 0; i < maxNumOfFragments; ++i) {
    const std::string fragname = "frag" + TlUtils::xtos(i);
    const int norbcut = FlFrag.getNumOfOrbitals(i);

    // transformed to QCLO based Kohn-Sham matrix
    {
      this->loggerStartTitle("DfTransFmatrix");

      DfTransFmatrix dfTransFmatrix(this->pPdfParam_, this->m_bExecDiis);
      dfTransFmatrix.DfTrsFmatQclo(fragname, norbcut);

      this->loggerEndTitle();
    }

    // add level shift to Kohn-Sham matrix
    {
      int start_iter =
          (*(this->pPdfParam_))["level_shift/start_iteration"].getInt();
      if (((*(this->pPdfParam_))["SCF"]["level_shift"].getBoolean()) &&
          (this->m_nIteration >= start_iter)) {
        this->loggerStartTitle("DfLevelshift");

        DfLevelshift LS(this->pPdfParam_, this->m_nIteration);
        LS.DfLshiftQclo(fragname, norbcut);

        this->loggerEndTitle();
      }
    }

    // Diagonalize Kohn-Sham matrix
    {
      this->loggerStartTitle("DfDiagonal");

      DfDiagonal dfDiagonal(this->pPdfParam_);

      switch (this->m_nMethodType) {
        case METHOD_RKS:
          dfDiagonal.DfDiagQclo(RUN_RKS, fragname, norbcut);
          break;

        case METHOD_UKS:
          dfDiagonal.DfDiagQclo(RUN_UKS_ALPHA, fragname, norbcut);
          dfDiagonal.DfDiagQclo(RUN_UKS_BETA, fragname, norbcut);
          break;

        case METHOD_ROKS:
          dfDiagonal.DfDiagQclo(RUN_ROKS, fragname, norbcut);
          break;

        default:
          CnErr.abort();
      }

      this->loggerEndTitle();
    }

    // transformed to original nonorth. A.O.based space
    {
      this->loggerStartTitle("DfTransatob");

      // DfTransatob dfTransatob(this->m_flGbi);
      DfTransatob dfTransatob(this->pPdfParam_);
      dfTransatob.DfTrsatobQclo(fragname, norbcut);

      this->loggerEndTitle();
    }
  }

  switch (this->m_nMethodType) {
    case METHOD_RKS:
      combineCqclo("rks", this->m_nIteration);
      break;

    case METHOD_UKS:
      combineCqclo("uks-alpha", this->m_nIteration);  // UKS alpha spin
      combineCqclo("uks-beta", this->m_nIteration);   // UKS beta spin
      break;

    case METHOD_ROKS:
      combineCqclo("roks", this->m_nIteration);
      break;

    default:
      CnErr.abort();
      break;
  }
}

void DfQclo::combineCqclo(const std::string& runtype, int iteration) {
  // read Cqclo matrix of each fragment
  int number_dimension_mo = 0;
  std::vector<int> number_dimension_qclo(this->number_fragment);
  for (int frag = 0; frag < this->number_fragment; frag++) {
    const int mo = FlFrag.getNumOfOrbitals(frag);
    number_dimension_qclo[frag] = mo;
    number_dimension_mo += mo;
  }

  // combine Cqclo matrix of each matrix
  TlDenseGeneralMatrix_Lapack A(this->m_nNumOfAOs, number_dimension_mo);
  for (int frag = 0; frag < this->number_fragment; frag++) {
    const std::string fname = "fl_Mtr_C.matrix.frag" + TlUtils::xtos(frag) +
                              "." + runtype + TlUtils::xtos(iteration);
    TlDenseGeneralMatrix_Lapack B;
    B.load(fname);

    if (B.getNumOfRows() != this->m_nNumOfAOs ||
        B.getNumOfCols() != number_dimension_qclo[frag]) {
      this->log_.warn(TlUtils::format("rowDim of %d = %d", fname.c_str(),
                                      B.getNumOfRows()));
      this->log_.warn(TlUtils::format("colDim of %d = %d", fname.c_str(),
                                      B.getNumOfCols()));
      this->log_.warn(
          TlUtils::format("number_dimension_basis = %d", this->m_nNumOfAOs));
      this->log_.warn(TlUtils::format("number_dimension_qclo = %d",
                                      number_dimension_qclo[frag]));
      this->log_.warn("DfQclo dimension is not consistency, but continue.");
    }

    Fl_Tbl_Fragment Tfrag(Fl_Geometry((*this->pPdfParam_)["coordinates"]));
    for (int qclo = 0; qclo < number_dimension_qclo[frag]; qclo++) {
      int mo = 0;
      if (runtype == "rks" || runtype == "roks") {
        mo = Tfrag.getQclo(frag, qclo);
      } else if (runtype == "uks-alpha") {
        mo = Tfrag.getQcloAlpha(frag, qclo);
      } else if (runtype == "uks-beta") {
        mo = Tfrag.getQcloBeta(frag, qclo);
      }
      for (int basis = 0; basis < this->m_nNumOfAOs; ++basis) {
        A.set(basis, mo, B.get(basis, qclo));
      }
    }
  }

  //     if (pout_debug <= -5) {
  //         Log << "@@@@ -5>=pout_debug\n";
  //         Log << "combined C matrix\n";
  //         A.print(Log);
  //     }

  // write combined C matrix
  A.save("fl_Mtr_C.matrix." + runtype + TlUtils::xtos(iteration));

  // read eigenvalue of each fragment
  std::vector<TlDenseVector_Lapack> eigval_frag(this->number_fragment);
  {
    for (int frag = 0; frag < this->number_fragment; frag++) {
      eigval_frag[frag].load("fl_Work/fl_Vct_Eigval.frag" +
                             TlUtils::xtos(frag) + "." + runtype +
                             TlUtils::xtos(iteration));

      assert(eigval_frag[frag].getSize() == number_dimension_qclo[frag]);
    }
  }

  // merge eigenvalue
  TlDenseVector_Lapack eigval(number_dimension_mo);
  {
    int mo = 0;
    for (int frag = 0; frag < this->number_fragment; frag++) {
      for (int qclo = 0; qclo < number_dimension_qclo[frag]; qclo++) {
        eigval.set(mo, eigval_frag[frag].get(qclo));
        mo++;
      }
    }
  }

  // sort eigenvalue, eigen vector
  for (int i = 0; i < number_dimension_mo - 1; i++) {
    for (int j = 0; j < number_dimension_mo - 1; j++) {
      if (eigval.get(j) > eigval.get(j + 1)) {
        // swap
        // std::swap(eigval[j], eigval[j + 1]);
        {
          const double tmp = eigval.get(j);
          eigval.set(j, eigval.get(j + 1));
          eigval.set(j + 1, tmp);
        }
        for (int k = 0; k < this->m_nNumOfAOs; ++k) {
          // swap
          double tmp = A.get(k, j);
          A.set(k, j, A.get(k, j + 1));
          A.set(k, j + 1, tmp);
        }
      }
    }
  }

  // for debug
  for (int i = 0; i < number_dimension_mo; i++) {
    this->log_.debug(TlUtils::format("%d th Eigval = %f", i, eigval.get(i)));
  }

  // write merged eigenvalue
  if ("rks" == runtype || "roks" == runtype) {
    eigval.save(this->getEigenvaluesPath(RUN_RKS, this->m_nIteration));
  } else if ("uks-alpha" == runtype) {
    eigval.save(this->getEigenvaluesPath(RUN_UKS_ALPHA, this->m_nIteration));
  } else if ("uks-beta" == runtype) {
    eigval.save(this->getEigenvaluesPath(RUN_UKS_BETA, this->m_nIteration));
  } else {
    CnErr.abort("DfQclo", "", "main", "runtype is illegal");
  }

  // write sorted LCAO coefficent
  const std::string fname = "result.guess.lcao." + runtype + ".2";
  std::ofstream fo;
  fo.open(fname.c_str(), std::ios::out);
  if (!fo) {
    std::cerr << "FILE OPEN ERROR" << std::endl;
    CnErr.abort("DfQclo", "", "combineCqclo()", "FILE OPEN ERROR");
  }

  fo << "TEXT"
     << "\n";
  fo << this->m_nNumOfAOs << "\n";
  fo << number_dimension_mo << "\n";

  for (int i = 0; i < this->m_nNumOfAOs; ++i) {
    for (int j = 0; j < number_dimension_mo; j++) {
      fo << TlUtils::format(" %10.6lf", A.get(i, j));
    }
    fo << "\n";
  }
  fo << "\n";
  fo.close();
}
