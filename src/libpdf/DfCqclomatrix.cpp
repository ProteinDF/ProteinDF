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
#include <fstream>
#include <string>

#include "DfCqclomatrix.h"
#include "Fl_Tbl_AtomFragment.h"
#include "Fl_Tbl_Fragment.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

DfCqclomatrix::DfCqclomatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {}

DfCqclomatrix::~DfCqclomatrix() {}

void DfCqclomatrix::main() {
  const int natom = this->m_nNumOfAtoms;
  Fl_Fragment FlFrag(Fl_Geometry((*this->pPdfParam_)["coordinates"]));
  this->number_fragment = FlFrag.getNumOfFragments();

  int* atom_fragment = new int[natom];
  int* basis_fragment = new int[this->m_nNumOfAOs];
  std::vector<int> numOfFragmentBasis(this->m_nNumOfAOs, 0);

  // read AtomFragmentTable
  {
    std::ifstream fi;
    int atom, fragment;

    fi.open("fl_Table/AtomFragmentTable");
    for (int i = 0; i < natom; i++) {
      fi >> atom >> fragment;
      atom_fragment[atom] = fragment;
    }
    fi.close();
  }

  // read OrbitalTable and make table which correlates basis and fragment
  {
    std::ifstream fi;
    int basis, atom, fragment, dummy;
    char linebuf[256];

    fi.open("fl_Table/OrbitalTable");
    for (int i = 0; i < m_nNumOfAOs; i++) {
      fi >> basis >> dummy >> atom;
      fragment = atom_fragment[atom];
      basis_fragment[basis] = fragment;
      numOfFragmentBasis[fragment]++;
      fi.getline(linebuf, 256);
    }
    fi.close();
  }

  // read S matrix
  TlDenseSymmetricMatrix_Lapack Spq =
      DfObject::getSpqMatrix<TlDenseSymmetricMatrix_Lapack>();

  // fragment loop
  int count_basis = 0;
  double* eigval = new double[m_nNumOfAOs];
  int* eigfrag = new int[m_nNumOfAOs];

  for (int frag = 0; frag < number_fragment; frag++) {
    // If number of qclo in the fragment is not set in fl_Fragment,
    // the value is calculated and written to fl_Fragment.

    // create partial S matrix
    TlDenseSymmetricMatrix_Lapack fragmentS(numOfFragmentBasis[frag]);
    int counti = 0;
    for (int i = 0; i < m_nNumOfAOs; i++) {
      if (basis_fragment[i] == frag) {
        int countj = 0;
        for (int j = 0; j <= i; j++) {
          if (basis_fragment[j] == frag) {
            fragmentS.set(counti, countj, Spq.get(i, j));
            countj++;
          }
        }
        counti++;
      }
    }
    this->log_.info(TlUtils::format("partial S matrix frag=%d", frag));

    TlDenseVector_Lapack eigVal;
    TlDenseGeneralMatrix_Lapack eigVec;
    fragmentS.eig(&eigVal, &eigVec);

    this->log_.info("eigenvalue of partial S matrix");
    {
      std::stringstream ss;
      ss << eigVal << std::endl;
      this->log_.info(ss.str());
    }

    for (int i = 0; i < numOfFragmentBasis[frag]; i++) {
      eigval[count_basis + i] = eigVal.get(i);
      eigfrag[count_basis + i] = frag;
    }
    count_basis += numOfFragmentBasis[frag];
  }

  // read guess.lcao (temporary),
  // and write the number of independent basis
  {
    std::string type;
    switch (this->m_nMethodType) {
      case METHOD_RKS:
        type = "rks";  // RKS
        break;

      case METHOD_UKS:
        type = "uks-alpha";  // UKS alpha spin
        // type = "uks-beta";   // UKS beta spin
        break;

      case METHOD_ROKS:
        type = "roks";  // ROKS
        break;

      default:
        CnErr.abort();
        break;
    }

    std::ifstream fi;
    fi.open(std::string("./guess.lcao." + type).c_str(), std::ios::in);
    if (fi.rdstate()) {
      std::cerr << "cannot open file " << ("./guess.lcao." + type) << std::endl;
      abort();
    }

    std::string dummy_line;
    int row_dimension, col_dimension;

    fi >> dummy_line;
    fi >> row_dimension >> col_dimension;
    if (row_dimension != m_nNumOfAOs) {
      CnErr.abort("DfCqclomatrix", "", "main",
                  "inputted guess lcao has illegal dimension");
    }
    this->m_nNumOfMOs = col_dimension;

    (*(this->pPdfParam_))["num_of_MOs"] = this->m_nNumOfMOs;
  }

  // count independant basis in each fragment
  // by the search of N smallest eigenvalues
  // ( N is the number of dependent basis )
  {
    for (int i = 0; i < this->m_nNumOfAOs - this->m_nNumOfMOs; i++) {
      for (int j = i + 1; j < m_nNumOfAOs; j++) {
        if (eigval[j] < eigval[i]) {
          double tmp1 = eigval[i];
          eigval[i] = eigval[j];
          eigval[j] = tmp1;
          int tmp2 = eigfrag[i];
          eigfrag[i] = eigfrag[j];
          eigfrag[j] = tmp2;
        }
      }
      numOfFragmentBasis[eigfrag[i]]--;
    }
  }

  // put data and rewrite fl_Fragment
  {
    FlFrag.calcDefault(numOfFragmentBasis);
    //     FlFrag.open( "fl_Input/fl_Fragment", "write" );
    //     FlFrag.write();
    //     FlFrag.close();
  }

  // write FragmentTable
  Fl_Tbl_Fragment Tfrag(Fl_Geometry((*this->pPdfParam_)["coordinates"]));

  // write matrices
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->main("rks");  // RKS
      break;

    case METHOD_UKS:
      this->main("uks-alpha");  // UKS alpha spin
      this->main("uks-beta");   // UKS beta spin
      break;

    case METHOD_ROKS:
      this->main("roks");  // ROKS
      break;

    default:
      CnErr.abort();
      break;
  }

  // release memory
  if (eigval) delete[] eigval;
  if (eigfrag) delete[] eigfrag;
  if (atom_fragment) delete[] atom_fragment;
  if (basis_fragment) delete[] basis_fragment;
  // if( number_fragmentbasis ) delete [] number_fragmentbasis;
}

void DfCqclomatrix::main(std::string type) {
  TlDenseGeneralMatrix_Lapack guess_lcao;
  Fl_Tbl_Fragment Tfrag(Fl_Geometry((*this->pPdfParam_)["coordinates"]));

  // read guess.lcao
  {
    std::ifstream fi;

    fi.open(std::string("./guess.lcao." + type).c_str(), std::ios::in);
    if (fi.rdstate()) {
      std::cerr << "cannot open file " << ("./guess.lcao." + type) << std::endl;
      abort();
    }

    std::string dummy_line;
    int row_dimension, col_dimension;

    fi >> dummy_line;
    fi >> row_dimension >> col_dimension;

    if (row_dimension != m_nNumOfAOs) {
      CnErr.abort("DfCqclomatrix", "", "main",
                  "inputted guess lcao has illegal dimension");
    }

    this->m_nNumOfMOs = col_dimension;

    guess_lcao.resize(this->m_nNumOfAOs, this->m_nNumOfMOs);

    for (int i = 0; i < this->m_nNumOfAOs; i++) {
      for (int j = 0; j < this->m_nNumOfMOs; j++) {
        double v;
        fi >> v;
        guess_lcao.set(i, j, v);
      }
    }

    this->log_.info("guess lcao " + type);
    {
      std::stringstream ss;
      ss << guess_lcao << std::endl;
      this->log_.info(ss.str());
    }
  }

  // write Cqclo of each fragment
  for (int frag = 0; frag < number_fragment; frag++) {
    int number_fragqclo = Tfrag.getNumberFragmentqclo(frag);
    TlDenseGeneralMatrix_Lapack Cqclo(m_nNumOfAOs, number_fragqclo);

    for (int fragqclo = 0; fragqclo < number_fragqclo; fragqclo++) {
      int qclo = 0;
      if (type == "rks" || type == "roks") {
        qclo = Tfrag.getQclo(frag, fragqclo);
      } else if (type == "uks-alpha") {
        qclo = Tfrag.getQcloAlpha(frag, fragqclo);
      } else if (type == "uks-beta") {
        qclo = Tfrag.getQcloBeta(frag, fragqclo);
      }
      for (int basis = 0; basis < this->m_nNumOfAOs; basis++) {
        Cqclo.set(basis, fragqclo, guess_lcao.get(basis, qclo));
      }
    }

    Cqclo.save("fl_Work/fl_Mtr_C.matrix.frag" + TlUtils::xtos(frag) + "." +
               type + "0");
  }

  // write Cprime.matrix.frag#.(type)0 for Level-Shift
  for (int frag = 0; frag < number_fragment; frag++) {
    int number_fragqclo = Tfrag.getNumberFragmentqclo(frag);
    TlDenseGeneralMatrix_Lapack Cprime_frag(number_fragqclo, number_fragqclo);

    // Cprime_frag = Cqclo^-1 * Cqclo
    for (int i = 0; i < number_fragqclo; i++) {
      Cprime_frag.set(i, i, 1.0);
    }

    std::string fname = "fl_Work/fl_Mtr_Cprime.matrix.frag" +
                        TlUtils::xtos(frag) + "." + type + "0";
    Cprime_frag.save(fname);
  }
}
