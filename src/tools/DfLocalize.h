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

#ifndef DFLOCALIZE_H
#define DFLOCALIZE_H

#include <string>
#include <vector>

#include "DfObject.h"
#include "TlOrbitalInfo.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"

class DfLocalize : public DfObject {
 protected:
  struct Orb_QA_Item {
   public:
    Orb_QA_Item(std::size_t o = 0, double q = 0.0) : orb(o), qa(q){};

   public:
    std::size_t orb;
    double qa;
  };

  struct do_OrbQAItem_sort_functor_cmp {
    bool operator()(const Orb_QA_Item& a, const Orb_QA_Item& b) const {
      return (a.qa > b.qa);
    };
  };

  struct JobItem {
   public:
    JobItem(std::size_t i = 0, std::size_t j = 0) : orb_i(i), orb_j(j) {}

   public:
    std::size_t orb_i;
    std::size_t orb_j;
  };

 public:
  DfLocalize(TlSerializeData* pPdfParam);
  virtual ~DfLocalize();

 public:
  virtual void localize(const std::string& inputCMatrixPath = "");

 protected:
  void makeGroup();
  void makeQATable();
  void calcQA(const std::size_t orb_i, const std::size_t orb_j, double* pA_ij,
              double* pB_ij);

  void getRotatingMatrix(const double A_ij, const double B_ij,
                         const double normAB,
                         TlDenseGeneralMatrix_Lapack* pRot);
  void rotateCmatrix(std::size_t orb_i, std::size_t orb_j,
                     const TlDenseGeneralMatrix_Lapack& rot);

 protected:
  double calcQA(const std::size_t orb_i);

  void makeJobList();
  virtual int getJobItem(DfLocalize::JobItem* pJob, bool isInitialized = false);

 protected:
  int maxIteration_;
  double threshold_;

  TlOrbitalInfo orbInfo_;
  int numOfOcc_;

  std::size_t startOrb_;
  std::size_t endOrb_;

  std::string SMatrixPath_;
  std::string CMatrixPath_;

  std::vector<std::vector<std::size_t> > group_;
  TlDenseSymmetricMatrix_Lapack S_;
  TlDenseGeneralMatrix_Lapack C_;

  std::vector<Orb_QA_Item> orb_QA_table_;
  std::vector<JobItem> jobList_;
};

#endif  // DFLOCALIZE_H
