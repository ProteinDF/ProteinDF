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
#include "TlLogging.h"
#include "TlOrbitalInfo.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"

class DfLocalize : public DfObject {
    //    protected:
    //     struct Orb_QA_Item {
    //        public:
    //         Orb_QA_Item(DfObject::index_type o = 0, double q = 0.0) : orb(o), qa(q){};

    //        public:
    //         DfObject::index_type orb;
    //         double qa;
    //     };

    //     struct do_OrbQAItem_sort_functor_cmp {
    //         bool operator()(const Orb_QA_Item& a, const Orb_QA_Item& b) const {
    //             return (a.qa > b.qa);
    //         };
    //     };

   public:
    DfLocalize(TlSerializeData* pPdfParam);
    virtual ~DfLocalize();

   public:
    void setRestart(const bool yn);
    void setCMatrixPath(const std::string& path);

    void exec();

    virtual double localize(TlDenseGeneralMatrix_Lapack* pC);

   protected:
    double localize_core(TlDenseGeneralMatrix_Lapack* pC, const index_type startMO1, const index_type endMO1,
                         const index_type startMO2, const index_type endMO2);

   protected:
    virtual void initialize();
    void checkOpenMP();

    bool getSMatrix(TlDenseSymmetricMatrix_Lapack* pS);
    std::string getCMatrixPath();
    void getCMatrix(TlDenseGeneralMatrix_Lapack* pC);

    void makeGroup();
    double calcG(const TlDenseGeneralMatrix_Lapack& C, const index_type startMO, const index_type endMO);

    void getRotatingMatrix(const double A_ij, const double B_ij, const double normAB,
                           TlDenseGeneralMatrix_Lapack* pRot);
    void rotateVectors(TlDenseVector_Lapack* pCpi, TlDenseVector_Lapack* pCpj, const TlDenseGeneralMatrix_Lapack& rot);

   protected:
    double calcQA_ii(const TlDenseGeneralMatrix_Lapack& C, const index_type orb_i);
    void calcQA_ij(const TlDenseVector_Lapack& Cpi, const TlDenseVector_Lapack& Cpj, double* pA_ij, double* pB_ij);

   protected:
    TlLogging& log_;

    bool isRestart_;
    std::string CMatrixPath_;

    int lo_iteration_;
    int maxIteration_;
    double threshold_;

    TlOrbitalInfo orbInfo_;
    int numOfOcc_;

    index_type startOrb_;
    index_type endOrb_;

    double G_;

    std::vector<TlDenseVector_Lapack> group_;

    TlDenseSymmetricMatrix_Lapack S_;

    // std::vector<Orb_QA_Item> orb_QA_table_;

   protected:
    typedef std::pair<index_type, index_type> TaskItem;
    std::vector<TaskItem> getTaskList(const index_type startMO, const index_type endMO);
    std::vector<TaskItem> getTaskList(const index_type startMO1, const index_type endMO1, const index_type startMO2,
                                      const index_type endMO2);  // for parallel

   protected:
    void initLockMO(const index_type numOfMOs);
    void lockMO(const index_type MO);
    void unlockMO(const index_type MO);
    bool isLockedMO(const index_type mo1, const index_type mo2) const;

    std::vector<char> lockMOs_;
};

#endif  // DFLOCALIZE_H
