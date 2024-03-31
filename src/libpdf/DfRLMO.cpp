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

#include "DfRLMO.h"

#include "TlMatrix.h"

DfRLMO::DfRLMO(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {}

DfRLMO::~DfRLMO() {}

void DfRLMO::exec(const std::vector<index_type>& startBlockAOs) {
    const int iteration = this->m_nIteration;

    const TlDenseSymmetricMatrix_Lapack D_AO =
        this->getPInMatrix<TlDenseSymmetricMatrix_Lapack>(RUN_RKS, iteration);
    const index_type numOfAOs = D_AO.getNumOfRows();
    // D_AO.save("D_AO.mat");

    // CD
    {
        TlDenseSymmetricMatrix_Lapack CD = 0.5 * D_AO;
        std::vector<int> pivot;
        int b = choleskyFactorization(&CD, &pivot);
        std::cerr << "CD=" << b << std::endl;
        CD.save("CD.mat");

        TlVector p(pivot.size());
        for (int i = 0; i < pivot.size(); ++i) {
            p[i] = pivot[i];
        }
        p.save("pivot.vtr");

        TlMatrix L(numOfAOs, numOfAOs);
        for (index_type i = 0; i < numOfAOs; ++i) {
            for (index_type j = 0; j <= i; ++j) {
                const double v = CD.get(i, j);
                L.set(i, j, v);
            }
        }
        L.save("L0.mat");

        // pivot
        TlMatrix M(numOfAOs, numOfAOs);
        for (index_type i = 0; i < numOfAOs; ++i) {
            index_type target = pivot[i] - 1;
            // if (i != target) {
            std::cerr << TlUtils::format("pivot: %d <=> %d", i, target)
                      << std::endl;
            TlMatrix row_part1 = L.getBlockMatrix(i, 0, 1, numOfAOs);
            TlMatrix row_part2 = L.getBlockMatrix(target, 0, 1, numOfAOs);
            M.setBlockMatrix(target, 0, row_part2);
            M.setBlockMatrix(i, 0, row_part1);

            TlMatrix col_part1 = L.getBlockMatrix(0, i, numOfAOs, 1);
            TlMatrix col_part2 = L.getBlockMatrix(0, target, numOfAOs, 1);
            M.setBlockMatrix(0, target, col_part2);
            M.setBlockMatrix(0, i, col_part1);
            // }
        }
        M.save("M.mat");

        TlMatrix tM = M;
        tM.transpose();
        TlMatrix MM = M * tM;
        MM.save("MM.mat");
    }

    // std::cerr << this->getCMatrixPath(RUN_RKS, this->m_nIteration) <<
    // std::endl;  const TlMatrix C_AO = this->getCMatrix<TlMatrix>(RUN_RKS,
    // this->m_nIteration);
    TlMatrix C_CMO_AO;
    C_CMO_AO.load(this->getCMatrixPath(RUN_RKS, this->m_nIteration));
    C_CMO_AO.save("C_CMO_AO.mat");

    // {
    //     const TlDenseSymmetricMatrix_Lapack S =
    //     this->getSpqMatrix<TlDenseSymmetricMatrix_Lapack>();
    //     TlMatrix tC_CMO_AO = C_CMO_AO;
    //     tC_CMO_AO.transpose();
    //     TlMatrix CSC = tC_CMO_AO * S * C_CMO_AO;
    //     CSC.save("CSC.mat");
    // }

    // make X
    const TlDenseSymmetricMatrix_Lapack X = this->getX();
    // check
    // {
    //     // XX == Spq
    //     TlMatrix XX = X * X;
    //     XX.save("XX.mtx");
    // }

    TlDenseSymmetricMatrix_Lapack Xinv = X;
    Xinv.inverse();
    Xinv.save("Xinv.mat");
    // check
    // {
    //     TlMatrix XX = X * Xinv;
    //     XX.save("XXinv.mtx");
    // }

    // TH
    // {
    //     TlMatrix XC = X * C_CMO_AO;
    //     XC.save("XC.mat");
    // }

    const TlDenseSymmetricMatrix_Lapack D_OAO = X * D_AO * X;
    // D_OAO.save("D_OAO.mtx");

    // check
    // {
    //     // D_OAO2 == _2D_OAO
    //     TlMatrix D_OAO2 = D_OAO * D_OAO;
    //     D_OAO2.save("D_OAO2.mtx");
    //     TlDenseSymmetricMatrix_Lapack _2D_OAO = 2.0 * D_OAO;
    //     _2D_OAO.save("_2D_OAO.mtx");
    // }

    const TlMatrix T = this->getT(D_OAO, startBlockAOs);
    T.save("T.mat");
    TlMatrix Tt = T;
    Tt.transpose();
    // Tt.save("Tt.mtx");

    // check
    // {
    //     TT == E
    //     TlMatrix TTt = T * Tt;
    //     TTt.save("TTt.mtx");

    //     TlMatrix TtT = Tt * T;
    //     TtT.save("TtT.mtx");
    // }

    const TlDenseSymmetricMatrix_Lapack D_RO = Tt * D_OAO * T;
    D_RO.save("D_RO.mat");

    const TlMatrix C_CMO_RO = Tt * X * C_CMO_AO;
    // C_CMO_RO.save("C_CMO_RO.mtx");

    // TH
    {
        TlMatrix C_CMO_RO = Tt * X * C_CMO_AO;
        C_CMO_RO.save("C_CMO_RO.mat");
    }

    TlMatrix U;
    {
        TlVector e;
        D_RO.diagonal(&e, &U);
        e.save("e.vct");
        U.save("U0.mat");

        const index_type numOfRows = U.getNumOfRows();
        const index_type numOfCols = U.getNumOfCols();
        TlMatrix tmp(numOfRows, numOfCols);
        for (int r = 0; r < numOfRows; ++r) {
            for (int c = 0; c < numOfCols; ++c) {
                tmp.set(r, c, U.get(r, numOfCols - c - 1));
            }
        }
        U = tmp;
        U.save("U.mat");
    }

    TlMatrix Ut = U;
    Ut.transpose();

    // check
    // {
    //     // UU = E
    //     TlMatrix UU = U * Ut;
    //     UU.save("UU.mtx");
    // }

    TlMatrix D_RLMO = Ut * D_RO * U;
    D_RLMO.save("D_RLMO.mat");

    TlMatrix C_CMO_RLMO = Ut * Tt * X * C_CMO_AO;
    C_CMO_RLMO.save("C_CMO_RLMO.mat");

    const TlMatrix XinvT = Xinv * T;
    const TlMatrix XinvTU = XinvT * U;
    const TlMatrix& C_RLMO_AO = XinvTU;
    XinvT.save("X_1T.mat");
    XinvTU.save("X_1TU.mat");

    // check
    // {
    //     TlMatrix C_RLMO_AO_t = C_RLMO_AO;
    //     TlMatrix tCC = C_RLMO_AO_t * C_CMO_RLMO;
    //     tCC.save("tCC.mtx");

    //     // int numOfRows = CC.getNumOfRows();
    //     // int numOfCols = CC.getNumOfCols();
    //     // TlMatrix CC_(numOfRows, numOfCols);
    //     // for (int c = 0; c < numOfCols; ++c) {
    //     //     for (int r = 0; r <numOfRows; ++r) {
    //     //         CC_.set(r, numOfCols - c -1, CC.get(r, c));
    //     //     }
    //     // }
    //     // CC_.save("CC_.mtx");
    // }

    // TlMatrix C_TH = C_AO * C_RLMO_AO;
    // C_TH.save("C_TH.mtx");
}

TlDenseSymmetricMatrix_Lapack DfRLMO::getX() {
    const TlDenseSymmetricMatrix_Lapack S =
        this->getSpqMatrix<TlDenseSymmetricMatrix_Lapack>();
    TlVector e;
    TlMatrix V;
    S.diagonal(&e, &V);

    TlMatrix V_dagger = V;
    V_dagger.transpose();

    const std::size_t size = e.getSize();
    TlMatrix E(size, size);
    for (std::size_t i = 0; i < size; ++i) {
        // E.set(i, i, 1.0 / std::sqrt(e[i]));
        E.set(i, i, std::sqrt(e[i]));
    }

    TlDenseSymmetricMatrix_Lapack X = V * E * V_dagger;
    X.save("X.mat");

    return X;
}

TlMatrix DfRLMO::getT(const TlDenseSymmetricMatrix_Lapack& D,
                      const std::vector<index_type>& startBlockAOs) {
    const index_type numOfAOs = this->m_nNumOfAOs;
    const std::size_t numOfBlocks = startBlockAOs.size();
    std::vector<index_type> startBlockAOs_mod = startBlockAOs;  // 番兵用
    startBlockAOs_mod.push_back(numOfAOs);

    std::vector<TlDenseSymmetricMatrix_Lapack> D_subblocks(numOfBlocks);

    index_type numOfOccupied = 0;
    index_type numOfVacant = 0;
    TlMatrix T_occ(numOfAOs, numOfAOs);
    TlMatrix T_vac(numOfAOs, numOfAOs);
    for (int blockID = 0; blockID < numOfBlocks; ++blockID) {
        const index_type startAO = startBlockAOs_mod[blockID];
        const index_type endAO = startBlockAOs_mod[blockID + 1];
        std::cerr << TlUtils::format("T[%d] %d -- %d", blockID, startAO, endAO)
                  << std::endl;
        const index_type distance = endAO - startAO;

        const TlDenseSymmetricMatrix_Lapack Dsub =
            D.getBlockMatrix(startAO, startAO, distance, distance);
        std::cerr << TlUtils::format("Dsub size = %d x %d", Dsub.getNumOfRows(),
                                     Dsub.getNumOfCols())
                  << std::endl;

        TlMatrix eigvec;
        TlVector eigval;
        Dsub.diagonal(&eigval, &eigvec);

        std::vector<index_type> orbTable(distance);
        for (index_type orb = 0; orb < distance; ++orb) {
            const double e = eigval[orb];

            bool isOccupied = false;
            if (std::fabs(e - 2.0) < 0.25) {
                // occupied
                isOccupied = true;
            } else if (std::fabs(e - 0.0) < 0.25) {
                // vacant
                isOccupied = false;
            } else {
                std::string check_str = "";
                if (e > 1.0) {
                    isOccupied = true;
                    check_str = "OCC";
                } else {
                    isOccupied = false;
                    check_str = "VAC";
                }

                std::cerr << TlUtils::format(
                                 "[WARN] blockID=%d, eigval=%f -> [%s]",
                                 blockID, e, check_str.c_str())
                          << std::endl;
            }

            if (isOccupied == true) {
                for (index_type i = 0; i < distance; ++i) {
                    T_occ.set(startAO + i, numOfOccupied, eigvec.get(i, orb));
                }
                ++numOfOccupied;
            } else {
                for (index_type i = 0; i < distance; ++i) {
                    T_vac.set(startAO + i, numOfVacant, eigvec.get(i, orb));
                }
                ++numOfVacant;
            }
        }
    }

    std::cerr << TlUtils::format("occ: %d/%d", numOfOccupied, numOfAOs)
              << std::endl;
    std::cerr << TlUtils::format("vac: %d/%d", numOfVacant, numOfAOs)
              << std::endl;
    if ((numOfOccupied + numOfVacant) != numOfAOs) {
        std::cerr << "[WARN] not consistent!" << std::endl;
    }
    T_occ.resize(numOfAOs, numOfOccupied);
    T_vac.resize(numOfAOs, numOfVacant);

    TlMatrix T(numOfAOs, numOfAOs);
    T.setBlockMatrix(0, 0, T_occ);
    T.setBlockMatrix(0, numOfOccupied, T_vac);

    return T;
}
