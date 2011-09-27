#include "DfRLMO.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

DfRLMO::DfRLMO(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
}


DfRLMO::~DfRLMO()
{
}


void DfRLMO::exec(const std::vector<index_type>& startBlockAOs)
{
    const int iteration = this->m_nIteration;
    const TlSymmetricMatrix D_AO = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, iteration);
    D_AO.save("D_AO.mtx");
    
    // make X
    const TlSymmetricMatrix X = this->getX();
    X.save("X.mtx");

    // check
    {
        TlMatrix XX = X * X;
        XX.save("XX.mtx");
    }
    
    const TlSymmetricMatrix D_OAO = X * D_AO * X;
    D_OAO.save("D_OAO.mtx");

    // check
    {
        TlMatrix D_OAO2 = D_OAO * D_OAO;
        D_OAO2.save("D_OAO2.mtx");
        TlSymmetricMatrix _2D_OAO = 2.0 * D_OAO;
        _2D_OAO.save("_2D_OAO.mtx");
    }
    
    const TlMatrix T = this->getT(D_OAO, startBlockAOs);
    T.save("T.mtx");
    TlMatrix Tt = T;
    Tt.transpose();
    
    const TlSymmetricMatrix D_RO = Tt * D_OAO * T;
    D_RO.save("D_RO.mtx");

    const TlMatrix C_AO = this->getCMatrix<TlMatrix>(RUN_RKS, this->m_nIteration);
    const TlMatrix C_RO = Tt * X * C_AO;
    C_RO.save("C_RO.mtx");

    TlMatrix U;
    {
        TlVector e;
        D_RO.diagonal(&e, &U);
    }
    TlMatrix Ut = U;
    Ut.transpose();
    
    TlMatrix D_RLMO = Ut * D_RO * U;
    D_RLMO.save("D_RLMO.mtx");

    TlMatrix C_RLMO = Ut * Tt * X * C_AO;
    C_RLMO.save("C_RLMO.mtx");

    TlSymmetricMatrix Xinv = X;
    Xinv.inverse();
    Xinv.save("Xinv.mtx");

    const TlMatrix C_RLMO_AO = Xinv * T * U;
    C_RLMO_AO.save("C_RLMO_AO.mtx");

    // check
    {
        TlMatrix C_RLMO_AO_t = C_RLMO_AO;
        TlMatrix CC = C_RLMO_AO_t * C_RLMO;
        CC.save("CC.mtx");

        // int numOfRows = CC.getNumOfRows();
        // int numOfCols = CC.getNumOfCols();
        // TlMatrix CC_(numOfRows, numOfCols);
        // for (int c = 0; c < numOfCols; ++c) {
        //     for (int r = 0; r <numOfRows; ++r) {
        //         CC_.set(r, numOfCols - c -1, CC.get(r, c));
        //     }
        // }
        // CC_.save("CC_.mtx");
    }
}


TlSymmetricMatrix DfRLMO::getX()
{
    const TlSymmetricMatrix S = this->getSpqMatrix<TlSymmetricMatrix>();
    TlVector e;
    TlMatrix V;
    S.diagonal(&e, &V);
        
    TlMatrix V_dagger = V;
    V_dagger.transpose();
        
    const std::size_t size = e.getSize();
    TlMatrix E(size, size);
    for (std::size_t i = 0; i < size; ++i) {
        //E.set(i, i, 1.0 / std::sqrt(e[i]));
        E.set(i, i, std::sqrt(e[i]));
    }
    
    TlSymmetricMatrix X = V * E * V_dagger;
    X.save("X.mtx");

    return X;
}


TlMatrix DfRLMO::getT(const TlSymmetricMatrix& D,
                      const std::vector<index_type>& startBlockAOs)
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    const std::size_t numOfBlocks = startBlockAOs.size();
    std::vector<index_type> startBlockAOs_mod = startBlockAOs; // 番兵用
    startBlockAOs_mod.push_back(numOfAOs);
    
    std::vector<TlSymmetricMatrix> D_subblocks(numOfBlocks);

    index_type numOfOccupied = 0;
    index_type numOfVacant = 0;
    TlMatrix T_occ(numOfAOs, numOfAOs);
    TlMatrix T_vac(numOfAOs, numOfAOs);
    for (int blockID = 0; blockID < numOfBlocks; ++blockID) {
        const index_type startAO = startBlockAOs_mod[blockID];
        const index_type endAO = startBlockAOs_mod[blockID +1];
        const index_type distance = endAO - startAO;

        const TlSymmetricMatrix Dsub = D.getBlockMatrix(startAO, startAO, distance, distance);

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

                std::cerr << TlUtils::format("[WARN] blockID=%d, eigval=%f -> [%s]",
                                             blockID, e, check_str.c_str())
                          << std::endl;
            }

            if (isOccupied == true) {
                for (index_type i = 0; i < distance; ++i) {
                    T_occ.set(startAO + i, numOfOccupied,
                              eigvec.get(i, orb));
                }
                ++numOfOccupied;
            } else {
                for (index_type i = 0; i < distance; ++i) {
                    T_vac.set(startAO + i, numOfVacant,
                              eigvec.get(i, orb));
                }
                ++numOfVacant;
            }
        }
    }

    std::cerr << TlUtils::format("occ: %d/%d", numOfOccupied, numOfAOs) << std::endl;
    std::cerr << TlUtils::format("vac: %d/%d", numOfVacant, numOfAOs) << std::endl;
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

