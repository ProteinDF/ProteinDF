#ifndef DFINITIALGUESSHUCKEL_H
#define DFINITIALGUESSHUCKEL_H

#include <string>

#include "DfInitialGuess.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlOrbitalInfo.h"

/// 拡張Huckel法による初期値を作成する
class DfInitialGuessHuckel : public DfInitialGuess {
public:
    DfInitialGuessHuckel(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuessHuckel();

public:
    virtual void createGuess();
    
protected:
    void initialize();

    template<typename SymmetricMatrixType, typename MatrixType>
    void createGuess();

    template<typename SymmetricMatrixType>
    SymmetricMatrixType getHuckelMatrix();

    double getHii(const std::string& sAtomName, int nOrbitalType);

    template<typename SymmetricMatrixType, typename MatrixType>
    void generatePMatrix(const MatrixType& C, const TlVector& occ,
                         SymmetricMatrixType& P1, SymmetricMatrixType& P2);

protected:
    static const double EPS;
    std::map<std::string, std::map<int, double> > m_Hii; // [atom_number][orbital_type] = value
};

template<typename SymmetricMatrixType, typename MatrixType>
void DfInitialGuessHuckel::createGuess()
{
    SymmetricMatrixType F;
    switch (initialGuessType_) {
    case GUESS_HUCKEL:
        F = this->getHuckelMatrix<SymmetricMatrixType>();
        break;
    case GUESS_CORE:
        F = DfObject::getHpqMatrix<SymmetricMatrixType>();
        break;
    default:
        this->log_.critical("unknown guess type");
        abort();
        break;
    }

    MatrixType X = DfObject::getXMatrix<MatrixType>();
    {
        MatrixType tX = X;
        tX.transpose();

        F = tX * F * X;
    }

    MatrixType C;
    {
        TlVector eigval;
        F.diagonal(&eigval, &C);
    }
    C = X * C;

    // P
    switch (this->m_nMethodType) {
    case METHOD_RKS: {
        const TlVector occ = this->createOccupation(RUN_RKS);
        this->saveCMatrix(RUN_RKS, 0, C);
        this->makeDensityMatrix();

        // SymmetricMatrixType P, P1, P2;
        // this->generatePMatrix(C, occ, P1, P2);
        // P = 2.0 * P1;
        // P += P2;
        // DfObject::savePpqMatrix(RUN_RKS, 0, P);
    }
    break;

    case METHOD_UKS: {
        const TlVector occA = this->createOccupation(RUN_UKS_ALPHA);
        const TlVector occB = this->createOccupation(RUN_UKS_BETA);
        this->saveCMatrix(RUN_UKS_ALPHA, 0, C);
        this->saveCMatrix(RUN_UKS_BETA, 0, C);
        this->makeDensityMatrix();

        // SymmetricMatrixType P, P1, P2;
        // this->generatePMatrix(C, occA, P1, P2);
        // P = 2.0 * P1;
        // P += P2;
        // DfObject::savePpqMatrix(RUN_UKS_ALPHA, 0, P);

        // this->generatePMatrix(C, occB, P1, P2);
        // P = 2.0 * P1;
        // P += P2;
        // DfObject::savePpqMatrix(RUN_UKS_BETA, 0, P);
    }
    break;

    case METHOD_ROKS:
        this->log_.critical("to implement");
        abort();
        break;

    default:
        this->log_.critical("unknown method");
        abort();
        break;
    }
}

template <typename SymmetricMatrixType>
SymmetricMatrixType DfInitialGuessHuckel::getHuckelMatrix()
{
    TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                          (*this->pPdfParam_)["basis_sets"]);

    const int nNumOfOrbs = orbInfo.getNumOfOrbitals();

    SymmetricMatrixType Hpq = DfObject::getHpqMatrix<SymmetricMatrixType>();
    assert(nNumOfOrbs == Hpq.getNumOfRows());

    SymmetricMatrixType Spq = DfObject::getSpqMatrix<SymmetricMatrixType>();
    assert(nNumOfOrbs == Spq.getNumOfRows());

    const double K = 1.75;

    // create Fock matrix using modified extended Huckel
    SymmetricMatrixType F(nNumOfOrbs);
    int nPrevAtomI = -1;
    int nPrevShellType = -1;
    for (int i=0; i < nNumOfOrbs; i++) {
        const int nAtomI = orbInfo.getAtomIndex(i);
        const int nShellTypeI = orbInfo.getShellType(i);
        const double Hii = this->getHii(orbInfo.getAtomName(i), nShellTypeI);

        // i != j
        for (int j=0; j < i; j++) {
            const int nAtomJ = orbInfo.getAtomIndex(j);
            if (nAtomI != nAtomJ) {
                const double Hjj = this->getHii(orbInfo.getAtomName(j), orbInfo.getShellType(j));
                F(i, j) = 0.5 * K * (Hii + Hjj) * Spq(i, j);
            }
        }

        // i == j
        if ((nPrevAtomI != nAtomI) || (nPrevShellType != nShellTypeI)) {
            F(i, i) = Hii;
            nPrevAtomI = nAtomI;
            nPrevShellType = nShellTypeI;
        }
    }

    // F.print(std::cout);
    return F;
}


template<typename SymmetricMatrixType, typename MatrixType>
void DfInitialGuessHuckel::generatePMatrix(const MatrixType& C, const TlVector& occ,
                                           SymmetricMatrixType& P1, SymmetricMatrixType& P2)
{
    // generate density matrix
    const int nNumOfRows = C.getNumOfRows();
    const int nNumOfCols = C.getNumOfCols();

    P1 = SymmetricMatrixType(nNumOfRows);
    P2 = SymmetricMatrixType(nNumOfRows);

    for (int p = 0; p < nNumOfRows; p++) {
        for (int q = 0; q < nNumOfCols; q++) {
            double p1 = 0.0;
            double p2 = 0.0;

            for (int i = 0; i < occ.getSize(); i++) {
                double occNum = occ[i];

                if (std::fabs(occNum - 2.0) < EPS) {
                    p1 += C(p, i) * C(q, i);
                } else if (std::fabs(occNum - 1.0) < EPS) {
                    p2 += C(p, i) * C(q, i);
                }
            }

            P1(p, q) = p1;
            P2(p, q) = p2;
        }
    }
}

#endif // DFINITIALGUESSHUCKEL_H
