#include <cassert>
#include "DfEri2.h"

#include "DfEri.h"
#include "FileX.h"

#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "TlLogX.h"
#include "TlTime.h"

DfEri2::DfEri2(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam)
{
    this->m_bUseDensityMatrix = true;
}

DfEri2::~DfEri2()
{
}

TlSymmetricMatrix DfEri2::getKMatrix()
{
    TlSymmetricMatrix K;

    if (this->m_bUseDensityMatrix == true) {
        // 密度行列を使うルーチン
        std::cout << "DfEri2::getKMatrix() use density." << std::endl;

        TlSymmetricMatrix PA = this->getPMatrix(this->m_nIteration -1);
        {
            TlSymmetricMatrix prevPA = this->getPMatrix(this->m_nIteration -2);
            PA -= prevPA;
        }
        TlSymmetricMatrix diffK = this->getKMatrixUsingDensityMatrix(PA);
        K = this->loadKMatrix(this->m_nIteration -1);
        K += diffK;
        this->saveKMatrix(K, this->m_nIteration);
    } else {
        // MO vectorを使うルーチン
        std::cout << "DfEri2::getKMatrix() use MO vector." << std::endl;

        const TlMatrix C = this->loadCutoffCMatrix();
        K = this->getKMatrixUsingCMatrix(C);
    }

    K *= -1.0;

    return K;
}

TlSymmetricMatrix DfEri2::getKMatrix(const TlSymmetricMatrix& P)
{
    TlSymmetricMatrix K = this->getKMatrixUsingDensityMatrix(P);
    K *= -1.0;

    // for debug
    this->saveKMatrix(K, this->m_nIteration -1);

    return K;
}

DfEri* DfEri2::getDfEri()
{
    DfEri* pDfEri = new DfEri(this->pPdfParam_);

    return pDfEri;
}

TlSymmetricMatrix DfEri2::generateInvSquareVMatrix()
{
    TlSymmetricMatrix V;
    V.load("fl_Work/fl_Mtr_Sab.matrix");
    assert(V.getNumOfRows() == this->m_nNumOfAux);

    TlVector s(this->m_nNumOfAux);
    TlMatrix U(this->m_nNumOfAux, this->m_nNumOfAux);
    //int nTrancatedSize = 0;
    {
        //     TlVector tmp_s;
        //     TlMatrix tmp_U;

        V.diagonal(&s, &U);
        assert(s.getSize() == this->m_nNumOfAux);
        assert(U.getNumOfRows() == this->m_nNumOfAux);
        assert(U.getNumOfCols() == this->m_nNumOfAux);

        //     // curoff
        //     {
        //       const int nSize = s.getSize();
        //       for (int i = 0; i < nSize; ++i){
        //  if (tmp_s[i] > 0.007){
        //    s[nTrancatedSize] = tmp_s[i];
        //    for (int x = 0; x < nSize; ++x){
        //      U(x, nTrancatedSize) = tmp_U(x, i);
        //    }
        //    ++nTrancatedSize;
        //  }
        //       }
        //     }
        //     s.resize(nTrancatedSize);
        //     U.resize(this->number_de_basis, nTrancatedSize);
    }

    TlMatrix tU = U;
    tU.transpose();

    // U * s^(-1/2)
    //const int nRows = U.getNumOfRows();
    //const int nCols = U.getNumOfCols();
    //assert(nRows == nTrancatedSize);
    //assert(nCols == this->number_de_basis);

    TlSymmetricMatrix r(s.getSize());
    TlVector::size_type max_i = s.getSize();
    for (TlVector::size_type i = 0; i < max_i; ++i) {
        r(i, i) = 1.0 / (std::sqrt(s[i]));
    }

    // v = Sab^(-1/2) = Ut * s^(-1/2) * U
    const TlSymmetricMatrix invSquareV = U * r * tU;

    //std::cerr << "[DfEri2] V^(-1/2)\n";
    //invSquareV.print(std::cerr);

    //invSquareV.save("fl_Work/fl_Mtr_invSquareV.matrix");
    return invSquareV;
}

TlSymmetricMatrix DfEri2::loadInvSquareVMatrix()
{
    TlSymmetricMatrix invSquareV;
    invSquareV.load("fl_Work/fl_Mtr_invSquareV.matrix");

    return invSquareV;
}

TlSymmetricMatrix DfEri2::loadKMatrix(const int nIteration)
{
    const std::string sFileName = "fl_Work/fl_Mtr_K.rks" + TlUtils::xtos(nIteration);

    TlSymmetricMatrix K(this->m_nNumOfAOs);

    if (FileX::isExist(sFileName) == true) {
        K.load(sFileName);
    }

    return K;
}

void DfEri2::saveKMatrix(const TlSymmetricMatrix& K, const int nIteration)
{
    K.save("fl_Work/fl_Mtr_K.rks" + TlUtils::xtos(nIteration));
}

TlMatrix DfEri2::loadCutoffCMatrix()
{
    TlMatrix C;
    {
        std::string type = "rks";
        TlMatrix LCAO;
        {
            // このファイルを読み込むことはRKS
            const std::string sFileName = "fl_Work/fl_Mtr_C.matrix." + type + TlUtils::xtos(this->m_nIteration -1);
            LCAO.load(sFileName);
        }

        TlVector Occ;
        {
            std::string sFileName = "fl_Work/fl_Occupation";
            if (type == "uks-alpha") {
                sFileName += "_Alpha";
            } else if (type == "uks-beta") {
                sFileName += "_Beta";
            }

            Occ.load(sFileName);
        }

        assert(static_cast<TlVector::size_type>(LCAO.getNumOfCols()) == Occ.getSize());

        int nLUMO = 0;
        for (int i = (Occ.getSize() -1); i >= 0; --i) {
            if (std::fabs(Occ[i] - 2.0) < 1.0E-10) {
                nLUMO = i +1;
                break;
            }
        }

        C = LCAO.getBlockMatrix(0, 0, LCAO.getNumOfRows(), nLUMO);
//     std::cout << "sigma = " << nLUMO << std::endl;
//     std::cout << TlUtils::format("C size = (%d, %d)", C.getNumOfRows(), C.getNumOfCols()) << std::endl;
    }

    return C;
}

TlSymmetricMatrix DfEri2::getKMatrixUsingCMatrix(const TlMatrix& C)
{
    const int nBas = this->m_nNumOfAOs;
    const int nAux = this->m_nNumOfAux;

    TlSymmetricMatrix K(nBas);

    TlSymmetricMatrix invSquareV = this->loadInvSquareVMatrix();

    // generate K
    // 密度行列と違い、C行列を読み込んでいるので、半分にする必要が無い
    {
        DfEri* pDfEri = this->getDfEri();

        for (int m = 0; m < nAux; ++m) {
            TlMatrix X;
            {
                const TlVector v = invSquareV.getColVector(m);
                TlSymmetricMatrix T(nBas);

                pDfEri->getDeltaHpqAForEri2(v, T);

                X = T * C;
            }

            TlMatrix tX(X);
            tX.transpose();

            K += X * tX;
        }

        delete pDfEri;
        pDfEri = NULL;
    }

    return K;
}

TlSymmetricMatrix DfEri2::getPMatrix(const int iteration)
{
    const std::string sFileName = this->getPpqMatrixPath(RUN_RKS, iteration);

    TlSymmetricMatrix P(this->m_nNumOfAOs);
    if (FileX::isExist(sFileName) == true) {
        P.load(sFileName);

        // rks の処理
        P *= 0.5;
    }

    return P;
}

// 密度行列を使った場合
TlSymmetricMatrix DfEri2::getKMatrixUsingDensityMatrix(const TlSymmetricMatrix& PA)
{
    const int nBas = this->m_nNumOfAOs;
    const int nAux = this->m_nNumOfAux;

    TlSymmetricMatrix K(nBas);
    const TlSymmetricMatrix invSquareV = this->loadInvSquareVMatrix();

    // for screening
    double dCutoffCoef = 1.0;
    double dAbsoluteMaxElement = PA.getMaxAbsoluteElement();
    if (dAbsoluteMaxElement < 1.0) {
        dCutoffCoef /= dAbsoluteMaxElement;
    }

    // for timing data
    double dEriTime = 0.0;
    double dMatrixOperationTime = 0.0;
    TlTime tlTime;

    // generate K
    {
        DfEri* pDfEri = this->getDfEri();

        for (int m = 0; m < nAux; ++m) {
            TlSymmetricMatrix W;
            {
                const TlVector v = invSquareV.getColVector(m);
                TlSymmetricMatrix T(nBas);

                tlTime.start();
                pDfEri->getDeltaHpqAForEri2(v, T, dCutoffCoef);
                dEriTime += tlTime.getElapseTime();

                //T.save("fl_Work/fl_Mtr_T." + TlUtils::xtos(m));

                tlTime.start();
                TlMatrix tT(T);
                tT.transpose();

                TlMatrix PAT = PA * T;

                W = tT * PAT;
                dMatrixOperationTime += tlTime.getElapseTime();
            }

            K += W;
        }

        delete pDfEri;
        pDfEri = NULL;
    }

    // for timing data
    TlLogX& Log = TlLogX::getInstance();
    Log << TlUtils::format("ERI time              = %+e\n", dEriTime)
    << TlUtils::format("matrix operation time = %+e\n", dMatrixOperationTime);

    return K;
}

