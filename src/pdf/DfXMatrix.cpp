#include <cassert>
#include <fstream>
#include <string>
#include <cmath>

#include "DfXMatrix.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "Fl_GlobalinputX.h"


DfXMatrix::DfXMatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
    assert(pPdfParam != NULL);
    const TlSerializeData& pdfParam = *pPdfParam;
    this->threshold_trancation = pdfParam["orbital-independence-threshold"].getDouble();
}


DfXMatrix::~DfXMatrix()
{
}


void DfXMatrix::main()
{
    this->exec(); // template function for parallel
}


void DfXMatrix::saveNumOfIndependentBasis()
{
    // parameterファイルに設定
    (*(this->pPdfParam_))["num_of_MOs"] = this->m_nNumOfMOs;
}


void DfXMatrix::exec()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    TlVector s;
    TlMatrix U;

    const int numOfAOs = this->m_nNumOfAOs;

    // MO dimension will be defined by follow code.
    this->loggerTime(" start");
    {
        this->loggerTime(" diagonalization of S matrix");
        TlVector EigVal;
        TlMatrix EigVec;
        {
            TlSymmetricMatrix Spq = this->getSpqMatrix<TlSymmetricMatrix>();
            Spq.diagonal(&EigVal, &EigVec);
        }
        assert(EigVal.getSize() == numOfAOs);


        // MO dimension is defined.
        this->loggerTime(" truncation of linear dependent");
        const int independent_orbital_number = pdfParam["num_of_MOs"].getInt();
        if (independent_orbital_number > 0) {
            this->m_nNumOfMOs = independent_orbital_number;
            this->logger(TlUtils::format(" set number_independent_basis (number) = %d\n", this->m_nNumOfMOs));
        } else {
            const double threshold = this->threshold_trancation;
            int cutoffCount = 0;
            for (int k = 0; k < numOfAOs; ++k) {
                if (EigVal.get(k) < threshold) {
                    ++cutoffCount;
                } else {
                    // 固有値は小さい方から格納されているはずである。
                    break;
                }
            }
            this->m_nNumOfMOs = numOfAOs - cutoffCount;

            this->saveNumOfIndependentBasis();
        }

        this->loggerTime(" generation of U matrix");
        const int numOfMOs = this->m_nNumOfMOs;
        const int cutoffBasis = numOfAOs - numOfMOs;
        s.resize(numOfMOs);
        U.resize(numOfAOs, numOfMOs);
#pragma omp parallel for
        for (int k = 0; k < numOfMOs; ++k) {
            const int index = cutoffBasis + k;
            s.set(k, EigVal.get(index));

            for (int y = 0; y < numOfAOs; ++y) {
                U.set(y, k, EigVec.get(y, index));
            }
        }
    }


    // 以下の計算で必要なため、sの各要素の平方根を求める
    const int numOfMOs = this->m_nNumOfMOs;
    TlVector sqrt_s(numOfMOs);
#pragma omp parallel for
    for (int i = 0; i < numOfMOs; ++i) {
        sqrt_s.set(i, std::sqrt(s.get(i)));
    }


    this->loggerTime(" generation of X matrix");
    {
        //TlSymmetricMatrix S_12(numOfMOs);
        TlMatrix S_12(numOfMOs, numOfMOs);
        for (int i = 0; i < numOfMOs; ++i) {
            S_12.set(i, i, (1.0 / sqrt_s.get(i)));
        }

        //const TlMatrix X = U * S_12;
        TlMatrix X = U;
        X *= S_12;

        DfObject::saveXMatrix(X);
    }


    this->loggerTime(" generation and invert of X matrix");
    {
        TlSymmetricMatrix S(numOfMOs);
        for (int i = 0; i < numOfMOs; ++i) {
            S.set(i, i, sqrt_s.get(i));
        }

        U.transpose();
        //const TlMatrix invX = S * U;
        TlMatrix invX = S;
        invX *= U;

        invX.save("fl_Work/fl_Mtr_InvX.matrix");
    }

    this->loggerTime(" finalize");

    // for a check
    if (TlUtils::toUpper(pdfParam["DfXmatrix/check-matrixes"].getStr()) == "YES") {
        this->checkMatrixes();
    }
}


void DfXMatrix::checkMatrixes()
{
    this->logger("@@@@ we get four matrices\n");

    TlMatrix X = DfObject::getXMatrix<TlMatrix>();
    TlMatrix invX;
    invX.load("fl_Work/fl_Mtr_InvX.matrix");

    {
        TlMatrix B = X * invX;
        this->log_.info("neary equal to be 1 ?, X matrix * inverce of X matrix");
        {
            std::stringstream ss;
            B.print(ss);
            this->log_.info(ss.str());
        }
    }

    {
        TlMatrix C = invX * X;
        this->log_.info("must be 1, inverce of X matrix * X matrix");
        {
            std::stringstream ss;
            C.print(ss);
            this->log_.info(ss.str());
        }
    }
}

