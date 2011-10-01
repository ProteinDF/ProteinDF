#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfXMatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"

#define FAST_TRANCATE

DfXMatrix_Parallel::DfXMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfXMatrix(pPdfParam)
{
    this->fastTrancate_ = (*pPdfParam)["fast_trancate"].getBoolean();
}


DfXMatrix_Parallel::~DfXMatrix_Parallel()
{
}


void DfXMatrix_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfXMatrix::logger(str);
    }
}


void DfXMatrix_Parallel::saveNumOfIndependentBasis()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfXMatrix::saveNumOfIndependentBasis();
    }
}


void DfXMatrix_Parallel::main()
{
#ifdef HAVE_SCALAPACK
    if (m_bUsingSCALAPACK == true) {
        this->exec_ScaLAPACK();
    } else {
        this->exec_LAPACK();
    }
#else
    {
        this->exec_LAPACK();
    }
#endif // HAVE_SCALAPACK
}


void DfXMatrix_Parallel::exec_LAPACK()
{
    this->logger("X matrix is created using LAPACK.\n");

    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfXMatrix::main();
    }
    rComm.barrier();
}


void DfXMatrix_Parallel::exec_ScaLAPACK()
{
    //TlCommunicate& rComm = TlCommunicate::getInstance();
    //rComm.barrier();
    
    this->logger("X matrix is created using ScaLAPACK.\n");

    const int numOfAOs = this->m_nNumOfAOs;
    TlVector s;
    TlDistributeMatrix U(this->m_nNumOfAOs, this->m_nNumOfAOs);

    // MO dimension will be defined by follow code.
    this->loggerTime(" start");
    {
        this->loggerTime(" diagonalization of S matrix");
        TlVector EigVal;
        TlDistributeMatrix EigVec;
        {
            TlDistributeSymmetricMatrix Spq = this->getSpqMatrix<TlDistributeSymmetricMatrix>();
            Spq.diagonal(&EigVal, &EigVec);
        }
        assert(EigVal.getSize() == numOfAOs);


        // MO dimension is defined.
        this->loggerTime(" truncation of linear dependent");
        const int independent_orbital_number = (*(this->pPdfParam_))["num_of_MOs"].getInt();
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

#ifdef FAST_TRANCATE
        {
            this->logger(" fast trancate routine applied.\n");
            TlDistributeMatrix trans(numOfAOs, numOfMOs);
            for (int i = 0; i < numOfMOs; ++i) {
                int index = cutoffBasis + i;
                trans.set(index, i, 1.0);
            }
            U = EigVec * trans;
            
            for (int i = 0; i < numOfMOs; ++i) {
                int index = cutoffBasis + i;
                s.set(i, EigVal.get(index));
            }
        }
#else
        {
            for (int i = 0; i < numOfMOs; ++i) {
                const int index = cutoffBasis + i;
                s.set(i, EigVal.get(index));
                
                for (int y = 0; y < numOfAOs; ++y) {
                    U.set(y, i, EigVec.get(y, index));
                }
            }
        }
#endif // FAST_TRANCATE        
    }

    // 以下の計算で必要なため、sの各要素の平方根を求める
    const int numOfMOs = this->m_nNumOfMOs;
    TlVector sqrt_s(numOfMOs);
    for (int i = 0; i < numOfMOs; ++i) {
        sqrt_s.set(i, std::sqrt(s.get(i)));
    }

    this->loggerTime(" generation of X matrix");
    {
        TlDistributeSymmetricMatrix S_12(numOfMOs);
        for (int i = 0; i < numOfMOs; ++i) {
            S_12.set(i, i, (1.0 / sqrt_s.get(i)));
        }

        const TlDistributeMatrix X = U * S_12;

        X.save(this->getXMatrixPath());
    }

    this->loggerTime(" generation and invert of X matrix");
    {
        TlDistributeSymmetricMatrix S(numOfMOs);
        for (int i = 0; i < numOfMOs; ++i) {
            S.set(i, i, sqrt_s.get(i));
        }

        U.transpose();

        const TlDistributeMatrix invX = S * U;

        invX.save(this->getInvXMatrixPath());
    }

    this->loggerTime(" finalize");
}



