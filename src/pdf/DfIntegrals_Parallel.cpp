#include "DfIntegrals_Parallel.h"
#include "DfHpq_Parallel.h"
#include "DfOverlap_Parallel.h"
#include "DfOverlapX_Parallel.h"
#include "DfEri_Parallel.h"
#include "DfEriX_Parallel.h"

#include "DfHpqX_Parallel.h"
#include "DfCD_Parallel.h"
#include "DfXMatrix_Parallel.h"
#include "DfInvMatrix_Parallel.h"
#include "DfGenerateGrid_Parallel.h"

#include "Fl_Geometry.h"
#include "TlCommunicate.h"
#include "TlSymmetricMatrix.h"

DfIntegrals_Parallel::DfIntegrals_Parallel(TlSerializeData* pParam)
    : DfIntegrals(pParam)
{
     TlCommunicate& rComm = TlCommunicate::getInstance();
     DfObject::rank_ = rComm.getRank();
//     std::cerr << TlUtils::format("[%d] DfIntegrals_Parallel constructor called.", rComm.getRank())
//               << std::endl;
}


DfIntegrals_Parallel::~DfIntegrals_Parallel()
{
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     std::cerr << TlUtils::format("[%d] DfIntegrals_Parallel destructor called.", rComm.getRank())
//               << std::endl;
}


void DfIntegrals_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfIntegrals::logger(str);
    }
}


DfCD* DfIntegrals_Parallel::getDfCDObject()
{
    DfCD *pDfCD = new DfCD_Parallel(this->pPdfParam_);
    return pDfCD;
}


DfXMatrix* DfIntegrals_Parallel::getDfXMatrixObject()
{
    DfXMatrix* pDfXMatrix = new DfXMatrix_Parallel(this->pPdfParam_);

    return pDfXMatrix;
}


DfInvMatrix* DfIntegrals_Parallel::getDfInvMatrixObject()
{
    DfInvMatrix* pDfInvMatrix = new DfInvMatrix_Parallel(this->pPdfParam_);

    return pDfInvMatrix;
}


DfGenerateGrid* DfIntegrals_Parallel::getDfGenerateGridObject()
{
    DfGenerateGrid* pDfGenerateGrid = new DfGenerateGrid_Parallel(this->pPdfParam_);

    return pDfGenerateGrid;
}


void saveInvSquareVMatrix(const TlSymmetricMatrix& v)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster()) {
        v.save("fl_Work/fl_Mtr_invSquareV.matrix");
    }
    rComm.barrier();
}


void DfIntegrals_Parallel::outputStartTitle(const std::string& stepName, const char lineChar)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfIntegrals::outputStartTitle(stepName, lineChar);
    }
}

void DfIntegrals_Parallel::outputEndTitle(const std::string& stepName, const char lineChar)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfIntegrals::outputEndTitle(stepName, lineChar);
    }
}


void DfIntegrals_Parallel::createHpqMatrix()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->createHpqMatrix_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK

    this->createHpqMatrix_LAPACK();
}


void DfIntegrals_Parallel::createHpqMatrix_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const std::size_t needMem = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) * sizeof(double);
    if ((this->isWorkOnDisk_ == true) || (this->procMaxMemSize_ < needMem)) {
        this->logger(" H_pq is build up on disk.\n");
        TlMatrix::useMemManager(true);
    } else {
        this->logger(" H_pq is build up on memory.\n");
        TlMatrix::useMemManager(false);
    }
    TlSymmetricMatrix Hpq(this->m_nNumOfAOs);
    TlSymmetricMatrix Hpq2(this->m_nNumOfAOs);

    if (this->isUseNewEngine_ == true) {
        this->logger(" use new engine\n");
        DfHpqX_Parallel dfHpq(this->pPdfParam_);
        dfHpq.getHpq(&Hpq, &Hpq2);
    } else {
        DfHpq_Parallel dfHpq(this->pPdfParam_);
        dfHpq.getHpq(&Hpq, &Hpq2);
    }

    if (rComm.isMaster() == true) {
        this->saveHpqMatrix(Hpq);
        this->saveHpq2Matrix(Hpq2);
    }
}


void DfIntegrals_Parallel::createHpqMatrix_ScaLAPACK()
{
    TlDistributeSymmetricMatrix Hpq(this->m_nNumOfAOs);
    TlDistributeSymmetricMatrix Hpq2(this->m_nNumOfAOs);
    
    this->logger(" Hpq build using distribute matrix.\n");
    if (this->isUseNewEngine_ == true) {
        this->logger(" use new engine\n");
        DfHpqX_Parallel dfHpq(this->pPdfParam_);
        dfHpq.getHpqD(&Hpq, &Hpq2);
    } else {
        DfHpq_Parallel dfHpq(this->pPdfParam_);
        dfHpq.getHpq(&Hpq, &Hpq2);
    }
    
    this->saveHpqMatrix(Hpq);
    this->saveHpq2Matrix(Hpq2);
}


void DfIntegrals_Parallel::createOverlapMatrix()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->createOverlapMatrix_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK

    this->createOverlapMatrix_LAPACK();
}


void DfIntegrals_Parallel::createOverlapMatrix_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();
    DfOverlap_Parallel dfOverlap(this->pPdfParam_);
    DfOverlapX_Parallel dfOverlapX(this->pPdfParam_);

    // Spq
    if ((calcState & DfIntegrals::Spq) == 0) {
        this->outputStartTitle("Spq");

        std::size_t needMem = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2 * sizeof(double);
        needMem += rComm.getWorkMemSize();
        if ((this->isWorkOnDisk_ == true) || (this->procMaxMemSize_ < needMem)) {
            this->logger(" S_(p q) is build on disk.\n");
            TlMatrix::useMemManager(true);
        } else {
            TlMatrix::useMemManager(false);
            this->logger(" S_(p q) is build on memory.\n");
        }
        
        TlSymmetricMatrix Spq(this->m_nNumOfAOs);
        if (this->isUseNewEngine_ == true) {
            this->logger(" use new engine.\n");
            dfOverlapX.getSpq(&Spq);
        } else {
            dfOverlap.getSpq(&Spq);
        }

        if (rComm.isMaster() == true) {
            this->saveSpqMatrix(Spq);
        }
        rComm.barrier();
        
        this->outputEndTitle();

        calcState |= DfIntegrals::Spq;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }

    // Sab2
    if ((calcState & DfIntegrals::Sab2) == 0) {
        this->outputStartTitle("Sab2");

        const std::size_t needMem = this->m_nNumOfAux * (this->m_nNumOfAux + 1) / 2 * sizeof(double);
        if ((this->isWorkOnDisk_ == true) || (this->procMaxMemSize_ < needMem)) {
            this->logger(" S_(alpha beta) is build on disk.\n");
            TlMatrix::useMemManager(true);
        } else {
            this->logger(" S_(alpha beta) is build on memory.\n");
            TlMatrix::useMemManager(false);
        }

        TlSymmetricMatrix Sab2(this->m_nNumOfAux);
        if (this->isUseNewEngine_ == true) {
            this->logger(" use new engine.\n");
            dfOverlapX.getSab(&Sab2);
        } else {
            dfOverlap.getSab2(&Sab2);
        }

        if (rComm.isMaster() == true) {
            this->saveSab2Matrix(Sab2);
        }
        rComm.barrier();
        
        this->outputEndTitle();

        calcState |= DfIntegrals::Sab2;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }

    // Sgd
    if ((calcState & DfIntegrals::Sgd) == 0) {
        if (this->m_bIsXCFitting == true) {
            this->outputStartTitle("Sgd");

            const std::size_t needMem = this->numOfAuxXC_ * (this->numOfAuxXC_ + 1) / 2 * sizeof(double);
            if ((isWorkOnDisk_ == true) ||
                (this->procMaxMemSize_ < needMem)) {
                this->logger(" S_(gamma delta) is build on disk.\n");
                TlMatrix::useMemManager(true);
            } else {
                this->logger(" S_(gamma delta) is build on memory.\n");
                TlMatrix::useMemManager(false);
            }

            TlSymmetricMatrix Sgd(this->numOfAuxXC_);
            dfOverlap.getSgd(&Sgd);

            if (rComm.isMaster() == true) {
                this->saveSgdMatrix(Sgd);
            }
            rComm.barrier();
        
            this->outputEndTitle();
        }

        calcState |= DfIntegrals::Sgd;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }

    // Na
    if ((calcState & DfIntegrals::Na) == 0) {
        this->outputStartTitle("N_alpha");

        const std::size_t needMem = this->m_nNumOfAux * sizeof(double);
        if ((this->isWorkOnDisk_ == true) ||
            (this->procMaxMemSize_ < needMem)) {
            this->logger(" [alpha] is build on disk.\n");
            TlMatrix::useMemManager(true);
        } else {
            this->logger(" [alpha] is build on memory.\n");
            TlMatrix::useMemManager(false);
        }

        TlVector Na(this->m_nNumOfAux);
        dfOverlap.getNa(&Na);
        this->saveNalpha(Na);
        
        this->outputEndTitle();

        calcState |= DfIntegrals::Na;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }
}


void DfIntegrals_Parallel::createOverlapMatrix_ScaLAPACK()
{
    // TlCommunicate& rComm = TlCommunicate::getInstance();
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();
    DfOverlap_Parallel dfOverlap(this->pPdfParam_);
    DfOverlapX_Parallel dfOverlapX(this->pPdfParam_);

    // Spq
    if ((calcState & DfIntegrals::Spq) == 0) {
        this->outputStartTitle("Spq");

        this->logger(" S_(p q) is build using distribute matrix.\n");
        TlDistributeSymmetricMatrix Spq(this->m_nNumOfAOs);
        if (this->isUseNewEngine_ == true) {
            this->logger(" use new engine.\n");
            dfOverlapX.getSpqD(&Spq);
        } else {
            dfOverlap.getSpq(&Spq);
        }
            this->saveSpqMatrix(Spq);
            
        
        this->outputEndTitle();

        calcState |= DfIntegrals::Spq;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }

    // Sab2
    if ((calcState & DfIntegrals::Sab2) == 0) {
        this->outputStartTitle("Sab2");

        this->logger(" S_(alpha beta) is build using on distribute matrix.\n");
        TlDistributeSymmetricMatrix Sab2(this->m_nNumOfAux);
        if (this->isUseNewEngine_ == true) {
            this->logger(" use new engine.\n");
            dfOverlapX.getSabD(&Sab2);
        } else {
            dfOverlap.getSab2(&Sab2);
        }
        this->saveSab2Matrix(Sab2);

        this->outputEndTitle();

        calcState |= DfIntegrals::Sab2;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }

    // Sgd
    if ((calcState & DfIntegrals::Sgd) == 0) {
        if (this->m_bIsXCFitting == true) {
            this->outputStartTitle("Sgd");

            TlDistributeSymmetricMatrix Sgd(this->numOfAuxXC_);
            dfOverlap.getSgd(&Sgd);
            this->saveSgdMatrix(Sgd);

            this->outputEndTitle();
        }

        calcState |= DfIntegrals::Sgd;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }

    // Na
    if ((calcState & DfIntegrals::Na) == 0) {
        this->outputStartTitle("N_alpha");

        const std::size_t needMem = this->m_nNumOfAux * sizeof(double);
        if ((this->isWorkOnDisk_ == true) ||
            (this->procMaxMemSize_ < needMem)) {
            this->logger(" [alpha] is build on disk.\n");
            TlMatrix::useMemManager(true);
        } else {
            this->logger(" [alpha] is build on memory.\n");
            TlMatrix::useMemManager(false);
        }

        TlVector Na(this->m_nNumOfAux);
        dfOverlap.getNa(&Na);
        this->saveNalpha(Na);
        
        this->outputEndTitle();

        calcState |= DfIntegrals::Na;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }
}


void DfIntegrals_Parallel::createERIMatrix()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->createERIMatrix_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK

    this->createERIMatrix_LAPACK();
}

    
void DfIntegrals_Parallel::createERIMatrix_LAPACK()
{
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

    if ((calcState & DfIntegrals::Sab) == 0) {
        this->outputStartTitle("Sab");

        const std::size_t needMem = this->m_nNumOfAux * (this->m_nNumOfAux + 1) / 2 * sizeof(double);
        if ((this->isWorkOnDisk_ == true) ||
            (this->procMaxMemSize_ < needMem)) {
            this->logger(" <alpha|beta> is build on disk.\n");
            TlMatrix::useMemManager(true);
        } else {
            this->logger(" <alpha|beta> is build on memory.\n");
            TlMatrix::useMemManager(false);
        }
        
        TlSymmetricMatrix Sab(this->m_nNumOfAux);

        if (this->isUseNewEngine_ == true) {
            this->logger(" use new engine.\n");
            DfEriX_Parallel dfEri(this->pPdfParam_);
            dfEri.getJab(&Sab);
        } else {
            DfEri_Parallel dfEri(this->pPdfParam_);
            dfEri.getSab(&Sab);
        }

        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            this->saveSabMatrix(Sab);
        }
        
        this->outputEndTitle();

        calcState |= DfIntegrals::Sab;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }
}


void DfIntegrals_Parallel::createERIMatrix_ScaLAPACK()
{
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

    if ((calcState & DfIntegrals::Sab) == 0) {
        this->outputStartTitle("Sab");

        DfEri_Parallel dfEri(this->pPdfParam_);
        TlDistributeSymmetricMatrix Sab(this->m_nNumOfAux);
        dfEri.getSab(&Sab);
        this->saveSabMatrix(Sab);
        
        this->outputEndTitle();

        calcState |= DfIntegrals::Sab;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }
}
