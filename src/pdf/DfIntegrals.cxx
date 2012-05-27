#include <cassert>

#include "DfIntegrals.h"
#include "DfHpq.h"
#include "DfHpqX.h"
#include "DfOverlap.h"
#include "DfOverlapX.h"
#include "DfEri.h"
#include "DfEriX.h"
#include "DfEri2.h"
#include "DfCD.h"

#include "DfXMatrix.h"
#include "DfInvMatrix.h"
#include "DfGenerateGrid.h"

#include "Fl_Geometry.h"
#include "Fl_GlobalinputX.h"
#include "TlTime.h"
#include "TlMsgPack.h"

DfIntegrals::DfIntegrals(TlSerializeData* pParam, const std::string& saveParamPath)
    : DfObject(pParam), saveParamPath_(saveParamPath) {
}


DfIntegrals::~DfIntegrals()
{
}


void DfIntegrals::saveParam()
{
    if (this->saveParamPath_.empty() != true) {
        TlMsgPack mpac(*(this->pPdfParam_));
        mpac.save(this->saveParamPath_);
    }
}


void DfIntegrals::main()
{
    // initialize --------------------------------------------------------
    if (this->isRestart_ == true) {
        this->logger(" restart calculation is enabled.\n");
    } else {
        (*this->pPdfParam_)["control"]["integrals_state"].set(0);
    }
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

    // Hpq --------------- -----------------------------------------------
    if ((this->isRestart_ == false) || ((calcState & DfIntegrals::Hpq) == 0)) {
        this->outputStartTitle("Hpq");
        this->createHpqMatrix();
        this->outputEndTitle();

        calcState |= DfIntegrals::Hpq;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }

    // Spq, Sab2, Sgd, Na ------------------------------------------------
    this->createOverlapMatrix();

    // Sab ---------------------------------------------------------------
    this->createERIMatrix();

//   if ((this->m_bIsXcFitting != true) &&
//       (this->isRI_K_ == true)) {
//     // V^(-1/2) ------------------------------------------------------
//     this->outputStartTitle("V^-1/2");
//     {
//       DfEri2 dfEri2(this->flGlobalinput_);
//       TlSymmetricMatrix V = dfEri2.generateInvSquareVMatrix();

//       this->saveInvSquareVMatrix(V);
//     }
//     this->outputEndTitle();
//   }

    // Cholesky Vector
    this->createCholeskyVectors();

    // X matrix
    this->createXMatrix();

    // inverse matrix
    this->createInverseMatrixes();

    // initialize grids
    this->createGrids();

    // flush
    this->matrixCache_.flush();
}


DfCD* DfIntegrals::getDfCDObject()
{
    DfCD *pDfCD = new DfCD(this->pPdfParam_);
    return pDfCD;
}

DfXMatrix* DfIntegrals::getDfXMatrixObject()
{
    DfXMatrix* pDfXMatrix = new DfXMatrix(this->pPdfParam_);

    return pDfXMatrix;
}


DfInvMatrix* DfIntegrals::getDfInvMatrixObject()
{
    DfInvMatrix* pDfInvMatrix = new DfInvMatrix(this->pPdfParam_);

    return pDfInvMatrix;
}


DfGenerateGrid* DfIntegrals::getDfGenerateGridObject()
{
    DfGenerateGrid* pDfGenerateGrid = new DfGenerateGrid(this->pPdfParam_);

    return pDfGenerateGrid;
}


void DfIntegrals::createHpqMatrix()
{
    TlSymmetricMatrix Hpq(this->m_nNumOfAOs);
    TlSymmetricMatrix Hpq2(this->m_nNumOfAOs);

    if (this->isUseNewEngine_ == true) {
        this->logger(" use new engine.\n");
        DfHpqX dfHpqX = DfHpqX(this->pPdfParam_);
        dfHpqX.getHpq(&Hpq, &Hpq2);
    } else {
        DfHpq dfHpq = DfHpq(this->pPdfParam_);
        dfHpq.getHpq(&Hpq, &Hpq2);
    }

    this->saveHpqMatrix(Hpq);
    this->saveHpq2Matrix(Hpq2);
}


void DfIntegrals::createOverlapMatrix()
{
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();
    DfOverlap dfOverlap(this->pPdfParam_);
    DfOverlapX dfOverlapX(this->pPdfParam_);

    // Spq
    if ((calcState & DfIntegrals::Spq) == 0) {
        this->outputStartTitle("Spq");

        TlSymmetricMatrix Spq(this->m_nNumOfAOs);
        if (this->isUseNewEngine_ == true) {
            this->logger(" use new engine.\n");
            dfOverlapX.getSpq(&Spq);
        } else {
            dfOverlap.getSpq(&Spq);
        }
        this->saveSpqMatrix(Spq);
        
        this->outputEndTitle();

        calcState |= DfIntegrals::Spq;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }

    if (this->K_engine_ == K_ENGINE_RI_K) {
        // Sgd
        if ((calcState & DfIntegrals::Sgd) == 0) {
            if (this->m_bIsXCFitting == true) {
                this->outputStartTitle("Sgd");
                
                TlSymmetricMatrix Sgd(this->numOfAuxXC_);
                dfOverlap.getSgd(&Sgd);
                this->saveSgdMatrix(Sgd);
                
                this->outputEndTitle();
            }
            
            calcState |= DfIntegrals::Sgd;
            (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
            this->saveParam();
        }
    }

    if (this->J_engine_ == J_ENGINE_RI_J) {
        // Sab2
        if ((calcState & DfIntegrals::Sab2) == 0) {
            this->outputStartTitle("Sab2");
            
            TlSymmetricMatrix Sab2(this->m_nNumOfAux);
            if (this->isUseNewEngine_ == true) {
                this->logger(" use new engine.\n");
                dfOverlapX.getSab(&Sab2);
            } else {
                dfOverlap.getSab2(&Sab2);
            }
            this->saveSab2Matrix(Sab2);
            
            this->outputEndTitle();
            
            calcState |= DfIntegrals::Sab2;
            (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
            this->saveParam();
        }

        // Na
        if ((calcState & DfIntegrals::Na) == 0) {
            this->outputStartTitle("N_alpha");
            
            TlVector Na(this->m_nNumOfAux);
            if (this->isUseNewEngine_ == true) {
                this->logger(" use new engine.\n");
                dfOverlapX.getNalpha(&Na);
            } else {
                dfOverlap.getNa(&Na);
            }
            this->saveNalpha(Na);
            
            this->outputEndTitle();
            
            calcState |= DfIntegrals::Na;
            (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
            this->saveParam();
        }
    }
}


void DfIntegrals::createERIMatrix()
{
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

    if (this->J_engine_ == J_ENGINE_RI_J) {
        if ((calcState & DfIntegrals::Sab) == 0) {
            this->outputStartTitle("Sab");
            
            TlSymmetricMatrix Sab(this->m_nNumOfAux);
            DfEri dfEri(this->pPdfParam_);
            DfEriX dfEriX(this->pPdfParam_);
            
            if (this->isUseNewEngine_ == true) {
                this->logger(" use new engine.\n");
                dfEriX.getJab(&Sab);
            } else {
                dfEri.getSab(&Sab);
            }
            this->saveSabMatrix(Sab);
            
            this->outputEndTitle();
            
            calcState |= DfIntegrals::Sab;
            (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
            this->saveParam();
        }
    }
}


void DfIntegrals::createCholeskyVectors()
{
    if ((this->J_engine_ == J_ENGINE_CD) || (this->K_engine_ == K_ENGINE_CD)) {
        unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();
        if ((calcState & DfIntegrals::CD) == 0) { 
            this->outputStartTitle("Cholesky Vectors");
            DfCD *pDfCD = this->getDfCDObject();
            pDfCD->calcCholeskyVectors();
            
            delete pDfCD;
            pDfCD = NULL;
        }

        this->outputEndTitle();

        calcState |= DfIntegrals::CD;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }
}


void DfIntegrals::createXMatrix()
{
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

    if ((calcState & DfIntegrals::X) == 0) {
        this->outputStartTitle("X matrix");
        DfXMatrix* pDfXMatrix = this->getDfXMatrixObject();
        pDfXMatrix->main();

        delete pDfXMatrix;
        pDfXMatrix = NULL;

        this->outputEndTitle();

        calcState |= DfIntegrals::X;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }
}


void DfIntegrals::createInverseMatrixes()
{
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

    if (this->J_engine_ == J_ENGINE_RI_J) {
        if ((calcState & DfIntegrals::INV) == 0) {
            this->outputStartTitle("inverse matrix");
            DfInvMatrix* pDfInvMatrix = this->getDfInvMatrixObject();
            pDfInvMatrix->DfInvMain();
            
            delete pDfInvMatrix;
            pDfInvMatrix = NULL;
            
            this->outputEndTitle();
            
            calcState |= DfIntegrals::INV;
            (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
            this->saveParam();
        }
    }
}


void DfIntegrals::createGrids()
{
    unsigned int calcState = (*this->pPdfParam_)["control"]["integrals_state"].getUInt();

    if ((calcState & DfIntegrals::GRID) == 0) {
        this->outputStartTitle("Grid generation");
        
        DfGenerateGrid* pDfGenerateGrid = this->getDfGenerateGridObject();
        pDfGenerateGrid->dfGrdMain();
        
        delete pDfGenerateGrid;
        pDfGenerateGrid = NULL;
        
        this->outputEndTitle();
        
        calcState |= DfIntegrals::GRID;
        (*this->pPdfParam_)["control"]["integrals_state"].set(calcState);
        this->saveParam();
    }
}

void DfIntegrals::saveInvSquareVMatrix(const TlSymmetricMatrix& v)
{
    v.save("fl_Work/fl_Mtr_invSquareV.matrix");
}


void DfIntegrals::outputStartTitle(const std::string& stepName, const char lineChar)
{
    const std::string title = ">>>> " + stepName;
    this->log_.info(title);
}

void DfIntegrals::outputEndTitle(const std::string& stepName, const char lineChar)
{
    const std::string title = "<<<< " + stepName + " ";
    this->log_.info(title);
}
