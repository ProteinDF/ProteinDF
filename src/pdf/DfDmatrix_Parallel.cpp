#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfDmatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlTime.h"
#include "TlFile.h"

DfDmatrix_Parallel::DfDmatrix_Parallel(TlSerializeData* pPdfParam)
    : DfDmatrix(pPdfParam)
{
}


DfDmatrix_Parallel::~DfDmatrix_Parallel()
{
}


void DfDmatrix_Parallel::main(const DfObject::RUN_TYPE runType)
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->main_SCALAPACK(runType);
    } else {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfDmatrix::main(runType);
        }
        rComm.barrier();
    }
#else
    {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        if (rComm.isMaster() == true) {
            DfDmatrix::main(runType);
        }
        rComm.barrier();
    }
#endif // HAVE_SCALAPACK
}

void DfDmatrix_Parallel::main_SCALAPACK(const DfObject::RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("build density matrix using ScaLAPACK.");

    // occupation
    TlVector currOcc;
    switch (this->orbitalCorrespondenceMethod_) {
    case OCM_OVERLAP:
        this->log_.info(" orbital correspondence method: MO-overlap");
        currOcc = this->getOccupationUsingOverlap<TlDistributeMatrix>(runType);
        if (rComm.isMaster() == true) {
            currOcc.save(this->getOccupationPath(runType));
        }
        break;

    case OCM_PROJECTION:
        this->log_.info(" orbital correspondence method: MO-projection");
        currOcc = this->getOccupationUsingProjection<TlDistributeMatrix, TlDistributeSymmetricMatrix>(runType);
        if (rComm.isMaster() == true) {
            currOcc.save(this->getOccupationPath(runType));
        }
        break;

    default:
        if (rComm.isMaster() == true) {
            this->log_.info(" orbital correspondence method: none");
            if (TlFile::isExist(this->getOccupationPath(runType)) == true) {
                currOcc.load(this->getOccupationPath(runType));
            } else {
                currOcc = this->createOccupation(runType);
                currOcc.save(this->getOccupationPath(runType));
            }
        }
        rComm.broadcast(currOcc);
        break;
    }
    
    rComm.barrier();
    this->generateDensityMatrix<TlDistributeMatrix, TlDistributeSymmetricMatrix>(runType, currOcc);
}


void DfDmatrix_Parallel::checkOccupation(const TlVector& prevOcc, const TlVector& currOcc)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDmatrix::checkOccupation(prevOcc, currOcc);
    }
}


void DfDmatrix_Parallel::printOccupation(const TlVector& occ)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDmatrix::printOccupation(occ);
    }
}


TlVector DfDmatrix_Parallel::getOccupation(const DfObject::RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlVector occ;
    if (rComm.isMaster() == true) {
        occ = DfDmatrix::getOccupation(runType);
    }
    rComm.broadcast(occ);

    return occ;
}


