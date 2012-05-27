#include "DfHpqX_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"

DfHpqX_Parallel::DfHpqX_Parallel(TlSerializeData* pPdfParam)
    : DfHpqX(pPdfParam) {
}


DfHpqX_Parallel::~DfHpqX_Parallel()
{
}


void DfHpqX_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfHpqX::logger(str);
    }
}


void DfHpqX_Parallel::getHpqD(TlDistributeSymmetricMatrix* pHpq,
                              TlDistributeSymmetricMatrix* pHpq2)
{
    const int numOfAOs = this->m_nNumOfAOs;
    
    // make coordinates
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();
    const int numOfDummyAtoms = flGeom.getDummyatom();
    const int numOfRealAtoms = numOfAtoms - numOfDummyAtoms;
    std::vector<TlAtom> Cs(numOfRealAtoms);
    std::vector<TlAtom> Xs(numOfDummyAtoms);
    std::size_t realAtomIndex = 0;
    std::size_t dummyAtomIndex = 0;
    for (int i = 0; i < numOfAtoms; ++i) {
        const std::string atomName = flGeom.getAtom(i);
        const TlPosition p = flGeom.getCoordinate(i);
        const double charge = flGeom.getCharge(i);
        const TlAtom atom(atomName, p, charge);
        if (atomName == "X") {
            Xs[dummyAtomIndex] = atom;
            ++dummyAtomIndex;
        } else {
            Cs[realAtomIndex] = atom;
            ++realAtomIndex;
        }
    }
    assert(realAtomIndex == Cs.size());
    assert(dummyAtomIndex == Xs.size());
    
    DfHpqEngine engine;

    pHpq->resize(numOfAOs);
    pHpq2->resize(numOfAOs);
    TlSparseSymmetricMatrix tmpHpq(numOfAOs);
    TlSparseSymmetricMatrix tmpHpq2(numOfAOs);

    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue2(this->orbitalInfo_,
                                        true,
                                        this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getHpq_part(this->orbitalInfo_,
                          taskList,
                          Cs, Xs,
                          &tmpHpq, &tmpHpq2);
        
        hasTask = pTaskCtrl->getQueue2(this->orbitalInfo_,
                                       true,
                                       this->grainSize_, &taskList);
    } 

    if (this->chargeExtrapolateNumber_ > 0) {
        tmpHpq2 /= static_cast<double>(this->chargeExtrapolateNumber_);
    }

    this->loggerTime(" finalize");
    pHpq->mergeSparseMatrix(tmpHpq);
    pHpq2->mergeSparseMatrix(tmpHpq2);

    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();
}


DfTaskCtrl* DfHpqX_Parallel::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl_Parallel(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfHpqX_Parallel::finalize(TlSymmetricMatrix* pHpq, TlSymmetricMatrix* pHpq2)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pHpq);
    rComm.allReduce_SUM(*pHpq2);
}


