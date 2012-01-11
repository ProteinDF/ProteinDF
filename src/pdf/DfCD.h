#ifndef DFCD_H
#define DFCD_H

#include "DfObject.h"
#include "DfTaskCtrl.h"
#include "TlSparseSymmetricMatrix.h"

class TlOrbitalInfo;
class DfEriEngine;

class DfCD : public DfObject
{
public:
    DfCD(TlSerializeData* pPdfParam);
    virtual ~DfCD();

public:
    void makeSuperMatrix();

protected:
    void createEngines();
    void destroyEngines();
    
    void makeSuperMatrix_kernel(const TlOrbitalInfo& orbitalInfo,
                                const std::vector<DfTaskCtrl::Task4>& taskList,
                                TlSymmetricMatrix* pT);
    void storeT(const index_type shellIndexP, const int maxStepsP,
                const index_type shellIndexQ, const int maxStepsQ,
                const index_type shellIndexR, const int maxStepsR,
                const index_type shellIndexS, const int maxStepsS,
                const DfEriEngine& engine,
                TlSymmetricMatrix* pT);

    TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);

    std::size_t index(index_type p, index_type q) const;
    DfTaskCtrl* getDfTaskCtrlObject() const;
    void finalize(TlMatrix* pMtx);

protected:
    DfEriEngine* pEriEngines_;
    double cutoffThreshold_;
    double cutoffEpsilon3_;
};

#endif // DFCD_H
