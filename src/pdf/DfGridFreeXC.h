#ifndef DFGRIDFREEXC_H
#define DFGRIDFREEXC_H

#include "DfObject.h"
#include "DfOverlapEngine.h"
#include "DfTaskCtrl.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"
#include "TlSymmetricMatrix.h"
#include "TlSparseSymmetricMatrix.h"

class DfGridFreeXC : public DfObject 
{
public:
    DfGridFreeXC(TlSerializeData* pPdfParam);
    virtual ~DfGridFreeXC();

public:
    void buildFxc();

protected:
    static const int MAX_SHELL_TYPE;
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;

    struct ShellPair {
    public:
        ShellPair(index_type index1 =0, index_type index2 =0) : shellIndex1(index1), shellIndex2(index2) {
        }
        
    public:
        index_type shellIndex1;
        index_type shellIndex2;
    };
    typedef std::vector<ShellPair> ShellPairArray;
    typedef std::vector<ShellPairArray> ShellPairArrayTable;
    
protected:
    virtual void getM(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);

    TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo);
    void getM_part(const TlOrbitalInfoObject& orbitalInfo,
                   const std::vector<DfTaskCtrl::Task4>& taskList,
                   const TlMatrixObject& P, TlMatrixObject* pM);
    void storeM(const index_type shellIndexP, const int maxStepsP,
                const index_type shellIndexQ, const int maxStepsQ,
                const index_type shellIndexR, const int maxStepsR,
                const index_type shellIndexS, const int maxStepsS,
                const DfOverlapEngine& engine,
                const TlMatrixObject& P,
                TlMatrixObject* pM);

    virtual void createEngines();
    virtual void destroyEngines();
    virtual DfTaskCtrl* getDfTaskCtrlObject() const;
    virtual void finalize(TlSymmetricMatrix* pMtx);

    TlSymmetricMatrix get_F_lamda(const TlVector lamda);


    void getM_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);
    ShellPairArrayTable getShellPairArrayTable(const ShellArrayTable& shellArrayTable);


protected:
    DfOverlapEngine* pOvpEngines_;
};

#endif // DFGRIDFREEXC_H

