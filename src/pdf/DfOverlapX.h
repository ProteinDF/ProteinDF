#ifndef DFOVERLAPX_H
#define DFOVERLAPX_H

#include <vector>
#include "DfObject.h"
#include "DfOverlapEngine.h"
#include "DfTaskCtrl.h"

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"

class DfOverlapX : public DfObject {
public:
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;
    
public:
    DfOverlapX(TlSerializeData* pPdfParam);
    virtual ~DfOverlapX();

public:
    void getSpq(TlSymmetricMatrix* pSpq);
    void getSab(TlSymmetricMatrix* pSab);
    void getForce(const TlSymmetricMatrix& W, TlMatrix* pForce);

protected:
    static const int MAX_SHELL_TYPE;
    
    enum {
        X = 0,
        Y = 1,
        Z = 2
    };

protected:
    /// DfOverlapEngineオブジェクトを作成する
    ///
    /// OpenMPスレッド数のオブジェクトを作成する。
    /// 富士通コンパイラではコンストラクタ中で
    /// オブジェクトを作成できないため。
    void createEngines();

    /// DfOverlapEngineオブジェクトを破棄する
    void destroyEngines();

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlMatrix* pMtx);
    virtual void finalize(TlSymmetricMatrix* pMtx);
    virtual void finalize(TlVector* pVct);
    
protected:
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo,
                     TlMatrixObject* pMatrix);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pMatrix);
    
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);
    DfOverlapEngine::PGTOs getPGTOs(const TlOrbitalInfoObject& orbitalInfo,
                                    const int shellIndex);

    void getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                           const int shellTypeP, const int shellTypeQ,
                           const index_type shellIndexP,
                           const ShellArray& shellArrayQ,
                           const TlSymmetricMatrix& W,
                           TlMatrix* pForce);

protected:
    DfOverlapEngine* pEngines_;
};

#endif // DFOVERLAPX_H
