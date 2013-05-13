#ifndef DFOVERLAPX_H
#define DFOVERLAPX_H

#include <vector>
#include "DfObject.h"
#include "DfOverlapEngine.h"
#include "DfTaskCtrl.h"

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

class TlOrbitalInfoObject;

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
    virtual void getSgd(TlSymmetricMatrix* pSgd);
    void getNalpha(TlVector* pNalpha);
    void getForce(const TlSymmetricMatrix& W, TlMatrix* pForce);

    /// 変換行列を作成する
    void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                     const TlOrbitalInfoObject& orbInfo2,
                     TlMatrix* pS);

public:
    /// calc <pq gamma>
    virtual void get_pqg(const TlVector& myu, TlSymmetricMatrix* pF);
    virtual void get_pqg(const TlVector& myu, const TlVector& epsilon,
                         TlSymmetricMatrix* pF,
                         TlSymmetricMatrix* pE);

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
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo1,
                     const TlOrbitalInfoObject& orbitalInfo2,
                     TlMatrixObject* pMatrix);


    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pMatrix);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo1,
                          const TlOrbitalInfoObject& orbitalInfo2,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pMatrix);
    // <pq gamma>
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo_XC,
                     const TlVector& myu,
                     const TlOrbitalInfoObject& orbitalInfo,
                     TlMatrixObject* pF);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo_XC,
                          const TlVector& myu,
                          const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pMatrix);

    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo_XC,
                     const TlVector& myu,
                     const TlVector& epsilon,
                     const TlOrbitalInfoObject& orbitalInfo,
                     TlMatrixObject* pF,
                     TlMatrixObject* pE);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo_XC,
                          const TlVector& myu,
                          const TlVector& epsilon,
                          const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task2>& taskList,
                          TlMatrixObject* pF,
                          TlMatrixObject* pE);

    // <alpha>
    void calcOverlap(const TlOrbitalInfoObject& orbitalInfo,
                     TlVectorObject* pVector);
    void calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                          const std::vector<DfTaskCtrl::Task>& taskList,
                          TlVectorObject* pVector);
    
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);

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
