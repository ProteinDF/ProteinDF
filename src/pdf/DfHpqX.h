#ifndef DFHPQX_H
#define DFHPQX_H

#include "DfObject.h"
#include "DfHpqEngine.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlOrbitalInfo.h"
#include "DfTaskCtrl.h"

class DfHpqX : public DfObject {
public:
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;

public:
    DfHpqX(TlSerializeData* pPdfParam);
    virtual ~DfHpqX();

public:
    void getHpq(TlSymmetricMatrix* pHpq, TlSymmetricMatrix* pHpq2);

    /// Hpq由来の力成分を求める
    ///
    /// @param [in] P 密度行列
    /// @param [out] pForce 原子数(dummy chargeを含む)×3(x, y, z成分)の行列
    void getForce(const TlSymmetricMatrix& P,
                  TlMatrix* pForce);

protected:
    /// DfHpqEngineオブジェクトを作成する
    ///
    /// OpenMPスレッド数のオブジェクトを作成する。
    /// 富士通コンパイラではコンストラクタ中で
    /// オブジェクトを作成できないため。
    void createEngines();

    /// DfHpqEngineオブジェクトを破棄する
    void destroyEngines();

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlSymmetricMatrix* pHpq, TlSymmetricMatrix* pHpq2);
    
protected:
    void getHpq_part(const TlOrbitalInfoObject& orbitalInfo,
                     const std::vector<DfTaskCtrl::Task2>& taskList,
                     const std::vector<TlAtom>& Cs,
                     const std::vector<TlAtom>& Xs,
                     TlMatrixObject* pHpq,
                     TlMatrixObject* pHpq2);
    void getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                           const int shellTypeP, const int shellTypeQ,
                           const index_type shellIndexP,
                           const ShellArray& shellArrayQ,
                           const TlMatrixObject& P,
                           TlMatrix* pForce);
    
protected:
    void makeShellArrayTable();
    DfHpqEngine::PGTOs getPGTOs(const index_type shellIndex);

protected:
    enum {
        X = 0,
        Y = 1,
        Z = 2
    };

protected:
    static const int MAX_SHELL_TYPE;

    DfHpqEngine* pEngines_;
    TlOrbitalInfo orbitalInfo_;
    ShellArrayTable shellArrayTable_;
};

#endif // DFHPQX_H
