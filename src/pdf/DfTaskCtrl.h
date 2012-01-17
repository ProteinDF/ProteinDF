#ifndef DFTASKCTRL_H
#define DFTASKCTRL_H

#include <vector>
#include "DfObject.h"
#include "DfEriEngine.h"
#include "TlOrbitalInfoObject.h"
#include "TlMatrixObject.h"
#include "TlSparseSymmetricMatrix.h"

class DfTaskCtrl : public DfObject {
public:
    struct Task2 {
    public:
        index_type shellIndex1;
        index_type shellIndex2;
    };

    struct Task4 {
    public:
        index_type shellIndex1;
        index_type shellIndex2;
        index_type shellIndex3;
        index_type shellIndex4;
    };

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
    
public:
    typedef std::vector<index_type> ShellArray;
    typedef std::vector<ShellArray> ShellArrayTable;

protected:
    typedef std::map<index_type, std::vector<std::vector<index_type> > > DistributedCutoffTable;
    
public:
    DfTaskCtrl(TlSerializeData* pPdfParam);
    virtual ~DfTaskCtrl();

public:
    void setCutoffThreshold(const double value);
    double getCutoffThreshold() const;
    
    virtual void cutoffReport();
    virtual bool getQueue(const TlOrbitalInfoObject& orbitalInfo,
                          const bool isCutoffByDistibution,
                          const int maxGrainSize,
                          std::vector<Task2>* pTask,
                          bool initialize = false);

    virtual bool getQueue4(const TlOrbitalInfoObject& orbitalInfo,
                            const TlSparseSymmetricMatrix& schwarzTable,
                            const int maxGrainSize,
                            std::vector<Task4>* pTaskList,
                            bool initialize = false);

    virtual bool getQueue4_K(const TlOrbitalInfoObject& orbitalInfo,
                             const TlSparseSymmetricMatrix& schwarzTable,
                             const TlMatrixObject& P,
                             const std::vector<TlMatrixObject::index_type>& rowIndexes,
                             const std::vector<TlMatrixObject::index_type>& colIndexes,
                             const int maxGrainSize,
                             std::vector<Task4>* pTaskList,
                             bool initialize = false);
    virtual bool getQueue4_K0(const TlOrbitalInfoObject& orbitalInfo,
                             const TlSparseSymmetricMatrix& schwarzTable,
                             const TlMatrixObject& P,
                             const std::vector<TlMatrixObject::index_type>& rowIndexes,
                             const std::vector<TlMatrixObject::index_type>& colIndexes,
                             const int maxGrainSize,
                             std::vector<Task4>* pTaskList,
                             bool initialize = false);

    virtual bool getQueue_Force4(const TlOrbitalInfoObject& orbitalInfo,
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const int maxGrainSize,
                                 std::vector<Task4>* pTaskList,
                                 bool initialize = false);
    
protected:
    void clearCutoffStats(const TlOrbitalInfoObject& orbitalInfo);

    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);

    ShellArray selectShellArrayByDistribution(const ShellArray& inShellArray,
                                              const index_type companionShellIndex,
                                              const TlOrbitalInfoObject& orbitalInfo);

    ShellPairArrayTable getShellPairArrayTable(const TlOrbitalInfoObject& orbitalInfo,
                                               const ShellArrayTable& shellArrayTable);

    ShellPairArrayTable selectShellPairArrayTableByDensity(const ShellPairArrayTable& inShellPairArrayTable,
                                                           const TlOrbitalInfoObject& orbitalInfo);
    
    TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo,
                                             DfEriEngine* pEngine);
    
    int getShellQuartetType(const TlOrbitalInfoObject& orbitalInfo,
                            const int shellTypeP,
                            const int shellTypeQ,
                            const int shellTypeR,
                            const int shellTypeS);
    bool isAliveBySchwarzCutoff(const index_type shellIndexP,
                                const index_type shellIndexQ,
                                const index_type shellIndexR,
                                const index_type shellIndexS,
                                const int shellQuartetType,
                                const TlSparseSymmetricMatrix& schwarzTable,
                                const double threshold);

    double getMaxValue(const TlMatrixObject& P,
                       const TlMatrixObject::index_type row, const TlMatrixObject::index_type d_row, 
                       const TlMatrixObject::index_type col, const TlMatrixObject::index_type d_col);

    /// distributed pair用のカットオフテーブルを作成する
    ///
    /// キーの軌道番号に対して有効な軌道番号の配列を返す
    DistributedCutoffTable makeDistributedCutoffTable(const TlOrbitalInfoObject& orbitalInfo);

    ShellArray selectShellArrayByDistribution(const ShellArray& inShellArray,
                                              const index_type companionShellIndex);
    
protected:
    int maxShellType_;
    
    double cutoffThreshold_;
    
    double lengthScaleParameter_;

    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.32
    double cutoffEpsilon1_;

    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.32
    double cutoffEpsilon2_;

    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.33
    double cutoffEpsilon3_;

    /// カットオフ情報保持変数
    mutable std::vector<unsigned long> cutoffAll_E1_;
    mutable std::vector<unsigned long> cutoffAlive_E1_;
    mutable std::vector<unsigned long> cutoffAll_E2_;
    mutable std::vector<unsigned long> cutoffAlive_E2_;
    mutable std::vector<unsigned long> cutoffAll_schwarz_;
    mutable std::vector<unsigned long> cutoffAlive_schwarz_;
};

#endif // DFTASKCTRL_H
