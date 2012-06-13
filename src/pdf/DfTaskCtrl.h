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
    struct Task {
    public:
        index_type shellIndex1;
    };

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
    void setCutoffEpsilon_density(const double value);
    double getCutoffEpsilon_density() const;
    void setCutoffEpsilon_distribution(const double value);
    double getCutoffEpsilon_distribution() const;
    // void setCutoffEpsilon_primitive(const double value);
    // double getCutoffEpsilon_primitive() const;

public:
    virtual bool getQueue(const TlOrbitalInfoObject& orbitalInfo,
                          const int maxGrainSize,
                          std::vector<Task>* pTask,
                          bool initialize = false);

    /// 2つのインデックスのタスクを返す(軌道情報が同じ場合)
    virtual bool getQueue2(const TlOrbitalInfoObject& orbitalInfo,
                           const bool isCutoffByDistibution,
                           const int maxGrainSize,
                           std::vector<Task2>* pTask,
                           bool initialize = false);

    /// 2つのインデックスのタスクを返す(軌道情報が異なる場合)
    virtual bool getQueue2(const TlOrbitalInfoObject& orbitalInfo1,
                           const TlOrbitalInfoObject& orbitalInfo2,
                           const bool isCutoffByDistibution,
                           const int maxGrainSize,
                           std::vector<Task2>* pTask,
                           bool initialize = false);

    virtual bool getQueue4(const TlOrbitalInfoObject& orbitalInfo,
                           const TlSparseSymmetricMatrix& schwarzTable,
                           const int maxGrainSize,
                           std::vector<Task4>* pTaskList,
                           bool initialize = false);
    
    // local matrix用?
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
    void clearCutoffStats(const int maxShellType);

    /// pre-screeningに関わるレポートを出力する
    virtual void prescreeningReport();

    /// Schwartzのカットオフレポートを出力する
    virtual void cutoffReport();

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

    /// distributed pair用のカットオフテーブルを作成する(軌道情報が異なる場合)
    ///
    /// キーの軌道番号に対して有効な軌道番号の配列を返す
    DistributedCutoffTable makeDistributedCutoffTable(const TlOrbitalInfoObject& orbitalInfo1,
                                                      const TlOrbitalInfoObject& orbitalInfo2);

    ShellArray selectShellArrayByDistribution(const ShellArray& inShellArray,
                                              const index_type companionShellIndex);

    std::size_t getTotalCalcAmount2(const TlOrbitalInfoObject& orbitalInfo,
                                     const ShellArrayTable& shellArrayTable);
    std::size_t getTotalCalcAmount2(const TlOrbitalInfoObject& orbitalInfo1,
                                    const TlOrbitalInfoObject& orbitalInfo2,
                                    const ShellArrayTable& shellArrayTable1,
                                    const ShellArrayTable& shellArrayTable2);
    std::size_t getTotalCalcAmount2(const TlOrbitalInfoObject& orbitalInfo,
                                    const ShellArrayTable& shellArrayTable,
                                    const DistributedCutoffTable& dct);
    std::size_t getTotalCalcAmount2(const TlOrbitalInfoObject& orbitalInfo1,
                                    const TlOrbitalInfoObject& orbitalInfo2,
                                    const ShellArrayTable& shellArrayTable1,
                                    const ShellArrayTable& shellArrayTable2,
                                    const DistributedCutoffTable& dct);
    std::size_t getTotalCalcAmount4(const TlOrbitalInfoObject& orbitalInfo,
                                    const ShellPairArrayTable& shellPairArrayTable,
                                    const DistributedCutoffTable& dct);
    
protected:
    int maxShellType_;
    
    double lengthScaleParameter_;

    /// カットオフ用閾値
    /// for Schwarz 
    double cutoffThreshold_;

    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.20
    double cutoffEpsilon_density_;

    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.32
    double cutoffEpsilon_distribution_;

    /// カットオフ用閾値
    /// J. Chem. Phys.,105,2726 (1996) : eq.33
    // double cutoffEpsilon_primitive_;

    /// カットオフ情報保持変数
    mutable std::vector<unsigned long> cutoffAll_density_;
    mutable std::vector<unsigned long> cutoffAlive_density_;
    mutable std::vector<unsigned long> cutoffAll_distribution_;
    mutable std::vector<unsigned long> cutoffAlive_distribution_;
    mutable std::vector<unsigned long> cutoffAll_schwarz_;
    mutable std::vector<unsigned long> cutoffAlive_schwarz_;
};

#endif // DFTASKCTRL_H
