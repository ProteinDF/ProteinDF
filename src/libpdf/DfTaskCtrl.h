// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

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
        Task4()
            : shellIndex1(0), shellIndex2(0),
              shellIndex3(0), shellIndex4(0) {
        }
            
        Task4(const Task4& rhs) 
            : shellIndex1(rhs.shellIndex1), shellIndex2(rhs.shellIndex2),
              shellIndex3(rhs.shellIndex3), shellIndex4(rhs.shellIndex4) {
        }

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
    
    virtual bool getQueue4(const TlOrbitalInfoObject& orbitalInfo_PQ,
                           const TlOrbitalInfoObject& orbitalInfo_RS,
                           const TlSparseSymmetricMatrix& schwarzTable_PQ,
                           const TlSparseSymmetricMatrix& schwarzTable_RS,
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

public:
    /// Schwartzのカットオフレポートを出力する
    virtual void cutoffReport();

protected:
    void clearCutoffStats(const int maxShellType);

    /// pre-screeningに関わるレポートを出力する
    virtual void prescreeningReport();

protected:
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);

    ShellArray selectShellArrayByDistribution(const ShellArray& inShellArray,
                                              const index_type companionShellIndex,
                                              const TlOrbitalInfoObject& orbitalInfo);

    // P v.s. R
    ShellPairArrayTable getShellPairArrayTable(const TlOrbitalInfoObject& orbitalInfo,
                                               const ShellArrayTable& shellArrayTable);
    ShellPairArrayTable getShellPairArrayTable(const TlOrbitalInfoObject& orbitalInfo1,
                                               const ShellArrayTable& shellArrayTable1,
                                               const TlOrbitalInfoObject& orbitalInfo2,
                                               const ShellArrayTable& shellArrayTable2);

    ShellPairArrayTable selectShellPairArrayTableByDensity(const ShellPairArrayTable& inShellPairArrayTable,
                                                           const TlOrbitalInfoObject& orbitalInfo);
    ShellPairArrayTable selectShellPairArrayTableByDensity(const ShellPairArrayTable& inShellPairArrayTable,
                                                           const TlOrbitalInfoObject& orbitalInfo1,
                                                           const TlOrbitalInfoObject& orbitalInfo2);
    
    TlSparseSymmetricMatrix makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo,
                                             DfEriEngine* pEngine);
    
    // int getShellQuartetType(const TlOrbitalInfoObject& orbitalInfo,
    //                         const int shellTypeP,
    //                         const int shellTypeQ,
    //                         const int shellTypeR,
    //                         const int shellTypeS);
    int getShellQuartetType(const TlOrbitalInfoObject& orbitalInfo_PQ,
                            const TlOrbitalInfoObject& orbitalInfo_RS,
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
    bool isAliveBySchwarzCutoff(const index_type shellIndexP,
                                const index_type shellIndexQ,
                                const index_type shellIndexR,
                                const index_type shellIndexS,
                                const int shellQuartetType,
                                const TlSparseSymmetricMatrix& schwarzTable_PQ,
                                const TlSparseSymmetricMatrix& schwarzTable_RS,
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
    std::size_t getTotalCalcAmount4(const TlOrbitalInfoObject& orbitalInfo1,
                                    const TlOrbitalInfoObject& orbitalInfo2,
                                    const ShellPairArrayTable& shellPairArrayTable,
                                    const DistributedCutoffTable& dct1,
                                    const DistributedCutoffTable& dct2);
    
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
