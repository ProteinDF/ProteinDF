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
    struct IndexPair2 {
    public:
        explicit IndexPair2(index_type i1 =0, index_type i2 =0) : index1_(i1), index2_(i2) {
            if (this->index1_ > this->index2_) {
                std::swap(this->index1_, this->index2_);
            }
        }
        
        std::size_t index() const {
            // 'U' format
            assert(this->index1_ <= this->index2_);
            return this->index1_ + this->index2_ * (this->index2_ +1) / 2;
        }
        
        bool operator<(const IndexPair2& rhs) const {
            return (this->index() < rhs.index());
        }
        
        index_type index1() const {
            return this->index1_;
        }
        
        index_type index2() const {
            return this->index2_;
        }
        
    private:
        index_type index1_;
        index_type index2_;
    };
    
    struct IndexPair4 {
    public:
        explicit IndexPair4(index_type i1 =0, index_type i2 =0,
                            index_type i3 =0, index_type i4 =0) 
            : indexPair1_(i1, i2), indexPair2_(i3, i4) {
            if (this->indexPair1_.index() > this->indexPair2_.index()) {
                std::swap(this->indexPair1_, this->indexPair2_);
            }
        }
        
        std::size_t index() const {
            std::size_t pair_index1 = this->indexPair1_.index();
            std::size_t pair_index2 = this->indexPair2_.index();
            assert(pair_index1 <= pair_index2);
            return pair_index1 + pair_index2 * (pair_index2 +1) / 2;
        }
        
        bool operator<(const IndexPair4& rhs) const {
            return (this->index() < rhs.index());
        }
        
        bool operator==(const IndexPair4& rhs) const {
            return (this->index() == rhs.index());
        }
        
        index_type index1() const {
            return this->indexPair1_.index1();
        }
        
        index_type index2() const {
            return this->indexPair1_.index2();
        }
        
        index_type index3() const {
            return this->indexPair2_.index1();
        }
        
        index_type index4() const {
            return this->indexPair2_.index2();
        }
    
    private:
        IndexPair2 indexPair1_;
        IndexPair2 indexPair2_;
    };

    typedef std::vector<IndexPair2> PQ_PairArray;
    
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

    void get_F_lamda(const TlVector lamda,
                     TlSymmetricMatrix* pF_lamda,
                     TlSymmetricMatrix* pE_lamda);

    void getM_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM);
    ShellArrayTable makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo);
    ShellPairArrayTable getShellPairArrayTable(const ShellArrayTable& shellArrayTable);

public:
    virtual void calcCholeskyVectors_onTheFly();

protected:
    void calcDiagonals(TlSparseSymmetricMatrix *pSchwartzTable,
                       PQ_PairArray *pI2PQ,
                       TlVector *pDiagonals);
    void calcDiagonals_kernel(const std::vector<DfTaskCtrl::Task2>& taskList,
                              TlSparseSymmetricMatrix *pSchwartzTable,
                              TlSparseSymmetricMatrix *pDiagonalMat,
                              PQ_PairArray *pI2PQ);
    void saveI2PQ(const PQ_PairArray& I2PQ);
    void saveL(const TlMatrix& L);

    std::vector<double>
    getSuperMatrixElements(const index_type G_row,
                           const std::vector<index_type>& G_col_list,
                           const PQ_PairArray& I2PQ,
                           const TlSparseSymmetricMatrix& schwartzTable);

    std::vector<DfGridFreeXC::IndexPair4> 
    getCalcList(const index_type G_row,
                const std::vector<index_type>& G_col_list,
                const PQ_PairArray& I2PQ);

    void calcElements(const std::vector<IndexPair4>& calcList,
                      const TlSparseSymmetricMatrix& schwartzTable);

    bool isAliveBySchwartzCutoff(const index_type shellIndexP,
                                 const index_type shellIndexQ,
                                 const index_type shellIndexR,
                                 const index_type shellIndexS,
                                 const int shellQuartetType,
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const double threshold);
    void initializeCutoffStats();
    void schwartzCutoffReport();
    mutable std::vector<unsigned long> cutoffAll_schwartz_;
    mutable std::vector<unsigned long> cutoffAlive_schwartz_;

    std::vector<double>
    setElements(const index_type G_row,
                const std::vector<index_type> G_col_list,
                const PQ_PairArray& I2PQ);

    void getM_byCD(TlSymmetricMatrix* pM);
    TlSymmetricMatrix getPMatrix();
    TlMatrix getL();
    PQ_PairArray getI2PQ();
    void divideCholeskyBasis(const index_type numOfCBs,
                             index_type *pStart, index_type *pEnd);
    TlSymmetricMatrix getCholeskyVector(const TlVector& L_col,
                                        const PQ_PairArray& I2PQ);

protected:
    DfOverlapEngine* pOvpEngines_;
    TlOrbitalInfo orbitalInfo_;

    index_type numOfPQs_;
    double tau_;
    double epsilon_;

    /// キャッシュ
    typedef std::map<IndexPair4, std::vector<double> > ElementsCacheType;
    ElementsCacheType elements_cache_;
};

#endif // DFGRIDFREEXC_H
