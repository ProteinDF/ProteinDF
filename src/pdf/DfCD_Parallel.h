#ifndef DFCD_PARALLEL_H
#define DFCD_PARALLEL_H

#include <cstdlib>
#include "DfCD.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlTime.h"

class DfCD_Parallel : public DfCD {
public:
    DfCD_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfCD_Parallel();

public:
    // virtual void calcCholeskyVectors();
    void getJ_distributed(TlDistributeSymmetricMatrix *pJ);
    void getK_distributed(const RUN_TYPE runType,
                          TlDistributeSymmetricMatrix *pK);

protected:
    void makeSuperMatrix_distribute();
    TlDistributeSymmetricMatrix 
    getGMatrix_distribute(const TlOrbitalInfoObject& orbitalInfo, 
                          const TlSparseSymmetricMatrix& schwarzTable,
                          const index_type numOfItilde,
                          const PQ2I_Type& PQ2I);
    void makeL(const TlDistributeSymmetricMatrix& G);

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;
    virtual void finalize(TlSymmetricMatrix *pMat);
    virtual void finalize(TlSparseSymmetricMatrix *pMat);
    virtual void finalize_I2PQ(I2PQ_Type* pI2PQ);

    virtual void saveI2PQ(const I2PQ_Type& I2PQ);
    virtual I2PQ_Type getI2PQ();

    virtual void saveL(const TlMatrix& L);
    virtual TlMatrix getL();

    virtual TlSymmetricMatrix getPMatrix();

    virtual void divideCholeskyBasis(const index_type numOfCBs,
                                     index_type *pStart, index_type *pEnd);

    TlDistributeSymmetricMatrix 
    getCholeskyVector_distribute(const TlVector& L_col,
                                 const I2PQ_Type& I2PQ);
    
    // NEW ---------------------------------------------------------------------
protected:
    class RowVectorMatrix {
    public:
        RowVectorMatrix(index_type row = 1, index_type col = 1);
        ~RowVectorMatrix();
        
    public:
        void resize(index_type row, index_type col);
        index_type getNumOfRows() const {
            return this->globalRows_;
        };
        index_type getNumOfCols() const {
            return this->globalCols_;
        };

        void set(index_type row, index_type col, double value);
        
        TlVector getRowVector(index_type row) const;
        
        int getPEinChargeByRow(index_type row) const;
        
        TlMatrix getTlMatrix() const;
        TlDistributeMatrix getTlDistributeMatrix() const;

    private:
        struct RowVector {
        public:
            explicit RowVector(index_type r =0, index_type c =1) 
                : row(r), cols(c) {
            };
            
            bool operator<(const RowVector& rhs) const {
                return (this->row < rhs.row);
            };
            
        public:
            index_type row;
            TlVector cols;
        };
        
    private:
        index_type globalRows_;
        index_type globalCols_;
        std::vector<RowVector> data_;

        std::vector<int> row_PE_table_;
    };

protected:
    virtual void calcCholeskyVectors_onTheFly();
    virtual std::vector<double>
    getSuperMatrixElements(const index_type G_row,
                           const std::vector<index_type>& G_col_list,
                           const I2PQ_Type& I2PQ,
                           const TlSparseSymmetricMatrix& schwartzTable);
    void saveL(const RowVectorMatrix& L);

private:
    TlTime CD_all_time_;
    TlTime CD_ERI_time_;
    TlTime CD_calc_time_;
    TlTime CD_calc_d1_time_;
    TlTime CD_calc_d2_time_;
};

#endif // DFCD_PARALLEL_H

