#ifndef DFTRANSATOB_H
#define DFTRANSATOB_H

#include <string>
#include <cassert>

#include "DfObject.h"
#include "TlUtils.h"

// function : "C = X * C'"  (with Bender scheme)
//            RKS and UKS are suported
//
// input    : "X"
//            "C'"
//
// output   : "C"
//

/** X 行列とC’行列を用いて、規格直交基底からCGTO基底に変換する変換行列を求め、ファイルに出力するクラス
 */
class DfTransatob : public DfObject {
public:
    DfTransatob(TlSerializeData* pPdfParam);
    virtual ~DfTransatob();

    virtual void DfTrsatobMain();
    virtual void DfTrsatobQclo(const std::string& fragname, int norbcut);

protected:
    template<typename MatrixType>
    void main(const RUN_TYPE runtype, const std::string& fragname ="", bool bPdfQcloMode =false);
};

// template
template<typename MatrixType>
void DfTransatob::main(const RUN_TYPE runType, const std::string& fragname, bool bPdfQcloMode)
{
    // "read X matrix"
    MatrixType X;
    {
        X = this->getXMatrix<MatrixType>();

        if (X.getNumOfRows() != this->m_nNumOfAOs ||
                X.getNumOfCols() != this->m_nNumOfMOs) {
            this->logger(TlUtils::format("rowDim of fl_Mtr_X.matrix = %d\n", X.getNumOfRows()));
            this->logger(TlUtils::format("colDim of fl_Mtr_X.matrix = %d\n", X.getNumOfCols()));
            this->logger(TlUtils::format("number_ao_basis = %d", this->m_nNumOfAOs));
            this->logger(TlUtils::format("number_independant_basis = %d\n", this->m_nNumOfMOs));
            this->logger("DfTransatob dimension is not consistency, but continue\n");
        }
    }

    // "read C' matrix"
    MatrixType Cprime;
    {
        const std::string fragment = "";
        Cprime = this->getCprimeMatrix<MatrixType>(runType, this->m_nIteration, fragment);

        if (Cprime.getNumOfRows() != this->m_nNumOfMOs ||
                Cprime.getNumOfCols() != this->m_nNumOfMOs) {
            this->logger(TlUtils::format("row of C' = %d\n", Cprime.getNumOfRows()));
            this->logger(TlUtils::format("col of C' = %d\n", Cprime.getNumOfRows()));
            this->logger(TlUtils::format("number_mo_basis = %d\n", this->m_nNumOfMOs));
            this->logger("DfTransatob dimension is not consistency, but continue\n");
        }
    }

    // calculate "C = X * C'"
    const MatrixType C = X * Cprime;
    this->saveCMatrix(runType, this->m_nIteration, C);
}


#endif // DFTRANSATOB_H
