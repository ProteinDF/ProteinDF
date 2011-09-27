#ifndef DFPRESCF_H
#define DFPRESCF_H

#include <string>
#include <vector>

#include "DfObject.h"
#include "TlMatrix.h"
#include "CnError.h"

/// DfScfクラスの前処理を行うクラス
class DfPreScf : public DfObject {
public:
    DfPreScf(TlSerializeData* pPdfParam);
    virtual ~DfPreScf();

    void prepareGuess();

protected:
    virtual void createInitialGuessUsingLCAO(const RUN_TYPE runType);
    virtual void createOccupation(const RUN_TYPE runType);

    std::vector<int> getLevel(std::string sLevel);

    /// LCAO行列を取得する
    template <typename MatrixType>
    MatrixType getLCAO(const RUN_TYPE runType);

    /// 行列C0を保存する
    template <typename MatrixType>
    void saveC0(const RUN_TYPE runType, const MatrixType& C0);

    /// 占有軌道情報を取得する
    virtual TlVector getOccupation(const RUN_TYPE runType);

    /// 占有軌道情報を保存する
    void saveOccupation(const RUN_TYPE runType, const TlVector& rOccupation);

    /// Cprime0 を作成し、保存する
    template <typename MatrixType>
    void buildCprime0(const RUN_TYPE runType, const MatrixType& C);
};


template <typename MatrixType>
MatrixType DfPreScf::getLCAO(const RUN_TYPE runType)
{
    MatrixType lcaoMatrix(this->m_nNumOfAOs, this->m_nNumOfMOs);

    {
        std::ifstream fi;
        const std::string sFile = std::string("./guess.lcao.") + this->m_sRunTypeSuffix[runType];
        fi.open(sFile.c_str(), std::ios::in);
        if (fi.rdstate()) {
            std::cerr << "cannot open file " << ("./guess.lcao." + sFile) << std::endl;
            CnErr.abort();
        }

        std::string dummy_line;
        fi >> dummy_line;

        int row_dimension, col_dimension;
        fi >> row_dimension >> col_dimension;
        if (row_dimension != this->m_nNumOfAOs) {
            CnErr.abort("DfPreScf", "", "prepare_occupation_and_or_mo", "inputted guess lcao has illegal dimension");
        }

        if (this->m_nNumOfMOs < col_dimension) {
            this->logger("The number of column dimension in inputed LCAO is larger than independent basis.\n");
            this->logger("Excess elements are discarded.\n");
        }

        const int maxRows = this->m_nNumOfAOs;
        const int maxCols = this->m_nNumOfMOs;
        const int excessCols = col_dimension - this->m_nNumOfMOs;
        for (int i = 0; i < maxRows; ++i) {
            for (int j = 0; j < maxCols; ++j) {
                fi >> lcaoMatrix(i, j);
            }

            for (int j = 0; j < excessCols; ++j) {
                double spoil;
                fi >> spoil;
            }
        }
    }

    return lcaoMatrix;
}


// C0 の保存
template <typename MatrixType>
void DfPreScf::saveC0(const RUN_TYPE runType, const MatrixType& C0)
{
    C0.save(this->getCMatrixPath(runType, 0));
}


template <typename MatrixType>
void DfPreScf::buildCprime0(const RUN_TYPE runType, const MatrixType& C)
{
    MatrixType Xinv(this->m_nNumOfAOs, this->m_nNumOfMOs);
    Xinv.load(this->getInvXMatrixPath());

    //  Xinv(RSFD) = Xinv(RSFD) * guess_lcao(RSFD)
    Xinv *= C;

    Xinv.save(this->getCprimeMatrixPath(runType, 0));
}


#endif // DFPRESCF_H
