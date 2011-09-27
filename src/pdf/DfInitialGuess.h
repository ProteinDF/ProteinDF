#ifndef DFINITIALGUESS_H
#define DFINITIALGUESS_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <fstream>
#include "DfObject.h"
#include "CnError.h"
#include "TlUtils.h"
class TlVector;

/// 初期電子密度を求めるクラス
class DfInitialGuess : public DfObject {
public:
    DfInitialGuess(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuess();

    void exec();

protected:
    virtual void createRho();
    virtual void saveRho1(RUN_TYPE runType);

    void createOccupation();
    virtual void createOccupation(const RUN_TYPE runType);

    
    virtual void createInitialGuessUsingHuckel();
    virtual void createInitialGuessUsingCore();
    virtual void createInitialGuessUsingHarris();

    void createInitialGuessUsingLCAO();
    virtual void createInitialGuessUsingLCAO(const RUN_TYPE runType);

    std::vector<int> getLevel(std::string sLevel);

    /// 占有軌道情報を取得する
    virtual TlVector getOccupation(const RUN_TYPE runType);

    /// 占有軌道情報を保存する
    virtual void saveOccupation(const RUN_TYPE runType, const TlVector& rOccupation);

protected:
    /// LCAO行列を取得する
    template <typename MatrixType>
    MatrixType getLCAO(const RUN_TYPE runType);

    /// 行列C0を保存する
    template <typename MatrixType>
    void saveC0(const RUN_TYPE runType, const MatrixType& C0);

    /// Cprime0 を作成し、保存する
    template <typename MatrixType>
    void buildCprime0(const RUN_TYPE runType, const MatrixType& C);
};


template <typename MatrixType>
MatrixType DfInitialGuess::getLCAO(const RUN_TYPE runType)
{
    //MatrixType lcaoMatrix(this->m_nNumOfAOs, this->m_nNumOfMOs);
    MatrixType lcaoMatrix;
    
    {
        std::ifstream fi;
        const std::string sFile = std::string("./guess.lcao.") + this->m_sRunTypeSuffix[runType];
        fi.open(sFile.c_str(), std::ios::in);
        if (fi.rdstate()) {
            CnErr.abort(TlUtils::format("cannot open file %s.\n", sFile.c_str()));
        }

        std::string dummy_line;
        fi >> dummy_line;

        int row_dimension, col_dimension;
        fi >> row_dimension >> col_dimension;
        if (row_dimension != this->m_nNumOfAOs) {
            CnErr.abort("DfPreScf", "", "prepare_occupation_and_or_mo", "inputted guess lcao has illegal dimension");
        }

//         if (this->m_nNumOfMOs < col_dimension) {
//             this->logger("The number of column dimension in inputed LCAO is larger than independent basis.\n");
//             this->logger("Excess elements are discarded.\n");
//         }
        lcaoMatrix.resize(row_dimension, col_dimension);
        
        //const int maxRows = this->m_nNumOfAOs;
        //const int maxCols = this->m_nNumOfMOs;
        //const int excessCols = col_dimension - this->m_nNumOfMOs;
        const int maxRows = row_dimension;
        const int maxCols = col_dimension;
        for (int i = 0; i < maxRows; ++i) {
            for (int j = 0; j < maxCols; ++j) {
                fi >> lcaoMatrix(i, j);
            }

//             for (int j = 0; j < excessCols; ++j) {
//                 double spoil;
//                 fi >> spoil;
//             }
        }
    }

    return lcaoMatrix;
}


template <typename MatrixType>
void DfInitialGuess::saveC0(const RUN_TYPE runType, const MatrixType& C0)
{
    C0.save(this->getCMatrixPath(runType, 0));
}

template <typename MatrixType>
void DfInitialGuess::buildCprime0(const RUN_TYPE runType, const MatrixType& C)
{
    MatrixType Xinv;
    Xinv.load(this->getInvXMatrixPath());

    //  Xinv(RSFD) = Xinv(RSFD) * guess_lcao(RSFD)
    Xinv *= C;

    Xinv.save(this->getCprimeMatrixPath(runType, 0));
}


#endif // DFINITIALGUESS_H


