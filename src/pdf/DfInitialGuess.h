#ifndef DFINITIALGUESS_H
#define DFINITIALGUESS_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <fstream>
#include "DfObject.h"
#include "CnError.h"
#include "TlUtils.h"

class DfDmatrix;
class TlVector;

/// 初期電子密度を求めるクラス
class DfInitialGuess : public DfObject {
public:
    DfInitialGuess(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuess();

    void exec();

protected:
    virtual void createRho();

    void createOccupation();
    virtual TlVector createOccupation(const RUN_TYPE runType);
    
    virtual void createInitialGuessUsingHuckel();
    virtual void createInitialGuessUsingCore();
    virtual void createInitialGuessUsingHarris();

    void createInitialGuessUsingLCAO();
    virtual void createInitialGuessUsingLCAO(const RUN_TYPE runType);

    std::vector<int> getLevel(const std::string& level);

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

    /// 密度行列を作成する
    ///
    /// DfDmatrixクラスを用いる
    void makeDensityMatrix();

    /// DfDmatrixクラスオブジェクトを作成して返す
    virtual DfDmatrix* getDfDmatrixObject(TlSerializeData* param);
    
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
        lcaoMatrix.resize(row_dimension, col_dimension);
        
        const int maxRows = row_dimension;
        const int maxCols = col_dimension;
        for (int i = 0; i < maxRows; ++i) {
            for (int j = 0; j < maxCols; ++j) {
                fi >> lcaoMatrix(i, j);
            }
        }
    }

    return lcaoMatrix;
}

template <typename MatrixType>
void DfInitialGuess::saveC0(const RUN_TYPE runType, const MatrixType& C0)
{
    DfObject::saveCMatrix(runType, 0, C0);
}

template <typename MatrixType>
void DfInitialGuess::buildCprime0(const RUN_TYPE runType, const MatrixType& C)
{
    MatrixType Xinv = DfObject::getXInvMatrix<MatrixType>();

    //  Xinv(RSFD) = Xinv(RSFD) * guess_lcao(RSFD)
    Xinv *= C;

    Xinv.save(this->getCprimeMatrixPath(runType, 0));
}


#endif // DFINITIALGUESS_H


