#ifndef DFDMATRIX_H
#define DFDMATRIX_H

#include <string>
#include <vector>

#include "DfObject.h"
#include "TlSymmetricMatrix.h"
#include "CnError.h"
#include "TlTime.h"

class TlVector;
class TlMatrix;

/// 分子軌道の占有数の決定と密度行列の作成を行うクラス
/// 
/// 収束させるための処理として、軌道の重なりの方法、射影演算子法を用いる
class DfDmatrix : public DfObject {
public:
    DfDmatrix(TlSerializeData* pPdfParam);
    virtual ~DfDmatrix();

public:
    void DfDmatrixMain();

protected:
    virtual void main(DfObject::RUN_TYPE runType);

    virtual TlVector getOccupation(DfObject::RUN_TYPE runType);

    template<typename MatrixType>
    TlVector getOccupationUsingOverlap(DfObject::RUN_TYPE runType);

    template<typename MatrixType, typename SymmetricMatrixType>
    TlVector getOccupationUsingProjection(DfObject::RUN_TYPE runType);

    virtual void checkOccupation(const TlVector& prevOcc, const TlVector& currOcc);
    virtual void printOccupation(const TlVector& occ);

    template<typename MatrixType, typename SymmetricMatrixType>
    void generateDensityMatrix(DfObject::RUN_TYPE runType, const TlVector& currOcc);

    template<typename MatrixType, typename SymmetricMatrixType>
    SymmetricMatrixType calcDensMatrix(const MatrixType& inputC, const TlVector& f, double base);

    void printTwoVectors(const TlVector& a, const TlVector& b, const std::string& title, int pnumcol);

    // for simple
    TlVector createOccupation(DfObject::RUN_TYPE runType);
    std::vector<int> getLevel(std::string sLevel);

protected:
    enum ORBITAL_CORRESPONDENCE_METHOD {
        OCM_NONE,
        OCM_OVERLAP,
        OCM_PROJECTION
    };

protected:
    ORBITAL_CORRESPONDENCE_METHOD orbitalCorrespondenceMethod_;

    
//     /// 軌道関連づけを行う(true)かどうか
//     bool isOrbitalCorrespondence_;

//     /// 1 回目のiteration から上記の方法を用いるかのキーワード
//     std::string orbital_overlap_first;

//     /// 軌道重なりの方法もしくは軌道射影法の指定を行うキーワード
//     std::string orbital_overlap_method;

//     /// 軌道重なりの方法をiteration 何回目まで行うかを指定するキーワード
//     int mo_overlap_iter;
};

//
// templates
//

// 軌道の重なりの対応のルーチン ====================================
template<typename MatrixType>
TlVector DfDmatrix::getOccupationUsingOverlap(DfObject::RUN_TYPE runType)
{
    const int nNumOfMOs = this->m_nNumOfMOs;

    this->log_.info(" MO overlap method is started.\n");

    this->log_.info(" load previous C' matrix");
    MatrixType prevCprime(nNumOfMOs, nNumOfMOs);
    {
        // read current orbital in orthonormal basis
        MatrixType Cprime;
        Cprime = DfObject::getCprimeMatrix<MatrixType>(runType, this->m_nIteration);

        // read previous orbital in orthonormal basis
        if ((this->m_nIteration != 1) ||
            (this->initialGuessType_ == GUESS_LCAO ||
             this->initialGuessType_ == GUESS_HUCKEL)) {
            prevCprime = DfObject::getCprimeMatrix<MatrixType>(runType, this->m_nIteration -1);
        } else if (this->m_nIteration == 1) {
            prevCprime = Cprime;
        }

        // calculate and print MO overlap between previous and current orbitals
        // C <-- C_pre^dagger * C_cur, in program  D <-- B^dagger * A
        prevCprime.transpose();
        prevCprime *= Cprime;
    }

    TlVector prevOcc = this->getOccupation(runType);

    // construct occupation number of current orbital with MO overlap matrix
    // 旧MO(pre)がどの新MO(crr)との重なりが一番大きいかを探す
    this->log_.info(" check overlap");
    TlVector currOcc(nNumOfMOs);
    {
        std::vector<bool> g(nNumOfMOs, false);
        bool bListHeaderOutput = false;
        for (int pre = 0; pre < nNumOfMOs; ++pre) {
            const double dPrevOcc = prevOcc[pre];
            if ((std::fabs(dPrevOcc - 1.00) < 1.0e-10) || // 占有数が 1 or 2 の時
                (std::fabs(dPrevOcc - 2.00) < 1.0e-10)) {

                const TlVector prevCprimeRowVector = prevCprime.getRowVector(pre);
                double xval = 0.0;
                int xord = -1;
                for (int crr = 0; crr < nNumOfMOs; ++crr) {
                    const double val = std::fabs(prevCprimeRowVector[crr]);
                    if ((g[crr] == false) && (val > xval)) {
                        xval = val; //fabs(prevCprime(pre, crr));
                        xord = crr;
                    }
                }

                if (xord == -1) {
                    this->log_.info(TlUtils::format(" crr MO %d th is not corresponded!\n", pre));
                    //CnErr.abort("DfDmatrix", "", "", "MO Overlap is not corresponded");
                }

                if (xval < 0.3) {
                    if (bListHeaderOutput == false) {
                        this->log_.info(" ##### MO Overlap is less than 0.3 #####\n");
                        bListHeaderOutput = true;
                    }
                    this->log_.info(TlUtils::format("pre MO %6d th -> crr MO %6d th %12.8lf\n",
                                                 (pre+1), (xord+1), xval));
                }
                currOcc[xord] = dPrevOcc;
                g[xord] = true;
            }
        }
    }

    this->log_.info(" check occupation vectors");
    this->checkOccupation(prevOcc, currOcc);

//     if (this->m_nIteration != 1) {
//         if (orbital_overlap != "on") {
//             this->logger(" orbital overlap correspondence is not carried out by the inputted keyword\n");
//             this->logger(" temporal occupation from orbital correspondence\n");
//             this->printOccupation(currOcc);

//             currOcc = prevOcc;
//         }
//     } else {
//         if (orbital_overlap_first != "on") {
//             this->logger(" orbital overlap correspondence is not carried out by the inputted keyword\n");
//             this->logger(" temporal occupation from orbital correspondence\n");
//             this->printOccupation(currOcc);

//             // 下から電子を詰めるようにする
//             this->logger(" the electrons are occupied from the 1 st orbital\n");

//             currOcc = prevOcc;

//             currOcc.sortByGrater();
//         }
//     }

    this->log_.info(" finish");
    return currOcc;
}


// 射影演算子法 ====================================================
template<typename MatrixType, typename SymmetricMatrixType>
TlVector DfDmatrix::getOccupationUsingProjection(const DfObject::RUN_TYPE runType)
{
    this->log_.info("orbital_overlap_method is mo-projection\n");

    const int nNumOfMOs = this->m_nNumOfMOs;
    //const int nNumOfAOs = this->m_nNumOfAOs;

    TlVector prevOcc = this->getOccupation(runType);
    TlVector currOcc(this->m_nNumOfMOs);

//     if (((this->m_nIteration == 1) &&
//          !((this->initialGuessType_ == GUESS_LCAO || this->initialGuessType_ == GUESS_HUCKEL) &&
//            orbital_overlap_first == "on"))) {

//         this->logger(" ... but ( iteration==1 || !(scf_start_guess==\"lcao\" && orbital_overlap_first==\"on\") ) is TRUE, skipping\n");
//         currOcc = prevOcc;
//     } else
    {
        // read molecular orbital occupation of (n-1) SCF iteration
        int num_mo_closed  = 0;
        int num_mo_open    = 0;
        int num_mo_virtual = 0;
        for (int k = 0; k < nNumOfMOs; ++k) {
            if (std::fabs(prevOcc[k] -2.00) < 1.0e-10) {
                num_mo_closed++;
            } else if (std::fabs(prevOcc[k] -1.00) < 1.0e-10) {
                num_mo_open++;
            } else if (std::fabs(prevOcc[k]) < 1.0e-10) {
                num_mo_virtual++;
            }
        }

        this->log_.info(TlUtils::format(" closed  orbital = %5ld", num_mo_closed));
        this->log_.info(TlUtils::format(" open    orbital = %5ld", num_mo_open));
        this->log_.info(TlUtils::format(" virtual orbital = %5ld", num_mo_virtual));

        // read molecular orbital ^(n)
        MatrixType C;
        C = DfObject::getCMatrix<MatrixType>(runType, this->m_nIteration);

        MatrixType Ct = C;
        Ct.transpose();
        
        SymmetricMatrixType S;
        S = DfObject::getSpqMatrix<SymmetricMatrixType>();

        // calculation of projection diagonal
        // projection diagonal for closed part
        TlVector pd_c(nNumOfMOs);
        if (num_mo_closed != 0) {
            // read density matrix of (n-1) SCF iteration
            SymmetricMatrixType D = DfObject::getPCMatrix<SymmetricMatrixType>(this->m_nIteration -1);

            // S * Dc * S
            const MatrixType SDS = S * D * S;

            // diagonal
            const MatrixType judge = Ct * SDS * C;
            for (int i = 0; i < nNumOfMOs; ++i) {
                pd_c[i] = judge.get(i, i);
            }
        }

        // projection diagonal for open part
        TlVector pd_o(this->m_nNumOfMOs);
        if (num_mo_open != 0) {
            // read density matrix of (n-1) SCF iteration
            SymmetricMatrixType D = DfObject::getPOMatrix<SymmetricMatrixType>(this->m_nIteration -1);

            // S * Do * S
            const MatrixType SDS = S * D * S;
           
            // diagonal
            const MatrixType judge = Ct * SDS * C;
            for (int i = 0; i < nNumOfMOs; ++i) {
                pd_o[i] = judge.get(i, i);
            }
        }

        for (int i=0; i < nNumOfMOs; ++i) {
            if (pd_c[i] < 0) {
                pd_c[i] *= -1;
            }
            if (pd_o[i] < 0) {
                pd_o[i] *= -1;
            }
        }

        // set order of MO (0から始まる軌道の番号)
        TlVector pd_c_ord(nNumOfMOs);
        for (int k=0; k < nNumOfMOs; ++k) {
            pd_c_ord[k] = k;
        }

        TlVector pd_o_ord(nNumOfMOs);
        for (int k = 0; k < nNumOfMOs; ++k) {
            pd_o_ord[k] = k;
        }

        // sort pd_c
        for (int i = 0; i < nNumOfMOs -1; ++i) {
            for (int j = i; j < nNumOfMOs; ++j) {
                if (pd_c[i] < pd_c[j]) {
                    std::swap(pd_c[i], pd_c[j]);
                    std::swap(pd_c_ord[i], pd_c_ord[j]);
                }
            }
        }
        // sort pd_o
        for (int i = 0; i < nNumOfMOs -1; ++i) {
            for (int j = i; j < nNumOfMOs; ++j) {
                if (pd_o[i] < pd_o[j]) {
                    std::swap(pd_o[i], pd_o[j]);
                    std::swap(pd_o_ord[i], pd_o_ord[j]);
                }
            }
        }

        if (num_mo_closed != 0) {
            this->printTwoVectors(pd_c_ord, pd_c, "projection diagonal of closed MO(sorted)", 10);
        } else {
            this->log_.info("projection diagonal of closed MO(sorted)");
            this->log_.info(" .. nothing of closed MO");
        }
        if (num_mo_open != 0) {
            this->printTwoVectors(pd_o_ord, pd_o, "projection diagonal of open   MO(sorted)", 10);
        } else {
            this->log_.info("projection diagonal of open   MO(sorted)");
            this->log_.info(" .. nothing of open MO");
        }

        // store electrons
        {
            currOcc = TlVector(nNumOfMOs);
            
            for (int k = 0; k < num_mo_open; ++k) {
                currOcc[(int)pd_o_ord[k] ] = 1.0;
            }
            for (int k = 0; k < num_mo_closed; ++k) {
                currOcc[(int)pd_c_ord[k] ] = 2.0;
            }
        }

        // check
        {
            double mustbe = num_mo_closed * 2 + num_mo_open;
            double ssum = 0.0;
            for (int k = 0; k < nNumOfMOs; ++k) {
                ssum += currOcc[k];
            }

            this->log_.info(TlUtils::format(" sum of electrons (new occupation) = %10.2lf", ssum));
            this->log_.info(TlUtils::format(" sum of electrons (must be)        = %10.2lf", mustbe));

            if (fabs(mustbe-ssum) > 1.0e-10) {
                this->log_.info("projection occupation was failed. change program code");
                CnErr.abort("DfDmatrix", "", "main", "projection occupation was failed");
            }
        }
    }

    return currOcc;
}

template<typename MatrixType, typename SymmetricMatrixType>
void DfDmatrix::generateDensityMatrix(const DfObject::RUN_TYPE runType, const TlVector& currOcc)
{
    // read current orbital in ao basis
    MatrixType C;
    C = DfObject::getCMatrix<MatrixType>(runType, this->m_nIteration);

    // generate densty matrix in terms of guess lcao and occupation
    switch (runType) {
    case RUN_RKS: {
        // closed shell Density Matrix
        SymmetricMatrixType PC = this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc, 2.0); // D1

        // open shell Density Matrix
        //SymmetricMatrixType PO = this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc, 1.0); // D2

        // projection法を書き換えれば不要(?)
        PC *= 2.0;
        //PC += PO;

        this->savePpqMatrix(runType, this->m_nIteration, PC);
    }
    break;

    case RUN_UKS_ALPHA:
        // go through
    case RUN_UKS_BETA: {
        // closed shell Density Matrix
        SymmetricMatrixType PC = this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc, 2.0); // D1
        //DfObject::savePCMatrix(this->m_nIteration, PC);

        // open shell Density Matrix
        SymmetricMatrixType PO = this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc, 1.0); // D2
        //DfObject::savePOMatrix(this->m_nIteration, PO);
        
        PC *= 2.0;
        PC += PO;

        this->savePpqMatrix(runType, this->m_nIteration, PC);
    }
    break;

    case RUN_ROKS: {
        // closed shell Density Matrix
        SymmetricMatrixType PC = this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc, 2.0); // D1
        DfObject::savePCMatrix(this->m_nIteration, PC);

        // open shell Density Matrix
        SymmetricMatrixType PO = this->calcDensMatrix<MatrixType, SymmetricMatrixType>(C, currOcc, 1.0); // D2
        DfObject::savePOMatrix(this->m_nIteration, PO);
    }
    break;

    default:
        std::cerr << " DfDmatrix::generateDensityMatrix() program error." << __FILE__ << __LINE__ << std::endl;
        break;
    }
}

template<typename MatrixType, typename SymmetricMatrixType>
SymmetricMatrixType DfDmatrix::calcDensMatrix(const MatrixType& inputC, const TlVector& f, double base)
{
    const std::size_t nNumOfAOs = inputC.getNumOfRows();
    const std::size_t nNumOfMOs = inputC.getNumOfCols();
    
    MatrixType C = inputC;
    {
        MatrixType E(nNumOfMOs, nNumOfAOs);

        std::size_t max_k = std::min<std::size_t>(nNumOfMOs, f.getSize());
        for (std::size_t k = 0; k < max_k; ++k) {
            if (std::fabs(f[k] - base) < 1.0e-10) {
                E(k, k) = 1.0;
            }
        }
        C = C * E;
    }

    MatrixType Ct = C;
    Ct.transpose();

    // DEBUG版ではここで対称性のチェックを行う
    const SymmetricMatrixType P = C * Ct;

    return P;
}

#endif // DFDMATRIX_H

