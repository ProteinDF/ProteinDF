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

#ifndef DFCONVERGE_DIIS
#define DFCONVERGE_DIIS

#include "DfConverge_Damping.h"

class DfConverge_DIIS : public DfConverge_Damping {
public:
    DfConverge_DIIS(TlSerializeData* pPdfParam);
    virtual ~DfConverge_DIIS();

protected:
    virtual void convergeRhoTilde();
    virtual void convergeKSMatrix();
    virtual void convergePMatrix();

protected:
    // template <class MatrixType, class SymmetricMatrixType, class VectorType>
    // void convergeRhoTilde(const DfObject::RUN_TYPE runType);

    template <class MatrixType, class SymmetricMatrixType>
    void convergeKSMatrix(const DfObject::RUN_TYPE runType);

    template <class MatrixType, class SymmetricMatrixType>
    void convergePMatrix(const DfObject::RUN_TYPE runType);

protected:
    template<class MatrixType, class SymmetricMatrixType>
    MatrixType getResidual(const RUN_TYPE runType,
                           const int itr);

    template<class MatrixType, class SymmetricMatrixType>
    TlMatrix buildBMatrix(const RUN_TYPE runType, 
                          const int itr, int last);

    template<class MatrixType, class SymmetricMatrixType>
    TlMatrix getCoef(const RUN_TYPE runType,
                     const int itr, int last);

protected:
    int startIterationOfDIIS_;
    int numOfLastItems_;
};


/// 残差を計算
///
/// 与えられたイテレーション番号(n)に対して
/// F(n)P(n-1)S - SP(n-1)F(n)を計算する
template<typename MatrixType, typename SymmetricMatrixType>
MatrixType DfConverge_DIIS::getResidual(const RUN_TYPE runType,
                                        const int itr)
{
    const std::string r_path = TlUtils::format("fl_Work/r_rks.%d.mat", itr);

    MatrixType r;    
    if (TlFile::isExist(r_path)) {
        r.load(r_path);
    } else {
        MatrixType FPS, SPF;
        {
            MatrixType FP, PF;
            {
                this->log_.info(TlUtils::format("residual: F(%d), P(%d)", itr, itr-1));
                const SymmetricMatrixType F = 
                    DfObject::getFpqMatrix<SymmetricMatrixType>(runType, itr);
                SymmetricMatrixType P = 
                    0.5 * DfObject::getPpqMatrix<SymmetricMatrixType>(runType, itr -1);
                FP = F * P;
                PF = P * F;
            }
            
            const SymmetricMatrixType S = DfObject::getSpqMatrix<SymmetricMatrixType>();
            FPS = FP * S;
            SPF = S * PF;
        }

        r = FPS - SPF;
        this->log_.info(TlUtils::format("save %s", r_path.c_str()));
        r.save(r_path);
    }

    return r;
}


/// B行列を作成する
template<class MatrixType, class SymmetricMatrixType>
TlMatrix DfConverge_DIIS::buildBMatrix(const RUN_TYPE runType,
                                       const int startItr,
                                       const int cycles)
{
    TlMatrix B(cycles +1, cycles +1);
    for (int i = 0; i < cycles; ++i) {
        this->log_.info(TlUtils::format("B[%d] -> itr %d", i, startItr +i));
        MatrixType r_i = this->getResidual<MatrixType, SymmetricMatrixType>(runType,
                                                                            startItr +i);

        for (int j = 0; j < i; ++j) {
            const MatrixType r_j = this->getResidual<MatrixType, SymmetricMatrixType>(runType,
                                                                                      startItr +j);
            
            MatrixType tmp = r_i;
            const double v = tmp.dot(r_j).sum();
            B.set(i, j, v); 
            B.set(j, i, v); 
        }
        B.set(i, i, r_i.dot(r_i).sum());
    }

    for (int i = 0; i < cycles; ++i) {
        B.set(i, cycles, -1.0);
        B.set(cycles, i, -1.0);
    }
    B.set(cycles, cycles, 0.0);
    
    B.save(TlUtils::format("fl_Work/B_rks.%d.mat", this->m_nIteration));

    return B;
}


/// 係数cを求める
template<class MatrixType, class SymmetricMatrixType>
TlMatrix DfConverge_DIIS::getCoef(const RUN_TYPE runType,
                                  const int startItr,
                                  int cycles)
{
    const TlMatrix B = this->buildBMatrix<MatrixType, SymmetricMatrixType>(runType, startItr, cycles);
    assert(B.getNumOfCols() == cycles +1);

    TlMatrix y(cycles +1, 1);
    for (int i = 0; i < cycles; ++i) {
        y.set(i, 0, 0.0);
    }
    y.set(cycles, 0, -1.0);

    const MatrixType c = B.solveLinearLeastSquaresProblem(y);
    c.save(TlUtils::format("fl_Work/c.%d.mat", this->m_nIteration));

    // debug
    MatrixType BC = B * c;
    BC.save(TlUtils::format("fl_Work/BC.%d.mat", this->m_nIteration));

    return c;
}

// template <class MatrixType, class SymmetricMatrixType, class VectorType>
// void DfConverge_DIIS::convergeRhoTilde(const DfObject::RUN_TYPE runType)
// {
//     this->log_.info(TlUtils::format("DIIS for rho: %d", this->m_nIteration));
//     const int itr = this->m_nIteration;

//     bool updateByDIIS = false;
//     if (itr >= this->startIterationOfDIIS_) {
//         // 残差ベクトルを求める
//         MatrixType r = this->getResidual<MatrixType, SymmetricMatrixType>(runType, itr -1); // 一つ前
        
//         // 残差ベクトルの要素が閾値(0.1)よりも大きい場合はDIISは行わない
//         const double maxError = r.getMaxAbsoluteElement();
//         this->log_.info(TlUtils::format("e_max = %e", maxError));
//         if (maxError < 0.1) {
//             int diis_e_max_pass_itr = (*this->pPdfParam_)["DIIS_e_max_pass_itr"].getInt();
//             if (diis_e_max_pass_itr == -1) {
//                 (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = itr;
//                 diis_e_max_pass_itr = itr;
//             }
//             this->log_.info(TlUtils::format("e_max pass itr = %d", diis_e_max_pass_itr));

//             int last = this->numOfLastItems_;
//             const int startItr = std::max(diis_e_max_pass_itr, itr - last +1); // 現itrを追加すると+1
//             last = std::min(last, itr - startItr +1);
//             this->log_.info(TlUtils::format("itr=%d, start=%d, last=%d",
//                                             itr, startItr, last));

//             // 係数行列の取得
//             const TlMatrix c = this->getCoef<MatrixType, SymmetricMatrixType>(runType, startItr -1, last);
//             assert(c.getNumOfRows() == last +1); // +1 はラグランジュ未定乗数法のため
            
//             this->log_.info("update rho using DIIS method.");
            
//             VectorType newRho(this->m_nNumOfAux);
//             VectorType tmpRho;
//             for (int i = 0; i < last; ++i) {
//                 const std::string inFilePath = DfObject::getRhoPath(runType, startItr +i);
//                 assert(TlFile::isExist(inFilePath));
//                 this->log_.info(inFilePath);
//                 tmpRho.load(inFilePath);
                
//                 const double coef = c.get(i, 0);
//                 this->log_.info(TlUtils::format("%d-th coef=% 8.3f",
//                                                 startItr +i, coef));
//                 newRho += coef * tmpRho;
//             }
            
//             this->log_.info(TlUtils::format("save new rho: %d", itr));
//             DfObject::saveRho(runType, itr, newRho);
//             updateByDIIS = true;
//         } else {
//             (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = -1;
//         }
//     } 

//     if (!updateByDIIS) {
//         DfConverge_Damping::convergeRhoTilde<VectorType>(runType);
//     }
// }

template <class MatrixType, class SymmetricMatrixType>
void DfConverge_DIIS::convergeKSMatrix(const DfObject::RUN_TYPE runType)
{
    this->log_.info(TlUtils::format("converge KS matrix using DIIS. #%d", this->m_nIteration));
    const int itr = this->m_nIteration;

    // backup
    {
        const std::string path = this->getFpqMatrixPath(runType, this->m_nIteration);
        SymmetricMatrixType F;
        F.load(path);
        F.save(path + ".diis");
    }

    bool updateByDIIS = false;
    if (itr >= this->startIterationOfDIIS_) {
        // 残差ベクトルを求める
        MatrixType r = this->getResidual<MatrixType, SymmetricMatrixType>(runType, itr);
        
        // 残差ベクトルの要素が閾値(0.1)よりも大きい場合はDIISは行わない
        const double maxError = r.getMaxAbsoluteElement();
        this->log_.info(TlUtils::format("e_max = %e", maxError));
        if (maxError < 0.1) {
            int diis_e_max_pass_itr = (*this->pPdfParam_)["DIIS_e_max_pass_itr"].getInt();
            if (diis_e_max_pass_itr == -1) {
                (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = itr;
                diis_e_max_pass_itr = itr;
            }
            this->log_.info(TlUtils::format("e_max pass itr = %d", diis_e_max_pass_itr));

            int last = this->numOfLastItems_;
            const int startItr = std::max(diis_e_max_pass_itr, itr - last +1);
            last = std::min(last, itr - startItr +1);
            this->log_.info(TlUtils::format("itr=%d, start=%d, last=%d",
                                            itr, startItr, last));

            // 係数行列の取得
            const TlMatrix c = this->getCoef<MatrixType, SymmetricMatrixType>(runType, startItr, last);
            assert(c.getNumOfRows() == last +1); // +1 はラグランジュ未定乗数法のため
            
            this->log_.info("update density matrix using DIIS method.");
            const index_type numOfAOs = this->m_nNumOfAOs;
            
            SymmetricMatrixType newF(numOfAOs);
            SymmetricMatrixType Fi;
            for (int i = 0; i < last; ++i) {
                const std::string inFilePath = DfObject::getFpqMatrixPath(runType, startItr +i) + ".diis";
                assert(TlFile::isExist(inFilePath));
                Fi.load(inFilePath);
                
                const double coef = c.get(i, 0);
                this->log_.info(TlUtils::format("%d-th coef=% 8.3f",
                                                startItr +i, coef));
                newF += coef * Fi;
            }
            
            this->log_.info(TlUtils::format("save new F matrix: %d", itr));
            DfObject::saveFpqMatrix(runType, itr, newF);
            updateByDIIS = true;
        } else {
            (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = -1;
        }
    } 

    if (!updateByDIIS) {
        DfConverge_Damping::convergeKSMatrix<SymmetricMatrixType>(runType);
    }
}

template <class MatrixType, class SymmetricMatrixType>
void DfConverge_DIIS::convergePMatrix(const DfObject::RUN_TYPE runType)
{
    this->log_.info(TlUtils::format("converge dennsity matrix using DIIS. #%d", this->m_nIteration));
    const int itr = this->m_nIteration;

    bool updateByDIIS = false;
    if (itr >= this->startIterationOfDIIS_) {
        // 残差ベクトルを求める
        MatrixType r = this->getResidual<MatrixType, SymmetricMatrixType>(runType, itr -1);

        // 残差ベクトルの要素が閾値(0.1)よりも大きい場合はDIISは行わない
        const double maxError = r.getMaxAbsoluteElement();
        this->log_.info(TlUtils::format("e_max = %e", maxError));
        if (maxError < 0.1) {
            int diis_e_max_pass_itr = (*this->pPdfParam_)["DIIS_e_max_pass_itr"].getInt();
            if (diis_e_max_pass_itr == -1) {
                (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = itr;
                diis_e_max_pass_itr = itr;
            }
            this->log_.info(TlUtils::format("e_max pass itr = %d", diis_e_max_pass_itr));

            int last = this->numOfLastItems_;
            const int startItr = std::max(diis_e_max_pass_itr, itr - last +1);
            last = std::min(last, itr - startItr +1);
            this->log_.info(TlUtils::format("itr=%d, start=%d, last=%d",
                                            itr, startItr, last));
            
            // 係数行列の取得
            const TlMatrix c = this->getCoef<MatrixType, SymmetricMatrixType>(runType, startItr -1, last);
            assert(c.getNumOfRows() == last +1); // +1 はラグランジュ未定乗数法のため
            
            this->log_.info("update density matrix using DIIS method.");
            
            SymmetricMatrixType newP(this->m_nNumOfAOs);
            SymmetricMatrixType Pi;
            for (int i = 0; i < last; ++i) {
                const std::string inFilePath = DfObject::getPpqMatrixPath(runType, startItr +i -1);
                assert(TlFile::isExist(inFilePath));
                Pi.load(inFilePath);
                
                const double coef = c.get(i, 0);
                this->log_.info(TlUtils::format("%d-th coef=% 8.3f",
                                                startItr +i -1, coef));
                newP += coef * Pi;
            }
            
            this->log_.info(TlUtils::format("save new P matrix: %d", itr -1));
            DfObject::savePpqMatrix(runType, itr -1, newP);
            updateByDIIS = true;
        } else {
            (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = -1;
        }
    }

    if (!updateByDIIS) {
        DfConverge_Damping::convergePMatrix<SymmetricMatrixType>(runType);
    }
}

#endif // DFCONVERGE_DIIS
