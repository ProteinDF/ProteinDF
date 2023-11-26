#ifndef DF_CONVERGE_DIIS_TEMPL_H
#define DF_CONVERGE_DIIS_TEMPL_H

#include <string>

#include "df_converge_damping_template.h"

template <class Matrix, class SymmetricMatrix, class Vector>
class DfConverge_Damping_Diis_Template : public DfConverge_Damping_Template<SymmetricMatrix, Vector> {
public:
    DfConverge_Damping_Diis_Template(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Diis_Template();

protected:
    virtual void convergeRhoTilde(const DfObject::RUN_TYPE runType);
    virtual void convergeKSMatrix(const DfObject::RUN_TYPE runType);
    virtual void convergePMatrix(const DfObject::RUN_TYPE runType);

protected:
    //     virtual void convergeRhoTilde();
    //     virtual void convergeKSMatrix();
    //     virtual void convergePMatrix();

protected:
    Matrix getResidual(const DfObject::RUN_TYPE runType, const int itr);
    Matrix buildBMatrix(const DfObject::RUN_TYPE runType, const int itr, int last);
    Matrix getCoef(const DfObject::RUN_TYPE runType, const int itr, int last);

protected:
    int startIterationOfDIIS_;
    int numOfLastItems_;
};

template <class Matrix, class SymmetricMatrix, class Vector>
DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::DfConverge_Damping_Diis_Template(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Template<SymmetricMatrix, Vector>(pPdfParam) {
}

template <class Matrix, class SymmetricMatrix, class Vector>
DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::~DfConverge_Damping_Diis_Template() {
}

/// 残差を計算
///
/// 与えられたイテレーション番号(n)に対して
/// F(n)P(n-1)S - SP(n-1)F(n)を計算する
template <class Matrix, class SymmetricMatrix, class Vector>
Matrix DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::getResidual(const DfObject::RUN_TYPE runType, const int itr) {
    const std::string r_path = DfObject::getDiisResidualMatrixPath(runType, itr);

    Matrix r;
    if (TlFile::isExistFile(r_path)) {
        r.load(r_path);
    } else {
        Matrix FPS, SPF;
        {
            Matrix FP, PF;
            {
                this->log_.info(TlUtils::format("residual: F(%d), P(%d)", itr, itr - 1));
                const SymmetricMatrix F = DfObject::getFpqMatrix<SymmetricMatrix>(runType, itr);
                SymmetricMatrix P = 0.5 * DfObject::getPInMatrix<SymmetricMatrix>(runType, itr);
                FP = F * P;
                PF = P * F;
            }

            const SymmetricMatrix S = DfObject::getSpqMatrix<SymmetricMatrix>();
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
template <class Matrix, class SymmetricMatrix, class Vector>
Matrix DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::buildBMatrix(const DfObject::RUN_TYPE runType, const int startItr, const int cycles) {
    Matrix B(cycles + 1, cycles + 1);
    for (int i = 0; i < cycles; ++i) {
        this->log_.info(TlUtils::format("B[%d] -> itr %d", i, startItr + i));
        Matrix r_i = this->getResidual(runType, startItr + i);

        for (int j = 0; j < i; ++j) {
            const Matrix r_j = this->getResidual(runType, startItr + j);

            Matrix tmp = r_i;
            const double v = tmp.dotInPlace(r_j).sum();
            B.set(i, j, v);
            B.set(j, i, v);
        }
        B.set(i, i, r_i.dotInPlace(r_i).sum());
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
template <class Matrix, class SymmetricMatrix, class Vector>
Matrix DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::getCoef(const DfObject::RUN_TYPE runType, const int startItr, int cycles) {
    const Matrix B = this->buildBMatrix(runType, startItr, cycles);
    assert(B.getNumOfCols() == cycles + 1);

    Matrix y(cycles + 1, 1);
    for (int i = 0; i < cycles; ++i) {
        y.set(i, 0, 0.0);
    }
    y.set(cycles, 0, -1.0);

    const Matrix c = B.getLeastSquaresSolution(y);
    c.save(TlUtils::format("fl_Work/c.%d.mat", this->m_nIteration));

    // debug
    Matrix BC = B * c;
    BC.save(TlUtils::format("fl_Work/BC.%d.mat", this->m_nIteration));

    return c;
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergeRhoTilde(const DfObject::RUN_TYPE runType) {
    this->log_.info(TlUtils::format("converge rho vector using DIIS. #%d", this->m_nIteration));
    const int itr = this->m_nIteration;

    bool updateByDIIS = false;
    if (itr >= this->startIterationOfDIIS_) {
        // 残差ベクトルを求める
        Matrix r = this->getResidual(runType, itr - 1);  // 一つ前

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
            const int startItr = std::max(diis_e_max_pass_itr, itr - last + 1);  // 現itrを追加すると+1
            last = std::min(last, itr - startItr + 1);
            this->log_.info(TlUtils::format("itr=%d, start=%d, last=%d", itr, startItr, last));

            // 係数行列の取得
            const Matrix c = this->getCoef(runType, startItr - 1, last);
            assert(c.getNumOfRows() == last + 1);  // +1 はラグランジュ未定乗数法のため

            this->log_.info("update rho using DIIS method.");

            Vector newRho(this->m_nNumOfAux);
            Vector tmpRho;
            for (int i = 0; i < last; ++i) {
                const std::string inFilePath = DfObject::getRhoPath(runType, startItr + i);
                assert(TlFile::isExistFile(inFilePath));
                this->log_.info(inFilePath);
                tmpRho.load(inFilePath);

                const double coef = c.get(i, 0);
                this->log_.info(TlUtils::format("%d-th coef=% 8.3f", startItr + i, coef));
                newRho += coef * tmpRho;
            }

            this->log_.info(TlUtils::format("save new rho: %d", itr));
            DfObject::saveRho(runType, itr, newRho);
            updateByDIIS = true;
        } else {
            (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = -1;
        }
    }

    if (!updateByDIIS) {
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeRhoTilde(runType);
    }
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergeKSMatrix(const DfObject::RUN_TYPE runType) {
    this->log_.info(TlUtils::format("converge KS matrix using DIIS. #%d", this->m_nIteration));
    const int itr = this->m_nIteration;

    // backup
    {
        const std::string path = this->getFpqMatrixPath(runType, this->m_nIteration);
        SymmetricMatrix F;
        F.load(path);
        F.save(path + ".diis");
    }

    bool updateByDIIS = false;
    if (itr >= this->startIterationOfDIIS_) {
        // 残差ベクトルを求める
        Matrix r = this->getResidual(runType, itr);

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
            const int startItr = std::max(diis_e_max_pass_itr, itr - last + 1);
            last = std::min(last, itr - startItr + 1);
            this->log_.info(TlUtils::format("itr=%d, start=%d, last=%d", itr, startItr, last));

            // 係数行列の取得
            const Matrix c = this->getCoef(runType, startItr, last);
            assert(c.getNumOfRows() == last + 1);  // +1 はラグランジュ未定乗数法のため

            this->log_.info("update density matrix using DIIS method.");
            const DfObject::index_type numOfAOs = this->m_nNumOfAOs;

            SymmetricMatrix newF(numOfAOs);
            SymmetricMatrix Fi;
            for (int i = 0; i < last; ++i) {
                const std::string inFilePath = DfObject::getFpqMatrixPath(runType, startItr + i) + ".diis";
                assert(TlFile::isExistFile(inFilePath));
                Fi.load(inFilePath);

                const double coef = c.get(i, 0);
                this->log_.info(TlUtils::format("%d-th coef=% 8.3f", startItr + i, coef));
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
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergeKSMatrix(runType);
    }
}

template <class Matrix, class SymmetricMatrix, class Vector>
void DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergePMatrix(const DfObject::RUN_TYPE runType) {
    this->log_.info(TlUtils::format("converge density matrix using DIIS. #%d/%d", this->m_nIteration, this->startIterationOfDIIS_));
    const int itr = this->m_nIteration;

    bool updateByDIIS = false;
    if (itr >= this->startIterationOfDIIS_) {
        // 残差ベクトルを求める
        Matrix r = this->getResidual(runType, itr - 1);

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
            const int startItr = std::max(diis_e_max_pass_itr, itr - last + 1);
            last = std::min(last, itr - startItr + 1);
            this->log_.info(TlUtils::format("itr=%d, start=%d, last=%d", itr, startItr, last));

            // 係数行列の取得
            const Matrix c = this->getCoef(runType, startItr - 1, last);
            assert(c.getNumOfRows() == last + 1);  // +1 はラグランジュ未定乗数法のため

            this->log_.info("update density matrix using DIIS method.");

            SymmetricMatrix newP(this->m_nNumOfAOs);
            SymmetricMatrix Pi;
            for (int i = 0; i < last; ++i) {
                const std::string inFilePath = DfObject::getPInMatrixPath(runType, startItr + i);
                assert(TlFile::isExistFile(inFilePath));
                Pi.load(inFilePath);

                const double coef = c.get(i, 0);
                this->log_.info(TlUtils::format("%d-th coef=% 8.3f", startItr + i - 1, coef));
                newP += coef * Pi;
            }

            this->log_.info(TlUtils::format("save new P matrix: %d", itr - 1));
            DfObject::savePInMatrix(runType, itr, newP);
            updateByDIIS = true;
        } else {
            this->log_.info("e_max has not been satisfied.");
            (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = -1;
        }
    } else {
        this->log_.info(TlUtils::format("the iteration, %d, is lower than start number, %d.", itr, this->startIterationOfDIIS_));
    }

    if (!updateByDIIS) {
        DfConverge_Damping_Template<SymmetricMatrix, Vector>::convergePMatrix(runType);
    }
}

// template <class Matrix, class SymmetricMatrix, class Vector>
// void DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergeRhoTilde() {
//     this->log_.critical("NOT support: DIIS for rho.");

//     // switch (this->m_nMethodType) {
//     // case METHOD_RKS:
//     //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
//     //     TlDenseSymmetricMatrix_Lapack,
//     //     TlDenseVector_Lapack>(DfObject::RUN_RKS); break;
//     // case METHOD_UKS:
//     //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
//     //     TlDenseSymmetricMatrix_Lapack,
//     //     TlDenseVector_Lapack>(DfObject::RUN_UKS_ALPHA);
//     //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
//     //     TlDenseSymmetricMatrix_Lapack,
//     //     TlDenseVector_Lapack>(DfObject::RUN_UKS_BETA);
//     //     break;
//     // case METHOD_ROKS:
//     //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
//     //     TlDenseSymmetricMatrix_Lapack,
//     //     TlDenseVector_Lapack>(DfObject::RUN_UKS_ALPHA);
//     //     this->convergeRhoTilde<TlDenseGeneralMatrix_Lapack,
//     //     TlDenseSymmetricMatrix_Lapack,
//     //     TlDenseVector_Lapack>(DfObject::RUN_UKS_BETA);
//     //     break;
//     // default:
//     //     std::cerr << "program error. @DfConverge_Damping::convergeRhoTilde()"
//     //     << std::endl; break;
//     // }
// }

// template <class Matrix, class SymmetricMatrix, class Vector>
// void DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergeKSMatrix() {
//     switch (this->m_nMethodType) {
//         case DfObject::METHOD_RKS:
//             this->convergeKSMatrix(DfObject::RUN_RKS);
//             break;
//         case DfObject::METHOD_UKS:
//             this->convergeKSMatrix(DfObject::RUN_UKS_ALPHA);
//             this->convergeKSMatrix(DfObject::RUN_UKS_BETA);
//             break;
//         case DfObject::METHOD_ROKS:
//             this->convergeKSMatrix(DfObject::RUN_ROKS_CLOSED);
//             this->convergeKSMatrix(DfObject::RUN_ROKS_OPEN);
//             break;
//         default:
//             std::cerr
//                 << "program error. @DfConverge_Damping::convergeKSMatrix()"
//                 << std::endl;
//             break;
//     }
// }

// template <class Matrix, class SymmetricMatrix, class Vector>
// void DfConverge_Damping_Diis_Template<Matrix, SymmetricMatrix, Vector>::convergePMatrix() {
//     switch (this->m_nMethodType) {
//         case DfObject::METHOD_RKS:
//             this->convergePMatrix(DfObject::RUN_RKS);
//             break;
//         case DfObject::METHOD_UKS:
//             this->convergePMatrix(DfObject::RUN_UKS_ALPHA);
//             this->convergePMatrix(DfObject::RUN_UKS_BETA);
//             break;
//         case DfObject::METHOD_ROKS:
//             this->convergePMatrix(DfObject::RUN_ROKS_CLOSED);
//             this->convergePMatrix(DfObject::RUN_ROKS_OPEN);
//             break;
//         default:
//             std::cerr << "program error. @DfConverge_Damping::convergePMatrix()"
//                       << std::endl;
//             break;
//     }
// }

#endif  // DF_CONVERGE_DIIS_TEMPL_H
