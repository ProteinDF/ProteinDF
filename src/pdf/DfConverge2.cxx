#include <fstream>
#include <cmath>

#include "DfConverge2.h"
#include "CnError.h"
#include "TlUtils.h"
#include "TlSymmetricMatrix.h"

#define LAMBDA  0.0

DfConverge2::DfConverge2(TlSerializeData* pPdfParam, int num_iter) : DfObject(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;

    // DIIS
    this->diis_dimension           = pdfParam["scf-acceleration/diis/number-of-diis"].getInt();
    this->diis_property            = pdfParam["scf-acceleration/diis/property"].getInt();
    this->diis_start_number        = pdfParam["scf-acceleration/diis/start-number"].getInt();
    this->diis_start_extrapolation = pdfParam["scf-acceleration/diis/start-extrapolation"].getInt();
    this->diis_type                = pdfParam["scf-acceleration/diis/type"].getStr();

    // DIIS actual dimension // Ｂ行列の次元とする
    this->diis_actual_dimension = this->m_nIteration - diis_start_number;
    this->diis_actual_dimension += 1;

    if (pdfParam["scf-acceleration"].getStr() == "mix") {
        this->log_.info("The mix acceleration method is employed.");
        if (this->diis_start_extrapolation != (this->diis_dimension + this->diis_start_number)) {
            this->log_.warn("In this method, start-extrapolation = number-of-diis + start-number.");
            this->log_.warn("Force to set the number of start-extrapolation.");
            this->diis_start_extrapolation = this->diis_dimension + this->diis_start_number;
        }

        int difference = this->m_nIteration - this->diis_start_number;
        int quotient   = (difference > 0 ? difference : 0) / (this->diis_dimension +1);

        this->diis_actual_dimension     = difference - quotient * (diis_dimension+1);
        this->diis_actual_dimension    += 1;
        this->diis_start_number        += quotient * (diis_dimension+1);
        this->diis_start_extrapolation += quotient * (diis_dimension+1);
    }

    if (this->diis_actual_dimension > (this->diis_dimension +1)) {
        this->diis_actual_dimension = this->diis_dimension +1;
    }

    // output informations
//     Log << "number_iteration    = " << m_nIteration << "\n";
//     Log << "number_mo_basis     = " << number_mo_basis  << "\n";
//     Log << "number_ao_basis     = " << number_ao_basis  << "\n";
//     Log << "//DIIS//\n";
//     Log << "diis_dimension           (scf-acceleration/diis/number-of-diis)        = " << diis_dimension           << " B行列の次元ー１\n";
//     Log << "diis_property            (scf-acceleration/diis/property)              = " << diis_property            << " 使ってない;FPS-SPFのみ" << "\n";
//     Log << "diis_start_number        (scf-acceleration/diis/start-number)          = " << diis_start_number        << "\n";
//     Log << "diis_start_extrapolation (scf-acceleration/diis/start-extrapolation)   = " << diis_start_extrapolation << "\n";
//     Log << "diis_type                (scf-acceleration/diis/type)                  = " << diis_type                << " Pulay版しか動かない"    << "\n";
//     Log << "diis_actual_dimension    ( dimension of B )                            = " << diis_actual_dimension    << "\n";
}


DfConverge2::~DfConverge2()
{
}

// 返値: DIIS を行うかどうかを返す
bool DfConverge2::DfConv2Main()
{
    bool bAnswer = true;

    if (!(this->diis_actual_dimension >= 2)) {
        this->log_.info("DIIS acceleralation is not started yet.");
        bAnswer = false;
    } else {
        this->log_.info("DIIS acceleralation is started.");

        TlMatrix Ev;
        TlMatrix Bm;
        TlVector c;

        switch (this->m_nMethodType) {
        case METHOD_RKS:
            {
                // RKS
                this->calculate_ErrorVector(RUN_RKS, m_nIteration, Ev);
                this->read_PreviousBmatrix("rks", m_nIteration, Bm);
                this->update_Bmatrix("rks", m_nIteration, Bm, Ev);
                if (diis_start_extrapolation <= m_nIteration) {
                    this->log_.info("solve DIIS equation and extrapolate khon-sham matrix.");
                    this->solve_DIIS(Bm, c);
                    this->interpolate_KhonSham("rks", m_nIteration, c);
                    
                    bAnswer = true;
                } else {
                    this->log_.info("only store updated B matrix, do not solve DIIS equation.");
                    bAnswer = false;
                }
            }
            break;

        case METHOD_UKS:
            {
                // alpha
                this->calculate_ErrorVector(RUN_UKS_ALPHA, m_nIteration, Ev);
                this->read_PreviousBmatrix("uks-alpha", m_nIteration, Bm);
                this->update_Bmatrix("uks-alpha", m_nIteration, Bm, Ev);
                if (diis_start_extrapolation <= m_nIteration) {
                    this->log_.info("solve DIIS equation and extrapolate khon-sham matrix.");
                    this->solve_DIIS(Bm, c);
                    this->interpolate_KhonSham("uks-alpha", m_nIteration, c);
                    
                    bAnswer = true;
                } else {
                    this->log_.info("only store updated B matrix, do not solve DIIS equation.");
                    
                    bAnswer = false;
                }
                // beta
                this->calculate_ErrorVector(RUN_UKS_BETA, m_nIteration, Ev);
                this->read_PreviousBmatrix("uks-beta", m_nIteration, Bm);
                this->update_Bmatrix("uks-beta", m_nIteration, Bm, Ev);
                if (diis_start_extrapolation <= m_nIteration) {
                    this->log_.info("solve DIIS equation and extrapolate khon-sham matrix.");
                    this->solve_DIIS(Bm, c);
                    this->interpolate_KhonSham("uks-beta",  m_nIteration, c);
                    
                    bAnswer = true;
                } else {
                    this->log_.info("only store updated B matrix, do not solve DIIS equation.");
                    bAnswer = false;
                }
            }
            break;

        case METHOD_ROKS:
            {
                this->calculate_ErrorVector(RUN_ROKS, m_nIteration, Ev);
                this->read_PreviousBmatrix("roks", m_nIteration, Bm);
                this->update_Bmatrix("roks", m_nIteration, Bm, Ev);
                if (diis_start_extrapolation <= m_nIteration) {
                    this->log_.info("solve DIIS equation and extrapolate khon-sham matrix.");
                    this->solve_DIIS(Bm, c);
                    this->interpolate_KhonSham("roks", m_nIteration, c);
                    
                    bAnswer = true;
                } else {
                    this->log_.info("only store updated B matrix, do not solve DIIS equation.");
                    bAnswer = false;
                }
            }
            break;

        default:
            CnErr.abort();
            break;
        }
    }

    return bAnswer;
}


// フル行列２個と、三角行列１個のメモリが必要である。
void DfConverge2::calculate_ErrorVector(const RUN_TYPE runType, int iteration, TlMatrix& A)
{
    TlSymmetricMatrix P(this->m_nNumOfAOs);
    TlSymmetricMatrix F(this->m_nNumOfAOs);
    TlSymmetricMatrix S(this->m_nNumOfAOs);

    if (runType == RUN_ROKS) {
        P.load(this->getP1pqMatrixPath(iteration -1));

        // read previous density matrix of open part
        F.load(this->getP2pqMatrixPath(iteration -1));

        // total density matrix
        P *= 2.0;
        P += F;
    } else {
        // for rks, uks-alpha or uks-beta case

        // read previous density matrix
        P.load(this->getPpqMatrixPath(runType, iteration -1));
    }

    // read S matrix
    DfObject::saveSpqMatrix(S);

    // "read F matrix"
    F = DfObject::getFpqMatrix<TlSymmetricMatrix>(runType, iteration);

    // generate FPS - SPF
    TlMatrix B(this->m_nNumOfAOs, this->m_nNumOfAOs);
    {
        // calc FPS
        {
            A = F * P;

            for (int i=0; i < this->m_nNumOfAOs; i++) {
                for (int j=0; j < this->m_nNumOfAOs; j++) {
                    for (int k=0; k < this->m_nNumOfAOs; k++) {
                        B(i, j) += A(i, k) * S(k, j);
                    }
                }
            }
        }

        // calc SPF
        for (int i=0; i < this->m_nNumOfAOs; i++) {
            for (int j=0; j < this->m_nNumOfAOs; j++) {
                A(i, j) = B(j, i);
            }
        }

        // calc FPS - SPF
        B -= A;
    }

    // generate error vector e_i ;  e_i = X^dagger (FPS - SPF) X
    {
        A = DfObject::getXMatrix<TlMatrix>();

        B *= A;

        A.transpose();

        A *= B;
    }

    // write error vector for next DIIS
    {
        this->log_.info("write Error Vector to a file.");

        assert(A.getNumOfCols() == this->m_nNumOfMOs);
        assert(A.getNumOfRows() == this->m_nNumOfMOs);

        // "write Error Vector"
        std::string suffix = DfObject::m_sRunTypeSuffix[runType];
        A.save("fl_Work/fl_Mtr_Errvct.matrix." + suffix + TlUtils::xtos(iteration));
    }

    // X^dagger (FPS - SPF) X スーパベクトルのノルムらしい----収束するとゼロとなる量である
    {
        double  norm  = 0.0;
        for (int i=0; i < this->m_nNumOfMOs; i++) {
            for (int j=0; j < this->m_nNumOfMOs; j++) {
                norm += A(i, j) * A(i, j);
            }
        }

        this->log_.info(TlUtils::format("norm of the error vector (%d iteration) = %f.",
                                        iteration, norm));
        if (diis_start_extrapolation > m_nIteration) {
            this->log_.info("do not extrapolation.");
        }
    }
}

void DfConverge2::read_PreviousBmatrix(const std::string& type, int iteration, TlMatrix& X)
{
    if (diis_actual_dimension >= 3) {
        // read previous B matrix
        X.load("fl_Work/fl_Mtr_Bmatrix.matrix." + type + TlUtils::xtos(iteration -1));
    } else {
        // create 1st B matrix
        X.resize(1, 1);
        X(0, 0) = 0.0;
    }
}

void DfConverge2::update_Bmatrix(const std::string& type, int iteration, TlMatrix& B, TlMatrix& E)
{
    TlMatrix G(this->m_nNumOfMOs, this->m_nNumOfMOs);  // for previous Error vectors

    // shift previous B matrix if need
    if (diis_actual_dimension <= m_nIteration - diis_start_number) {
        //Log << "B shift DEBUG, before" << "\n";
        for (int i = 1; i < diis_actual_dimension; i++) {
            for (int j=1; j < diis_actual_dimension; j++) {
                B(i-1, j-1) = B(i, j);
            }
        }

        for (int i=1; i < diis_actual_dimension; i++) {
            B(diis_actual_dimension -1, i) = 0;
            B(i, diis_actual_dimension -1) = 0;
        }
    } else {
        B.resize(diis_actual_dimension, diis_actual_dimension);
    }

    // fill 1st row and 1st col elements
    B(0, 0) = 0.0;
    for (int i=1; i < diis_actual_dimension; i++) {
        B(i, 0) = -1.0;
        B(0, i) = -1.0;
    }

    // B行列の追加要素
    for (int bm_k = 1; bm_k < diis_actual_dimension; bm_k++) {
        G.load("fl_Work/fl_Mtr_Errvct.matrix." + type + TlUtils::xtos(iteration -bm_k +1));

        // inner product between error vectors, <G|E>
        for (int pk = 0; pk < this->m_nNumOfMOs; pk++) {
            for (int qk = 0; qk < this->m_nNumOfMOs; qk++) {
                B((diis_actual_dimension - bm_k +1) -1, diis_actual_dimension -1) += G(pk, qk) * E(pk, qk);
            }
        }
        B(diis_actual_dimension -1, (diis_actual_dimension - bm_k +1) -1)
        = B((diis_actual_dimension - bm_k +1) -1, diis_actual_dimension -1);
    }

    // write B matrix
    B.save("fl_Work/fl_Mtr_Bmatrix.matrix." + type + TlUtils::xtos(iteration));
}

void DfConverge2::interpolate_KhonSham(const std::string& type, int iteration, TlVector& c)
{
    TlSymmetricMatrix F(this->m_nNumOfAOs);

    for (int k = 1; k < diis_actual_dimension; k++) {
        TlSymmetricMatrix T;
        T.load("fl_Work/fl_Mtr_Fpq.matrix." + type + TlUtils::xtos(k + (iteration - diis_actual_dimension) +1));

        for (int i=0; i < this->m_nNumOfAOs; i++) {
            for (int j=0; j<=i; j++) {
                F(i, j) += c[k] * T(i, j);
            }
        }
    }

    // write F matrix to a file
    F.save("fl_Work/fl_Mtr_Fpq.matrix." + type + TlUtils::xtos(-10));
}

void DfConverge2::solve_DIIS(TlMatrix& B, TlVector& c)
{
    TlVector L(B.getNumOfRows());

    c = TlVector(B.getNumOfCols());
    c[0] = -LAMBDA;

    L[0] = -1.0;
    for (int i = 1; i < L.getSize(); i++) {
        L[i] = 0.0;
    }

    // 連立一次方程式を解く
    // 元は LUdecomp となっていたが、これってガウスの消去法と同一か？
    this->gauss(B, c, L);
}

// B * X = L の連立一次方程式を解く
void DfConverge2::gauss(TlMatrix& B, TlVector& X, TlVector& L)
{
    // check
    assert(B.getNumOfCols() == X.getSize());
    assert(B.getNumOfRows() == L.getSize());
    assert(B.getNumOfRows() == B.getNumOfCols());

    int dimension = X.getSize();
    double*  w  = new double [ dimension ];
    int*    ip = new int   [ dimension ];

    //initialize ip
    for (int i = 0; i < dimension; i++) {
        ip[i] = i;
        w[i]  = 0.0;
    }

    const double eps = 1e-30; // threshhold for pivot exchange
    // gauss elemination
    for (int k = 0; k < dimension; k++) {
        int l = k;
        double al = fabs(B(ip[l], k));
        for (int i = k +1; i < dimension; i++) {
            if (fabs(B(ip[i], k)) > al) {
                l = i;
                al = fabs(B(ip[l], k));
            }
        }

        // row echange
        if (l != k) {
            int tmp = ip[k];
            ip[k] = ip[l];
            ip[l] = tmp;
        }

        // check singularity
        if (fabs(B(ip[k], k)) < eps) {
            std::cerr << "singular matrix in LUdecomp: ";
            std::cerr << TlUtils::format("B(%ld,%ld)=%le < %le (eps)", ip[k], k, B(ip[k], k), eps) << std::endl;
            CnErr.abort();
        }

        // gauss elemination start
        B(ip[k], k) = 1.0 / B(ip[k], k);
        for (int i = k+1; i < dimension; i++) {
            B(ip[i], k) *= B(ip[k], k);
            for (int j = k+1; j < dimension; j++) {
                w[j] = B(ip[i], j) - B(ip[i], k) * B(ip[k], j);
            }
            for (int j = k +1; j < dimension; j++) {
                B(ip[i], j) = w[j];
            }
        }
    }

    // solve linear equations by forward, backward substitution
    // forward substitution
    for (int i=0; i < dimension; i++) {
        double t = L[ip[i] ];
        for (int j = 0; j < i; j++) {
            t -= B(ip[i], j) * X[j];
        }
        X[i] = t;
    }

    // backward substitution
    for (int i = dimension -1; i >= 0; i--) {
        double t = X[i];
        for (int j=i +1; j < dimension; j++) {
            t -= B(ip[i], j) * X[j];
        }
        X[i] = t * B(ip[i], i);
    }

    delete[] w;
    delete[] ip;
}
