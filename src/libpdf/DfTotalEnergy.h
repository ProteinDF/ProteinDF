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

#ifndef DFTOTALENERGY_H
#define DFTOTALENERGY_H

#include "DfObject.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"

#include "Fl_Geometry.h"
#include "CnError.h"
#include "DfXCFunctional.h"

class DfEriX;
class DfOverlapX;

/// 得られた密度行列、交換相関ポテンシャルをもとに全エネルギーの計算を行う
class DfTotalEnergy : public DfObject {
public:
    DfTotalEnergy(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy();

    /// nsp、sp、roks の判断を行い、各エネルギーの計算ルーチンを呼び、全エネルギーの計算を行う
    virtual void exec();

public:
    /// ダミー原子を除いた時の全エネルギー計算を行う
    virtual void calculate_real_energy();

protected:
    /// 核—核反発エネルギーの計算を行う
    double calculate_energy_nuclear_repulsion();

    /// J[rho~, rho~]の計算を行う
    template<typename SymmetricMatrixType, typename VectorType>
    double calcJRhoTildeRhoTilde(const VectorType& rho);

    /// 全エネルギーを出力する
    void write_total_energy(double E_Total) const;

protected:
    template<typename DfOverlapType, typename DfEriType, typename SymmetricMatrixType, typename VectorType>
    void exec_template();

    /// 密度行列を用意する
    template<typename SymmetricMatrixType>
    SymmetricMatrixType getPpq(METHOD_TYPE methodType);

    /// 密度行列の展開係数を用意する
    template<typename VectorType>
    VectorType getRho(METHOD_TYPE methodType);

    /// 交換相関ポテンシャル展開係数を用意する
    template<typename VectorType>
    VectorType getEps(RUN_TYPE runType);

    /// (一電子部分 + J[rho, rho~] + 交換相関ポテンシャ部分)の計算を行う
    template<typename SymmetricMatrixType>
    double calculate_E_WITH_DIRECT(const SymmetricMatrixType& D);

    /// 一電子部分の計算を行う
    template<typename SymmetricMatrixType>
    double calcOneElectronPart(const SymmetricMatrixType& D);

    /// Jの計算を行う
    template<typename SymmetricMatrixType>
    double calcJ(const SymmetricMatrixType& P);

    /// J[rho, rho~]の計算を行う（積分を直接計算）
    template<typename DfEriType, typename SymmetricMatrixType, typename VectorType>
    double calcJRhoRhoTilde_DIRECT(const SymmetricMatrixType& D);
    
    /// 交換相関ポテンシャル部分の計算を行う
    template<typename DfOverlapType, typename SymmetricMatrixType, typename VectorType>
    double calcExc_DIRECT(const SymmetricMatrixType& D, const VectorType& E);
    
    template<typename SymmetricMatrixType>
    double calcExc(const RUN_TYPE runType,
                   const SymmetricMatrixType& P);

    /// Kの計算を行う
    template<typename SymmetricMatrixType>
    double calcK(const RUN_TYPE runType, const SymmetricMatrixType& P);

    template <class SymmetricMatrixType>
    void calcRealEnergy();

    template<class SymmetricMatrixType>
    double calcEnergyFromDummy();

protected:
    virtual void output();

    virtual DfEriX* getDfEriX() const;
    virtual DfOverlapX* getDfOverlapX() const;

protected:
    double m_dE_OneElectronPart;
    double J_term_;
    double m_dE_J_Rho_RhoTilde;
    double m_dE_J_RhoTilde_RhoTilde;
    double m_dExc;

    double K_term_;
    double E_KA_;
    double E_KB_;

    double m_dE_NuclearRepulsion;
    double m_dE_OEP_JRR_Exc;
    double E_disp_; /// Grimme dispersion Energy

    const static double TOO_SMALL;
    static double m_dNuclearRepulsion;
};


template<typename DfOverlapType, typename DfEriType, typename SymmetricMatrixType, typename VectorType>
void DfTotalEnergy::exec_template()
{
    // setup density matrix
    SymmetricMatrixType Ppq, PpqA, PpqB;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        Ppq = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration);
        break;

    case METHOD_UKS:
        PpqA = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA, this->m_nIteration);
        PpqB = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_UKS_BETA,  this->m_nIteration);
        Ppq = PpqA + PpqB;
        break;

    case METHOD_ROKS:
        PpqB = 0.5 * DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSED, this->m_nIteration);
        PpqA = PpqB + DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN,  this->m_nIteration);
        Ppq = PpqA + PpqB;
        break;

    default:
        this->log_.critical("program error");
    }

    // 核-核反発
    this->m_dE_NuclearRepulsion = this->calculate_energy_nuclear_repulsion();
 
    // クーロン項
    switch (this->J_engine_) {
    case J_ENGINE_RI_J:
        {
            // use Threeindexintegrals
            const VectorType rho = this->getRho<VectorType>(this->m_nMethodType);
            this->m_dE_J_RhoTilde_RhoTilde =
                this->calcJRhoTildeRhoTilde<SymmetricMatrixType, VectorType>(rho);

            if (! ((this->m_bMemorySave == false) && (this->m_bDiskUtilization == false))) {
                // NOT use Threeindexintegrals
                this->m_dE_J_Rho_RhoTilde =
                    this->calcJRhoRhoTilde_DIRECT<DfEriType, SymmetricMatrixType, VectorType>(Ppq);
            }
        }
        break;

    case J_ENGINE_CONVENTIONAL:
    case J_ENGINE_CD:
        this->J_term_ = this->calcJ(Ppq);
        break;

    default:
        break;
    }
    
    if ((this->m_bMemorySave == false) && (this->m_bDiskUtilization == false)) {
        // use Threeindexintegrals
        assert(this->m_bIsXCFitting == true);
        this->m_dE_OEP_JRR_Exc = this->calculate_E_WITH_DIRECT(Ppq);
    } else {
        // NOT use Threeindexintegrals
        this->m_dE_OneElectronPart = this->calcOneElectronPart(Ppq);

        DfXCFunctional dfXCFunctional(this->pPdfParam_);
        switch (this->m_nMethodType) {
        case METHOD_RKS:
            {
                if (dfXCFunctional.getXcType() != DfXCFunctional::HF) {
                    if (this->m_bIsXCFitting == true) {
                        const VectorType Eps = this->getEps<VectorType>(RUN_RKS);
                        this->m_dExc = this->calcExc_DIRECT<DfOverlapType, SymmetricMatrixType, VectorType>(Ppq, Eps);
                    } else {
                        if (this->XC_engine_ != XC_ENGINE_GRID) {
                            this->m_dExc = this->calcExc(RUN_RKS, 0.5 * Ppq) * 2.0;
                        } else {
                            DfXCFunctional dfXCFunctional(this->pPdfParam_);
                            this->m_dExc = dfXCFunctional.getEnergy();
                            if (this->enableGrimmeDispersion_ == true) {
                                this->E_disp_ = dfXCFunctional.getGrimmeDispersionEnergy();
                            }
                        }
                    }
                }
                if (dfXCFunctional.isHybridFunctional() == true) {
                    this->K_term_ = this->calcK(RUN_RKS, 0.5 * Ppq) * 2.0;
                }
            }
            break;

        case METHOD_UKS:
            {
                if (dfXCFunctional.getXcType() != DfXCFunctional::HF) {
                    if (this->m_bIsXCFitting == true) {
                        const VectorType EpsA = this->getEps<VectorType>(RUN_UKS_ALPHA);
                        const VectorType EpsB = this->getEps<VectorType>(RUN_UKS_BETA);
                        this->m_dExc  = this->calcExc_DIRECT<DfOverlapType, SymmetricMatrixType, VectorType>(PpqA, EpsA);
                        this->m_dExc += this->calcExc_DIRECT<DfOverlapType, SymmetricMatrixType, VectorType>(PpqB, EpsB);
                    } else {
                        if (this->XC_engine_ != XC_ENGINE_GRID) {
                            this->m_dExc  = this->calcExc(RUN_UKS_ALPHA, PpqA);
                            this->m_dExc += this->calcExc(RUN_UKS_BETA,  PpqB);
                        } else {
                            this->m_dExc = dfXCFunctional.getEnergy();
                            if (this->enableGrimmeDispersion_ == true) {
                                this->E_disp_ = dfXCFunctional.getGrimmeDispersionEnergy();
                            }
                        }
                        
                        this->K_term_ = this->E_KA_ + this->E_KB_;
                    }
                }

                if (dfXCFunctional.isHybridFunctional() == true) {
                    this->E_KA_ = this->calcK(RUN_UKS_ALPHA, PpqA);
                    this->E_KB_ = this->calcK(RUN_UKS_BETA,  PpqB);
                    this->K_term_ = this->E_KA_ + this->E_KB_;
                }
            }
            break;
            
        case METHOD_ROKS:
            {
                if (dfXCFunctional.getXcType() != DfXCFunctional::HF) {
                    if (this->m_bIsXCFitting == true) {
                        const VectorType Eps_closed = 0.5 * this->getEps<VectorType>(RUN_ROKS_CLOSED);
                        const VectorType Eps_open = this->getEps<VectorType>(RUN_ROKS_CLOSED);
                        const VectorType EpsA = Eps_closed + Eps_open;
                        const VectorType EpsB = Eps_closed;
                        this->m_dExc  = this->calcExc_DIRECT<DfOverlapType, SymmetricMatrixType, VectorType>(PpqA, EpsA);
                        this->m_dExc += this->calcExc_DIRECT<DfOverlapType, SymmetricMatrixType, VectorType>(PpqB, EpsB);
                    } else {
                        if (this->XC_engine_ != XC_ENGINE_GRID) {
                            this->m_dExc  = this->calcExc(RUN_ROKS_ALPHA, PpqA);
                            this->m_dExc += this->calcExc(RUN_ROKS_BETA,  PpqB);
                        } else {
                            this->m_dExc = dfXCFunctional.getEnergy();
                            if (this->enableGrimmeDispersion_ == true) {
                                this->E_disp_ = dfXCFunctional.getGrimmeDispersionEnergy();
                            }
                        }
                    }
                }

                if (dfXCFunctional.isHybridFunctional() == true) {
                    this->E_KA_ = this->calcK(RUN_UKS_ALPHA, PpqA);
                    this->E_KB_ = this->calcK(RUN_UKS_BETA,  PpqB);
                    this->K_term_ = this->E_KA_ + this->E_KB_;
                }
            }
            break;

        default:
            CnErr.abort("DfTotalEnergy::exec(): unknown type.");
            break;
        }
    }

    // 表示
    this->output();
}


// 全密度行列を用意する
template<typename SymmetricMatrixType>
SymmetricMatrixType DfTotalEnergy::getPpq(const METHOD_TYPE methodType)
{
    SymmetricMatrixType Ppq;

    switch (methodType) {
    case METHOD_RKS:
        {
            Ppq = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration);
            if (Ppq.getNumOfRows() != this->m_nNumOfAOs || Ppq.getNumOfCols() != this->m_nNumOfAOs) {
                CnErr.abort("DfTotalEnergy", "getPpq", "", "program error");
            }
        }
        break;

    case METHOD_UKS:
        {
            // Ppq matrix for alpha spin
            Ppq = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA, this->m_nIteration);
            if (Ppq.getNumOfRows() != this->m_nNumOfAOs || Ppq.getNumOfCols() != this->m_nNumOfAOs) {
                CnErr.abort("DfTotalEnergy", "getPpq", "", "program error");
            }
            // Ppq matrix for beta spin
            SymmetricMatrixType Ppq_b;
            Ppq_b = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_UKS_BETA, this->m_nIteration);
            if (Ppq_b.getNumOfRows() != this->m_nNumOfAOs || Ppq_b.getNumOfCols() != this->m_nNumOfAOs) {
                CnErr.abort("DfTotalEnergy", "getPpq", "", "program error");
            }

            // Ppq matrix for alpha + beta spin
            Ppq += Ppq_b;
        }
        break;

    case METHOD_ROKS:
        {
            SymmetricMatrixType P2pq;
            {
                Ppq = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSED, this->m_nIteration);
                if (Ppq.getNumOfRows() != this->m_nNumOfAOs || Ppq.getNumOfCols() != this->m_nNumOfAOs) {
                    CnErr.abort("DfTotalEnergy", "DfTotalEnergyMain", "", "program error");
                }
            }
            {
                P2pq = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN, this->m_nIteration);
                if (P2pq.getNumOfRows() != this->m_nNumOfAOs || P2pq.getNumOfCols() != this->m_nNumOfAOs) {
                    CnErr.abort("DfTotalEnergy", "DfTotalEnergyMain", "", "program error");
                }
            }
            Ppq += P2pq;
        }
        break;
        
    default:
        CnErr.abort("DfTotalEnergy::getPpq(): unknown method");
        break;
    }

    return Ppq;
}

// 全電子密度分布を用意する
template<typename VectorType>
VectorType DfTotalEnergy::getRho(const METHOD_TYPE methodType)
{
    VectorType rho;

    switch (methodType) {
    case METHOD_RKS:
        rho.load(this->getRhoPath(RUN_RKS, this->m_nIteration));
        break;
        
    case METHOD_UKS: // go down
    case METHOD_ROKS:
        {
            VectorType rho_a, rho_b;
            rho_a.load(this->getRhoPath(RUN_UKS_ALPHA, this->m_nIteration));
            rho_b.load(this->getRhoPath(RUN_UKS_BETA, this->m_nIteration));
            rho = rho_a + rho_b;
        }
        break;
        
    default:
        CnErr.abort(" DfTotalEnergy::getRho error.\n");
        break;
    }

    return rho;
}

template<typename VectorType>
VectorType DfTotalEnergy::getEps(const RUN_TYPE runType)
{
    VectorType E;

    switch (runType) {
    case RUN_RKS:
        E.load("fl_Work/fl_Vct_Epsilon");
        if (static_cast<TlVector::size_type>(this->numOfAuxXC_) != E.getSize()) {
            this->logger(TlUtils::format("dimension_eps = %d\n", E.getSize()));
            this->logger(TlUtils::format("numOfAuxXC_ = %d\n", this->numOfAuxXC_));
            this->logger("DfTotalEnergy dimension is not consistency, but continue\n");
        }
        break;

    case RUN_UKS_ALPHA:
        E.load("fl_Work/fl_Vct_Epsilona");
        if (static_cast<TlVector::size_type>(this->numOfAuxXC_) != E.getSize()) {
            this->logger(TlUtils::format("dimension of epsa = %d\n", E.getSize()));
            this->logger(TlUtils::format("numOfAuxXC_   = %d\n", this->numOfAuxXC_));
            this->logger("DfTotalEnergy dimension is not consistency, but continue\n");
        }
        break;

    case RUN_UKS_BETA:
        E.load("fl_Work/fl_Vct_Epsilonb");
        if (static_cast<TlVector::size_type>(this->numOfAuxXC_) != E.getSize()) {
            this->logger(TlUtils::format("dimension of epsb = %d\n", E.getSize()));
            this->logger(TlUtils::format("numOfAuxXC_   = %d\n", this->numOfAuxXC_));
            this->logger("DfTotalEnergy dimension is not consistency, but continue\n");
        }
        break;

    default:
        CnErr.abort("DfTotalEnergy::getEps() unknown type.");
        break;
    }
    
    return E;
}

template<typename SymmetricMatrixType>
double DfTotalEnergy::calculate_E_WITH_DIRECT(const SymmetricMatrixType& D)
{
    SymmetricMatrixType tmpEpq;
    if (this->m_nIteration > 1) {
        const std::string fname = "fl_Work/fl_Mtr_Epqtmp" + TlUtils::xtos(this->m_nIteration -1);
        tmpEpq.load(fname);
    } else {
        tmpEpq.load(this->getHpqMatrixPath());
        if (tmpEpq.getNumOfRows() != this->m_nNumOfAOs || tmpEpq.getNumOfCols() != this->m_nNumOfAOs) {
            CnErr.abort("DfTotalEnergy", "calculate_oep+JRR+Exc", "", "program error");
        }
    }

    // add dummy charge
    //{
        //Fl_Geometry geom(Fl_Geometry::getDefaultFileName());
        //const int dNumOfDummyAtoms = geom.getDummyatom();
        //if (dNumOfDummyAtoms != 0) {
    if (this->m_nNumOfDummyAtoms > 0) {
        const int chgextra_number = (*(this->pPdfParam_))["charge-extrapolate-number"].getInt();
        if (((this->m_nIteration == 1) && (chgextra_number == 0)) ||
            ((this->m_nIteration != 1) && (this->m_nIteration <= (chgextra_number +1)))) {
            // for the case of gradually addition of dummy-charge term
            SymmetricMatrixType Hpq2;
            Hpq2.load(this->getHpq2MatrixPath());
            tmpEpq += Hpq2;
        }
    }
    //}

    {
        SymmetricMatrixType E;
        E.load("fl_Work/fl_Mtr_Epqtmp" + TlUtils::xtos(m_nIteration));

        tmpEpq += E;
    }
    tmpEpq.save("fl_Work/fl_Mtr_Epqtmp" +  TlUtils::xtos(m_nIteration));

    return tmpEpq.dot(D).sum();
}

// energy for one electron part
template<typename SymmetricMatrixType>
double DfTotalEnergy::calcOneElectronPart(const SymmetricMatrixType& D)
{
    SymmetricMatrixType Hpq;
    Hpq.load(this->getHpqMatrixPath());
    if (Hpq.getNumOfRows() != this->m_nNumOfAOs || Hpq.getNumOfCols() != this->m_nNumOfAOs) {
        CnErr.abort("DfTotalEnergy", "calculate_one_electron_part", "", "program error");
    }

    // add dummy charge
    if (this->m_nNumOfDummyAtoms > 0) {
        const int chgextra_number = (*(this->pPdfParam_))["charge-extrapolate-number"].getInt();
        
        SymmetricMatrixType Hpq2;
        Hpq2.load(this->getHpq2MatrixPath());

        if (chgextra_number > 1) {
            const int coef = std::min(this->m_nIteration, chgextra_number);
            Hpq += Hpq2 * double(coef);
        } else {
            Hpq += Hpq2;
        }
    }

    return Hpq.dot(D).sum();
}

template<typename SymmetricMatrixType>
double DfTotalEnergy::calcJ(const SymmetricMatrixType& P)
{
    SymmetricMatrixType J = DfObject::getJMatrix<SymmetricMatrixType>(this->m_nIteration);
    return 0.5 * J.dot(P).sum();
}


/// J[rho~, rho~] (( "1/2*rho*rho'*Sab" ))の計算を行う
template<typename SymmetricMatrixType, typename VectorType>
double DfTotalEnergy::calcJRhoTildeRhoTilde(const VectorType& rho)
{
    double answer = 0.0;
    if (this->J_engine_ == J_ENGINE_RI_J) {
        SymmetricMatrixType Sab;
        Sab.load(this->getSabMatrixPath());
        answer = -0.5 * (rho * (Sab * rho));
    }
    
    return answer;
}


// energy for a part of coulomb term ( rou*Ppq*Pqa )
template<typename DfEriType, typename SymmetricMatrixType, typename VectorType>
double DfTotalEnergy::calcJRhoRhoTilde_DIRECT(const SymmetricMatrixType& D)
{
    SymmetricMatrixType J = DfObject::getJMatrix<SymmetricMatrixType>(this->m_nIteration);

    const double coef = (this->J_engine_ == J_ENGINE_RI_J) ? 1.0 : 0.5;
    
    return coef * (J.dot(D).sum());
}


// energy for xc energy term (compare myu*Ppq*Pqa with 4/3*Ex1)
template<typename DfOverlapType, typename SymmetricMatrixType, typename VectorType>
double DfTotalEnergy::calcExc_DIRECT(const SymmetricMatrixType& D, const VectorType& E)
{
    // put eps*[pqs] into B
    SymmetricMatrixType B(this->m_nNumOfAOs);

    DfOverlapType dfOverlap(this->pPdfParam_);
    dfOverlap.get_pqg(E, &B);

    return B.dot(D).sum();
}


template<typename SymmetricMatrixType>
double DfTotalEnergy::calcExc(const RUN_TYPE runType,
                              const SymmetricMatrixType& P)
{
    SymmetricMatrixType Exc = DfObject::getExcMatrix<SymmetricMatrixType>(runType, this->m_nIteration);
    const double answer = Exc.dot(P).sum();
    
    return answer;
}


template<typename SymmetricMatrixType>
double DfTotalEnergy::calcK(const RUN_TYPE runType,
                            const SymmetricMatrixType& P)
{
    double answer = 0.0;
    DfXCFunctional dfXCFunctional(this->pPdfParam_);
    if (dfXCFunctional.isHybridFunctional() == true) {
        SymmetricMatrixType K = DfObject::getHFxMatrix<SymmetricMatrixType>(runType, this->m_nIteration);
        answer = 0.5 * dfXCFunctional.getFockExchangeCoefficient() * K.dot(P).sum();
    }

    return answer;
}


template <class SymmetricMatrixType>
void DfTotalEnergy::calcRealEnergy()
{
    this->logger(" total energy --- Information related with dummy atom:\n\n");

    const double current_energy = (*this->pPdfParam_)["TEs"][this->m_nIteration].getDouble();
    
    //check whether dummy atom exists or not
    if (this->m_nNumOfDummyAtoms == 0) {
        this->logger(" no dummy atoms are found.\n");
        return;
    }

    //enegy derived from dummy charge and real energy
    double energy_from_dummy = this->calcEnergyFromDummy<SymmetricMatrixType>();
    double real_energy = current_energy - energy_from_dummy;

    this->logger(TlUtils::format(" total energy excluding dummy charge = %18.16lf\n",real_energy));
    return;
}


template<class SymmetricMatrixType>
double DfTotalEnergy::calcEnergyFromDummy()
{
    // read Density Matrix
    SymmetricMatrixType P;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        P = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration);
        break;

    default:
        abort();// program error
    }

    // part of Hpq2
    const int chargeExtrapolateNumber = (*(this->pPdfParam_))["charge-extrapolate-number"].getInt();

    SymmetricMatrixType Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrixType>();

    if (chargeExtrapolateNumber > 0) {
        const int coef = std::min(this->m_nIteration, chargeExtrapolateNumber);
        Hpq2 *= (double)coef;
    }

    const double energy_from_Hpq2 = Hpq2.dot(P).sum();
    this->logger(TlUtils::format(" energy derived from core H = %18.16lf\n", energy_from_Hpq2));

    //part of nuclear-nuclear repulsion
    const int start_number_dummy = this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;

    Fl_Geometry geom((*this->pPdfParam_)["coordinates"]);
    double energy_from_nuclear = 0.0;
    for (int i = start_number_dummy; i < this->m_nNumOfAtoms; ++i) {
        const double ci = geom.getCharge(i);
        const TlPosition pi = geom.getCoordinate(i);

        for (int j = 0; j < i; ++j) {
            const double cj = geom.getCharge(j);
            const TlPosition pj = geom.getCoordinate(j);
            const double dist = pi.distanceFrom(pj);

            energy_from_nuclear += ci * cj / dist;
        }
    }
    this->logger(TlUtils::format(" energy derived from nuclear part = %18.16lf\n", energy_from_nuclear));

    // sum
    const double energy_from_dummy = energy_from_Hpq2 + energy_from_nuclear;

    this->logger(TlUtils::format(" total energy derived from dummy atom = %18.16lf\n", energy_from_dummy));

    return energy_from_dummy;
}


#endif // DFTOTALENERGY_H


