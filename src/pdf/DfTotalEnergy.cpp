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

#include <cmath>

#include "DfTotalEnergy.h"
#include "DfEriX.h"
#include "DfOverlapX.h"
#include "DfXCFunctional.h"

#include "Fl_Geometry.h"
#include "TlUtils.h"
#include "CnError.h"

const double DfTotalEnergy::TOO_SMALL = 1.0E-16;
double DfTotalEnergy::m_dNuclearRepulsion = 0.0; // 核-核反発

DfTotalEnergy::DfTotalEnergy(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam),
      m_dE_OneElectronPart(0.0),
      J_term_(0.0), m_dE_J_Rho_RhoTilde(0.0), m_dE_J_RhoTilde_RhoTilde(0.0),
      m_dExc(0.0), K_term_(0.0),
      E_KA_(0.0), E_KB_(0.0),
      m_dE_NuclearRepulsion(0.0), m_dE_OEP_JRR_Exc(0.0),
      E_disp_(0.0)
{
}


DfTotalEnergy::~DfTotalEnergy()
{
}


void DfTotalEnergy::exec()
{
    this->exec_template<DfOverlapX, DfEriX, TlSymmetricMatrix, TlVector>();
}


void DfTotalEnergy::output()
{
    double E_Total  = 0.0;

    // 表示
    if ((this->m_bMemorySave == false) && (this->m_bDiskUtilization == false)) {
        E_Total += this->m_dE_OEP_JRR_Exc;
        E_Total += this->m_dE_J_RhoTilde_RhoTilde;
        E_Total += this->m_dE_NuclearRepulsion;

        this->logger("------------------------------------------------\n");
        this->logger(TlUtils::format(" Ts+Vn+J[Rho~,Rho~]+Exc = %28.16lf\n", this->m_dE_OEP_JRR_Exc));
        this->logger(TlUtils::format(" J[Rho~,Rho~]           = %28.16lf\n", this->m_dE_J_RhoTilde_RhoTilde));
        this->logger(TlUtils::format(" Enuclei                = %28.16lf\n", this->m_dE_NuclearRepulsion));
        this->logger(TlUtils::format(" TE                     = %28.16lf\n", E_Total));
        this->logger("------------------------------------------------\n");
    } else {
        this->log_.info("------------------------------------------------");

        this->log_.info(TlUtils::format(" Ts+Vn          = %28.16lf\n", this->m_dE_OneElectronPart));
        E_Total += this->m_dE_OneElectronPart;

        switch (this->J_engine_) {
        case J_ENGINE_RI_J:
            this->log_.info(TlUtils::format(" E_J[Rho, Rho~] = %28.16lf\n", this->m_dE_J_Rho_RhoTilde));
            E_Total += this->m_dE_J_Rho_RhoTilde;
            this->log_.info(TlUtils::format(" E_J[Rho~,Rho~] = %28.16lf\n", this->m_dE_J_RhoTilde_RhoTilde));
            E_Total += this->m_dE_J_RhoTilde_RhoTilde;
            break;
            
        case J_ENGINE_CONVENTIONAL:
        case J_ENGINE_CD:
            this->log_.info(TlUtils::format(" E_J            = %28.16lf\n", this->J_term_));
            E_Total += this->J_term_;
            break;

        default:
            break;
        }

        this->log_.info(TlUtils::format(" E_xc(pure)     = %28.16lf\n", this->m_dExc));
        E_Total += this->m_dExc;

        if (this->enableGrimmeDispersion_ == true) {
            this->log_.info(TlUtils::format(" E_cx(+disp.)   = %28.16lf\n", E_Total + this->E_disp_));
        }

        switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->log_.info(TlUtils::format(" E_K            = %28.16lf\n", this->K_term_));
            break;

        case METHOD_UKS:
            this->log_.info(TlUtils::format(" E_K            = %28.16lf\n", this->K_term_));
            this->log_.info(TlUtils::format("   E_K(alpha)   = %28.16lf\n", this->E_KA_));
            this->log_.info(TlUtils::format("   E_K(beta)    = %28.16lf\n", this->E_KB_));
            break;

        case METHOD_ROKS:
            break;

        default:
            this->log_.critical("program error");
            break;
        }
        E_Total += this->K_term_;

        this->log_.info(TlUtils::format(" E_nuclei       = %28.16lf\n", this->m_dE_NuclearRepulsion));
        E_Total += this->m_dE_NuclearRepulsion;

        this->log_.info(TlUtils::format(" TE             = %28.16lf\n", E_Total));
        this->log_.info("------------------------------------------------");
        // this->logger("------------------------------------------------\n");
        // this->logger(TlUtils::format(" Ts+Vn        = %28.16lf\n", this->m_dE_OneElectronPart));
        // this->logger(TlUtils::format(" J[Rho, Rho~] = %28.16lf\n", this->m_dE_J_Rho_RhoTilde));
        // if (this->J_engine_ == J_ENGINE_RI_J) {
        //     this->logger(TlUtils::format(" J[Rho~,Rho~] = %28.16lf\n", this->m_dE_J_RhoTilde_RhoTilde));
        // }
        // this->logger(TlUtils::format(" Exc          = %28.16lf\n", this->m_dExc));
        // this->logger(TlUtils::format(" Enuclei      = %28.16lf\n", this->m_dE_NuclearRepulsion));
        // this->logger(TlUtils::format(" TE           = %28.16lf\n", E_Total));
        // this->logger("------------------------------------------------\n");
    }

    std::cout << TlUtils::format(" %3d th TE = %18.16lf", this->m_nIteration, E_Total) << std::endl;
    this->write_total_energy(E_Total);
}


// energy for nuclear repulsion
double DfTotalEnergy::calculate_energy_nuclear_repulsion()
{
    double E_nuclear_repulsion = 0.0;

    if (std::fabs(this->m_dNuclearRepulsion) > TOO_SMALL) {
        this->logger(" energy of nuclear repulstion has already been calculated.\n");
        E_nuclear_repulsion = this->m_dNuclearRepulsion;
    } else {
        // read nuclear charge
        const Fl_Geometry geom((*this->pPdfParam_)["coordinates"]);

        //calculate nuclear repulsion
        for (int i = 0; i < this->m_nNumOfAtoms; ++i) {
            const double ci = geom.getCharge(i);
            const TlPosition pi = geom.getCoordinate(i);

            for (int j = i + 1; j < this->m_nNumOfAtoms; ++j) {
                const double cj = geom.getCharge(j);
                const TlPosition pj = geom.getCoordinate(j);

                //double distance =  sqrt((x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) + (z_i - z_j)*(z_i - z_j));
                const double distance = pi.distanceFrom(pj);
                E_nuclear_repulsion += ci * cj / distance;
            }
        }

        // stored to static variable
        this->m_dNuclearRepulsion = E_nuclear_repulsion;
    }

    return E_nuclear_repulsion;
}


void DfTotalEnergy::write_total_energy(const double E_Total) const
{
    (*this->pPdfParam_)["TEs"][this->m_nIteration] = E_Total;
}

// total energy including dummy charge
void DfTotalEnergy::calculate_real_energy()
{
    this->calcRealEnergy<TlSymmetricMatrix>();
}


DfEriX* DfTotalEnergy::getDfEriX() const
{
    DfEriX* pDfEri = new DfEriX(this->pPdfParam_);

    return pDfEri;
}


DfOverlapX* DfTotalEnergy::getDfOverlapX() const
{
    DfOverlapX* pDfOverlap = new DfOverlapX(this->pPdfParam_);

    return pDfOverlap;
}

