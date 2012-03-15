#include <cmath>

#include "DfTotalEnergy.h"
#include "DfEri.h"
#include "DfOverlap.h"
#include "DfEri2.h"
#include "DfXCFunctional.h"

#include "Fl_Geometry.h"
#include "Fl_Int_Pqa.h"
#include "TlUtils.h"
#include "CnError.h"

const double DfTotalEnergy::TOO_SMALL = 1.0E-16;
double DfTotalEnergy::m_dNuclearRepulsion = 0.0; // 核-核反発

DfTotalEnergy::DfTotalEnergy(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
}


DfTotalEnergy::~DfTotalEnergy()
{
}


void DfTotalEnergy::exec()
{
    this->exec_template<DfOverlap, DfEri, TlSymmetricMatrix, TlVector>();
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
        E_Total += this->m_dE_OneElectronPart;
        E_Total += this->m_dE_J_Rho_RhoTilde;
        if (this->isRI_J_ == true) {
            E_Total += this->m_dE_J_RhoTilde_RhoTilde;
        }
        E_Total += this->m_dExc;
        E_Total += this->m_dE_NuclearRepulsion;

        this->logger("------------------------------------------------\n");
        this->logger(TlUtils::format(" Ts+Vn        = %28.16lf\n", this->m_dE_OneElectronPart));
        this->logger(TlUtils::format(" J[Rho, Rho~] = %28.16lf\n", this->m_dE_J_Rho_RhoTilde));
        if (this->isRI_J_ == true) {
            this->logger(TlUtils::format(" J[Rho~,Rho~] = %28.16lf\n", this->m_dE_J_RhoTilde_RhoTilde));
        }
        this->logger(TlUtils::format(" Exc          = %28.16lf\n", this->m_dExc));
        this->logger(TlUtils::format(" Enuclei      = %28.16lf\n", this->m_dE_NuclearRepulsion));
        this->logger(TlUtils::format(" TE           = %28.16lf\n", E_Total));
        if (this->enableGrimmeDispersion_ == true) {
            this->logger(TlUtils::format(" TE(+disp.)   = %28.16lf\n", E_Total + this->E_disp_));
        }
        this->logger("------------------------------------------------\n");
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
        Fl_Geometry geom(Fl_Geometry::getDefaultFileName());

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
    (*this->pPdfParam_)["TE"][this->m_nIteration] = E_Total;
}

// total energy including dummy charge
void DfTotalEnergy::calculate_real_energy()
{
    this->calcRealEnergy<TlSymmetricMatrix>();
}


DfEri* DfTotalEnergy::getDfEri() const
{
    DfEri* pDfEri = new DfEri(this->pPdfParam_);

    return pDfEri;
}


DfOverlap* DfTotalEnergy::getDfOverlap() const
{
    DfOverlap* pDfOverlap = new DfOverlap(this->pPdfParam_);

    return pDfOverlap;
}

