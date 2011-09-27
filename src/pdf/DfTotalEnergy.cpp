#include <cmath>

#include "DfTotalEnergy.h"
#include "Common.h"
#include "DfEri.h"
#include "DfOverlap.h"
#include "DfEri2.h"
#include "DfTwoElectronIntegral.h"
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
    TlVector enevec;

    if (this->m_nIteration != 1) {
        enevec.load("fl_Work/fl_Vct_Energy");
    }

    enevec.push_back(E_Total);
    assert(static_cast<int>(enevec.getSize()) == this->m_nIteration);

    enevec.save("fl_Work/fl_Vct_Energy");
}


// energy for a part of coulomb term ( rho*Ppq*Pqa )
double DfTotalEnergy::calculate_J_Rho_Rhotilda_WITH_FILE(const TlSymmetricMatrix& D, const TlVector& R)
{
    Fl_Int_Pqa IntPqa;  // coulomb three index integrals(tmp)
    //IntPqa.open( "fl_Temp", IntPqa.getfilename(), "read" );

    TlSymmetricMatrix B(this->m_nNumOfAOs);
    {
        int icount, npqalp, nalp, alpha;
        int*   indalp = new int   [ MAXNA   ];
        int*   npq    = new int   [ MAXNA   ];
        int*   indp   = new int   [ MAXNPQA ];
        int*   indq   = new int   [ MAXNPQA ];
        int*   indpq  = new int   [ MAXNPQA ];
        int*   indgam = new int   [ MAXNA   ];
        double* pqA    = new double [ MAXNPQA ];
        double* pqG    = new double [ MAXNPQA ];

        do {
            IntPqa.read(&icount, &npqalp, &nalp, indalp, npq, indp, indq, pqA);
            int ind = 0;
            for (int i=0; i<nalp; i++) {
                alpha = indalp[i];
                for (int j=0; j<npq[i]; j++) {
                    indpq[ind] = indp[ind] * (indp[ind]+1)/2 + indq[ind];
                    B(indp[ind], indq[ind]) += R[alpha] * pqA[ind];
                    ind++;
                }
            }
        } while (icount != 0);

        delete [] indalp;
        delete [] npq;
        delete [] indp;
        delete [] indq;
        delete [] indpq;
        delete [] indgam;
        delete [] pqA;
        delete [] pqG;
    }

    double E_J_Rou_Routilda = 0.0;
    for (int i = 0; i < this->m_nNumOfAOs; ++i) {
        // case: i != j
        for (int j = 0; j < i; ++j) {
            E_J_Rou_Routilda += 2.0 * D(i,j) * B(i,j);
        }
        // case: i == j
        E_J_Rou_Routilda += D(i,i) * B(i,i);
    }

    //this->logger(TlUtils::format("    Energy of a part of Coulomb Energy Part = %28.10lf\n", E_J_Rou_Routilda ));

    return E_J_Rou_Routilda;
}


// energy for xc energy term (compaire myu*Ppq*Pqa with 4/3*Ex1)
double DfTotalEnergy::calculate_Exc_WITH_FILE(const TlSymmetricMatrix& D, const TlVector& E)
{
    //LogX& Log = TlLogX::getInstance();

    Fl_Int_Pqg IntPqg;  // xc three index integrals(tmp)

    //IntPqg.open("fl_Temp", IntPqg.getfilename(), "read");

    TlSymmetricMatrix B(this->m_nNumOfAOs);
    {
        int icount, npqgam, ngam, gamma;
        int*   indalp = new int   [ MAXNA   ];
        int*   npq    = new int   [ MAXNA   ];
        int*   indp   = new int   [ MAXNPQA ];
        int*   indq   = new int   [ MAXNPQA ];
        int*   indpq  = new int   [ MAXNPQA ];
        int*   indgam = new int   [ MAXNA   ];
        double* pqA    = new double [ MAXNPQA ];
        double* pqG    = new double [ MAXNPQA ];
        do {
            IntPqg.read(&icount, &npqgam, &ngam, indgam, npq, indp, indq, pqG);

            int ind=0;
            for (int i=0; i<ngam; i++) {
                gamma = indgam[i];
                for (int j=0; j<npq[i]; j++) {
                    indpq[ind] = indp[ind] * (indp[ind]+1)/2 + indq[ind];
                    B(indp[ind], indq[ind]) += E[gamma] * pqG[ind];
                    ind++;
                }
            }
        } while (icount != 0);

        delete [] indalp;
        delete [] npq;
        delete [] indp;
        delete [] indq;
        delete [] indpq;
        delete [] indgam;
        delete [] pqA;
        delete [] pqG;
    }

    double E_Exc = 0.0;
    for (int i = 0; i < this->m_nNumOfAOs; ++i) {
        // case: i != j
        for (int j = 0; j < i; ++j) {
            E_Exc += 2.0 * D(i,j) * B(i,j);
        }
        // case: i == j
        E_Exc += D(i,i) * B(i,i);
    }

    //this->logger(TlUtils::format("    Energy of Exc Part                      = %28.10lf\n", E_Exc));

    return E_Exc;
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

