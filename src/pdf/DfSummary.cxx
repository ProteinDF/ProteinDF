#include <cmath>
#include "DfSummary.h"
#include "CnError.h"
#include "Fl_Geometry.h"

#include "TlUtils.h"
#include "TlMath.h"
#include "TlSymmetricMatrix.h"

DfSummary::DfSummary(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
}


DfSummary::~DfSummary()
{
}


void DfSummary::exec()
{
    this->exec<TlMatrix>();
}


void DfSummary::printEigen(DfObject::RUN_TYPE runType)
{
    TlVector eigval;
    eigval.load(this->getEigenvaluesPath(runType, this->m_nIteration));

    std::stringstream ss;
    eigval.print(ss);
    this->log_.info(ss.str());
}


void DfSummary::printAux(const TlVector& rho, const TlVector& myu, const TlVector& nyu)
{
    const int rhoSize = rho.getSize();
    const int myuSize = myu.getSize();
    const int nyuSize = nyu.getSize();
    //TlLogX& Log = TlLogX::getInstance();

    // header
    this->logger("    GTO    ATOM       SHELL                 RHO");
    if (myuSize > 0) {
        this->logger("                MYU");
        if (nyuSize > 0) {
            this->logger("                NYU");
        }
    }
    this->logger("\n");

    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_sets_k"]);
    const int dim = std::max(rhoSize, std::max(myuSize, nyuSize));
    for (int i = 0; i < dim; ++i) {
        this->logger(TlUtils::format(" %6d %-2s %6s",
                                     i +1,
                                     orbInfo.getAtomName(i).c_str(),
                                     TlOrbitalInfoObject::basisTypeNameTbl_[orbInfo.getBasisType(i)]));

        if (i < rhoSize) {
            this->logger(TlUtils::format(" %18.6lf", rho[i]));
        } else {
            this->logger("                   ");
        }

        if (i < myuSize) {
            this->logger(TlUtils::format(" %18.6lf", myu[i]));
        } else {
            this->logger("                   ");
        }

        if (i < nyuSize) {
            this->logger(TlUtils::format(" %18.6lf", nyu[i]));
        } else {
            this->logger("                   ");
        }
        this->logger("\n");
    }
}


void DfSummary::printRhoPop(const TlVector& rho)
{
    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"],
                                (*this->pPdfParam_)["basis_sets_j"]);
    const int numOfAux = this->m_nNumOfAux;

    this->logger("    GTO    ATOM       SHELL                 RHO\n\n");
    for (int i = 0; i < numOfAux; ++i) {
        this->logger(TlUtils::format(" %6d %-2s %6s",
                                     i+1, 
                                     orbInfo.getAtomName(i).c_str(),
                                     TlOrbitalInfoObject::basisTypeNameTbl_[orbInfo.getBasisType(i)]));

        //const double popu = std::pow((M_PI / (2.0 * RGTO.getExponent(Tdens.getcgtonum(k),0))), 0.25) * rho[i];
        const double popu = std::pow((M_PI / (2.0 * orbInfo.getExponent(i, 0))), 0.25) * rho[i];
        this->logger(TlUtils::format(" %18.6lf\n", popu));
    }
    this->logger("\n");
}
