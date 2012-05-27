#include <cmath>
#include "DfSummary.h"
#include "CnError.h"
#include "Fl_Geometry.h"
#include "Fl_Gto_Density.h"
#include "Fl_Tbl_Density.h"

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

    const Fl_Tbl_Density Tdens(Fl_Geometry((*this->pPdfParam_)["coordinates"]));
    int previous_blflag = -1;
    const int dim = std::max(rhoSize, std::max(myuSize, nyuSize));
    for (int k = 0, blflag = Tdens.getInpAtomnum(k) +1; k < dim; ++k) {
        if (blflag != Tdens.getInpAtomnum(k)+1) {
            blflag = Tdens.getInpAtomnum(k)+1;
        }

        if (previous_blflag != blflag) {
            this->logger(TlUtils::format("  %5ld   %-2s   %3ld   %-9s ",
                                         k+1, Tdens.getInpAtomtype(k).c_str(),
                                         Tdens.getInpAtomnum(k)+1,
                                         Tdens.getNote1(k).c_str()));
            previous_blflag = blflag;
        } else {
            this->logger(TlUtils::format("  %5ld              %-9s ",
                                         k+1, Tdens.getNote1(k).c_str()));
        }

        if (k < rhoSize) {
            this->logger(TlUtils::format(" %18.6lf", rho[k]));
        } else {
            this->logger("                   ");
        }

        if (k < myuSize) {
            this->logger(TlUtils::format(" %18.6lf", myu[k]));
        } else {
            this->logger("                   ");
        }

        if (k < nyuSize) {
            this->logger(TlUtils::format(" %18.6lf", nyu[k]));
        } else {
            this->logger("                   ");
        }
        this->logger("\n");
    }
}


void DfSummary::printRhoPop(const TlVector& rho)
{
    const Fl_Tbl_Density Tdens((*this->pPdfParam_)["coordinates"]);
    const Fl_Gto_Density RGTO;
    const int numOfAux = this->m_nNumOfAux;

    this->logger("    GTO    ATOM       SHELL                 RHO\n\n");
    int previous_blflag = -1;
    for (int  k = 0, blflag = Tdens.getInpAtomnum(k)+1; k < numOfAux; ++k) {
        if (blflag != Tdens.getInpAtomnum(k)+1) {
            blflag = Tdens.getInpAtomnum(k)+1;
        }

        if ('s' == Tdens.getNote1(k)[0]) {
            if (previous_blflag != blflag) {
                this->logger(TlUtils::format("  %5ld   %-2s   %3ld   %-9s ",
                                             k+1, Tdens.getInpAtomtype(k).c_str(),
                                             Tdens.getInpAtomnum(k)+1,
                                             Tdens.getNote1(k).c_str()));
                previous_blflag = blflag;
            } else {
                this->logger(TlUtils::format("  %5ld              %-9s ",
                                             k+1, Tdens.getNote1(k).c_str()));
            }
            const double popu = std::pow((M_PI / (2.0 * RGTO.getExponent(Tdens.getcgtonum(k),0))), 0.25) * rho[k];
            this->logger(TlUtils::format(" %18.6lf\n", popu));
        }
    }
    this->logger("\n");
}
