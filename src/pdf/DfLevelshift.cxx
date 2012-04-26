#include "DfLevelshift.h"
#include "CnError.h"
#include "Fl_Tbl_Fragment.h"
#include "TlUtils.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"

DfLevelshift::DfLevelshift(TlSerializeData* pPdfParam, int num_iter)
    : DfObject(pPdfParam)
{
}

DfLevelshift::~DfLevelshift()
{
}

void DfLevelshift::DfLshiftMain()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main(RUN_RKS, this->m_nIteration);
        break;

    case METHOD_UKS:
        this->main(RUN_UKS_ALPHA, this->m_nIteration);
        this->main(RUN_UKS_BETA,  this->m_nIteration);
        break;

    case METHOD_ROKS:
        this->main(RUN_ROKS, this->m_nIteration);
        break;

    default:
        CnErr.abort();
        break;
    }
}

// for extended QCLO method
void DfLevelshift::DfLshiftQclo(const std::string& fragname, int norbcut)
{
    //this->number_mo_basis = norbcut;
    this->m_nNumOfMOs = norbcut;

    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main(RUN_RKS, this->m_nIteration, fragname);
        break;

    case METHOD_UKS:
        this->main(RUN_UKS_ALPHA, this->m_nIteration, fragname);
        this->main(RUN_UKS_BETA,  this->m_nIteration, fragname);
        break;

    case METHOD_ROKS:
        this->main(RUN_ROKS, this->m_nIteration, fragname);
        break;

    default:
        CnErr.abort();
        break;
    }
}

void DfLevelshift::main(const RUN_TYPE runType, int iteration, const std::string& fragname, bool bPdfQcloMode)
{
    TlSerializeData& pdfParam = *(this->pPdfParam_);
    
    // construct level shift matrix for "F' matrix"
    const double ls_closed_mo = pdfParam["level-shift/ls-closed-mo"].getDouble();
    const double ls_open_mo = pdfParam["level-shift/ls-open-mo"].getDouble();
    const double ls_virtual_mo = pdfParam["level-shift/ls-virtual-mo"].getDouble();

    const double delta_group_closed = pdfParam["level-shift/delta-group-closed"].getDouble();
    const double delta_group_open = pdfParam["level-shift/delta-group-open"].getDouble();
    const double delta_group_virtual = pdfParam["level-shift/delta-group-virtual"].getDouble();

    this->log_.info("construct Level Shift Operator.");

    this->log_.info("RKS,UKS calculation");
    this->log_.info(TlUtils::format(" level shift for closed  MO = %8.2lf", ls_closed_mo));
    this->log_.info(TlUtils::format(" level shift for open    MO = %8.2lf", ls_open_mo));
    this->log_.info(TlUtils::format(" level shift for virtual MO = %8.2lf\n", ls_virtual_mo));
    this->log_.info(TlUtils::format(" level shift between closed  MO = %8.2lf", delta_group_closed));
    this->log_.info(TlUtils::format("    level shift between open    MO = %8.2lf", delta_group_open));
    this->log_.info(TlUtils::format("    level shift between virtual MO = %8.2lf", delta_group_virtual));

    // prepear "beta" vector (level-shift values)
    const int norbcut = this->m_nNumOfMOs;
    TlVector beta(norbcut); // b
    {
        TlVector vOcc;
        vOcc.load(DfObject::getOccupationPath(runType));
        assert(vOcc.getSize() == norbcut);

        double shift_closed  = 0.0;
        double shift_open    = 0.0;
        double shift_virtual = 0.0;

        int orb_id = 0;
        for (int k = 0; k < norbcut; k++) {
            if (bPdfQcloMode == true) {
                int frag_id = -1;
                Fl_Tbl_Fragment Tfrag(Fl_Geometry((*this->pPdfParam_)["coordinates"]));
                if (Tfrag.getFragment(k) != frag_id) {
                    continue;
                }
            }

            if (std::fabs(vOcc[k]) < 1.0E-10) {
                beta[orb_id] += ls_virtual_mo + shift_virtual;
                shift_virtual += delta_group_virtual;
            } else if (fabs(vOcc[k] - 2.0) < 1.0E-10) {
                beta[orb_id] += ls_closed_mo + shift_closed;
                shift_closed  += delta_group_closed;
            } else if (fabs(vOcc[k] - 1.0) < 1.0E-10) {
                beta[orb_id] += ls_open_mo + shift_open;
                shift_open  += delta_group_open;
            } else {
                CnErr.abort("DfLevelshift", "", "", "occ is illegal");
            }

            orb_id++;
        }
    }

    // prepar "C'" matrix
    TlMatrix Cprime(this->m_nNumOfMOs, this->m_nNumOfMOs);
    {
        if (iteration == 1 &&
            ((this->initialGuessType_ != GUESS_LCAO) && (this->initialGuessType_ != GUESS_HUCKEL))) {
            this->log_.info("construct level shift operator with fukue's Rou at iteration==1");
            this->log_.info("level shift value is simply add to F' matrix (guess Rou is solution of DFT)");

            // construct "C'" matrix for initial Rou
            const int numOfMOs = this->m_nNumOfMOs;
            for (int k = 0; k < numOfMOs; ++k) {
                Cprime(k, k) = 1.0;
            }
        } else {
            // "read previous C' matrix"
            Cprime = DfObject::getCprimeMatrix<TlMatrix>(runType, iteration -1);

            if (Cprime.getNumOfRows() != this->m_nNumOfMOs || Cprime.getNumOfCols() != this->m_nNumOfMOs) {
                this->log_.warn(TlUtils::format("rowDim of previous C' matrix = %d", Cprime.getNumOfRows()));
                this->log_.warn(TlUtils::format("colDIm of previous C' matrix = %d", Cprime.getNumOfCols()));
                this->log_.warn(TlUtils::format("number_mo_basis = %d", this->m_nNumOfMOs));
                this->log_.warn("DfLevelshift dimension is not consistency, but continue");
            }
        }
    }

    // calculate "C' * beta * C'^dagger" = "C' * (C' * beta)^dagger"
    {
        // calc. "C' * beta"
        const int numOfMOs = this->m_nNumOfMOs;
        TlMatrix Cprime_beta(numOfMOs, numOfMOs);
        for (int i = 0; i < numOfMOs; ++i) {
            for (int j = 0; j < numOfMOs; ++j) {
                Cprime_beta(i, j) = Cprime(i, j) * beta[j];
            }
        }

        // (C' * beta)^dagger
        Cprime_beta.transpose();

        // calc. "C' * (C'*beta)^dagger"
        Cprime *= Cprime_beta;
    }

    // "read F' matrix"
    TlSymmetricMatrix Fprime;
    {
        Fprime = DfObject::getFprimeMatrix<TlSymmetricMatrix>(runType, iteration);

        if (Fprime.getNumOfRows() != this->m_nNumOfMOs || Fprime.getNumOfCols() != this->m_nNumOfMOs) {
            this->log_.warn(TlUtils::format("rowDim of F' matrix = %d", Fprime.getNumOfRows()));
            this->log_.warn(TlUtils::format("colDIm of F' matrix = %d", Fprime.getNumOfCols()));
            this->log_.warn(TlUtils::format("number_mo_basis = %d", this->m_nNumOfMOs));
            this->log_.warn("DfLevelshift dimension is not consistency, but continue");
        }
    }

    // construct "F' + Level shift matrix "
    Fprime += Cprime;

    // write "shifted F' matrix"
    DfObject::saveFprimeMatrix(runType, iteration, Fprime);
}
