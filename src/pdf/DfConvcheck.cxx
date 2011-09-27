#include <cstdlib>
#include <fstream>
#include <cmath>

#include "CnError.h"
#include "DfConvcheck.h"
#include "TlUtils.h"
#include "TlSymmetricMatrix.h"

DfConvcheck::DfConvcheck(TlSerializeData* pPdfParam, int num_iter)
    : DfObject(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;
    this->convergence_type = pdfParam["model"]["convergence/type"].getStr();

    this->threshold_cri     = pdfParam["model"]["convergence/threshold"].getDouble();
    this->threshold_cri_ene = pdfParam["model"]["convergence/threshold-energy"].getDouble();

    this->dev_sd = 0.0;
    this->dev_dm = 0.0;
    this->dev_te = 0.0;
    this->dev_ks = 0.0;
    this->dev_cd = 0.0;
    this->dev_xc = 0.0;
    this->dev_xa = 0.0;
}


DfConvcheck::~DfConvcheck()
{
}

void DfConvcheck::DfConvcheckMain()
{
    this->converged_flag = 0;

    // no judgement for the first iteration
    if (this->m_nIteration == 1) {
        return;
    }
    this->main<TlSymmetricMatrix>(this->m_nIteration);

    // check for convergence
    double return_threshold = 0.0;
    if (this->convergence_type == "fock") {
        return_threshold = this->dev_ks;
    } else if (this->convergence_type == "density") {
        return_threshold = this->dev_dm;
    } else if (this->convergence_type == "energy") {
        return_threshold = this->dev_te;
    } else {
        return_threshold = this->dev_cd;
    }

    if ((return_threshold < this->threshold_cri) && (dev_te < this->threshold_cri_ene)) {
        this->converged_flag = 1;
    }

    this->showResults();
}


void DfConvcheck::showResults()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->logger(TlUtils::format("+ convergence informations\n"));
        this->logger(TlUtils::format("+     standard deviation of cd               = %14.4le\n", this->dev_sd));
        this->logger(TlUtils::format("+     maximum  deviation of cd coefficient   = %14.4le\n", this->dev_cd));
        this->logger(TlUtils::format("+     maximum  deviation of density matrix   = %14.4le\n", this->dev_dm));
        this->logger(TlUtils::format("+              deviation of total energy     = %14.4le\n", this->dev_te));
        this->logger(TlUtils::format("+     maximum  deviation of kohn-sham matrix = %14.4le\n", this->dev_ks));
        if (this->m_bIsXCFitting == true) {
            this->logger(TlUtils::format("+     maximum  deviation of xc coefficient   = %14.4le\n", this->dev_xc));
            this->logger(TlUtils::format("+     maximum  deviation of xa coefficient   = %14.4le", this->dev_xa));
        }
        break;
        
    case METHOD_UKS:
        this->logger("+ convergence informations\n");
        this->logger(TlUtils::format("+     standard deviation of alpha cd               = %14.4le\n", this->dev_sd_a));
        this->logger(TlUtils::format("+     standard deviation of beta  cd               = %14.4le\n", this->dev_sd_b));
        this->logger(TlUtils::format("+     maximum  deviation of alpha density matrix   = %14.4le\n", this->dev_dm_a));
        this->logger(TlUtils::format("+     maximum  deviation of beta  density matrix   = %14.4le\n", this->dev_dm_b));
        this->logger(TlUtils::format("+              deviation of total energy           = %14.4le\n", this->dev_te));
        this->logger(TlUtils::format("+     maximum  deviation of alpha kohn-sham matrix = %14.4le\n", this->dev_ks_a));
        this->logger(TlUtils::format("+     maximum  deviation of beta  kohn-sham matrix = %14.4le\n", this->dev_ks_b));
        this->logger(TlUtils::format("+     maximum  deviation of cd alpha coefficient   = %14.4le\n", this->dev_cd_a));
        this->logger(TlUtils::format("+     maximum  deviation of cd beta  coefficient   = %14.4le\n", this->dev_cd_b));
        if (this->m_bIsXCFitting == true) {
            this->logger(TlUtils::format("+     maximum  deviation of xc alpha coefficient   = %14.4le\n", this->dev_xc_a));
            this->logger(TlUtils::format("+     maximum  deviation of xc beta  coefficient   = %14.4le\n", this->dev_xc_b));
            this->logger(TlUtils::format("+     maximum  deviation of xa alpha coefficient   = %14.4le\n", this->dev_xa_a));
            this->logger(TlUtils::format("+     maximum  deviation of xa beta  coefficient   = %14.4le\n", this->dev_xa_b));
        }
        break;

    case METHOD_ROKS:
        this->logger("+ convergence informations\n");
        this->logger(TlUtils::format("+     standard deviation of alpha cd                     = %14.4le\n", this->dev_sd_a));
        this->logger(TlUtils::format("+     standard deviation of beta  cd                     = %14.4le\n", this->dev_sd_b));
        this->logger(TlUtils::format("+     maximum  deviation of closed part density matrix   = %14.4le\n", this->dev_dm_c));
        this->logger(TlUtils::format("+     maximum  deviation of open   part density matrix   = %14.4le\n", this->dev_dm_o));
        this->logger(TlUtils::format("+              deviation of total energy                 = %14.4le\n", this->dev_te));
        this->logger(TlUtils::format("+     maximum  deviation of kohn-sham matrix             = %14.4le\n", this->dev_ks));
        this->logger(TlUtils::format("+     maximum  deviation of cd alpha coefficient         = %14.4le\n", this->dev_cd_a));
        this->logger(TlUtils::format("+     maximum  deviation of cd beta  coefficient         = %14.4le\n", this->dev_cd_b));
        if (this->m_bIsXCFitting == true) {
            this->logger(TlUtils::format("+     maximum  deviation of xc alpha coefficient         = %14.4le\n", this->dev_xc_a));
            this->logger(TlUtils::format("+     maximum  deviation of xc beta  coefficient         = %14.4le\n", this->dev_xc_b));
            this->logger(TlUtils::format("+     maximum  deviation of xa alpha coefficient         = %14.4le\n", this->dev_xa_a));
            this->logger(TlUtils::format("+     maximum  deviation of xa beta  coefficient         = %14.4le\n", this->dev_xa_b));
        }
        break;

    default:
        std::cerr << "program error:" << __FILE__ << __LINE__ << std::endl;
        break;
    }
    
}

// total energy convergence
double DfConvcheck::dev_total_energy(int iteration)
{
    assert(iteration >= 2);

    TlVector Ene;
    Ene.load("fl_Work/fl_Vct_Energy");

    const double deviation_value = std::fabs(Ene[iteration -1] - Ene[iteration -1-1]);

    return deviation_value;
}


// cd coefficient convergence
double DfConvcheck::dev_cd_coefficient(const RUN_TYPE runType, const int iteration)
{
    assert(iteration >= 2);

    TlVector prevRho = DfObject::getRho<TlVector>(runType, iteration -1);
    TlVector rho = DfObject::getRho<TlVector>(runType, iteration);

    assert(prevRho.getSize() == this->m_nNumOfAux);
    assert(rho.getSize() == this->m_nNumOfAux);

    // get maximum deviation
    TlVector v = rho - prevRho;
    const double deviation_value = v.getMaxAbsoluteElement();

    return deviation_value;
}


// xc coefficient convergence
double DfConvcheck::dev_xc_coefficient(const RUN_TYPE runType, const int iteration)
{
    assert(iteration >= 2);

    double deviation_value = 0.0;

    if (this->m_bIsXCFitting == true) {
        TlVector prevMyu = DfObject::getMyu<TlVector>(runType, iteration -1);
        TlVector Myu = DfObject::getMyu<TlVector>(runType, iteration -1);

        assert(prevMyu.getSize() == this->m_nNumOfAux);
        assert(Myu.getSize() == this->m_nNumOfAux);

        // get maximum deviation
        TlVector v = Myu - prevMyu;
        deviation_value = v.getMaxAbsoluteElement();
    }

    return deviation_value;
}


double DfConvcheck::dev_xa_coefficient(const RUN_TYPE runType, const int iteration)
{
    assert(iteration >= 2);

    double deviation_value = 0.0;

    // xa coefficient convergence if NEED
    if (this->m_sXCFunctional == "xalpha") {
        TlVector prevNyu = DfObject::getNyu<TlVector>(runType, iteration -1);
        TlVector Nyu = DfObject::getNyu<TlVector>(runType, iteration);

        // get maximum deviation
        TlVector v = Nyu - prevNyu;
        deviation_value = v.getMaxAbsoluteElement();
    }

    return deviation_value;
}


