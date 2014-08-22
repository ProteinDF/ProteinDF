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
    this->convergence_type = pdfParam["convergence/type"].getStr();

    this->threshold_cri     = pdfParam["convergence/threshold"].getDouble();
    this->threshold_cri_ene = pdfParam["convergence/threshold-energy"].getDouble();

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
        this->log_.info(TlUtils::format("+ convergence informations\n"));
        this->log_.info(TlUtils::format("+     standard deviation of cd               = %14.4le", this->dev_sd));
        this->log_.info(TlUtils::format("+     maximum  deviation of cd coefficient   = %14.4le", this->dev_cd));
        this->log_.info(TlUtils::format("+     maximum  deviation of density matrix   = %14.4le", this->dev_dm));
        this->log_.info(TlUtils::format("+              deviation of total energy     = %14.4le", this->dev_te));
        this->log_.info(TlUtils::format("+     maximum  deviation of kohn-sham matrix = %14.4le", this->dev_ks));
        if (this->m_bIsXCFitting == true) {
            this->log_.info(TlUtils::format("+     maximum  deviation of xc coefficient   = %14.4le", this->dev_xc));
            this->log_.info(TlUtils::format("+     maximum  deviation of xa coefficient   = %14.4le", this->dev_xa));
        }
        break;
        
    case METHOD_UKS:
        this->log_.info("+ convergence informations");
        this->log_.info(TlUtils::format("+     standard deviation of alpha cd               = %14.4le", this->dev_sd_a));
        this->log_.info(TlUtils::format("+     standard deviation of beta  cd               = %14.4le", this->dev_sd_b));
        this->log_.info(TlUtils::format("+     maximum  deviation of alpha density matrix   = %14.4le", this->dev_dm_a));
        this->log_.info(TlUtils::format("+     maximum  deviation of beta  density matrix   = %14.4le", this->dev_dm_b));
        this->log_.info(TlUtils::format("+              deviation of total energy           = %14.4le", this->dev_te));
        this->log_.info(TlUtils::format("+     maximum  deviation of alpha kohn-sham matrix = %14.4le", this->dev_ks_a));
        this->log_.info(TlUtils::format("+     maximum  deviation of beta  kohn-sham matrix = %14.4le", this->dev_ks_b));
        this->log_.info(TlUtils::format("+     maximum  deviation of cd alpha coefficient   = %14.4le", this->dev_cd_a));
        this->log_.info(TlUtils::format("+     maximum  deviation of cd beta  coefficient   = %14.4le", this->dev_cd_b));
        if (this->m_bIsXCFitting == true) {
            this->log_.info(TlUtils::format("+     maximum  deviation of xc alpha coefficient   = %14.4le", this->dev_xc_a));
            this->log_.info(TlUtils::format("+     maximum  deviation of xc beta  coefficient   = %14.4le", this->dev_xc_b));
            this->log_.info(TlUtils::format("+     maximum  deviation of xa alpha coefficient   = %14.4le", this->dev_xa_a));
            this->log_.info(TlUtils::format("+     maximum  deviation of xa beta  coefficient   = %14.4le", this->dev_xa_b));
        }
        break;

    case METHOD_ROKS:
        this->log_.info("+ convergence informations");
        this->log_.info(TlUtils::format("+     standard deviation of alpha cd                     = %14.4le", this->dev_sd_a));
        this->log_.info(TlUtils::format("+     standard deviation of beta  cd                     = %14.4le", this->dev_sd_b));
        this->log_.info(TlUtils::format("+     maximum  deviation of closed part density matrix   = %14.4le", this->dev_dm_c));
        this->log_.info(TlUtils::format("+     maximum  deviation of open   part density matrix   = %14.4le", this->dev_dm_o));
        this->log_.info(TlUtils::format("+              deviation of total energy                 = %14.4le", this->dev_te));
        this->log_.info(TlUtils::format("+     maximum  deviation of kohn-sham matrix             = %14.4le", this->dev_ks));
        this->log_.info(TlUtils::format("+     maximum  deviation of cd alpha coefficient         = %14.4le", this->dev_cd_a));
        this->log_.info(TlUtils::format("+     maximum  deviation of cd beta  coefficient         = %14.4le", this->dev_cd_b));
        if (this->m_bIsXCFitting == true) {
            this->log_.info(TlUtils::format("+     maximum  deviation of xc alpha coefficient         = %14.4le", this->dev_xc_a));
            this->log_.info(TlUtils::format("+     maximum  deviation of xc beta  coefficient         = %14.4le", this->dev_xc_b));
            this->log_.info(TlUtils::format("+     maximum  deviation of xa alpha coefficient         = %14.4le", this->dev_xa_a));
            this->log_.info(TlUtils::format("+     maximum  deviation of xa beta  coefficient         = %14.4le", this->dev_xa_b));
        }
        break;

    default:
        std::cerr << "program error:" << __FILE__ << __LINE__ << std::endl;
        break;
    }
    
}

// total energy convergence
double DfConvcheck::dev_total_energy(const int iteration)
{
    assert(iteration >= 2);

    const double TE = (*this->pPdfParam_)["TE"][iteration].getDouble();
    const double prevTE = (*this->pPdfParam_)["TE"][iteration -1].getDouble();
    const double deviation_value = std::fabs(TE - prevTE);
    
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


