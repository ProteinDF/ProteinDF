#include "df_total_energy_tmpl.h"
#include "DfOverlapX.h"
#include "DfXCFunctional.h"
#include "Fl_Geometry.h"

template <class SymmetricMatrix, class Vector, class DfOverlap>
const double DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::TOO_SMALL =
    1.0E-16;
template <class SymmetricMatrix, class Vector, class DfOverlap>
double DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::E_nuc_ = 0.0;
template <class SymmetricMatrix, class Vector, class DfOverlap>
double DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::E_nuc_X_ = 0.0;

template <class SymmetricMatrix, class Vector, class DfOverlap>
DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::DfTotalEnergy_tmpl(
    TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      E_H_(0.0),
      E_H2_(0.0),
      E_J_(0, 0),
      E_K_(0.0),
      E_xc_(0.0),
      pDfXCFunctionalObject_(NULL),
      pDfOverlapObject_(NULL),
      pPA_(NULL),
      pPB_(NULL) {}

template <class SymmetricMatrix, class Vector, class DfOverlap>
DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::~DfTotalEnergy_tmpl() {
  delete this->pDfXcFuncionalObject_;
  this->pDfXcFunctionalObject_ = NULL;
  delete this->pDfOverlapObject_;
  this->pDfOverlapObject_;

  delete this->pPA_;
  this->pPA_ = NULL;
  delete this->PB_;
  this->pPB_ = NULL;
  delete this->pPAB_;
  this->pPAB_ = NULL;
}

// ----------------------------------------------------------------------------
// initialize
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector,
                        DfOverlap>::prepareXCFuncionalObject() {
  this->pDfXCFunctionalObject_ = new DfXCFunctional(this->pPdfParam_);
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector,
                        DfOverlap>::prepareDfOverlapObject() {
  this->pDfOverlapObject_ = new DfOverlapX(this->pPdfParam_);
}

// ----------------------------------------------------------------------------
// MAIN
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calc(
    const int iteration) {
  this->calcScfIteration_ = iteration;

  this->calcE_nuc();

  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      this->calcE_K(*(this->pPA_));
      this->E_K_ *= 2.0;

      this->calcE_xc(RUN_RKS, *(this->pPA_));
      this->E_xc_ *= 2.0;

      // PA <- (PA+PB)
      *(this->pPA_) *= 2.0;
    } break;

    case METHOD_UKS: {
      this->calcE_K(*(this->PA));
      this->calcE_K(*(this->PB));

      this->calcE_xc(RUN_UKS_ALPHA, *(this->pPA_));
      this->calcE_xc(RUN_UKS_BETA, *(this->pPB_));

      // PA <- (PA+PB)
      *(this->pPA_) = *(this->pPA_) + *(this->pPB_);
    } break;

    case METHOD_ROKS: {
      this->calcE_K(*(this->PA));
      this->calcE_K(*(this->PB));

      this->calcE_xc(RUN_UKS_ALPHA, *(this->pPA_));
      this->calcE_xc(RUN_UKS_BETA, *(this->pPB_));

      // PA <- (PA+PB)
      *(this->pPA_) = *(this->pPA_) + *(this->pPB_);
    } break;

    default:
      break;
  }

  this->calcE_H();
  this->calcE_J();

  if (this->enableGrimmeDispersion_ == true) {
    this->E_disp_ = this->pDfXCFunctionalObject_->getGrimmeDispersionEnergy();
  }
}

// ----------------------------------------------------------------------------
// energy for nuclear repulsion
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_nuc() {
  if (std::fabs(this->E_nuc_) > TOO_SMALL) {
    this->logger(" energy of nuclear repulstion has already been calculated.");
  } else {
    // read nuclear charge
    const Fl_Geometry geom((*this->pPdfParam_)["coordinates"]);
    const int numOfAtoms = this->m_nNumOfAtoms;

    // calculate nuclear repulsion
    double E_nuc = 0.0;
    double E_nuc_X = 0.0;
    for (int i = 0; i < numOfAtoms; ++i) {
      const double ci = geom.getCharge(i);
      const TlPosition pi = geom.getCoordinate(i);
      const bool isDummy = ("X" == this->getAtomSymbol(i));

      for (int j = 0; j < i; ++j) {
        const double cj = geom.getCharge(j);
        const TlPosition pj = geom.getCoordinate(j);
        const double distance = pi.distanceFrom(pj);

        const double E = ci * cj / distance;
        E_nuc += E;
        if (isDummy || ("X" == this->getAtomSymbol(j))) {
          E_nuc_X += E;
        }
      }
    }

    this->E_nuc_ = E_nuc;
    this->E_nuc_X_ = E_nuc_X;
  }
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::output() {
  double E_Total = 0.0;

  // 表示
  this->log_.info("------------------------------------------------");

  // 1e term
  this->log_.info(TlUtils::format(" Ts+Vn          = %28.16lf", this->E_H_));
  E_Total += this->E_H;

  // J term
  switch (this->J_engine_) {
    case J_ENGINE_RI_J: {
      if ((this->m_bMemorySave == false) &&
          (this->m_bDiskUtilization == false)) {
        E_Total += this->m_dE_OEP_JRR_Exc;

        this->logger(TlUtils::format(" Ts+Vn+J[Rho~,Rho~]+Exc = %28.16lf",
                                     this->m_dE_OEP_JRR_Exc));
      }
      this->log_.info(TlUtils::format(" E_J[Rho, Rho~] = %28.16lf\n",
                                      this->m_dE_J_Rho_RhoTilde));
      E_Total += this->m_dE_J_Rho_RhoTilde;
      this->log_.info(TlUtils::format(" E_J[Rho~,Rho~] = %28.16lf\n",
                                      this->m_dE_J_RhoTilde_RhoTilde));
      E_Total += this->m_dE_J_RhoTilde_RhoTilde;
    } break;

    case J_ENGINE_CONVENTIONAL:
    case J_ENGINE_CD: {
      this->log_.info(
          TlUtils::format(" E_J            = %28.16lf\n", this->E_J_));
      E_Total += this->E_J_;
    } break;

    default:
      break;
  }

  // 
  this->log_.info(
      TlUtils::format(" E_xc(pure)     = %28.16lf\n", this->m_dExc));
  E_Total += this->m_dExc;

  if (this->enableGrimmeDispersion_ == true) {
    this->log_.info(TlUtils::format(" E_cx(+disp.)   = %28.16lf\n",
                                    E_Total + this->E_disp_));
  }

  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->log_.info(
          TlUtils::format(" E_K            = %28.16lf\n", this->K_term_));
      break;

    case METHOD_UKS:
      this->log_.info(
          TlUtils::format(" E_K            = %28.16lf\n", this->K_term_));
      this->log_.info(
          TlUtils::format("   E_K(alpha)   = %28.16lf\n", this->E_KA_));
      this->log_.info(
          TlUtils::format("   E_K(beta)    = %28.16lf\n", this->E_KB_));
      break;

    case METHOD_ROKS:
      this->log_.info(
          TlUtils::format(" E_K            = %28.16lf\n", this->K_term_));
      this->log_.info(
          TlUtils::format("   E_K(alpha)   = %28.16lf\n", this->E_KA_));
      this->log_.info(
          TlUtils::format("   E_K(beta)    = %28.16lf\n", this->E_KB_));
      break;

    default:
      this->log_.critical("program error");
      break;
  }
  E_Total += this->K_term_;

  this->log_.info(TlUtils::format(" E_nuclei       = %28.16lf\n",
                                  this->E_nuc_));
  E_Total += this->E_nuc_;

  this->log_.info(TlUtils::format(" TE             = %28.16lf\n", E_Total));
  // this->logger("------------------------------------------------\n");
  // this->logger(TlUtils::format(" Ts+Vn        = %28.16lf\n",
  // this->m_dE_OneElectronPart)); this->logger(TlUtils::format(" J[Rho, Rho~]
  // = %28.16lf\n", this->m_dE_J_Rho_RhoTilde)); if (this->J_engine_ ==
  // J_ENGINE_RI_J) {
  //     this->logger(TlUtils::format(" J[Rho~,Rho~] = %28.16lf\n",
  //     this->m_dE_J_RhoTilde_RhoTilde));
  // }
  // this->logger(TlUtils::format(" Exc          = %28.16lf\n",
  // this->m_dExc)); this->logger(TlUtils::format(" Enuclei      =
  // %28.16lf\n", this->m_dE_NuclearRepulsion));
  // this->logger(TlUtils::format(" TE           = %28.16lf\n", E_Total));
  // this->logger("------------------------------------------------\n");

  std::cout << TlUtils::format(" %3d th TE = %18.16lf", this->m_nIteration,
                               E_Total)
            << std::endl;
  this->log_.info("------------------------------------------------");
  this->write_total_energy(E_Total);
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::updateParam() {
  (*this->pPdfParam_)["TEs"][this->m_nIteration] = this->E_total_;
}
