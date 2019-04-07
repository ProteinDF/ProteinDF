#ifndef DF_TOTALENERGY_TMPL_H
#define DF_TOTALENERGY_TMPL_H

#include "CnError.h"
#include "DfObject.h"
#include "DfOverlapX.h"
#include "DfXCFunctional.h"

// ----------------------------------------------------------------------------
// template class
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlapType>
class DfTotalEnergy_tmpl : public DfObject {
 public:
  DfTotalEnergy_tmpl(TlSerializeData* pPdfParam);
  virtual ~DfTotalEnergy_tmpl();

 public:
  void calc(const int iteration);
  void output();
  void updateParam();

 protected:
  virtual void prepareXCFuncionalObject();
  virtual void prepareDfOverlapObject();

  void loadDensityMatrix();
  Vector getRho();
  Vector getEps(const RUN_TYPE runType);

  void calcE_nuc();
  void calcE_H();

  void calcE_J();
  void calcE_J_exact();
  void calcE_J_RI();
  void calcE_J_rho_rhoTilde_DIRECT();
  void calcE_J_rhoTilde_rhoTilde();

  void calcE_K(const RUN_TYPE runType, const SymmetricMatrix& PA);
  void calcE_K_exact(const RUN_TYPE runType, const SymmetricMatrix& PA);

  void calcE_xc(const RUN_TYPE runType, const SymmetricMatrix& PA);
  void calcE_xc_gridfree(const RUN_TYPE runType, const SymmetricMatrix& PA);
  void calcE_xc_DIRECT(const SymmetricMatrix& PA, const Vector& eps);

 protected:
  DfXCFunctional* pDfXCFunctionalObject_;
  DfOverlapX* pDfOverlapObject_;

  static const double TOO_SMALL;
  static double E_nuc_;
  static double E_nuc_X_;
  double E_H_;
  double E_H2_;
  double E_J_;
  double E_J_rho_rhoTilde_;
  double E_J_rhoTilde_rhoTilde_;
  double E_K_;
  double E_xc_;
  double E_disp_;

  int calcScfIteration_;
  SymmetricMatrix E_;

  SymmetricMatrix* pPA_;
  SymmetricMatrix* pPB_;
  SymmetricMatrix* pPAB_;
};

// ----------------------------------------------------------------------------
// template
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
const double DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::TOO_SMALL =
    1.0E-16;

template <class SymmetricMatrix, class Vector, class DfOverlap>
double DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::E_nuc_ = 0.0;

template <class SymmetricMatrix, class Vector, class DfOverlap>
double DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::E_nuc_X_ = 0.0;

// constructor & destructor ---------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::DfTotalEnergy_tmpl(
    TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      pDfXCFunctionalObject_(NULL),
      pDfOverlapObject_(NULL),
      E_H_(0.0),
      E_H2_(0.0),
      E_J_(0.0),
      E_J_rho_rhoTilde_(0.0),
      E_J_rhoTilde_rhoTilde_(0.0),
      E_K_(0.0),
      E_xc_(0.0),
      E_disp_(0.0),
      pPA_(NULL),
      pPB_(NULL),
      pPAB_(NULL) {
  //
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::~DfTotalEnergy_tmpl() {
  delete this->pDfXCFunctionalObject_;
  this->pDfXCFunctionalObject_ = NULL;
  delete this->pDfOverlapObject_;
  this->pDfOverlapObject_ = NULL;

  delete this->pPA_;
  this->pPA_ = NULL;
  delete this->pPB_;
  this->pPB_ = NULL;
  delete this->pPAB_;
  this->pPAB_ = NULL;
}

// initialize -----------------------------------------------------------------
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
// density matrix
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector,
                        DfOverlap>::loadDensityMatrix() {
  const int iteration = this->calcScfIteration_;

  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->pPA_ = new SymmetricMatrix(
          DfObject::getPpqMatrix<SymmetricMatrix>(RUN_RKS, iteration));
      this->pPAB_ = new SymmetricMatrix(2.0 * (*(this->pPA_)));
      break;

    case METHOD_UKS:
      this->pPA_ = new SymmetricMatrix(
          DfObject::getPpqMatrix<SymmetricMatrix>(RUN_UKS_ALPHA, iteration));
      this->pPB_ = new SymmetricMatrix(
          DfObject::getPpqMatrix<SymmetricMatrix>(RUN_UKS_BETA, iteration));
      this->pPAB_ = new SymmetricMatrix(*(this->pPA_) + *(this->pPB_));
      break;

    case METHOD_ROKS:
      this->pPB_ =
          new SymmetricMatrix(0.5 * DfObject::getPpqMatrix<SymmetricMatrix>(
                                        RUN_ROKS_CLOSED, iteration));
      this->pPA_ = new SymmetricMatrix(
          *(this->pPB_) +
          DfObject::getPpqMatrix<SymmetricMatrix>(RUN_ROKS_OPEN, iteration));
      this->pPAB_ = new SymmetricMatrix(*(this->pPA_) + *(this->pPB_));
      break;

    default:
      this->log_.critical("program error");
  }
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
Vector DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::getRho() {
  Vector rho;

  switch (this->m_nMethodType) {
    case METHOD_RKS:
      rho.load(this->getRhoPath(RUN_RKS, this->m_nIteration));
      break;

    case METHOD_UKS:  // go down
    case METHOD_ROKS: {
      Vector rho_a, rho_b;
      rho_a.load(this->getRhoPath(RUN_UKS_ALPHA, this->m_nIteration));
      rho_b.load(this->getRhoPath(RUN_UKS_BETA, this->m_nIteration));
      rho = rho_a + rho_b;
    } break;

    default:
      CnErr.abort(" DfTotalEnergy::getRho error.");
      break;
  }

  return rho;
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
Vector DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::getEps(
    const RUN_TYPE runType) {
  Vector eps;

  switch (runType) {
    case RUN_RKS:
      eps.load("fl_Work/fl_Vct_Epsilon");
      assert(static_cast<TlVectorAbstract::size_type>(this->numOfAuxXC_) ==
             eps.getSize());
      break;

    case RUN_UKS_ALPHA:
      eps.load("fl_Work/fl_Vct_Epsilona");
      assert(static_cast<TlVectorAbstract::size_type>(this->numOfAuxXC_) ==
             eps.getSize());
      break;

    case RUN_UKS_BETA:
      eps.load("fl_Work/fl_Vct_Epsilonb");
      assert(static_cast<TlVectorAbstract::size_type>(this->numOfAuxXC_) ==
             eps.getSize());
      break;

    default:
      CnErr.abort("DfTotalEnergy_tmpl::getEps() unknown type.");
      break;
  }

  return eps;
}

// ----------------------------------------------------------------------------
// MAIN
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calc(
    const int iteration) {
  this->calcE_nuc();

  this->calcScfIteration_ = iteration;
  this->loadDensityMatrix();
  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      this->calcE_K(RUN_RKS, *(this->pPA_));
      this->E_K_ *= 2.0;

      this->calcE_xc(RUN_RKS, *(this->pPA_));
      this->E_xc_ *= 2.0;
    } break;

    case METHOD_UKS: {
      this->calcE_K(RUN_UKS_ALPHA, *(this->pPA_));
      this->calcE_K(RUN_UKS_BETA, *(this->pPB_));

      this->calcE_xc(RUN_UKS_ALPHA, *(this->pPA_));
      this->calcE_xc(RUN_UKS_BETA, *(this->pPB_));
    } break;

    case METHOD_ROKS: {
      this->calcE_K(RUN_ROKS_ALPHA, *(this->pPA_));
      this->calcE_K(RUN_ROKS_BETA, *(this->pPB_));

      this->calcE_xc(RUN_ROKS_ALPHA, *(this->pPA_));
      this->calcE_xc(RUN_ROKS_BETA, *(this->pPB_));
    } break;

    default:
      break;
  }

  this->calcE_H();
  this->calcE_J();

  // if (this->enableGrimmeDispersion_ == true) {
  //   this->E_disp_ =
  //   this->pDfXCFunctionalObject_->getGrimmeDispersionEnergy();
  // }
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
      const bool isDummy = ("X" == geom.getAtomSymbol(i));

      for (int j = 0; j < i; ++j) {
        const double cj = geom.getCharge(j);
        const TlPosition pj = geom.getCoordinate(j);
        const double distance = pi.distanceFrom(pj);

        const double E = ci * cj / distance;
        E_nuc += E;
        if (isDummy || ("X" == geom.getAtomSymbol(j))) {
          E_nuc_X += E;
        }
      }
    }

    this->E_nuc_ = E_nuc;
    this->E_nuc_X_ = E_nuc_X;
  }
}

// ----------------------------------------------------------------------------
// H
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_H() {
  const SymmetricMatrix& P_AB = *(this->pPA_);

  SymmetricMatrix Hpq = DfObject::getHpqMatrix<SymmetricMatrix>();
  if (Hpq.getNumOfRows() != this->m_nNumOfAOs ||
      Hpq.getNumOfCols() != this->m_nNumOfAOs) {
    CnErr.abort("DfTotalEnergy", "calculate_one_electron_part", "",
                "program error");
  }

  const SymmetricMatrix E_H = Hpq.dotInPlace(P_AB);
  this->E_H_ = E_H.sum();

  // add dummy charge
  if (this->m_nNumOfDummyAtoms > 0) {
    const int chgextra_number =
        (*(this->pPdfParam_))["charge-extrapolate-number"].getInt();

    SymmetricMatrix Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrix>();

    if (chgextra_number > 1) {
      const int coef = std::min(this->m_nIteration, chgextra_number);
      Hpq2 *= coef;
    }

    const SymmetricMatrix E_H2 = Hpq2.dotInPlace(P_AB);
    this->E_H2_ = E_H2.sum();
  }
}

// ----------------------------------------------------------------------------
// J term
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_J() {
  switch (this->J_engine_) {
    case J_ENGINE_RI_J: {
      this->calcE_J_RI();
    } break;

    case J_ENGINE_CONVENTIONAL:
    case J_ENGINE_CD:
      this->calcE_J_exact();
      break;

    default:
      break;
  }
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_J_exact() {
  const SymmetricMatrix& P_AB = *(this->pPAB_);
  SymmetricMatrix J =
      DfObject::getJMatrix<SymmetricMatrix>(this->calcScfIteration_);

  SymmetricMatrix E_J = J.dotInPlace(P_AB);
  this->E_J_ = 0.25 * E_J.sum();
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_J_RI() {
  this->calcE_J_rhoTilde_rhoTilde();

  if (!((this->m_bMemorySave == false) &&
        (this->m_bDiskUtilization == false))) {
    // NOT use Threeindexintegrals
    this->calcE_J_rho_rhoTilde_DIRECT();
  }
}

// energy for a part of coulomb term ( rho*Ppq*Pqa )
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector,
                        DfOverlap>::calcE_J_rho_rhoTilde_DIRECT() {
  const SymmetricMatrix& P_AB = *(this->pPA_);
  SymmetricMatrix J =
      DfObject::getJMatrix<SymmetricMatrix>(this->calcScfIteration_);

  const double coef = (this->J_engine_ == J_ENGINE_RI_J) ? 1.0 : 0.5;
  SymmetricMatrix E_J = coef * J.dotInPlace(P_AB);
  this->E_J_rho_rhoTilde_ = E_J.sum();
}

/// J[rho~, rho~] (( "1/2*rho*rho'*Sab" ))の計算を行う
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector,
                        DfOverlap>::calcE_J_rhoTilde_rhoTilde() {
  SymmetricMatrix Sab = DfObject::getSabMatrix<SymmetricMatrix>();
  const Vector rho = this->getRho();
  this->E_J_rhoTilde_rhoTilde_ = -0.5 * (rho.dot(Sab * rho));
}

// ----------------------------------------------------------------------------
// K term
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_K(
    const RUN_TYPE runType, const SymmetricMatrix& PA) {
  this->calcE_K_exact(runType, PA);
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_K_exact(
    const RUN_TYPE runType, const SymmetricMatrix& PA) {
  DfXCFunctional dfXCFunctional(this->pPdfParam_);
  if (dfXCFunctional.isHybridFunctional() == true) {
    SymmetricMatrix K = DfObject::getHFxMatrix<SymmetricMatrix>(
        runType, this->calcScfIteration_);
    this->E_K_ += 0.5 * 0.5 * dfXCFunctional.getFockExchangeCoefficient() *
                  K.dotInPlace(PA).sum();
  }
}

// ----------------------------------------------------------------------------
// XC term
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_xc(
    const RUN_TYPE runType, const SymmetricMatrix& PA) {
  DfXCFunctional dfXCFunctional(this->pPdfParam_);
  if (dfXCFunctional.getXcType() != DfXCFunctional::HF) {
    if (this->m_bIsXCFitting == true) {
      const Vector Eps = this->getEps(runType);
      this->calcE_xc_DIRECT(PA, Eps);
    } else {
      if (this->XC_engine_ != XC_ENGINE_GRID) {
        this->calcE_xc_gridfree(runType, PA);
      } else {
        this->E_xc_ = dfXCFunctional.getEnergy();
      }
    }
  }
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_xc_gridfree(
    const RUN_TYPE runType, const SymmetricMatrix& PA) {
  SymmetricMatrix Exc =
      DfObject::getExcMatrix<SymmetricMatrix>(runType, this->calcScfIteration_);
  this->E_xc_ += Exc.dotInPlace(PA).sum();
}

// energy for xc energy term (compare myu*Ppq*Pqa with 4/3*Ex1)
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_xc_DIRECT(
    const SymmetricMatrix& PA, const Vector& eps) {
  // put eps*[pqs] into B
  SymmetricMatrix pqg(this->m_nNumOfAOs);

  // TODO: lapack only
  // this->pDfOverlapObject_->get_pqg(eps, &pqg);

  this->E_xc_ = pqg.dotInPlace(PA).sum();
}

// ----------------------------------------------------------------------------
// output
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::output() {
  double E_Total = 0.0;

  // 表示
  this->log_.info("------------------------------------------------");

  // 1e term
  this->log_.info(TlUtils::format(" Ts+Vn          = %28.16lf", this->E_H_));
  E_Total += this->E_H_;

  // J term
  switch (this->J_engine_) {
    case J_ENGINE_RI_J: {
      // if ((this->m_bMemorySave == false) &&
      //     (this->m_bDiskUtilization == false)) {
      //   E_Total += this->m_dE_OEP_JRR_Exc;

      //   this->logger(TlUtils::format(" Ts+Vn+J[Rho~,Rho~]+Exc = %28.16lf",
      //                                this->m_dE_OEP_JRR_Exc));
      // }
      this->log_.info(TlUtils::format(" E_J[Rho, Rho~] = %28.16lf\n",
                                      this->E_J_rho_rhoTilde_));
      // E_Total += this->m_dE_J_Rho_RhoTilde;
      this->log_.info(TlUtils::format(" E_J[Rho~,Rho~] = %28.16lf\n",
                                      this->E_J_rhoTilde_rhoTilde_));
      // E_Total += this->m_dE_J_RhoTilde_RhoTilde;
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
  this->log_.info(TlUtils::format(" E_xc(pure)     = %28.16lf\n", this->E_xc_));
  E_Total += this->E_xc_;

  if (this->enableGrimmeDispersion_ == true) {
    this->log_.info(TlUtils::format(" E_cx(+disp.)   = %28.16lf\n",
                                    E_Total + this->E_disp_));
  }

  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->log_.info(
          TlUtils::format(" E_K            = %28.16lf\n", this->E_K_));
      break;

    case METHOD_UKS:
      this->log_.info(
          TlUtils::format(" E_K            = %28.16lf\n", this->E_K_));
      // this->log_.info(
      //     TlUtils::format("   E_K(alpha)   = %28.16lf\n", this->E_KA_));
      // this->log_.info(
      //     TlUtils::format("   E_K(beta)    = %28.16lf\n", this->E_KB_));
      break;

    case METHOD_ROKS:
      this->log_.info(
          TlUtils::format(" E_K            = %28.16lf\n", this->E_K_));
      // this->log_.info(
      //     TlUtils::format("   E_K(alpha)   = %28.16lf\n", this->E_KA_));
      // this->log_.info(
      //     TlUtils::format("   E_K(beta)    = %28.16lf\n", this->E_KB_));
      break;

    default:
      this->log_.critical("program error");
      break;
  }
  E_Total += this->E_K_;

  this->log_.info(
      TlUtils::format(" E_nuclei       = %28.16lf\n", this->E_nuc_));
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
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::updateParam() {
  (*this->pPdfParam_)["TEs"][this->m_nIteration] = this->E_total_;
}

#endif  // DF_TOTALENERGY_TMPL_H
