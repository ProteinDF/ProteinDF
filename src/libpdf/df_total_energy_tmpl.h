#ifndef DF_TOTALENERGY_TMPL_H
#define DF_TOTALENERGY_TMPL_H

#include "DfObject.h"
#include "CnError.h"
#include "DfXCFunctional.h"
#include "DfOverlapX.h"

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

  void loadDensityMatrix(const int iteration);
  Vector getRho();

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
  double E_K_;
  double E_xc_;

  int calcScfIteration_;
  SymmetricMatrix E_;

  SymmetricMatrix* pPA_;
  SymmetricMatrix* pPB_;
};

// ----------------------------------------------------------------------------
// template
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::loadDensityMatrix(
    const int iteration) {
  switch (this->m_nMethodType) {
    case METHOD_RKS:
      this->pPA_ = new SymmetricMatrix(
          DfObject::getPpqMatrix<SymmetricMatrix>(RUN_RKS, iteration));
      break;

    case METHOD_UKS:
      this->pPA_ = new SymmetricMatrix(
          DfObject::getPpqMatrix<SymmetricMatrix>(RUN_UKS_ALPHA, iteration));
      this->pPB_ = new SymmetricMatrix(
          DfObject::getPpqMatrix<SymmetricMatrix>(RUN_UKS_BETA, iteration));
      break;

    case METHOD_ROKS:
      this->pPB_ =
          new SymmetricMatrix(0.5 * DfObject::getPpqMatrix<SymmetricMatrix>(
                                        RUN_ROKS_CLOSED, iteration));
      this->pPA_ = new SymmetricMatrix(
          this->PB_ +
          DfObject::getPpqMatrix<SymmetricMatrix>(RUN_ROKS_OPEN, iteration));
      break;

    default:
      this->log_.critical("program error");
  }
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
Vector DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::getRho() {
  Vector rho;

  switch (this->methodType_) {
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
      this->J_term_ = this->calcE_J_RI();
    } break;

    case J_ENGINE_CONVENTIONAL:
    case J_ENGINE_CD:
      this->J_term_ = this->calcE_J_exact();
      break;

    default:
      break;
  }
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_J_exact() {
  const SymmetricMatrix& P_AB = *(this->pPA_);
  SymmetricMatrix J =
      DfObject::getJMatrix<SymmetricMatrix>(this->calcScfIteration_);

  SymmetricMatrix E_J = J.dotInPlace(P_AB);
  this->E_J_ = E_J.sum();
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
void DfTotalEnergy_tmpl<SymmetricMatrix,
                        Vector, DfOverlap>::calcE_J_rho_rhoTilde_DIRECT() {
  const SymmetricMatrix& P_AB = *(this->pPA_);
  SymmetricMatrix J =
      DfObject::getJMatrix<SymmetricMatrix>(this->calcScfIteration_);

  const double coef = (this->J_engine_ == J_ENGINE_RI_J) ? 1.0 : 0.5;
  SymmetricMatrix E_J = coef * J.dotInPlace(P_AB);
  this->E_J_rho_rhoTilde_ = E_J.sum();
}

/// J[rho~, rho~] (( "1/2*rho*rho'*Sab" ))の計算を行う
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_J_rhoTilde_rhoTilde() {
  SymmetricMatrix Sab = DfObject::getSabMatrix<SymmetricMatrix>();
  const Vector rho = this->getRho(this->m_nMethodType);
  this->E_J_rhoTilde_rhoTilde_ = -0.5 * (rho * (Sab * rho));
}

// ----------------------------------------------------------------------------
// K term
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_K(
    const RUN_TYPE runType, const SymmetricMatrix& PA) {
  this->calcK_exact(runType, PA);
}

template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_K_exact(
    const RUN_TYPE runType, const SymmetricMatrix& PA) {
  DfXCFunctional dfXCFunctional(this->pPdfParam_);
  if (dfXCFunctional.isHybridFunctional() == true) {
    SymmetricMatrix K = DfObject::getHFxMatrix<SymmetricMatrix>(
        runType, this->calcScfIteration_);
    this->E_K_ += 0.5 * dfXCFunctional.getFockExchangeCoefficient() *
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
      this->E_xc_ =
          this->calcExc_DIRECT(PA, Eps);
    } else {
      if (this->XC_engine_ != XC_ENGINE_GRID) {
        this->calcE_xc_gridfree(PA);
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
      DfObject::getExcMatrix<SymmetricMatrix>(runType, this->calcScfteration_);
  this->E_xc_ += Exc.dotInPlace(PA).sum();
}

// energy for xc energy term (compare myu*Ppq*Pqa with 4/3*Ex1)
template <class SymmetricMatrix, class Vector, class DfOverlap>
void DfTotalEnergy_tmpl<SymmetricMatrix, Vector, DfOverlap>::calcE_xc_DIRECT(const SymmetricMatrix& PA, const Vector& eps) {
  // put eps*[pqs] into B
  SymmetricMatrix pqg(this->m_nNumOfAOs);
  this->pDfOverlapObject_->get_pqg(eps, &pqg);

  this->E_xc = pqg.dotInPlace(PA).sum();
}

#endif  // DF_TOTALENERGY_TMPL_H
