#ifndef DF_TOTALENERGY_TMPL_H
#define DF_TOTALENERGY_TMPL_H

#include <vector>

#include "CnError.h"
#include "DfObject.h"
#include "DfOverlapX.h"
#include "DfXCFunctional.h"
#include "TlOrbitalInfo.h"

// ----------------------------------------------------------------------------
// template class
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
class DfTotalEnergy_tmpl : public DfObject {
   public:
    DfTotalEnergy_tmpl(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy_tmpl();

   public:
    void run();
    void calc(const int iteration);
    void output();

    void get_IE(const std::vector<std::vector<DfObject::index_type> >&
                    atomIndexGroup) const;

    // double get_IE(const std::vector<int>& indeces1) const;
    // double get_IE(const std::vector<int>& indeces1,
    //               const std::vector<int>& indeces2) const;

   protected:
    virtual void saveParams();

   protected:
    virtual void prepareXCFunctionalObject();
    virtual void prepareDfOverlapObject();

   protected:
    virtual SymmetricMatrix getSpinDensityMatrix(const RUN_TYPE runType,
                                                 const int iteration);

    virtual Vector getRho();
    Vector getEps(const RUN_TYPE runType);

    void calcE_nuc();

    // calc energy from the core hamiltonian
    void calcE_h(const SymmetricMatrix& P_AB);
    // calc energy from the core hamiltonian (corresponding to dummy charges)
    void calcE_h_fromX(const SymmetricMatrix& P_AB);

    void calcE_J(const SymmetricMatrix& P_AB);
    void calcE_J_exact(const SymmetricMatrix& P_AB);
    void calcE_J_RI(const SymmetricMatrix& P_AB);
    void calcE_J_rho_rhoTilde_DIRECT(const SymmetricMatrix& P_AB);
    void calcE_J_rhoTilde_rhoTilde();

    double calcE_K(const RUN_TYPE runType, const SymmetricMatrix& PA,
                   SymmetricMatrix* pE_K = NULL);
    double calcE_K_exact(const RUN_TYPE runType, const SymmetricMatrix& PA,
                         SymmetricMatrix* pE_K = NULL);

    void calcE_XC(const RUN_TYPE runType, const SymmetricMatrix& PA,
                  SymmetricMatrix* pE_XC = NULL);
    void calcE_XC_gridfree(const RUN_TYPE runType, const SymmetricMatrix& PA,
                           SymmetricMatrix* pE_XC = NULL);
    void calcE_XC_DIRECT(const SymmetricMatrix& PA, const Vector& eps);

   protected:
    void calcS2();

   protected:
    SymmetricMatrix get_IE_nuc(
        const SymmetricMatrix& matrix,
        const std::vector<std::vector<DfObject::index_type> >& atomIndexGroup)
        const;
    SymmetricMatrix get_IE(
        const SymmetricMatrix& matrix,
        const std::vector<std::vector<DfObject::index_type> >& atomIndexGroup)
        const;

    std::vector<int> atomArray2AoArray(const std::vector<int>& atomArray) const;
    double get_IE_matrix(const SymmetricMatrix& matrix,
                         const std::vector<int>& atomIndeces1,
                         const std::vector<int>& atomIndeces2) const;

   protected:
    TlOrbitalInfo orbInfo_;
    DfXCFunctional* pDfXCFunctionalObject_;
    DfOverlapX* pDfOverlapObject_;

    static const double TOO_SMALL;
    static double E_nuc_;
    static double E_nuc_X_;
    double E_total_;
    double E_h_;
    double E_h_X_;
    double E_J_;
    double E_J_rho_rhoTilde_;
    double E_J_rhoTilde_rhoTilde_;
    double E_K_;
    double E_KA_;
    double E_KB_;
    double E_XC_;
    double E_disp_;

    int calcScfIteration_;
    SymmetricMatrix E_;

    SymmetricMatrix* pPA_;
    SymmetricMatrix* pPB_;

    // for interaction energy
    // SymmetricMatrix* pIE_;
};

// ----------------------------------------------------------------------------
// template
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
const double DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                                DfOverlapType>::TOO_SMALL = 1.0E-16;

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
double DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                          DfOverlapType>::E_nuc_ = 0.0;

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
double DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                          DfOverlapType>::E_nuc_X_ = 0.0;

// constructor & destructor ---------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    DfTotalEnergy_tmpl(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      orbInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]),
      pDfXCFunctionalObject_(NULL),
      pDfOverlapObject_(NULL),
      E_total_(0.0),
      E_h_(0.0),
      E_h_X_(0.0),
      E_J_(0.0),
      E_J_rho_rhoTilde_(0.0),
      E_J_rhoTilde_rhoTilde_(0.0),
      E_K_(0.0),
      E_KA_(0.0),
      E_KB_(0.0),
      E_XC_(0.0),
      E_disp_(0.0),
      pPA_(NULL),
      pPB_(NULL) {}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                   DfOverlapType>::~DfTotalEnergy_tmpl() {
    delete this->pDfXCFunctionalObject_;
    this->pDfXCFunctionalObject_ = NULL;
    delete this->pDfOverlapObject_;
    this->pDfOverlapObject_ = NULL;

    delete this->pPA_;
    this->pPA_ = NULL;
    delete this->pPB_;
    this->pPB_ = NULL;
}

// initialize -----------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::prepareXCFunctionalObject() {
    this->pDfXCFunctionalObject_ = new DfXCFunctional(this->pPdfParam_);
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::prepareDfOverlapObject() {
    this->pDfOverlapObject_ = new DfOverlapX(this->pPdfParam_);
}

// ----------------------------------------------------------------------------
// density matrix
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
SymmetricMatrix
DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                   DfOverlapType>::getSpinDensityMatrix(const RUN_TYPE runType,
                                                        const int iteration) {
    SymmetricMatrix matrix =
        DfObject::getSpinDensityMatrix<SymmetricMatrix>(runType, iteration);
    return matrix;
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
Vector DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                          DfOverlapType>::getRho() {
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

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
Vector DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                          DfOverlapType>::getEps(const RUN_TYPE runType) {
    Vector eps;
    switch (runType) {
        case RUN_RKS:
            eps.load("fl_Work/fl_Vct_Epsilon");
            assert(static_cast<TlVectorAbstract::size_type>(
                       this->numOfAuxXC_) == eps.getSize());
            break;

        case RUN_UKS_ALPHA:
            eps.load("fl_Work/fl_Vct_Epsilona");
            assert(static_cast<TlVectorAbstract::size_type>(
                       this->numOfAuxXC_) == eps.getSize());
            break;

        case RUN_UKS_BETA:
            eps.load("fl_Work/fl_Vct_Epsilonb");
            assert(static_cast<TlVectorAbstract::size_type>(
                       this->numOfAuxXC_) == eps.getSize());
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
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::run() {
    this->calc(this->m_nIteration);
    this->saveParams();
    this->output();
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::calc(const int iteration) {
    // calc parts
    this->calcE_nuc();

    this->calcScfIteration_ = iteration;

    SymmetricMatrix P_AB;
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            const SymmetricMatrix PA =
                this->getSpinDensityMatrix(RUN_RKS, iteration);

            // K
            {
                SymmetricMatrix E_K(this->m_nNumOfAOs);
                this->E_K_ = this->calcE_K(RUN_RKS, PA, &E_K);

                this->E_K_ *= 2.0;
                E_K *= 2.0;
                E_K.save(DfObject::m_sWorkDirPath + "/" + "e_K.rks.mat");
            }
            // XC
            {
                SymmetricMatrix E_XC(this->m_nNumOfAOs);
                this->calcE_XC(RUN_RKS, PA, &E_XC);

                E_XC *= 2.0;
                E_XC.save(DfObject::m_sWorkDirPath + "/" + "e_XC.rks.mat");
            }

            P_AB = 2.0 * PA;
        } break;

        case METHOD_UKS: {
            const SymmetricMatrix PA =
                this->getSpinDensityMatrix(RUN_UKS_ALPHA, iteration);
            const SymmetricMatrix PB =
                this->getSpinDensityMatrix(RUN_UKS_BETA, iteration);

            // K
            {
                SymmetricMatrix E_KA(this->m_nNumOfAOs),
                    E_KB(this->m_nNumOfAOs);
                this->E_KA_ = this->calcE_K(RUN_UKS_ALPHA, PA, &E_KA);
                this->E_KB_ = this->calcE_K(RUN_UKS_BETA, PB, &E_KB);
                this->E_K_ = this->E_KA_ + this->E_KB_;
                E_KA.save(DfObject::m_sWorkDirPath + "/" + "e_K.uks-alpha.mat");
                E_KA.save(DfObject::m_sWorkDirPath + "/" + "e_K.uks-beta.mat");
            }
            // XC
            {
                SymmetricMatrix E_XC_A(this->m_nNumOfAOs),
                    E_XC_B(this->m_nNumOfAOs);
                this->calcE_XC(RUN_UKS_ALPHA, PA, &E_XC_A);
                this->calcE_XC(RUN_UKS_BETA, PB, &E_XC_B);
                E_XC_A.save(DfObject::m_sWorkDirPath + "/" +
                            "e_XC.uks-alpha.mat");
                E_XC_B.save(DfObject::m_sWorkDirPath + "/" +
                            "e_XC.uks-beta.mat");
            }

            P_AB = PA + PB;
        } break;

        case METHOD_ROKS: {
            SymmetricMatrix PB =
                this->getSpinDensityMatrix(RUN_ROKS_CLOSED, iteration);
            PB *= 0.5;
            SymmetricMatrix PA =
                this->getSpinDensityMatrix(RUN_ROKS_OPEN, iteration);
            PA += PB;

            {
                SymmetricMatrix E_KA(this->m_nNumOfAOs),
                    E_KB(this->m_nNumOfAOs);
                this->E_KA_ = this->calcE_K(RUN_ROKS_ALPHA, PA, &E_KA);
                this->E_KB_ = this->calcE_K(RUN_ROKS_BETA, PB, &E_KB);
                this->E_K_ = this->E_KA_ + this->E_KB_;
                E_KA.save(DfObject::m_sWorkDirPath + "/" +
                          "e_K.roks-alpha.mat");
                E_KA.save(DfObject::m_sWorkDirPath + "/" + "e_K.roks-beta.mat");
            }
            {
                SymmetricMatrix E_XC_A(this->m_nNumOfAOs),
                    E_XC_B(this->m_nNumOfAOs);
                this->calcE_XC(RUN_ROKS_ALPHA, PA, &E_XC_A);
                this->calcE_XC(RUN_ROKS_BETA, PB, &E_XC_B);
                E_XC_A.save(DfObject::m_sWorkDirPath + "/" +
                            "e_XC.roks-alpha.mat");
                E_XC_B.save(DfObject::m_sWorkDirPath + "/" +
                            "e_XC.roks-beta.mat");
            }
            P_AB = PA + PB;
        } break;

        default:
            break;
    }

    this->calcE_h(P_AB);
    this->calcE_h_fromX(P_AB);
    this->calcE_J(P_AB);

    if (this->enableGrimmeDispersion_ == true) {
        this->E_disp_ =
            this->pDfXCFunctionalObject_->getGrimmeDispersionEnergy();
    }

    // calc total energy
    double E_total = 0.0;
    E_total += this->E_h_;
    E_total += this->E_h_X_;
    if (this->J_engine_ == J_ENGINE_RI_J) {
        E_total += this->E_J_rho_rhoTilde_;
        E_total += this->E_J_rhoTilde_rhoTilde_;
    } else {
        E_total += this->E_J_;
    }
    E_total += this->E_K_;
    E_total += this->E_XC_;
    E_total += this->E_nuc_;
    E_total += this->E_nuc_X_;
    this->E_total_ = E_total;
}

// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::saveParams() {
    (*this->pPdfParam_)["TEs"][this->m_nIteration] = this->E_total_;

    (*this->pPdfParam_)["total_energy"][this->m_nIteration]["h"] = this->E_h_;
    (*this->pPdfParam_)["total_energy"][this->m_nIteration]["h_X"] =
        this->E_h_X_;
    (*this->pPdfParam_)["total_energy"][this->m_nIteration]["J"] = this->E_J_;
    (*this->pPdfParam_)["total_energy"][this->m_nIteration]["J_rr~"] =
        this->E_J_rho_rhoTilde_;
    (*this->pPdfParam_)["total_energy"][this->m_nIteration]["J_r~r~"] =
        this->E_J_rhoTilde_rhoTilde_;
    (*this->pPdfParam_)["total_energy"][this->m_nIteration]["K"] = this->E_K_;
    (*this->pPdfParam_)["total_energy"][this->m_nIteration]["XC"] = this->E_XC_;
    (*this->pPdfParam_)["total_energy"]["nuc"] = this->E_nuc_;
    (*this->pPdfParam_)["total_energy"]["nuc_X"] = this->E_nuc_X_;
    (*this->pPdfParam_)["total_energy"]["disp"] = this->E_disp_;
}

// ----------------------------------------------------------------------------
// IE
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    get_IE(const std::vector<std::vector<DfObject::index_type> >&
               atomIndexGroup) const {
    // nuc
    SymmetricMatrix IE_nuc;
    {
        SymmetricMatrix E_nuc;
        E_nuc.load(DfObject::m_sWorkDirPath + "/" + "e_nuc.mat");

        IE_nuc = this->get_IE_nuc(E_nuc, atomIndexGroup);
        IE_nuc.save(DfObject::m_sWorkDirPath + "/" + "IE_nuc.mat");
    }

    SymmetricMatrix IE_e;
    SymmetricMatrix E_e(this->m_nNumOfAOs);
    // K
    {
        SymmetricMatrix E_K;
        E_K.load(DfObject::m_sWorkDirPath + "/" + "e_K.rks.mat");
        E_e += E_K;

        SymmetricMatrix IE_K = this->get_IE(E_K, atomIndexGroup);
        IE_K.save(DfObject::m_sWorkDirPath + "/" + "IE_K.rks.mat");
    }
    // XC
    {
        SymmetricMatrix E_XC;
        E_XC.load(DfObject::m_sWorkDirPath + "/" + "e_XC.rks.mat");
        E_e += E_XC;

        SymmetricMatrix IE_K = this->get_IE(E_XC, atomIndexGroup);
        IE_K.save(DfObject::m_sWorkDirPath + "/" + "IE_XC.rks.mat");
    }

    // h
    {
        SymmetricMatrix E_h;
        E_h.load(DfObject::m_sWorkDirPath + "/" + "e_h.mat");
        E_e += E_h;

        SymmetricMatrix IE_h = this->get_IE(E_h, atomIndexGroup);
        IE_h.save(DfObject::m_sWorkDirPath + "/" + "IE_h.mat");
    }
    if (this->m_nNumOfDummyAtoms > 0) {
        SymmetricMatrix E_h_X;
        E_h_X.load(DfObject::m_sWorkDirPath + "/" + "e_h_X.mat");
        E_e += E_h_X;

        SymmetricMatrix IE_h_X = this->get_IE(E_h_X, atomIndexGroup);
        IE_h_X.save(DfObject::m_sWorkDirPath + "/" + "IE_h_X.mat");
    }

    // J
    {
        SymmetricMatrix E_J;
        E_J.load(DfObject::m_sWorkDirPath + "/" + "e_J.mat");
        E_e += E_J;

        SymmetricMatrix IE_J = this->get_IE(E_J, atomIndexGroup);
        IE_J.save(DfObject::m_sWorkDirPath + "/" + "IE_J.mat");
    }

    // electron-term
    {
        IE_e = this->get_IE(E_e, atomIndexGroup);
        IE_e.save(DfObject::m_sWorkDirPath + "/" + "IE_e.mat");
    }

    const SymmetricMatrix IE = IE_nuc + IE_e;
    IE.save(DfObject::m_sWorkDirPath + "/" + "IE.mat");
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
SymmetricMatrix
DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    get_IE_nuc(const SymmetricMatrix& matrix,
               const std::vector<std::vector<DfObject::index_type> >&
                   atomIndexGroup) const {
    const int numOfGroups = atomIndexGroup.size();
    SymmetricMatrix IE(numOfGroups);
    for (int groupID1 = 0; groupID1 < numOfGroups; ++groupID1) {
        // non-diagonal
        for (int groupID2 = 0; groupID2 < groupID1; ++groupID2) {
            double e = this->get_IE_matrix(matrix, atomIndexGroup[groupID1],
                                           atomIndexGroup[groupID2]);
            e *= 2.0;
            IE.set(groupID1, groupID2, e);
        }

        // diagonal
        {
            const double e = this->get_IE_matrix(
                matrix, atomIndexGroup[groupID1], atomIndexGroup[groupID1]);
            IE.set(groupID1, groupID1, e);
        }
    }

    return IE;
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
SymmetricMatrix
DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    get_IE(const SymmetricMatrix& matrix,
           const std::vector<std::vector<DfObject::index_type> >&
               atomIndexGroup) const {
    const int numOfGroups = atomIndexGroup.size();
    SymmetricMatrix IE(numOfGroups);
    for (int groupID1 = 0; groupID1 < numOfGroups; ++groupID1) {
        const std::vector<int> aoIndeces1 =
            this->atomArray2AoArray(atomIndexGroup[groupID1]);

        // non-diagonal
        for (int groupID2 = 0; groupID2 < groupID1; ++groupID2) {
            const std::vector<int> aoIndeces2 =
                this->atomArray2AoArray(atomIndexGroup[groupID2]);
            double e = this->get_IE_matrix(matrix, aoIndeces1, aoIndeces2);
            e *= 2.0;
            IE.set(groupID1, groupID2, e);
        }

        // diagonal
        {
            const double e =
                this->get_IE_matrix(matrix, aoIndeces1, aoIndeces1);
            IE.set(groupID1, groupID1, e);
        }
    }

    return IE;
}

// template <class GeneralMatrix, class SymmetricMatrix, class Vector,
//           class DfOverlapType>
// double DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
//                           DfOverlapType>::get_IE(const std::vector<int>&
//                                                      atomIndeces1) const {
//     double ie_nuc =
//         this->get_IE_matrix(*(this->pIE_nuc_), atomIndeces1, atomIndeces1);

//     const std::vector<int> aoIndeces1 =
//     this->atomArray2AoArray(atomIndeces1); double ie =
//     this->get_IE_matrix(*(this->pIE_), aoIndeces1, aoIndeces1);

//     double answer = ie + ie_nuc;
//     std::cout << TlUtils::format("IE(nuc): %28.16lf", ie_nuc) << std::endl;
//     std::cout << TlUtils::format("IE(nie): %28.16lf", ie) << std::endl;
//     std::cout << TlUtils::format("IE     : %28.16lf", answer) << std::endl;

//     return answer;
// }

// template <class GeneralMatrix, class SymmetricMatrix, class Vector,
//           class DfOverlapType>
// double DfTotalEnergy_tmpl<
//     GeneralMatrix, SymmetricMatrix, Vector,
//     DfOverlapType>::get_IE(const std::vector<int>& atomIndeces1,
//                            const std::vector<int>& atomIndeces2) const {
//     double ie_nuc =
//         this->get_IE_matrix(*(this->pIE_nuc_), atomIndeces1, atomIndeces2);

//     const std::vector<int> aoIndeces1 =
//     this->atomArray2AoArray(atomIndeces1); const std::vector<int> aoIndeces2
//     = this->atomArray2AoArray(atomIndeces2); double ie =
//     this->get_IE_matrix(*(this->pIE_), aoIndeces1, aoIndeces2);

//     ie_nuc *= 2.0;
//     ie *= 2.0;
//     double answer = ie + ie_nuc;
//     std::cout << TlUtils::format("IE(nuc): %28.16lf", ie_nuc) << std::endl;
//     std::cout << TlUtils::format("IE(nie): %28.16lf", ie) << std::endl;
//     std::cout << TlUtils::format("IE     : %28.16lf", answer) << std::endl;

//     return answer;
// }

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
double DfTotalEnergy_tmpl<
    GeneralMatrix, SymmetricMatrix, Vector,
    DfOverlapType>::get_IE_matrix(const SymmetricMatrix& matrix,
                                  const std::vector<int>& indeces1,
                                  const std::vector<int>& indeces2) const {
    double IE = 0.0;
    std::vector<int>::const_iterator itEnd1 = indeces1.end();
    std::vector<int>::const_iterator itEnd2 = indeces2.end();
    for (std::vector<int>::const_iterator it1 = indeces1.begin(); it1 != itEnd1;
         ++it1) {
        for (std::vector<int>::const_iterator it2 = indeces2.begin();
             it2 != itEnd2; ++it2) {
            IE += matrix.get(*it1, *it2);
        }
    }

    return IE;
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
std::vector<int> DfTotalEnergy_tmpl<
    GeneralMatrix, SymmetricMatrix, Vector,
    DfOverlapType>::atomArray2AoArray(const std::vector<int>& atomArray) const {
    std::vector<int> aoArray;
    const int numOfAOs = this->orbInfo_.getNumOfOrbitals();
    for (int ao = 0; ao < numOfAOs; ++ao) {
        const int atomIndex = this->orbInfo_.getAtomIndex(ao);

        std::vector<int>::const_iterator itEnd = atomArray.end();
        for (std::vector<int>::const_iterator it = atomArray.begin();
             it != itEnd; ++it) {
            if (*it == atomIndex) {
                aoArray.push_back(ao);
                break;
            }
        }
    }

    return aoArray;
}

// ----------------------------------------------------------------------------
// energy for nuclear repulsion
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::calcE_nuc() {
    if (std::fabs(this->E_nuc_) > TOO_SMALL) {
        this->logger(
            " energy of nuclear repulstion has already been calculated.");
    } else {
        // read nuclear charge
        const Fl_Geometry geom((*this->pPdfParam_)["coordinates"]);
        const int numOfAtoms = this->m_nNumOfAtoms;

        SymmetricMatrix E_nuc(numOfAtoms);

        // calculate nuclear repulsion
        double e_nuc = 0.0;
        double e_nuc_X = 0.0;
        for (int i = 0; i < numOfAtoms; ++i) {
            const double ci = geom.getCharge(i);
            const TlPosition pi = geom.getCoordinate(i);
            const bool isDummy = ("X" == geom.getAtomSymbol(i));

            for (int j = 0; j < i; ++j) {
                const double cj = geom.getCharge(j);
                const TlPosition pj = geom.getCoordinate(j);
                const double distance = pi.distanceFrom(pj);

                const double e = ci * cj / distance;
                if (isDummy || ("X" == geom.getAtomSymbol(j))) {
                    e_nuc_X += e;
                } else {
                    e_nuc += e;
                }
                E_nuc.set(i, j, 0.5 * e);
            }
        }

        this->E_nuc_ = e_nuc;
        this->E_nuc_X_ = e_nuc_X;
        E_nuc.save(DfObject::m_sWorkDirPath + "/" + "e_nuc.mat");
    }
}

// ----------------------------------------------------------------------------
// H
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::calcE_h(const SymmetricMatrix& P_AB) {
    SymmetricMatrix Hpq = DfObject::getHpqMatrix<SymmetricMatrix>();
    if (Hpq.getNumOfRows() != this->m_nNumOfAOs ||
        Hpq.getNumOfCols() != this->m_nNumOfAOs) {
        CnErr.abort("DfTotalEnergy", "calculate_one_electron_part", "",
                    "program error");
    }

    const SymmetricMatrix E_h = Hpq.dotInPlace(P_AB);

    this->E_h_ = E_h.sum();
    E_h.save(DfObject::m_sWorkDirPath + "/" + "e_h.mat");
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    calcE_h_fromX(const SymmetricMatrix& P_AB) {
    // add dummy charge
    if (this->m_nNumOfDummyAtoms > 0) {
        const int chargeExtrapolateNumber =
            (*(this->pPdfParam_))["charge-extrapolate-number"].getInt();

        SymmetricMatrix Hpq2 = DfObject::getHpq2Matrix<SymmetricMatrix>();

        if (chargeExtrapolateNumber > 1) {
            const int coef =
                std::min(this->m_nIteration, chargeExtrapolateNumber);
            Hpq2 *= coef;
        }

        const SymmetricMatrix E_h_X = Hpq2.dotInPlace(P_AB);

        this->E_h_X_ = E_h_X.sum();
        E_h_X.save(DfObject::m_sWorkDirPath + "/" + "e_h_X.mat");
    }
}

// ----------------------------------------------------------------------------
// J term
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::calcE_J(const SymmetricMatrix& P_AB) {
    switch (this->J_engine_) {
        case J_ENGINE_RI_J: {
            this->calcE_J_RI(P_AB);
        } break;

        case J_ENGINE_CONVENTIONAL:
        case J_ENGINE_CD:
            this->calcE_J_exact(P_AB);
            break;

        default:
            break;
    }
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    calcE_J_exact(const SymmetricMatrix& P_AB) {
    SymmetricMatrix J =
        DfObject::getJMatrix<SymmetricMatrix>(this->calcScfIteration_);

    const SymmetricMatrix E_J = 0.5 * J.dotInPlace(P_AB);

    this->E_J_ = E_J.sum();
    E_J.save(DfObject::m_sWorkDirPath + "/" + "e_J.mat");
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    calcE_J_RI(const SymmetricMatrix& P_AB) {
    this->calcE_J_rhoTilde_rhoTilde();

    if (!((this->m_bMemorySave == false) &&
          (this->m_bDiskUtilization == false))) {
        // NOT use Threeindexintegrals
        this->calcE_J_rho_rhoTilde_DIRECT(P_AB);
    }
}

// energy for a part of coulomb term ( rho*Ppq*Pqa )
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    calcE_J_rho_rhoTilde_DIRECT(const SymmetricMatrix& P_AB) {
    SymmetricMatrix J =
        DfObject::getJMatrix<SymmetricMatrix>(this->calcScfIteration_);

    const double coef = (this->J_engine_ == J_ENGINE_RI_J) ? 1.0 : 0.5;
    SymmetricMatrix E_J = coef * J.dotInPlace(P_AB);
    this->E_J_rho_rhoTilde_ = E_J.sum();
}

/// J[rho~, rho~] (( "1/2*rho*rho'*Sab" ))の計算を行う
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::calcE_J_rhoTilde_rhoTilde() {
    SymmetricMatrix Sab = DfObject::getSabMatrix<SymmetricMatrix>();
    const Vector rho = this->getRho();
    this->E_J_rhoTilde_rhoTilde_ = -0.5 * (rho.dot(Sab * rho));
}

// ----------------------------------------------------------------------------
// K term
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
double DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                          DfOverlapType>::calcE_K(const RUN_TYPE runType,
                                                  const SymmetricMatrix& PA,
                                                  SymmetricMatrix* pE_K) {
    return this->calcE_K_exact(runType, PA, pE_K);
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
double
DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                   DfOverlapType>::calcE_K_exact(const RUN_TYPE runType,
                                                 const SymmetricMatrix& PA,
                                                 SymmetricMatrix* pE_K) {
    double answer = 0.0;

    DfXCFunctional dfXCFunctional(this->pPdfParam_);
    if (dfXCFunctional.isHybridFunctional() == true) {
        const SymmetricMatrix E_K =
            0.5 * dfXCFunctional.getFockExchangeCoefficient() *
            DfObject::getHFxMatrix<SymmetricMatrix>(runType,
                                                    this->calcScfIteration_)
                .dotInPlace(PA);

        answer = E_K.sum();
        if (pE_K != NULL) {
            *pE_K = E_K;
        }
    }

    return answer;
}

// ----------------------------------------------------------------------------
// XC term
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::calcE_XC(const RUN_TYPE runType,
                                                 const SymmetricMatrix& PA,
                                                 SymmetricMatrix* pE_XC) {
    DfXCFunctional dfXCFunctional(this->pPdfParam_);
    if (dfXCFunctional.getXcType() != DfXCFunctional::HF) {
        if (this->m_bIsXCFitting == true) {
            const Vector Eps = this->getEps(runType);
            this->calcE_XC_DIRECT(PA, Eps);
        } else {
            if (this->XC_engine_ != XC_ENGINE_GRID) {
                this->calcE_XC_gridfree(runType, PA, pE_XC);
            } else {
                this->E_XC_ = dfXCFunctional.getEnergy();
            }
        }
    }
}

template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    calcE_XC_gridfree(const RUN_TYPE runType, const SymmetricMatrix& PA,
                      SymmetricMatrix* pE_XC) {
    SymmetricMatrix E_XC = DfObject::getExcMatrix<SymmetricMatrix>(
        runType, this->calcScfIteration_);
    E_XC.dotInPlace(PA);

    if (runType == DfObject::RUN_RKS) {
        E_XC *= 2.0;
    }

    this->E_XC_ += E_XC.sum();
    if (pE_XC != NULL) {
        *pE_XC = E_XC;
    }
}

// energy for xc energy term (compare myu*Ppq*Pqa with 4/3*Ex1)
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector, DfOverlapType>::
    calcE_XC_DIRECT(const SymmetricMatrix& PA, const Vector& eps) {
    // put eps*[pqs] into B
    SymmetricMatrix pqg(this->m_nNumOfAOs);

    // TODO: lapack only
    // this->pDfOverlapObject_->get_pqg(eps, &pqg);

    this->E_XC_ = pqg.dotInPlace(PA).sum();
}

// ----------------------------------------------------------------------------
// S2
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::calcS2() {
    // calc (N_alpha - N_beta / 2)
    const double N_alpha = this->m_nNumOfAlphaElectrons;
    const double N_beta = this->m_nNumOfBetaElectrons;
    const double Nab2 = 0.5 * (N_alpha - N_beta);

    // S^2_exact
    double s2_exact = Nab2 * (Nab2 + 1.0);

    // S^2 (UKS)
    double csc2 = 0.0;
    {
        const int itr = this->m_nIteration;
        GeneralMatrix C_alpha =
            DfObject::getCMatrix<GeneralMatrix>(RUN_UKS_ALPHA, itr);
        GeneralMatrix C_beta =
            DfObject::getCMatrix<GeneralMatrix>(RUN_UKS_BETA, itr);
        const SymmetricMatrix S = DfObject::getSpqMatrix<SymmetricMatrix>();

        C_alpha.resize(this->m_nNumOfAOs, N_alpha);
        C_beta.resize(this->m_nNumOfAOs, N_beta);

        const GeneralMatrix CSC = C_alpha.transpose() * S * C_beta;
        const GeneralMatrix CSC2 = CSC.dot(CSC);
        csc2 = CSC2.sum();
    }

    const double s2 = s2_exact + N_beta - csc2;
    this->log_.info("calculate S^2:");
    this->log_.info(TlUtils::format(" S^2            = %28.16lf", s2));
    this->log_.info(TlUtils::format(" S^2 (expected) = %28.16lf", s2_exact));
}

// ----------------------------------------------------------------------------
// output
// ----------------------------------------------------------------------------
template <class GeneralMatrix, class SymmetricMatrix, class Vector,
          class DfOverlapType>
void DfTotalEnergy_tmpl<GeneralMatrix, SymmetricMatrix, Vector,
                        DfOverlapType>::output() {
    double E_Total = 0.0;

    // 表示
    this->log_.info("------------------------------------------------");

    // 1e term
    this->log_.info(TlUtils::format(" Ts+Vn          = %28.16lf", this->E_h_));
    this->log_.info(
        TlUtils::format(" Ts+Vn (X)      = %28.16lf", this->E_h_X_));
    E_Total += this->E_h_ + this->E_h_X_;

    // J term
    switch (this->J_engine_) {
        case J_ENGINE_RI_J: {
            // if ((this->m_bMemorySave == false) &&
            //     (this->m_bDiskUtilization == false)) {
            //   E_Total += this->m_dE_OEP_JRR_Exc;

            //   this->logger(TlUtils::format(" Ts+Vn+J[Rho~,Rho~]+Exc =
            //   %28.16lf",
            //                                this->m_dE_OEP_JRR_Exc));
            // }
            this->log_.info(TlUtils::format(" E_J[Rho, Rho~] = %28.16lf",
                                            this->E_J_rho_rhoTilde_));
            // E_Total += this->m_dE_J_Rho_RhoTilde;
            this->log_.info(TlUtils::format(" E_J[Rho~,Rho~] = %28.16lf",
                                            this->E_J_rhoTilde_rhoTilde_));
            // E_Total += this->m_dE_J_RhoTilde_RhoTilde;
        } break;

        case J_ENGINE_CONVENTIONAL:
        case J_ENGINE_CD: {
            this->log_.info(
                TlUtils::format(" E_J            = %28.16lf", this->E_J_));
            E_Total += this->E_J_;
        } break;

        default:
            break;
    }

    //
    this->log_.info(TlUtils::format(" E_xc(pure)     = %28.16lf", this->E_XC_));
    E_Total += this->E_XC_;

    if (this->enableGrimmeDispersion_ == true) {
        this->log_.info(TlUtils::format(" E_cx(+disp.)   = %28.16lf",
                                        E_Total + this->E_disp_));
    }

    switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->log_.info(
                TlUtils::format(" E_K            = %28.16lf", this->E_K_));
            break;

        case METHOD_UKS:
            this->log_.info(
                TlUtils::format(" E_K            = %28.16lf", this->E_K_));
            this->log_.info(
                TlUtils::format("   E_K(alpha)   = %28.16lf", this->E_KA_));
            this->log_.info(
                TlUtils::format("   E_K(beta)    = %28.16lf", this->E_KB_));
            break;

        case METHOD_ROKS:
            this->log_.info(
                TlUtils::format(" E_K            = %28.16lf", this->E_K_));
            this->log_.info(
                TlUtils::format("   E_K(alpha)   = %28.16lf", this->E_KA_));
            this->log_.info(
                TlUtils::format("   E_K(beta)    = %28.16lf", this->E_KB_));
            break;

        default:
            this->log_.critical("program error");
            break;
    }
    E_Total += this->E_K_;

    this->log_.info(
        TlUtils::format(" E_nuclei       = %28.16lf", this->E_nuc_));
    this->log_.info(
        TlUtils::format(" E_nuclei(X)    = %28.16lf", this->E_nuc_X_));
    E_Total += this->E_nuc_;
    E_Total += this->E_nuc_X_;

    this->log_.info(TlUtils::format(" TE             = %28.16lf", E_Total));
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

    if (this->m_nMethodType != METHOD_RKS) {
        this->calcS2();
    }

    this->log_.info("------------------------------------------------");
}

#endif  // DF_TOTALENERGY_TMPL_H
