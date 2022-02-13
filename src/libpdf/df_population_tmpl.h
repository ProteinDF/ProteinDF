#ifndef DF_POPULATION_TMPL_H
#define DF_POPULATION_TMPL_H

#include <numeric>
#include <valarray>

#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlOrbitalInfo.h"
#include "TlUtils.h"

/// Mulliken Population の計算を行い、電荷情報、スピン情報を出力するクラス
template <class SymmetricMatrix, class Vector>
class DfPopulation_tmpl : public DfObject {
public:
    DfPopulation_tmpl(TlSerializeData* pPdfParam);
    virtual ~DfPopulation_tmpl();

public:
    double getSumOfElectrons(const SymmetricMatrix& P);

    void getAtomPopulation(const int iteration);

protected:
    std::valarray<double> getAtomPopulation_runtype(const RUN_TYPE runType,
                                                    const int iteration);
    std::valarray<double> getGrossOrbPop(const DfObject::RUN_TYPE runType,
                                         const int iteration);
    std::valarray<double> getPS(const SymmetricMatrix& P);

    template <typename T>
    void getReport(const int iteration, T& out);

    template <typename T>
    void getAtomPopStr(const std::valarray<double>& trPS,
                       const bool isCalcNetCharge, T& out) const;

    template <typename T>
    void getOrbPopStr(const std::valarray<double>& trPS, T& out) const;

protected:
    TlOrbitalInfo orbitalInfo_;

    std::valarray<double> nucleiCharges_;

    std::valarray<double> grossOrbPopA_;
    std::valarray<double> grossOrbPopB_;

    std::valarray<double> grossAtomPopA_;
    std::valarray<double> grossAtomPopB_;
};

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector>
DfPopulation_tmpl<SymmetricMatrix, Vector>::DfPopulation_tmpl(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]) {
}

template <class SymmetricMatrix, class Vector>
DfPopulation_tmpl<SymmetricMatrix, Vector>::~DfPopulation_tmpl() {
}

template <class SymmetricMatrix, class Vector>
double DfPopulation_tmpl<SymmetricMatrix, Vector>::getSumOfElectrons(
    const SymmetricMatrix& P) {
    const std::valarray<double> trPS = this->getPS(P);
    return trPS.sum();
}

template <class SymmetricMatrix, class Vector>
void DfPopulation_tmpl<SymmetricMatrix, Vector>::getAtomPopulation(
    const int iteration) {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            const Vector atomPop =
                this->getAtomPopulation_runtype(RUN_RKS, iteration);
        } break;

        case METHOD_UKS: {
            const Vector atomPopA =
                this->getAtomPopulation_runtype(RUN_UKS_ALPHA, iteration);

            const Vector atomPopB =
                this->getAtomPopulation_runtype(RUN_UKS_BETA, iteration);
        } break;

        case METHOD_ROKS: {
            const Vector atomPopA =
                this->getAtomPopulation_runtype(RUN_ROKS_CLOSED, iteration);

            const Vector atomPopB =
                this->getAtomPopulation_runtype(RUN_ROKS_OPEN, iteration);
        } break;

        default:
            break;
    }
}

template <class SymmetricMatrix, class Vector>
std::valarray<double>
DfPopulation_tmpl<SymmetricMatrix, Vector>::getAtomPopulation_runtype(
    const RUN_TYPE runType, const int iteration) {
    const std::valarray<double> grossOrbPop =
        this->getGrossOrbPop(runType, iteration);
    Vector(grossOrbPop).save(this->getPopGrossOrbPath(runType, iteration));

    const std::valarray<double> grossAtomPop =
        this->getGrossAtomPop(grossOrbPop);
    Vector(grossAtomPop).save(this->getPopGrossAtomPath(runType, iteration));

    const std::valarray<double> atomPop = this->nucleiCharges_ - grossAtomPop;
    Vector(atomPop).save(this->getPopMullikenPath(runType, iteration));

    return atomPop;
}

template <class SymmetricMatrix, class Vector>
std::valarray<double>
DfPopulation_tmpl<SymmetricMatrix, Vector>::getGrossOrbPop(
    const DfObject::RUN_TYPE runType, const int iteration) {
    SymmetricMatrix P;
    assert((runType == RUN_RKS) || (runType == RUN_UKS_ALPHA) ||
           (runType == RUN_UKS_BETA) || (runType == RUN_ROKS_CLOSED) ||
           (runType == RUN_ROKS_OPEN));
    P = DfObject::getPpqMatrix<SymmetricMatrix>(runType, iteration);

    return this->getPS(P);
}

template <class SymmetricMatrix, class Vector>
std::valarray<double> DfPopulation_tmpl<SymmetricMatrix, Vector>::getPS(
    const SymmetricMatrix& P) {
    SymmetricMatrix PS = DfObject::getSpqMatrix<SymmetricMatrix>();

    PS.dotInPlace(P);

    const DfObject::index_type numOfAOs = this->m_nNumOfAOs;
    std::valarray<double> answer(numOfAOs);
    for (DfObject::index_type i = 0; i < numOfAOs; ++i) {
        const std::vector<double> v = PS.getColVector(i);
        answer[i] = std::accumulate(v.begin(), v.end(), 0.0);
    }

    return answer;
}

// -----------------------------------------------------------------------------
// OUTPUT
// -----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector>
template <typename T>
void DfPopulation_tmpl<SymmetricMatrix, Vector>::getReport(const int iteration,
                                                           T& out) {
    this->calcPop(iteration);

    const double nucleiCharge = this->getSumOfNucleiCharges();

    out << " Mulliken Population\n";
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            const double elecCharge = this->grossOrbPopA_.sum();
            const double netCharge = nucleiCharge - elecCharge;

            out << TlUtils::format(" total atomic charge = %7.2lf\n",
                                   nucleiCharge);
            out << TlUtils::format("          Net charge = %7.2lf\n",
                                   netCharge);
            out << "\n";

            out << " Gross atom population\n";
            this->getAtomPopStr(this->grossAtomPopA_, true, out);
            out << "\n";

            out << " Mulliken population analysis (orbial population)\n";
            this->getOrbPopStr(this->grossOrbPopA_, out);
        } break;

        case METHOD_UKS: {
            const double elecChargeA = this->grossOrbPopA_.sum();
            const double elecChargeB = this->grossOrbPopB_.sum();

            const double netCharge = nucleiCharge - (elecChargeA + elecChargeB);
            out << TlUtils::format(" total atomic charge = %7.2lf\n",
                                   nucleiCharge);
            out << TlUtils::format("          Net charge = %7.2lf\n",
                                   netCharge);
            out << "\n";

            out << " Mulliken Population(alpha)\n";
            this->getAtomPopStr(this->grossAtomPopA_, false, out);
            this->getOrbPopStr(this->grossOrbPopA_, out);

            out << " Mulliken Population(beta)\n";
            this->getAtomPopStr(this->grossAtomPopB_, false, out);
            this->getOrbPopStr(this->grossOrbPopB_, out);

            {
                const std::valarray<double> grossOrbPop =
                    this->grossOrbPopA_ + this->grossOrbPopB_;
                out << " Gross atom population (UKS total)\n";
                this->getAtomPopStr(grossOrbPop, false, out);
                out << " Mulliken population analysis (orbial population; UKS "
                       "total)\n";
                this->getOrbPopStr(grossOrbPop, out);
            }

            {
                const std::valarray<double> grossOrbPopSpin =
                    this->grossOrbPopA_ - this->grossOrbPopB_;
                out << " Spin density analysis for atom (alpha - beta)\n";
                this->getAtomPopStr(grossOrbPopSpin, false, out);
                out << " Spin density analysis for orbital (alpha - beta)\n";
                this->getOrbPopStr(grossOrbPopSpin, out);
            }
        } break;

        case METHOD_ROKS: {
            const double elecChargeA = this->grossOrbPopA_.sum();
            const double elecChargeB = this->grossOrbPopB_.sum();

            const double netCharge = nucleiCharge - (elecChargeA + elecChargeB);
            out << TlUtils::format(" total atomic charge = %7.2lf\n",
                                   nucleiCharge);
            out << TlUtils::format("          Net charge = %7.2lf\n",
                                   netCharge);
            out << "\n";

            out << " Mulliken Population(closed)\n";
            this->getAtomPopStr(this->grossAtomPopA_, false, out);
            this->getOrbPopStr(this->grossOrbPopA_, out);

            out << " Mulliken Population(open)\n";
            this->getAtomPopStr(this->grossAtomPopB_, false, out);
            this->getOrbPopStr(this->grossOrbPopB_, out);

            {
                const std::valarray<double> grossOrbPop =
                    this->grossOrbPopA_ + this->grossOrbPopB_;
                out << " Gross atom population (ROKS total)\n";
                this->getAtomPopStr(grossOrbPop, false, out);
                out << " Mulliken population analysis (orbial population; ROKS "
                       "total)\n";
                this->getOrbPopStr(grossOrbPop, out);
            }
        } break;

        default:
            // programer error
            std::abort();
            break;
    }
}

template <class SymmetricMatrix, class Vector>
template <typename T>
void DfPopulation_tmpl<SymmetricMatrix, Vector>::getAtomPopStr(
    const std::valarray<double>& trPS, const bool isCalcNetCharge,
    T& out) const {
    const index_type numOfAtoms = this->m_nNumOfAtoms;
    assert(trPS.size() == numOfAtoms);

    if (isCalcNetCharge == true) {
        out << "  index Atom    gross         net\n";
    } else {
        out << "  index Atom    gross\n";
    }

    const Fl_Geometry flGeom((*(this->pPdfParam_))["coordinates"]);
    for (index_type atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const double nucCharge = flGeom.getCharge(atomIndex);
        out << TlUtils::format(" %6d   %-2s ", atomIndex,
                               flGeom.getAtomSymbol(atomIndex).c_str());

        const double grossAtomPop = trPS[atomIndex];
        if (isCalcNetCharge == true) {
            out << TlUtils::format("% 12.6lf % 12.6lf\n", grossAtomPop,
                                   nucCharge - grossAtomPop);
        } else {
            out << TlUtils::format("% 12.6lf\n", grossAtomPop);
        }
    }
}

// print out Mulliken Analysis (Orbital Population)
template <class SymmetricMatrix, class Vector>
template <typename T>
void DfPopulation_tmpl<SymmetricMatrix, Vector>::getOrbPopStr(
    const std::valarray<double>& trPS, T& out) const {
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(trPS.size() == numOfAOs);

    out << "  index  Atom        Shell     gross\n";

    int prevAtomIndex = -1;
    for (index_type aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
        int atomIndex = this->orbitalInfo_.getAtomIndex(aoIndex);

        if (prevAtomIndex != atomIndex) {
            out << TlUtils::format(
                "  %6d   %-2s %6d %-10s ", aoIndex + 1,
                this->orbitalInfo_.getAtomName(aoIndex).c_str(),
                this->orbitalInfo_.getAtomIndex(aoIndex),
                this->orbitalInfo_.getBasisTypeName(aoIndex).c_str());
            prevAtomIndex = atomIndex;
        } else {
            out << TlUtils::format(
                "  %6d             %-10s ", aoIndex + 1,
                this->orbitalInfo_.getBasisTypeName(aoIndex).c_str());
        }
        out << TlUtils::format(" %12.6lf\n", trPS[aoIndex]);
    }
}

#endif  // DF_POPULATION_TMPL_H
