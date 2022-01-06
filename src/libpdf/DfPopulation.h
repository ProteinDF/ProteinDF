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

#ifndef DFPOPULATION_H
#define DFPOPULATION_H

#include <numeric>
#include <valarray>
#include <vector>

#include "CnError.h"
#include "DfObject.h"
#include "TlOrbitalInfo.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class TlDenseSymmetricMatrix_Lapack;

/// Mulliken Population の計算を行い、電荷情報、スピン情報を出力するクラス
class DfPopulation : public DfObject {
public:
    DfPopulation(TlSerializeData* pPdfParam);
    virtual ~DfPopulation();

public:
    void exec(const int iteration);

    template <class SymmetricMatrixType, class VectorType>
    void getAtomPopulation(const int iteration);

    // TlDenseGeneralMatrix_Lapack getAtomPopData(const int iteration);

    double getCharge(int atomIndex);

public:
    template <class SymmetricMatrixType>
    double getSumOfElectrons(const SymmetricMatrixType& P);

public:
    template <typename T>
    void getReport(const int iteration, T& out);

protected:
    // init
    // ----------------------------------------------------------------------
    void setNucleiCharges();

    // calc
    // ----------------------------------------------------------------------
    template <class SymmetricMatrixType, class VectorType>
    std::valarray<double> getAtomPopulation_runtype(const RUN_TYPE runType, const int iteration);

    virtual void calcPop(const int iteration);

    template <class SymmetricMatrixType>
    void calcPop(const int iteration);

    template <class SymmetricMatrixType>
    std::valarray<double> getGrossOrbPop(DfObject::RUN_TYPE runType, int iteration);

    template <class SymmetricMatrixType>
    std::valarray<double> getPS(const SymmetricMatrixType& P);

    std::valarray<double> getGrossAtomPop(const std::valarray<double>& grossOrbPop);

    // report
    // --------------------------------------------------------------------
    double getSumOfNucleiCharges() const;

    template <typename T>
    void getAtomPopStr(const std::valarray<double>& trPS, const bool isCalcNetCharge, T& out) const;

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

template <class SymmetricMatrixType>
double DfPopulation::getSumOfElectrons(const SymmetricMatrixType& P) {
    // std::cerr << "DfPopulation::getSumOfElectrons()" << std::endl;
    const std::valarray<double> trPS = this->getPS<SymmetricMatrixType>(P);
    return trPS.sum();
}

template <class SymmetricMatrixType, class VectorType>
void DfPopulation::getAtomPopulation(const int iteration) {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            const VectorType atomPop =
                this->getAtomPopulation_runtype<SymmetricMatrixType, VectorType>(RUN_RKS, iteration);
        } break;

        case METHOD_UKS: {
            const VectorType atomPopA =
                this->getAtomPopulation_runtype<SymmetricMatrixType, VectorType>(RUN_UKS_ALPHA, iteration);

            const VectorType atomPopB =
                this->getAtomPopulation_runtype<SymmetricMatrixType, VectorType>(RUN_UKS_BETA, iteration);
        } break;

        case METHOD_ROKS: {
            const VectorType atomPopA =
                this->getAtomPopulation_runtype<SymmetricMatrixType, VectorType>(RUN_ROKS_CLOSED, iteration);

            const VectorType atomPopB =
                this->getAtomPopulation_runtype<SymmetricMatrixType, VectorType>(RUN_ROKS_OPEN, iteration);
        } break;

        default:
            break;
    }
}

template <class SymmetricMatrixType, class VectorType>
std::valarray<double> DfPopulation::getAtomPopulation_runtype(const RUN_TYPE runType, const int iteration) {
    const std::valarray<double> grossOrbPop = this->getGrossOrbPop<SymmetricMatrixType>(runType, iteration);
    VectorType(grossOrbPop).save(this->getPopGrossOrbPath(runType, iteration));

    const std::valarray<double> grossAtomPop = this->getGrossAtomPop(grossOrbPop);
    VectorType(grossAtomPop).save(this->getPopGrossAtomPath(runType, iteration));

    const std::valarray<double> atomPop = this->nucleiCharges_ - grossAtomPop;
    VectorType(atomPop).save(this->getPopMullikenPath(runType, iteration));

    return atomPop;
}

template <class SymmetricMatrixType>
void DfPopulation::calcPop(const int iteration) {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            this->grossOrbPopA_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_RKS, iteration);
            this->grossAtomPopA_ = this->getGrossAtomPop(this->grossOrbPopA_);
        } break;

        case METHOD_UKS: {
            this->grossOrbPopA_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_UKS_ALPHA, iteration);
            this->grossAtomPopA_ = this->getGrossAtomPop(this->grossOrbPopA_);

            this->grossOrbPopB_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_UKS_BETA, iteration);
            this->grossAtomPopB_ = this->getGrossAtomPop(this->grossOrbPopB_);
        } break;

        case METHOD_ROKS: {
            this->grossOrbPopA_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_ROKS_CLOSED, iteration);
            this->grossAtomPopA_ = this->getGrossAtomPop(this->grossOrbPopA_);

            this->grossOrbPopB_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_ROKS_OPEN, iteration);
            this->grossAtomPopB_ = this->getGrossAtomPop(this->grossOrbPopB_);
        } break;

        default:
            // programer error
            std::abort();
            break;
    }
}

template <class SymmetricMatrixType>
std::valarray<double> DfPopulation::getGrossOrbPop(const DfObject::RUN_TYPE runType, const int iteration) {
    SymmetricMatrixType P;
    assert((runType == RUN_RKS) || (runType == RUN_UKS_ALPHA) || (runType == RUN_UKS_BETA) ||
           (runType == RUN_ROKS_CLOSED) || (runType == RUN_ROKS_OPEN));
    P = DfObject::getPpqMatrix<SymmetricMatrixType>(runType, iteration);

    return this->getPS<SymmetricMatrixType>(P);
}

template <class SymmetricMatrixType>
std::valarray<double> DfPopulation::getPS(const SymmetricMatrixType& P) {
    SymmetricMatrixType PS = DfObject::getSpqMatrix<SymmetricMatrixType>();

    PS.dotInPlace(P);

    const DfObject::index_type numOfAOs = this->m_nNumOfAOs;
    std::valarray<double> answer(numOfAOs);
    for (DfObject::index_type i = 0; i < numOfAOs; ++i) {
        const std::vector<double> v = PS.getColVector(i);
        answer[i] = std::accumulate(v.begin(), v.end(), 0.0);
        // answer[i] = v.sum();
    }

    return answer;
}

template <typename T>
void DfPopulation::getReport(const int iteration, T& out) {
    this->calcPop(iteration);

    const double nucleiCharge = this->getSumOfNucleiCharges();

    out << " Mulliken Population\n";
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            const double elecCharge = this->grossOrbPopA_.sum();
            const double netCharge = nucleiCharge - elecCharge;

            out << TlUtils::format(" total atomic charge = %7.2lf\n", nucleiCharge);
            out << TlUtils::format("          Net charge = %7.2lf\n", netCharge);
            out << "\n";

            out << " Gross atom population\n";
            this->getAtomPopStr(this->grossAtomPopA_, true, out);
            out << "\n";

            out << " Mulliken population analysis (orbial population)\n";
            this->getOrbPopStr(grossOrbPopA_, out);
        } break;

        case METHOD_UKS: {
            const double elecChargeA = this->grossOrbPopA_.sum();
            const double elecChargeB = this->grossOrbPopB_.sum();

            const double netCharge = nucleiCharge - (elecChargeA + elecChargeB);
            out << TlUtils::format(" total atomic charge = %7.2lf\n", nucleiCharge);
            out << TlUtils::format("          Net charge = %7.2lf\n", netCharge);
            out << "\n";

            out << " Mulliken Population(alpha)\n";
            this->getAtomPopStr(this->grossAtomPopA_, false, out);
            this->getOrbPopStr(this->grossOrbPopA_, out);

            out << " Mulliken Population(beta)\n";
            this->getAtomPopStr(this->grossAtomPopB_, false, out);
            this->getOrbPopStr(this->grossOrbPopB_, out);

            {
                const std::valarray<double> grossOrbPop = this->grossOrbPopA_ + this->grossOrbPopB_;
                out << " Gross atom population (UKS total)\n";
                this->getAtomPopStr(grossOrbPop, false, out);
                out << " Mulliken population analysis (orbial population; UKS "
                       "total)\n";
                this->getOrbPopStr(grossOrbPop, out);
            }

            {
                const std::valarray<double> grossOrbPopSpin = this->grossOrbPopA_ - this->grossOrbPopB_;
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
            out << TlUtils::format(" total atomic charge = %7.2lf\n", nucleiCharge);
            out << TlUtils::format("          Net charge = %7.2lf\n", netCharge);
            out << "\n";

            out << " Mulliken Population(closed)\n";
            this->getAtomPopStr(this->grossAtomPopA_, false, out);
            this->getOrbPopStr(this->grossOrbPopA_, out);

            out << " Mulliken Population(open)\n";
            this->getAtomPopStr(this->grossAtomPopB_, false, out);
            this->getOrbPopStr(this->grossOrbPopB_, out);

            {
                const std::valarray<double> grossOrbPop = this->grossOrbPopA_ + this->grossOrbPopB_;
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

template <typename T>
void DfPopulation::getAtomPopStr(const std::valarray<double>& trPS, const bool isCalcNetCharge, T& out) const {
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
        out << TlUtils::format(" %6d   %-2s ", atomIndex, flGeom.getAtomSymbol(atomIndex).c_str());

        const double grossAtomPop = trPS[atomIndex];
        if (isCalcNetCharge == true) {
            out << TlUtils::format("% 12.6lf % 12.6lf\n", grossAtomPop, nucCharge - grossAtomPop);
        } else {
            out << TlUtils::format("% 12.6lf\n", grossAtomPop);
        }
    }
}

// print out Mulliken Analysis (Orbital Population)
template <typename T>
void DfPopulation::getOrbPopStr(const std::valarray<double>& trPS, T& out) const {
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(trPS.size() == numOfAOs);

    out << "  index  Atom        Shell     gross\n";

    int prevAtomIndex = -1;
    for (index_type aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
        int atomIndex = this->orbitalInfo_.getAtomIndex(aoIndex);

        if (prevAtomIndex != atomIndex) {
            out << TlUtils::format(
                "  %6d   %-2s %6d %-10s ", aoIndex + 1, this->orbitalInfo_.getAtomName(aoIndex).c_str(),
                this->orbitalInfo_.getAtomIndex(aoIndex), this->orbitalInfo_.getBasisTypeName(aoIndex).c_str());
            prevAtomIndex = atomIndex;
        } else {
            out << TlUtils::format("  %6d             %-10s ", aoIndex + 1,
                                   this->orbitalInfo_.getBasisTypeName(aoIndex).c_str());
        }
        out << TlUtils::format(" %12.6lf\n", trPS[aoIndex]);
    }
}

#endif  // DFPOPULATION_H
