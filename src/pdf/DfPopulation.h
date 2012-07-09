#ifndef DFPOPULATION_H
#define DFPOPULATION_H

#include "DfObject.h"
#include "TlOrbitalInfo.h"
#include "TlVector.h"
#include "TlMatrix.h"

/** Mulliken Population の計算を行い、電荷情報、スピン情報を出力するクラス
 */
class DfPopulation : public DfObject {
public:
    DfPopulation(TlSerializeData* pPdfParam);
    virtual ~DfPopulation();

public:
    template <typename T>
    void getReport(const int iteration, T& out);

    TlMatrix getAtomPopData(int iteration);
    
protected:
    virtual void calcPop(const int iteration);

    template<class SymmetricMatrixType>
    void calcPop(const int iteration);

    double getNucleiCharge();

    template<class SymmetricMatrixType>
    TlVector getGrossOrbPop(DfObject::RUN_TYPE runType, int iteration);

    TlVector getGrossAtomPop(const TlVector& grossOrbPop);
    
    template <typename T>
    void getAtomPopStr(const TlVector& trPS, const bool isCalcNetCharge,
                       T& out) const;
    
    template <typename T>
    void getOrbPopStr(const TlVector& trPS, T& out) const;
    
protected:
    TlOrbitalInfo orbitalInfo_;

    TlVector grossOrbPopA_;
    TlVector grossOrbPopB_;

    TlVector grossAtomPopA_;
    TlVector grossAtomPopB_;
};


template <typename T>
void DfPopulation::getReport(const int iteration, T& out)
{
    this->calcPop(iteration);

    const double nucleiCharge = this->getNucleiCharge();

    out << " Mulliken Population\n";
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        {
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
        }
        break;

    case METHOD_ROKS:
        {
            const double elecCharge = this->grossOrbPopA_.sum();
            const double netCharge = nucleiCharge - elecCharge;

            out << TlUtils::format(" total atomic charge = %7.2lf\n", nucleiCharge);
            out << TlUtils::format("          Net charge = %7.2lf\n", netCharge);
            out << "\n";
            
            out << " Gross atom population\n";
            this->getAtomPopStr(this->grossAtomPopA_, true, out);
            out << "\n";
            
            out << " Mulliken population analysis (orbial population)\n";
            this->getOrbPopStr(this->grossOrbPopA_, out);
        }
        break;
        
    case METHOD_UKS:
        {
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
            this->getAtomPopStr(this->grossAtomPopA_, false, out);
            this->getOrbPopStr(this->grossOrbPopA_, out);

            {
                const TlVector grossOrbPop = this->grossOrbPopA_ + this->grossOrbPopB_;
                out << " Gross atom population (UKS total)\n";
                this->getAtomPopStr(grossOrbPop, false, out);
                out << " Mulliken population analysis (orbial population; UKS total)\n";
                this->getOrbPopStr(grossOrbPop, out);
            }
            
            {
                const TlVector grossOrbPopSpin = this->grossOrbPopA_ - this->grossOrbPopB_;
                out << " Spin density analysis for atom (alpha - beta)\n";
                this->getAtomPopStr(grossOrbPopSpin, false, out);
                out << " Spin density analysis for orbital (alpha - beta)\n";
                this->getOrbPopStr(grossOrbPopSpin, out);
            }
        }
        break;

    default:
        // programer error
        std::abort();
        break;
    }
}


template<class SymmetricMatrixType>
void DfPopulation::calcPop(const int iteration)
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        {
            this->grossOrbPopA_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_RKS, iteration);
            this->grossAtomPopA_ = this->getGrossAtomPop(this->grossOrbPopA_);
        }
        break;

    case METHOD_ROKS:
        {
            this->grossOrbPopA_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_ROKS, iteration);
            this->grossAtomPopA_ = this->getGrossAtomPop(this->grossOrbPopA_);
        }
        break;
        
    case METHOD_UKS:
        {
            this->grossOrbPopA_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_UKS_ALPHA, iteration);
            this->grossAtomPopA_ = this->getGrossAtomPop(this->grossOrbPopA_);
            //const double elecChargeA = this->grossOrbPopA_.sum();
            
            this->grossOrbPopB_ = this->getGrossOrbPop<SymmetricMatrixType>(DfObject::RUN_UKS_BETA, iteration);
            this->grossAtomPopB_ = this->getGrossAtomPop(this->grossOrbPopB_);
            //const double elecChargeB = this->grossOrbPopB_.sum();
        }
        break;

    default:
        // programer error
        std::abort();
        break;
    }
}


template<class SymmetricMatrixType>
TlVector DfPopulation::getGrossOrbPop(const DfObject::RUN_TYPE runType, const int iteration)
{
    SymmetricMatrixType P;
    if ((runType == DfObject::RUN_RKS) ||
        (runType == DfObject::RUN_UKS_ALPHA) ||
        (runType == DfObject::RUN_UKS_BETA)) {
        P = DfObject::getPpqMatrix<SymmetricMatrixType>(runType, iteration);
    } else if (runType == DfObject::RUN_ROKS) {
        P = DfObject::getPCMatrix<SymmetricMatrixType>(iteration); // close
        const SymmetricMatrixType Popen  = DfObject::getPOMatrix<SymmetricMatrixType>(iteration); // open

        P *= 2.0;
        P += Popen;
    }

    // read AO overlap matrix
    const SymmetricMatrixType S = DfObject::getSpqMatrix<SymmetricMatrixType>();

    // calculate Mulliken Analysis ( P * S )
    const SymmetricMatrixType PS = P * S;
    
    return PS.getDiagonalElements();
}


template <typename T>
void DfPopulation::getAtomPopStr(const TlVector& trPS, const bool isCalcNetCharge,
                                 T& out) const
{
    const index_type numOfAtoms = this->m_nNumOfAtoms;
    assert(trPS.getSize() == numOfAtoms);
    
    if (isCalcNetCharge == true) {
        out << "  index Atom    gross         net\n";
    } else {
        out << "  index Atom    gross\n";
    }

    const Fl_Geometry flGeom((*(this->pPdfParam_))["coordinates"]);
    for (index_type atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const double nucCharge = flGeom.getCharge(atomIndex);
        out << TlUtils::format(" %6d   %-2s ",
                               atomIndex,
                               flGeom.getAtom(atomIndex).c_str());
        
        const double grossAtomPop = trPS.get(atomIndex);
        if (isCalcNetCharge == true) {
            out << TlUtils::format("% 12.6lf % 12.6lf\n",
                                   grossAtomPop,
                                   nucCharge - grossAtomPop);
        } else {
            out << TlUtils::format("% 12.6lf\n", grossAtomPop);
        }
    }
}


// print out Mulliken Analysis (Orbital Population)
template <typename T>
void DfPopulation::getOrbPopStr(const TlVector& trPS, T& out) const
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(trPS.getSize() == numOfAOs);
    
    out << "  index  Atom        Shell     gross\n";

    int prevAtomIndex = -1;
    for (index_type aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
        int atomIndex = this->orbitalInfo_.getAtomIndex(aoIndex);
        
        if (prevAtomIndex != atomIndex) {
            out << TlUtils::format("  %6d   %-2s %6d %-10s ",
                                   aoIndex +1,
                                   this->orbitalInfo_.getAtomName(aoIndex).c_str(),
                                   this->orbitalInfo_.getAtomIndex(aoIndex),
                                   this->orbitalInfo_.getBasisTypeName(aoIndex).c_str());
            prevAtomIndex = atomIndex;
        } else {
            out << TlUtils::format("  %6d             %-10s ",
                                   aoIndex +1,
                                   this->orbitalInfo_.getBasisTypeName(aoIndex).c_str());
        }
        out << TlUtils::format(" %12.6lf\n", trPS.get(aoIndex));
    }
}


#endif // DFPOPULATION_H

