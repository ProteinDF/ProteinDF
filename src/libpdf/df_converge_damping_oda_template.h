#ifndef DF_CONVERGE_OdaTemplate_H
#define DF_CONVERGE_OdaTemplate_H

#include <cassert>

#include "df_converge_damping_template.h"

// ----------------------------------------------------------------------------
// template class
// ----------------------------------------------------------------------------
template <class SymmetricMatrix, class Vector>
class DfConverge_Damping_OdaTemplate : public DfConverge_Damping_Template<SymmetricMatrix, Vector> {
public:
    DfConverge_Damping_OdaTemplate(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_OdaTemplate();

protected:
    virtual double getDampingFactor();
    void getParameters();

protected:
    int odaStartIteration_;

    // ODA parameters
    double oda_a_;
    double oda_b_;
    double oda_c_;
    double oda_d_;
};

///////////////////////////////////////////////////////////////////////////////
// implementation
///////////////////////////////////////////////////////////////////////////////

template <class SymmetricMatrix, class Vector>
DfConverge_Damping_OdaTemplate<SymmetricMatrix, Vector>::DfConverge_Damping_OdaTemplate(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Template<SymmetricMatrix, Vector>(pPdfParam) {
    this->odaStartIteration_ = std::max((*pPdfParam)["scf_acceleration/oda/start"].getInt(), 10);
}

template <class SymmetricMatrix, class Vector>
DfConverge_Damping_OdaTemplate<SymmetricMatrix, Vector>::~DfConverge_Damping_OdaTemplate() {}

template <class SymmetricMatrix, class Vector>
double DfConverge_Damping_OdaTemplate<SymmetricMatrix, Vector>::getDampingFactor() {
    this->log_.info(TlUtils::format("ODA start iteration: %d", this->odaStartIteration_));
    if (this->m_nIteration <= this->odaStartIteration_) {
        return DfConverge_Damping_Template<SymmetricMatrix, Vector>::getDampingFactor();
    }

    const double stdDampingFactor = (*this->pPdfParam_)["scf_acceleration/damping/damping_factor"].getDouble();
    this->log_.info(TlUtils::format("scf_acceleration/damping/damping_factor = %f", stdDampingFactor));

    this->getParameters();
    const double a = this->oda_a_;
    const double b = this->oda_b_;
    const double c = this->oda_c_;
    const double d = this->oda_d_;
    this->log_.debug(TlUtils::format("ODA parameters: a=%f, b=%f, c=%f, d=%f", a, b, c, d));

    // lambda_m: l_m
    // f(l_m) = a*l_m^3 + b*l_m^2 + c*l_m + d
    // f'(l_m) = 3.0*a*l_m^2 + 2.0*b*l_m + c
    // l_m = (-2.0*b +- sqrt((2.0*b)^2 - 4*3.0*ac)) / 2*(3.0*a)
    double l_m = 1.0;
    {
        double min_f = a + b + c + d;  // in case of "l_m = 1.0"
        std::vector<double> points;
        // points.push_back(1.0);
        // points.push_back(0.0);
        points.push_back(1.0 - stdDampingFactor);

        // 3a^2 + 2b +c
        double a_ = 3.0 * a;
        double b_ = 2.0 * b;
        const double D = b_ * b_ - 4.0 * a_ * c;
        if (D > 0.0) {
            const double root_D = sqrt(D);
            double inv_a2 = 1.0 / (2.0 * a_);

            const double l1 = (-b_ + root_D) * inv_a2;
            const double l2 = (-b_ - root_D) * inv_a2;
            this->log_.info(TlUtils::format("ODA candidates: % 8.5f, % 8.5f", l1, l2));
            if ((0.0 < l1) && (l1 < 1.0)) {
                points.push_back(l1);
            }
            if ((0.0 < l2) && (l2 < 1.0)) {
                points.push_back(l2);
            }
        }

        std::vector<double>::const_iterator pEnd = points.end();
        for (std::vector<double>::const_iterator p = points.begin(); p != pEnd; ++p) {
            const double l1 = *p;
            const double l2 = l1 * l1;
            const double l3 = l2 * l1;
            const double f = a * l3 + b * l2 + c * l1 + d;
            this->log_.info(TlUtils::format("ODA arginf info: lambda=% 8.5f, f(lambda)=% 8.5f", l1, f));
            if (f < min_f) {
                min_f = f;
                l_m = l1;
            }
        }
    }

    this->dampingFactor_ = 1.0 - l_m;  // the damping parameter of the paper and implementation are reversed
    this->log_.info(TlUtils::format("ODA damping factor = %f", this->dampingFactor_));

    return this->dampingFactor_;
}

template <class SymmetricMatrix, class Vector>
void DfConverge_Damping_OdaTemplate<SymmetricMatrix, Vector>::getParameters() {
    std::vector<double> params(4);

    const int prevIter = this->m_nIteration - 2;
    const int currIter = this->m_nIteration - 1;
    assert(currIter >= 2);

    const double E0 = this->getTotalEnergy_elec(prevIter);
    const double E1 = this->getTotalEnergy_elec(currIter);

    this->oda_d_ = E0;

    switch (this->m_nMethodType) {
        case DfObject::METHOD_RKS: {
            SymmetricMatrix dP;
            {
                const SymmetricMatrix P0 = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_RKS, prevIter);
                const SymmetricMatrix P1 = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_RKS, currIter);

                dP = P1 - P0;
            }

            {
                SymmetricMatrix F0 = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_RKS, prevIter);
                this->oda_c_ = 2.0 * F0.dotInPlace(dP).sum();
            }

            {
                SymmetricMatrix F1 = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_RKS, currIter);
                const double a1 = 2.0 * F1.dotInPlace(dP).sum();
                this->oda_a_ = a1 - 2.0 * E1 + this->oda_c_ + 2.0 * this->oda_d_;
            }
        } break;

        case DfObject::METHOD_UKS: {
            SymmetricMatrix dPa, dPb;
            {
                const SymmetricMatrix P0 = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_UKS_ALPHA, prevIter);
                const SymmetricMatrix P1 = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_UKS_ALPHA, currIter);

                dPa = P1 - P0;
            }
            {
                const SymmetricMatrix P0 = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_UKS_BETA, prevIter);
                const SymmetricMatrix P1 = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_UKS_BETA, currIter);

                dPb = P1 - P0;
            }

            {
                SymmetricMatrix Fa0 = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_UKS_ALPHA, prevIter);
                SymmetricMatrix Fb0 = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_UKS_BETA, prevIter);
                this->oda_c_ = Fa0.dotInPlace(dPa).sum() + Fb0.dotInPlace(dPb).sum();
            }

            {
                SymmetricMatrix Fa1 = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_UKS_ALPHA, currIter);
                SymmetricMatrix Fb1 = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_UKS_BETA, currIter);
                const double a1 = Fa1.dotInPlace(dPa).sum() + Fb1.dotInPlace(dPb).sum();
                this->oda_a_ = a1 - 2.0 * E1 + this->oda_c_ + 2.0 * this->oda_d_;
            }
        } break;

        case DfObject::METHOD_ROKS: {
            SymmetricMatrix dPa, dPb;
            {
                const SymmetricMatrix P0o = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_ROKS_OPEN, prevIter);
                const SymmetricMatrix P1o = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_ROKS_OPEN, currIter);
                const SymmetricMatrix P0c = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_ROKS_CLOSED, prevIter);
                const SymmetricMatrix P1c = 0.5 * DfObject::getPOutMatrix<SymmetricMatrix>(DfObject::RUN_ROKS_CLOSED, currIter);

                const SymmetricMatrix P0a = 0.5 * P0c + P0o;
                const SymmetricMatrix P0b = 0.5 * P0c;
                const SymmetricMatrix P1a = 0.5 * P1c + P1o;
                const SymmetricMatrix P1b = 0.5 * P1c;

                dPa = P1a - P0a;
                dPb = P1b - P0b;
            }

            {
                SymmetricMatrix Fa0 = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_ROKS_ALPHA, prevIter);
                SymmetricMatrix Fb0 = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_ROKS_BETA, prevIter);
                this->oda_c_ = Fa0.dotInPlace(dPa).sum() + Fb0.dotInPlace(dPb).sum();
            }

            {
                SymmetricMatrix F1a = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_UKS_ALPHA, currIter);
                SymmetricMatrix F1b = DfObject::getFpqMatrix<SymmetricMatrix>(DfObject::RUN_UKS_BETA, currIter);
                const double a1 = F1a.dotInPlace(dPa).sum() + F1b.dotInPlace(dPb).sum();
                this->oda_a_ = a1 - 2.0 * E1 + this->oda_c_ + 2.0 * this->oda_d_;
            }
        } break;

        default: {
        } break;
    }

    this->oda_b_ = E1 - this->oda_a_ - this->oda_c_ - this->oda_d_;
}

#endif  // DF_CONVERGE_OdaTemplate_H
