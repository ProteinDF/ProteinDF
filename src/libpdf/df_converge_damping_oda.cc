#include "df_converge_damping_oda.h"

#include "tl_dense_symmetric_matrix_lapack.h"

DfConverge_Damping_Oda::DfConverge_Damping_Oda(TlSerializeData* pPdfParam)
    : DfConverge_Damping(pPdfParam) {
}

DfConverge_Damping_Oda::~DfConverge_Damping_Oda() {}

double DfConverge_Damping_Oda::getDampingFactor() const {
    return this->getDampingFactor_templ<TlDenseSymmetricMatrix_Lapack>();
}

template <typename SymmetricMatrixType>
double DfConverge_Damping_Oda::getDampingFactor_templ() const {
    DfObject::RUN_TYPE runType = DfObject::RUN_RKS;  // TODO

    const int prevIter = this->m_nIteration - 2;
    const int currIter = this->m_nIteration - 1;

    if (currIter <= 1) {
        return DfConverge_Damping::getDampingFactor();
    }

    const double E0 = this->getTotalEnergy_elec(prevIter);
    const double E1 = this->getTotalEnergy_elec(currIter);

    const double d = E0;

    SymmetricMatrixType dP;
    {
        const SymmetricMatrixType P0 = 0.5 * this->getPInMatrix<SymmetricMatrixType>(runType, prevIter);
        const SymmetricMatrixType P1 = 0.5 * this->getPInMatrix<SymmetricMatrixType>(runType, currIter);

        dP = P1 - P0;
    }

    double c = 0.0;
    {
        SymmetricMatrixType F0 =
            this->getFpqMatrix<SymmetricMatrixType>(runType, prevIter);
        c = 2.0 * F0.dotInPlace(dP).sum();
    }

    double a = 0.0;
    {
        SymmetricMatrixType F1 = this->getFpqMatrix<SymmetricMatrixType>(runType, currIter);
        double a1 = 2.0 * F1.dotInPlace(dP).sum();
        a = a1 - 2.0 * E1 + c + 2.0 * d;
    }

    const double b = E1 - a - c - d;

    this->log_.info(TlUtils::format("ODA parameters: %f, %f, %f, %f", a, b, c, d));

    // lambda_m: l_m
    // f(l_m) = a*l_m^3 + b*l_m^2 + c*l_m + d
    // f'(l_m) = a*l_m^2 + b*l_m + c
    // l_m = (-b +- sqrt(b^2 - 4ac)) / 2a
    double l_m = 1.0;
    {
        double min_f = a + b + c + d;  // in case of "l_m = 1.0"
        std::vector<double> points;
        points.push_back(1.0);
        points.push_back(0.0);
        // points.push_back(1.0 - DfConverge_Damping::getDampingFactor());

        // 3a^2 + 2b +c
        double a_ = 3.0 * a;
        double b_ = 2.0 * b;
        const double D = b_ * b_ - 4.0 * a_ * c;
        if (D > 0.0) {
            const double root_D = sqrt(D);
            double inv_a2 = 1.0 / (2.0 * a_);

            const double l1 = (-b_ + root_D) * inv_a2;
            const double l2 = (-b_ - root_D) * inv_a2;
            // this->log_.info(TlUtils::format("ODA candidate: %f", l1));
            // this->log_.info(TlUtils::format("ODA candidate: %f", l2));
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
            this->log_.info(TlUtils::format("ODA arginf info: lambda=%f, f(lambda)=%f", l1, f));
            if (f < min_f) {
                min_f = f;
                l_m = l1;
            }
        }
    }

    const double answer = 1.0 - l_m;  // the damping parameter of the paper and implementation are reversed
    this->log_.info(TlUtils::format("ODA lambda = %f", answer));

    return answer;
}
