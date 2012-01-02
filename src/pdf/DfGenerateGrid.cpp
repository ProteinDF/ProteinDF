#include <cmath>

#include "Common.h"
#include "CnError.h"
#include "DfGenerateGrid.h"
#include "DfXCFunctional.h"
#include "Fl_Gto_Orbital.h"
#include "Fl_Gto_Density.h"
#include "Fl_Gto_Xcpot.h"
#include "Fl_Gto_Xcpot2.h"
#include "TlPosition.h"
#include "GridDataManager.h"
#include "TlOrbitalInfo_Density.h"
#include "TlPrdctbl.h"
#include "TlUtils.h"
#include "TlTime.h"
#include "TlMath.h"
#include "TlSymmetricMatrix.h"

#define SQ2             1.414213562373095049
#define SQ1_2           0.707106781186547524
#define SQ1_3           0.577350269189625765
#define SMALL           1
#define TOOBIG          30.0
#define BOHR            0.52917706


DfGenerateGrid::DfGenerateGrid(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), flGeometry_(Fl_Geometry::getDefaultFileName()),
      gridDataFilePath_("fl_Work/grids.dat") {
    const TlSerializeData& pdfParam = *pPdfParam;
    
    this->xctype   = pdfParam["xc-potential"].getStr();

    this->weightCutoff_ = 1.0E-16;
    if (pdfParam["grid-weight-cutoff"].getStr() != "") {
        this->weightCutoff_ = pdfParam["grid-weight-cutoff"].getDouble();
    }

    const std::string sGridType = pdfParam["xc-potential/grid-type"].getStr();
    if (TlUtils::toUpper(sGridType) == "COARSE") {
        this->m_gridType = COARSE;
        this->nrgrid =  32;
        this->nOgrid =  72;
    } else if (TlUtils::toUpper(sGridType) == "MEDIUM") {
        this->m_gridType = MEDIUM;
        this->nrgrid =  32;
        this->nOgrid = 146;
    } else if (TlUtils::toUpper(sGridType) == "MEDIUM-FINE") {
        this->m_gridType = MEDIUM_FINE;
        this->nrgrid =  64;
        this->nOgrid = 146;
    } else if (TlUtils::toUpper(sGridType) == "FINE") {
        this->m_gridType = FINE;
        this->nrgrid =  64;
        this->nOgrid = 302;
    } else if (TlUtils::toUpper(sGridType) == "SG-1") {
        this->m_gridType = SG_1;
        this->nrgrid =  50;
        this->nOgrid = 196;   // 適当
    } else {
        this->logger("Selection of Grid Type is Wrong.\n");
        this->logger("You type in [" + sGridType + "].\n");
        CnErr.abort();
    }

    const int dNumOfAtoms = this->numOfRealAtoms_;
    this->coord_.resize(dNumOfAtoms);
    for (int i = 0; i < dNumOfAtoms; ++i) {
        this->coord_[i] = this->flGeometry_.getCoordinate(i);
    }

    this->distanceMatrix_.resize(dNumOfAtoms);
    for (int mc = 0; mc < dNumOfAtoms; ++mc) {
        const TlPosition pos_mc = this->coord_[mc];
        for (int nc = 0; nc < mc; ++nc) {
            const double dist = pos_mc.distanceFrom(this->coord_[nc]);
            if (dist < 1E-10) {
                std::cerr << " distance < 1E-10: " << mc << " th atom and " << nc << " th atom." << std::endl;
            }
            this->distanceMatrix_(mc, nc) = dist;
        }
    }

    this->numOfColsOfGrdMat_ = 5; // x, y, z, weight, atom_index
    {
        const int coef = (this->m_nMethodType == METHOD_RKS) ? 1 : 2;
        DfXCFunctional dfXcFunctional(this->pPdfParam_);
        if (dfXcFunctional.getFunctionalType() == DfXCFunctional::LDA) {
            this->numOfColsOfGrdMat_ += coef * 1; // rho only
        } else if (dfXcFunctional.getFunctionalType() == DfXCFunctional::GGA) {
            this->numOfColsOfGrdMat_ += coef * 4; // rho, gradRhoX, gradRhoY, gradRhoZ
        }
    }
}

DfGenerateGrid::~DfGenerateGrid()
{
}

int DfGenerateGrid::dfGrdMain()
{
    this->logger("start        : " + TlTime::getNow() + "\n");

    this->makeTable();
    this->logger("makeTable    : " + TlTime::getNow() + "\n");

    this->setCellPara();
    this->logger("setCellPara  : " + TlTime::getNow() + "\n");

    this->generateGrid();
    this->logger("generateGrid : " + TlTime::getNow() + "\n");

    this->logger("end          : " + TlTime::getNow() + "\n");

    return 0;
}

void DfGenerateGrid::makeTable()
{
    TlOrbitalInfo_Density orbInfoAuxCD((*this->pPdfParam_)["coordinates"],
                                       (*this->pPdfParam_)["basis_sets_j"]);
    const int maxNumOfAuxCDs = orbInfoAuxCD.getNumOfOrbitals();
    
    double dMaxExpAlpha = 0.0;
    for (int i = 0; i < maxNumOfAuxCDs; ++i) {
        dMaxExpAlpha = std::max(dMaxExpAlpha, orbInfoAuxCD.getExponent(i, 0));
    }
    
    const double r = std::max((TOOBIG / dMaxExpAlpha), 200.0);
    this->maxRadii_ = std::sqrt(r);
}

void DfGenerateGrid::setCellPara()
{
    // Set Atomic Radii
    // Array number is equal to the atomic number, then this begins at 1.
    // Note that Atomic Radius of H is twiced.
    const int maxKindOfAtoms = 110;
    this->radiusList_.resize(maxKindOfAtoms); //
    this->xGL_.resize(100);
    this->wGL_.resize(100);

    if ((this->m_gridType == COARSE) ||
        (this->m_gridType == MEDIUM) ||
        (this->m_gridType == MEDIUM_FINE) ||
        (this->m_gridType == FINE)) {

        this->radiusList_[  1] = 0.50;
        this->radiusList_[  2] = 2.00;
        this->radiusList_[  3] = 1.45;
        this->radiusList_[  4] = 1.05;
        this->radiusList_[  5] = 0.85;
        this->radiusList_[  6] = 0.70;
        this->radiusList_[  7] = 0.65;
        this->radiusList_[  8] = 0.60;
        this->radiusList_[  9] = 0.50;
        this->radiusList_[ 10] = 2.25;

        this->radiusList_[ 11] = 1.80;
        this->radiusList_[ 12] = 1.50;
        this->radiusList_[ 13] = 1.25;
        this->radiusList_[ 14] = 1.10;
        this->radiusList_[ 15] = 1.00;
        this->radiusList_[ 16] = 1.00;
        this->radiusList_[ 17] = 1.00;
        this->radiusList_[ 18] = 2.50;
        this->radiusList_[ 19] = 2.20;
        this->radiusList_[ 20] = 1.80;

        this->radiusList_[ 21] = 1.60;
        this->radiusList_[ 22] = 1.40;
        this->radiusList_[ 23] = 1.35;
        this->radiusList_[ 24] = 1.40;
        this->radiusList_[ 25] = 1.40;
        this->radiusList_[ 26] = 1.40;
        this->radiusList_[ 27] = 1.35;
        this->radiusList_[ 28] = 1.35;
        this->radiusList_[ 29] = 1.35;
        this->radiusList_[ 30] = 1.35;

        this->radiusList_[ 31] = 1.30;
        this->radiusList_[ 32] = 1.25;
        this->radiusList_[ 33] = 1.15;
        this->radiusList_[ 34] = 1.15;
        this->radiusList_[ 35] = 1.15;
        this->radiusList_[ 36] = 2.75;
        this->radiusList_[ 37] = 2.35;
        this->radiusList_[ 38] = 2.00;
        this->radiusList_[ 39] = 1.80;
        this->radiusList_[ 40] = 1.55;

        this->radiusList_[ 41] = 1.45;
        this->radiusList_[ 42] = 1.45;
        this->radiusList_[ 43] = 1.35;
        this->radiusList_[ 44] = 1.30;
        this->radiusList_[ 45] = 1.35;
        this->radiusList_[ 46] = 1.40;
        this->radiusList_[ 47] = 1.60;
        this->radiusList_[ 48] = 1.55;
        this->radiusList_[ 49] = 1.55;
        this->radiusList_[ 50] = 1.45;

        this->radiusList_[ 51] = 1.45;
        this->radiusList_[ 52] = 1.40;
        this->radiusList_[ 53] = 1.40;
        this->radiusList_[ 54] = 3.00;
        this->radiusList_[ 55] = 2.60;
        this->radiusList_[ 56] = 2.15;
        this->radiusList_[ 57] = 1.95;
        this->radiusList_[ 58] = 1.85;
        this->radiusList_[ 59] = 1.85;
        this->radiusList_[ 60] = 1.85;

        this->radiusList_[ 61] = 1.85;
        this->radiusList_[ 62] = 1.85;
        this->radiusList_[ 63] = 1.85;
        this->radiusList_[ 64] = 1.80;
        this->radiusList_[ 65] = 1.75;
        this->radiusList_[ 66] = 1.75;
        this->radiusList_[ 67] = 1.75;
        this->radiusList_[ 68] = 1.75;
        this->radiusList_[ 69] = 1.75;
        this->radiusList_[ 70] = 1.75;

        this->radiusList_[ 71] = 1.75;
        this->radiusList_[ 72] = 1.55;
        this->radiusList_[ 73] = 1.45;
        this->radiusList_[ 74] = 1.35;
        this->radiusList_[ 75] = 1.35;
        this->radiusList_[ 76] = 1.30;
        this->radiusList_[ 77] = 1.35;
        this->radiusList_[ 78] = 1.35;
        this->radiusList_[ 79] = 1.35;
        this->radiusList_[ 80] = 1.50;

        this->radiusList_[ 81] = 1.90;
        this->radiusList_[ 82] = 1.80;
        this->radiusList_[ 83] = 1.60;
        this->radiusList_[ 84] = 1.90;
        this->radiusList_[ 85] = 1.65;
        this->radiusList_[ 86] = 3.25;
        this->radiusList_[ 87] = 2.80;
        this->radiusList_[ 88] = 2.15;
        this->radiusList_[ 89] = 1.95;
        this->radiusList_[ 90] = 1.80;

        this->radiusList_[ 91] = 1.80;
        this->radiusList_[ 92] = 1.75;
        this->radiusList_[ 93] = 1.75;
        this->radiusList_[ 94] = 1.75;
        this->radiusList_[ 95] = 1.75;
        this->radiusList_[ 96] = 1.75;
        this->radiusList_[ 97] = 1.75;
        this->radiusList_[ 98] = 1.75;
        this->radiusList_[ 99] = 1.75;
        this->radiusList_[100] = 1.75;

        this->radiusList_[101] = 1.75;
        this->radiusList_[102] = 1.75;
        this->radiusList_[103] = 1.75;
        this->radiusList_[104] = 1.55;
        this->radiusList_[105] = 1.55;

    } else if (this->m_gridType == SG_1) {
        // for SG-1(Slater Atomic Radii)
        this->radiusList_[  1] = 1.0000;   //0.50;
        this->radiusList_[  2] = 0.5882;   //2.00;
        this->radiusList_[  3] = 3.0769;   //1.45;
        this->radiusList_[  4] = 2.0513;   //1.05;
        this->radiusList_[  5] = 1.5385;   //0.85;
        this->radiusList_[  6] = 1.2308;   //0.70;
        this->radiusList_[  7] = 1.0256;   //0.65;
        this->radiusList_[  8] = 0.8791;   //0.60;
        this->radiusList_[  9] = 0.7692;   //0.50;
        this->radiusList_[ 10] = 0.6838;   //2.25;

        this->radiusList_[ 11] = 4.0909;   //1.80;
        this->radiusList_[ 12] = 3.1579;   //1.50;
        this->radiusList_[ 13] = 2.5714;   //1.25;
        this->radiusList_[ 14] = 2.1687;   //1.10;
        this->radiusList_[ 15] = 1.8750;   //1.00;
        this->radiusList_[ 16] = 1.6514;   //1.00;
        this->radiusList_[ 17] = 1.4754;   //1.00;
        this->radiusList_[ 18] = 1.3333;   //2.50;
        this->radiusList_[ 19] = 2.20;     //2.20;
        this->radiusList_[ 20] = 2.03;     //1.80;

        this->radiusList_[ 21] = 3.023563;  //1.60;
        this->radiusList_[ 22] = 2.645618;  //1.40;
        this->radiusList_[ 23] = 2.551131;  //1.35;
        this->radiusList_[ 24] = 2.645618;
        this->radiusList_[ 25] = 2.645618;
        this->radiusList_[ 26] = 2.645618;
        this->radiusList_[ 27] = 2.551131;
        this->radiusList_[ 28] = 2.551131;
        this->radiusList_[ 29] = 2.551131;
        this->radiusList_[ 30] = 2.551131;

        this->radiusList_[ 31] = 2.456645;  //1.30;
        this->radiusList_[ 32] = 2.362159;  //1.25;
        this->radiusList_[ 33] = 2.173186;  //1.15;
        this->radiusList_[ 34] = 2.173186;
        this->radiusList_[ 35] = 2.173186;
        this->radiusList_[ 36] = 5.196749;  //2.75;??
        this->radiusList_[ 37] = 4.440858;  //2.35;
        this->radiusList_[ 38] = 3.779454;  //2.00;
        this->radiusList_[ 39] = 3.401508;  //1.80;
        this->radiusList_[ 40] = 2.551131;  //1.55;

        this->radiusList_[ 41] = 2.702309;  //1.45;
        this->radiusList_[ 42] = 2.740104;  //1.45;
        this->radiusList_[ 43] = 2.551131;  //1.35;
        this->radiusList_[ 44] = 2.456645;  //1.30;
        this->radiusList_[ 45] = 2.551131;  //1.35;
        this->radiusList_[ 46] = 2.645618;  //1.40;
        this->radiusList_[ 47] = 3.023563;  //1.60;
        this->radiusList_[ 48] = 2.929077;  //1.55;
        this->radiusList_[ 49] = 2.929077;  //1.55;
        this->radiusList_[ 50] = 2.740104;  //1.45;

        this->radiusList_[ 51] = 2.740104;  //1.45;
        this->radiusList_[ 52] = 2.645618;  //1.40;
        this->radiusList_[ 53] = 2.645618;  //1.40;
        this->radiusList_[ 54] = 5.669181;  //3.00;
        this->radiusList_[ 55] = 4.913290;  //2.60;
        this->radiusList_[ 56] = 4.062913;  //2.15;
        this->radiusList_[ 57] = 3.684967;  //1.95;
        this->radiusList_[ 58] = 3.495995;  //1.85;
        this->radiusList_[ 59] = 3.495995;  //1.85;
        this->radiusList_[ 60] = 3.495995;  //1.85;

        this->radiusList_[ 61] = 3.458200;  //1.83;
        this->radiusList_[ 62] = 3.495995;  //1.85;
        this->radiusList_[ 63] = 3.495995;  //1.85;
        this->radiusList_[ 64] = 3.401508;   //1.80;
        this->radiusList_[ 65] = 3.307022;   //1.75;
        this->radiusList_[ 66] = 3.307022;   //1.75;
        this->radiusList_[ 67] = 3.307022;   //1.75;
        this->radiusList_[ 68] = 3.307022;   //1.75;
        this->radiusList_[ 69] = 3.307022;   //1.75;
        this->radiusList_[ 70] = 3.307022;   //1.75;

        this->radiusList_[ 71] = 3.307022;   //1.75;
        this->radiusList_[ 72] = 2.929077;   //1.55;
        this->radiusList_[ 73] = 2.740104;   //1.45;
        this->radiusList_[ 74] = 2.551131;   //1.35;
        this->radiusList_[ 75] = 2.551131;   //1.35;
        this->radiusList_[ 76] = 2.456645;   //1.30;
        this->radiusList_[ 77] = 2.551131;   //1.35;
        this->radiusList_[ 78] = 2.551131;   //1.35;
        this->radiusList_[ 79] = 2.551131;   //1.35;
        this->radiusList_[ 80] = 2.834590;   //1.50;

        this->radiusList_[ 81] = 3.590481;   //1.90;
        this->radiusList_[ 82] = 3.401508;   //1.80;
        this->radiusList_[ 83] = 3.023563;   //1.60;
        this->radiusList_[ 84] = 3.590481;   //1.90;
        this->radiusList_[ 85] = 3.118049;   //1.65;
        this->radiusList_[ 86] = 6.141612;   //3.25;
        this->radiusList_[ 87] = 5.291235;   //2.80;
        this->radiusList_[ 88] = 4.062913;   //2.15;
        this->radiusList_[ 89] = 3.684967;   //1.95;
        this->radiusList_[ 90] = 3.401508;   //1.80;

        this->radiusList_[ 91] = 3.401508;   //1.80;
        this->radiusList_[ 92] = 3.307022;   //1.75;
        this->radiusList_[ 93] = 3.307022;   //1.75;
        this->radiusList_[ 94] = 3.307022;   //1.75;
        this->radiusList_[ 95] = 3.307022;   //1.75;
        this->radiusList_[ 96] = 3.307022;   //1.75;
        this->radiusList_[ 97] = 3.307022;   //1.75;
        this->radiusList_[ 98] = 3.307022;   //1.75;
        this->radiusList_[ 99] = 3.307022;   //1.75;
        this->radiusList_[100] = 3.307022;   //1.75;

        this->radiusList_[101] = 3.307022;   //1.75;
        this->radiusList_[102] = 3.307022;   //1.75;
        this->radiusList_[103] = 3.307022;   //1.75;
        this->radiusList_[104] = 2.929077;   //1.55;
        this->radiusList_[105] = 2.929077;   //1.55;

    }

    // Set Gauss-Legendre abscissas (xGL) and weights (wGL)
    if ((this->m_gridType == COARSE) ||
            (this->m_gridType == MEDIUM)) {

        this->xGL_[ 0] = -9.9726386184948157e-01;
        this->xGL_[ 1] = -9.8561151154526838e-01;
        this->xGL_[ 2] = -9.6476225558750639e-01;
        this->xGL_[ 3] = -9.3490607593773967e-01;
        this->xGL_[ 4] = -8.9632115576605209e-01;
        this->xGL_[ 5] = -8.4936761373256997e-01;
        this->xGL_[ 6] = -7.9448379596794239e-01;
        this->xGL_[ 7] = -7.3218211874028971e-01;
        this->xGL_[ 8] = -6.6304426693021523e-01;
        this->xGL_[ 9] = -5.8771575724076230e-01;
        this->xGL_[10] = -5.0689990893222936e-01;
        this->xGL_[11] = -4.2135127613063533e-01;
        this->xGL_[12] = -3.3186860228212767e-01;
        this->xGL_[13] = -2.3928736225213706e-01;
        this->xGL_[14] = -1.4447196158279649e-01;
        this->xGL_[15] = -4.8307665687738310e-02;
        this->xGL_[16] =  4.8307665687738310e-02;
        this->xGL_[17] =  1.4447196158279649e-01;
        this->xGL_[18] =  2.3928736225213706e-01;
        this->xGL_[19] =  3.3186860228212767e-01;
        this->xGL_[20] =  4.2135127613063533e-01;
        this->xGL_[21] =  5.0689990893222936e-01;
        this->xGL_[22] =  5.8771575724076230e-01;
        this->xGL_[23] =  6.6304426693021523e-01;
        this->xGL_[24] =  7.3218211874028971e-01;
        this->xGL_[25] =  7.9448379596794239e-01;
        this->xGL_[26] =  8.4936761373256997e-01;
        this->xGL_[27] =  8.9632115576605209e-01;
        this->xGL_[28] =  9.3490607593773967e-01;
        this->xGL_[29] =  9.6476225558750639e-01;
        this->xGL_[30] =  9.8561151154526838e-01;
        this->xGL_[31] =  9.9726386184948157e-01;

        this->wGL_[ 0] =  7.0186100094695213e-03;
        this->wGL_[ 1] =  1.6274394730905709e-02;
        this->wGL_[ 2] =  2.5392065309262139e-02;
        this->wGL_[ 3] =  3.4273862913021411e-02;
        this->wGL_[ 4] =  4.2835898022226704e-02;
        this->wGL_[ 5] =  5.0998059262376154e-02;
        this->wGL_[ 6] =  5.8684093478535128e-02;
        this->wGL_[ 7] =  6.5822222776361808e-02;
        this->wGL_[ 8] =  7.2345794108848616e-02;
        this->wGL_[ 9] =  7.8193895787070436e-02;
        this->wGL_[10] =  8.3311924226946721e-02;
        this->wGL_[11] =  8.7652093004403742e-02;
        this->wGL_[12] =  9.1173878695763905e-02;
        this->wGL_[13] =  9.3844399080804414e-02;
        this->wGL_[14] =  9.5638720079274847e-02;
        this->wGL_[15] =  9.6540088514727854e-02;
        this->wGL_[16] =  9.6540088514727854e-02;
        this->wGL_[17] =  9.5638720079274847e-02;
        this->wGL_[18] =  9.3844399080804414e-02;
        this->wGL_[19] =  9.1173878695763905e-02;
        this->wGL_[20] =  8.7652093004403742e-02;
        this->wGL_[21] =  8.3311924226946721e-02;
        this->wGL_[22] =  7.8193895787070436e-02;
        this->wGL_[23] =  7.2345794108848616e-02;
        this->wGL_[24] =  6.5822222776361808e-02;
        this->wGL_[25] =  5.8684093478535128e-02;
        this->wGL_[26] =  5.0998059262376154e-02;
        this->wGL_[27] =  4.2835898022226704e-02;
        this->wGL_[28] =  3.4273862913021411e-02;
        this->wGL_[29] =  2.5392065309262139e-02;
        this->wGL_[30] =  1.6274394730905709e-02;
        this->wGL_[31] =  7.0186100094695213e-03;

    } else if ((this->m_gridType == MEDIUM_FINE) ||
               (this->m_gridType == FINE) ||
               (this->m_gridType == SG_1)) {

        this->xGL_[ 0] = -9.9930504173577217e-01;
        this->xGL_[ 1] = -9.9634011677195522e-01;
        this->xGL_[ 2] = -9.9101337147674429e-01;
        this->xGL_[ 3] = -9.8333625388462598e-01;
        this->xGL_[ 4] = -9.7332682778991098e-01;
        this->xGL_[ 5] = -9.6100879965205377e-01;
        this->xGL_[ 6] = -9.4641137485840277e-01;
        this->xGL_[ 7] = -9.2956917213193957e-01;
        this->xGL_[ 8] = -9.1052213707850282e-01;
        this->xGL_[ 9] = -8.8931544599511414e-01;
        this->xGL_[10] = -8.6599939815409288e-01;
        this->xGL_[11] = -8.4062929625258043e-01;
        this->xGL_[12] = -8.1326531512279754e-01;
        this->xGL_[13] = -7.8397235894334139e-01;
        this->xGL_[14] = -7.5281990726053194e-01;
        this->xGL_[15] = -7.1988185017161088e-01;
        this->xGL_[16] = -6.8523631305423327e-01;
        this->xGL_[17] = -6.4896547125465731e-01;
        this->xGL_[18] = -6.1115535517239328e-01;
        this->xGL_[19] = -5.7189564620263400e-01;
        this->xGL_[20] = -5.3127946401989457e-01;
        this->xGL_[21] = -4.8940314570705296e-01;
        this->xGL_[22] = -4.4636601725346409e-01;
        this->xGL_[23] = -4.0227015796399163e-01;
        this->xGL_[24] = -3.5722015833766813e-01;
        this->xGL_[25] = -3.1132287199021097e-01;
        this->xGL_[26] = -2.6468716220876742e-01;
        this->xGL_[27] = -2.1742364374000708e-01;
        this->xGL_[28] = -1.6964442042399280e-01;
        this->xGL_[29] = -1.2146281929612056e-01;
        this->xGL_[30] = -7.2993121787799042e-02;
        this->xGL_[31] = -2.4350292663424429e-02;
        this->xGL_[32] =  2.4350292663424429e-02;
        this->xGL_[33] =  7.2993121787799042e-02;
        this->xGL_[34] =  1.2146281929612056e-01;
        this->xGL_[35] =  1.6964442042399280e-01;
        this->xGL_[36] =  2.1742364374000708e-01;
        this->xGL_[37] =  2.6468716220876742e-01;
        this->xGL_[38] =  3.1132287199021097e-01;
        this->xGL_[39] =  3.5722015833766813e-01;
        this->xGL_[40] =  4.0227015796399163e-01;
        this->xGL_[41] =  4.4636601725346409e-01;
        this->xGL_[42] =  4.8940314570705296e-01;
        this->xGL_[43] =  5.3127946401989457e-01;
        this->xGL_[44] =  5.7189564620263400e-01;
        this->xGL_[45] =  6.1115535517239328e-01;
        this->xGL_[46] =  6.4896547125465731e-01;
        this->xGL_[47] =  6.8523631305423327e-01;
        this->xGL_[48] =  7.1988185017161088e-01;
        this->xGL_[49] =  7.5281990726053194e-01;
        this->xGL_[50] =  7.8397235894334139e-01;
        this->xGL_[51] =  8.1326531512279754e-01;
        this->xGL_[52] =  8.4062929625258043e-01;
        this->xGL_[53] =  8.6599939815409288e-01;
        this->xGL_[54] =  8.8931544599511414e-01;
        this->xGL_[55] =  9.1052213707850282e-01;
        this->xGL_[56] =  9.2956917213193957e-01;
        this->xGL_[57] =  9.4641137485840277e-01;
        this->xGL_[58] =  9.6100879965205377e-01;
        this->xGL_[59] =  9.7332682778991098e-01;
        this->xGL_[60] =  9.8333625388462598e-01;
        this->xGL_[61] =  9.9101337147674429e-01;
        this->xGL_[62] =  9.9634011677195522e-01;
        this->xGL_[63] =  9.9930504173577217e-01;

        this->wGL_[ 0] =  1.7832807216962678e-03;
        this->wGL_[ 1] =  4.1470332605625208e-03;
        this->wGL_[ 2] =  6.5044579689783680e-03;
        this->wGL_[ 3] =  8.8467598263639348e-03;
        this->wGL_[ 4] =  1.1168139460131076e-02;
        this->wGL_[ 5] =  1.3463047896718736e-02;
        this->wGL_[ 6] =  1.5726030476024472e-02;
        this->wGL_[ 7] =  1.7951715775697225e-02;
        this->wGL_[ 8] =  2.0134823153530181e-02;
        this->wGL_[ 9] =  2.2270173808383212e-02;
        this->wGL_[10] =  2.4352702568710933e-02;
        this->wGL_[11] =  2.6377469715054645e-02;
        this->wGL_[12] =  2.8339672614259459e-02;
        this->wGL_[13] =  3.0234657072402426e-02;
        this->wGL_[14] =  3.2057928354851606e-02;
        this->wGL_[15] =  3.3805161837141613e-02;
        this->wGL_[16] =  3.5472213256882358e-02;
        this->wGL_[17] =  3.7055128540240068e-02;
        this->wGL_[18] =  3.8550153178615598e-02;
        this->wGL_[19] =  3.9953741132720363e-02;
        this->wGL_[20] =  4.1262563242623548e-02;
        this->wGL_[21] =  4.2473515123653625e-02;
        this->wGL_[22] =  4.3583724529323471e-02;
        this->wGL_[23] =  4.4590558163756580e-02;
        this->wGL_[24] =  4.5491627927418121e-02;
        this->wGL_[25] =  4.6284796581314396e-02;
        this->wGL_[26] =  4.6968182816209993e-02;
        this->wGL_[27] =  4.7540165714830315e-02;
        this->wGL_[28] =  4.7999388596458331e-02;
        this->wGL_[29] =  4.8344762234802899e-02;
        this->wGL_[30] =  4.8575467441503394e-02;
        this->wGL_[31] =  4.8690957009139731e-02;
        this->wGL_[32] =  4.8690957009139731e-02;
        this->wGL_[33] =  4.8575467441503394e-02;
        this->wGL_[34] =  4.8344762234802899e-02;
        this->wGL_[35] =  4.7999388596458331e-02;
        this->wGL_[36] =  4.7540165714830315e-02;
        this->wGL_[37] =  4.6968182816209993e-02;
        this->wGL_[38] =  4.6284796581314396e-02;
        this->wGL_[39] =  4.5491627927418121e-02;
        this->wGL_[40] =  4.4590558163756580e-02;
        this->wGL_[41] =  4.3583724529323471e-02;
        this->wGL_[42] =  4.2473515123653625e-02;
        this->wGL_[43] =  4.1262563242623548e-02;
        this->wGL_[44] =  3.9953741132720363e-02;
        this->wGL_[45] =  3.8550153178615598e-02;
        this->wGL_[46] =  3.7055128540240068e-02;
        this->wGL_[47] =  3.5472213256882358e-02;
        this->wGL_[48] =  3.3805161837141613e-02;
        this->wGL_[49] =  3.2057928354851606e-02;
        this->wGL_[50] =  3.0234657072402426e-02;
        this->wGL_[51] =  2.8339672614259459e-02;
        this->wGL_[52] =  2.6377469715054645e-02;
        this->wGL_[53] =  2.4352702568710933e-02;
        this->wGL_[54] =  2.2270173808383212e-02;
        this->wGL_[55] =  2.0134823153530181e-02;
        this->wGL_[56] =  1.7951715775697225e-02;
        this->wGL_[57] =  1.5726030476024472e-02;
        this->wGL_[58] =  1.3463047896718736e-02;
        this->wGL_[59] =  1.1168139460131076e-02;
        this->wGL_[60] =  8.8467598263639348e-03;
        this->wGL_[61] =  6.5044579689783680e-03;
        this->wGL_[62] =  4.1470332605625208e-03;
        this->wGL_[63] =  1.7832807216962678e-03;
    }
}


void DfGenerateGrid::generateGrid()
{
    std::size_t numOfGrids = 0;
    const int endAtom = this->numOfRealAtoms_;
    for (int atom = 0; atom < endAtom; ++atom) {
        std::vector<double> coordX;
        std::vector<double> coordY;
        std::vector<double> coordZ;
        std::vector<double> weight;

        if (this->m_gridType == SG_1) {
            this->generateGrid_SG1(atom, &coordX, &coordY, &coordZ, &weight);
        } else {
            this->generateGrid(atom, &coordX, &coordY, &coordZ, &weight);
        }

        // store grid matrix
        const std::size_t numOfAtomGrids = weight.size();
        this->grdMat_.resize(numOfGrids + numOfAtomGrids, this->numOfColsOfGrdMat_);
        for (std::size_t i = 0; i < numOfAtomGrids; ++i) {
            this->grdMat_.set(numOfGrids, 0, coordX[i]);
            this->grdMat_.set(numOfGrids, 1, coordY[i]);
            this->grdMat_.set(numOfGrids, 2, coordZ[i]);
            this->grdMat_.set(numOfGrids, 3, weight[i]);
            this->grdMat_.set(numOfGrids, 4, atom);
            ++numOfGrids;
        }
    }

    this->saveGridMatrix(0, this->grdMat_);
}


void DfGenerateGrid::generateGrid(const int iatom,
                                  std::vector<double>* pCoordX,
                                  std::vector<double>* pCoordY,
                                  std::vector<double>* pCoordZ,
                                  std::vector<double>* pWeight)
{
    assert(pCoordX != NULL);
    assert(pCoordY != NULL);
    assert(pCoordZ != NULL);
    assert(pWeight != NULL);
    assert((this->m_gridType == COARSE) || (this->m_gridType == MEDIUM) ||
           (this->m_gridType == MEDIUM_FINE) || (this->m_gridType == FINE));

    // 2. generate radial grid
    // 3. generate Omega  grid
    // 4. calculate weight by Fuzzy cell method

    std::vector<TlPosition> crdpoint(20000);
    std::vector<double> weightvec(20000);
    std::vector<double> Ps(this->numOfRealAtoms_);
    const int nNumOfAtoms = this->numOfRealAtoms_; //this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;
    int GPthrnum = 0;
    const double rM = 0.5 * this->radiusList_[TlPrdctbl::getAtomicNumber(this->flGeometry_.getAtom(iatom))] / BOHR;

    for (int radvec = 0; radvec < this->nrgrid; ++radvec) {
        const double r0  = rM * (1.0 + this->xGL_[radvec]) / (1.0 - this->xGL_[radvec]);
        const double dr0 = 2.0 * rM / ((1.0 - this->xGL_[radvec]) * (1.0 - this->xGL_[radvec]));
        const double weightpoint = 1.0 / (double) nOgrid;
        double weight = weightpoint * r0 * r0 * dr0 * this->wGL_[radvec] * 4.0 * M_PI;

        std::vector<TlPosition> grid(nOgrid);
        std::vector<double> lebWeight(nOgrid);
        this->points2(nOgrid, r0, this->coord_[iatom], weight, grid, lebWeight);

        // Loop for the grid number of Omega vector for normal Grid(Beck\'s Method)
        for (int Omega = 0; Omega < nOgrid; ++Omega) {
            weight = lebWeight[Omega];
            const TlPosition pos_O = grid[Omega];

            // グリッド省略判定
            int checkp = 1;
            for (int p = 0; p < this->numOfRealAtoms_; ++p) {
                const double rp = pos_O.distanceFrom(this->coord_[p]);

                if (rp < this->maxRadii_) {
                    checkp = 0;
                    break;
                }
            }

            if (checkp == 1) {
                goto JUMP1;
            }

            // Loop for the m-center
            for (int mc = 0; mc < nNumOfAtoms; ++mc) {
                double Psuij = 1.0;

                const int atomnum1 = TlPrdctbl::getAtomicNumber(this->flGeometry_.getAtom(mc));
                double radius1  = this->radiusList_[atomnum1];
                if (atomnum1 == 1) {
                    radius1 = 0.35; // for H
                }

                const double ri = pos_O.distanceFrom(this->coord_[mc]);

                // Loop for the n-center
                for (int nc = 0; nc < nNumOfAtoms; ++nc) {
                    const double rj = pos_O.distanceFrom(this->coord_[nc]);

                    if (nc != mc) {
                        const double Rij = this->coord_[mc].distanceFrom(this->coord_[nc]);

                        const int atomnum2 = TlPrdctbl::getAtomicNumber(this->flGeometry_.getAtom(nc));
                        double radius2  = this->radiusList_[atomnum2];
                        if (atomnum2 == 1) {
                            radius2 = 0.35; // for H
                        }
                        const double Kai = radius1 / radius2;
                        const double uijtmp = (Kai -1.0) / (Kai +1.0);
                        double aij = uijtmp / (uijtmp * uijtmp - 1.0);
                        if (aij >= 0.5) {
                            aij =  0.5;
                        }
                        if (aij <= -0.5) {
                            aij = -0.5;
                        }

                        // Calculate the Fuzzy Cell
                        double uij = (ri - rj) / Rij;
                        uij += aij * (1.0 - uij * uij);

                        double suij = 0.0;
                        if (uij > 0.9) {
                            Ps[mc] = 0.0;
                            goto JUMPMC1;
                        } else if (uij < -0.9) {
                            suij = 1.0;
                        } else {
                            const double f1u    = 1.5 * uij - 0.5 * uij * uij * uij;
                            const double f2u    = 1.5 * f1u - 0.5 * f1u * f1u * f1u;
                            const double f3u    = 1.5 * f2u - 0.5 * f2u * f2u * f2u;
                            suij = 0.5 * (1.0 - f3u);
                        }
                        // Pai suij
                        Psuij *= suij;
                    }
                }

                Ps[mc] = Psuij;
JUMPMC1:
                continue;
            }

            // Normalization
            {
                double Pstotal = 0.0;
                for (int mc = 0; mc < this->numOfRealAtoms_; ++mc) {
                    Pstotal += Ps[mc];
                }

                weight *= Ps[iatom] / Pstotal;
            }

            weightvec[GPthrnum] = weight;
            crdpoint[GPthrnum] = grid[Omega];
            ++GPthrnum;

JUMP1:
            continue; // グリッド省略でここにJUMPする
        }
    }

    pCoordX->resize(GPthrnum);
    pCoordY->resize(GPthrnum);
    pCoordZ->resize(GPthrnum);
    for (int i = 0; i < GPthrnum; ++i) {
        (*pCoordX)[i] = crdpoint[i].x();
        (*pCoordY)[i] = crdpoint[i].y();
        (*pCoordZ)[i] = crdpoint[i].z();
    }
    weightvec.resize(GPthrnum);
    *pWeight = weightvec;
}


void DfGenerateGrid::generateGrid_SG1(const int iAtom,
                                      std::vector<double>* pCoordX,
                                      std::vector<double>* pCoordY,
                                      std::vector<double>* pCoordZ,
                                      std::vector<double>* pWeight)
{
    assert(pCoordX != NULL);
    assert(pCoordY != NULL);
    assert(pCoordZ != NULL);
    assert(pWeight != NULL);

    std::vector<TlPosition> crdpoint(20000);
    std::vector<double> weightvec(20000);
    static const size_t maxGridSize = 20000;
    pCoordX->resize(maxGridSize);
    pCoordY->resize(maxGridSize);
    pCoordZ->resize(maxGridSize);
    pWeight->resize(maxGridSize);
    std::vector<double> Ps(this->numOfRealAtoms_);
    int GPthrnum = 0;

    const int atomnum = TlPrdctbl::getAtomicNumber(this->flGeometry_.getAtom(iAtom));
    const double rM = this->radiusList_[atomnum];
    const double inv_rM = 1.0 / rM;

    // Set the partitioning parameters alpha used in the SG-1 grid
    // set default value; for after atom #18
    double alpha0 = 0.0;
    double alpha1 = 0.0;
    double alpha2 = 0.0;
    double alpha3 = 10000.0;
    if ((atomnum == 1) || (atomnum == 2)) {
        alpha0 = 0.2500;
        alpha1 = 0.5000;
        alpha2 = 1.0000;
        alpha3 = 4.5000;
    } else if ((3 <= atomnum) && (atomnum <= 10)) {
        alpha0 = 0.1667;
        alpha1 = 0.5000;
        alpha2 = 0.9000;
        alpha3 = 3.5000;
    } else if ((11 <= atomnum) && (atomnum <= 18)) {
        alpha0 = 0.1000;
        alpha1 = 0.4000;
        alpha2 = 0.8000;
        alpha3 = 2.5000;
    }

    // Loop for the grid number of radial vector
    const int radvec_max = this->nrgrid;
    for (int radvec = 0; radvec < radvec_max; ++radvec) {
        const double r0 = rM * TlMath::pow(radvec + 1.0, 2) * TlMath::pow(radvec_max - radvec, -2);
        const double weight = 2.0 * rM * rM * rM * (radvec_max + 1.0)
                              * TlMath::pow(radvec + 1.0, 5)
                              * TlMath::pow(radvec_max - radvec, -7)
                              * 4.0 * M_PI;

        // Set the partitioning area
        int Ogrid = 86;
        const double judge = r0 * inv_rM;
        if (judge <= alpha0) {
            // (0, alpha0)
            Ogrid = 6;
        } else if (judge <= alpha1) {
            // (alpha0, alpha1]
            Ogrid = 38;
        } else if (judge <= alpha2) {
            // (alpha1, alpha2]
            Ogrid = 86;
        } else if (judge <= alpha3) {
            // (alpha2, alpha3]
            Ogrid = 194;
        }

        // The coordinates of the atom on which the present grid is centred are put into the array "Ogridr"
        std::vector<TlPosition> grid(Ogrid);
        std::vector<double> lebWeight(Ogrid);
        this->points2(Ogrid, r0, this->coord_[iAtom], weight,
                      grid, lebWeight);

        // Loop for the grid number of Omega vector for SG-1
        for (int Omega = 0; Omega < Ogrid; ++Omega) {
            double weight_omega = lebWeight[Omega];
            const TlPosition pos_O = grid[Omega];

            // グリッド省略判定
            const int numOfAtoms = this->numOfRealAtoms_;
            std::vector<double> rr(numOfAtoms);
            bool bPass = true;
            for (int p = 0; p < numOfAtoms; ++p) {
                rr[p] = pos_O.distanceFrom(this->coord_[p]);

                if (rr[p] < this->maxRadii_) {
                    bPass = false;
                }
            }
            if (bPass == true) {
                continue;
            }

            // Loop for the m-center
#pragma omp parallel for
            for (int mc = 0; mc < numOfAtoms; ++mc) {
                double Psuij = 1.0;
                for (int n = 0; n < numOfAtoms; ++n) {
                    // 隣接番号の原子からサーチ
                    int nc = iAtom + n;
                    if (nc >= numOfAtoms) {
                        nc -= numOfAtoms;
                    }

                    if (nc != mc) {
                        const double uij = (rr[mc] - rr[nc]) / this->distanceMatrix_(mc, nc);

                        double suij = 0.0;
                        if (uij > 0.9) {
                            Psuij = 0.0;
                            break;
                        } else if (uij < -0.9) {
                            suij = 1.0;
                        } else {
                            const double f1u = 1.5 * uij - 0.5 * uij * uij * uij;
                            const double f2u = 1.5 * f1u - 0.5 * f1u * f1u * f1u;
                            const double f3u = 1.5 * f2u - 0.5 * f2u * f2u * f2u;
                            suij = 0.5 * (1.0 - f3u);
                        }

                        Psuij *= suij;
                    }
                }

                Ps[mc] = Psuij;
            }

            // Normalization
            double Pstotal = 0.0;
#pragma omp parallel for reduction(+:Pstotal)
            for (int mc = 0; mc < numOfAtoms; ++mc) {
                Pstotal += Ps[mc];
            }
            weight_omega *= (Ps[iAtom] / Pstotal);

            // weight cutoff
            //#pragma omp critical (generateGrid_SG1)
            {
                if (std::fabs(weight_omega) > this->weightCutoff_) {
                    weightvec[GPthrnum] = weight_omega;
                    crdpoint[GPthrnum] = grid[Omega];
                    ++GPthrnum;
                }
            }
        }
    }
    //++GPthrnum;

    // save
    pCoordX->resize(GPthrnum);
    pCoordY->resize(GPthrnum);
    pCoordZ->resize(GPthrnum);
    for (int i = 0; i < GPthrnum; ++i) {
        (*pCoordX)[i] = crdpoint[i].x();
        (*pCoordY)[i] = crdpoint[i].y();
        (*pCoordZ)[i] = crdpoint[i].z();
    }
    weightvec.resize(GPthrnum);
    *pWeight = weightvec;
}


// Set the grid coordinates in "Ogridr" and the Lebedev weight in "w"
// Get atom core coordinates
void DfGenerateGrid::points2(const int nOgrid, const double r0, const TlPosition& core,
                             const double weight, std::vector<TlPosition>& Ogrid, std::vector<double>& w)
{
    double bm[4];
    double bl[4];

    // 初期化
    Ogrid.resize(nOgrid);
    w.resize(nOgrid);
    for (int i = 0; i < nOgrid; ++i) {
        w[i] = weight;
        Ogrid[i] = TlPosition(0.0, 0.0, 0.0);
    }

    if ((this->m_gridType == COARSE) ||
        (this->m_gridType == FINE)) {
        std::cout << "Sorry, do not support now." << std::endl;
    } else if (this->m_gridType == SG_1) {
        //  SG-1 Metod
        if (nOgrid == 6) {
            //*******Set the Lebedev Grid of 6 points**********************
            //Set the Lebedev weight
            for (int i = 0; i < nOgrid; ++i) {
                w[i] *= 0.166666666666667;
            }

            //Set the grid coordinates
            Ogrid[ 0] = TlPosition(1.0,  0.0,  0.0);
            Ogrid[ 1] = TlPosition(0.0,  1.0,  0.0);
            Ogrid[ 2] = TlPosition(0.0,  0.0,  1.0);
            Ogrid[ 3] = TlPosition(-1.0,  0.0,  0.0);
            Ogrid[ 4] = TlPosition(0.0, -1.0,  0.0);
            Ogrid[ 5] = TlPosition(0.0,  0.0, -1.0);

            // Scale and Shift
            for (int i = 0; i < nOgrid; ++i) {
                Ogrid[i] *= r0;
                Ogrid[i] += core;
            }

        } else if (nOgrid == 38) {
            //******Set the Lebedev Grid of 38 points*******************
            // Set the Lebedev parameter
            const double A1 = 0.00952380952387;
            const double A3 = 0.0321428571429;
            const double C1 = 0.0285714285714;
            const double p1 = 0.888073833977;
            const double q1 = 0.459700843381;

            // Set the Lebedev weight
            for (int i = 0; i < 6; ++i) {
                w[i] *=  A1;
            }
            for (int i = 6; i < 14; ++i) {
                w[i] *=  A3;
            }
            for (int i = 14; i < 38; ++i) {
                w[i] *= C1;
            }

            // Set the grid coordinates
            Ogrid[ 0] = TlPosition(1.0,  0.0,  0.0);
            Ogrid[ 1] = TlPosition(0.0,  1.0,  0.0);
            Ogrid[ 2] = TlPosition(0.0,  0.0,  1.0);
            Ogrid[ 3] = TlPosition(-1.0,  0.0,  0.0);
            Ogrid[ 4] = TlPosition(0.0, -1.0,  0.0);
            Ogrid[ 5] = TlPosition(0.0,  0.0, -1.0);

            Ogrid[ 6] = TlPosition(SQ1_3,  SQ1_3,  SQ1_3);
            Ogrid[ 7] = TlPosition(-SQ1_3,  SQ1_3,  SQ1_3);
            Ogrid[ 8] = TlPosition(SQ1_3, -SQ1_3,  SQ1_3);
            Ogrid[ 9] = TlPosition(SQ1_3,  SQ1_3, -SQ1_3);
            Ogrid[10] = TlPosition(-SQ1_3, -SQ1_3,  SQ1_3);
            Ogrid[11] = TlPosition(-SQ1_3,  SQ1_3, -SQ1_3);
            Ogrid[12] = TlPosition(SQ1_3, -SQ1_3, -SQ1_3);
            Ogrid[13] = TlPosition(-SQ1_3, -SQ1_3, -SQ1_3);

            // set ci1
            Ogrid[14] = TlPosition(p1,  q1, 0.0);
            Ogrid[15] = TlPosition(p1, -q1, 0.0);
            Ogrid[16] = TlPosition(-p1,  q1, 0.0);
            Ogrid[17] = TlPosition(-p1, -q1, 0.0);
            Ogrid[18] = TlPosition(p1, 0.0,  q1);
            Ogrid[19] = TlPosition(p1, 0.0, -q1);
            Ogrid[20] = TlPosition(-p1, 0.0,  q1);
            Ogrid[21] = TlPosition(-p1, 0.0, -q1);
            Ogrid[22] = TlPosition(0.0,  p1,  q1);
            Ogrid[23] = TlPosition(0.0,  p1, -q1);
            Ogrid[24] = TlPosition(0.0, -p1,  q1);
            Ogrid[25] = TlPosition(0.0, -p1, -q1);

            Ogrid[26] = TlPosition(q1,  p1, 0.0);
            Ogrid[27] = TlPosition(q1, -p1, 0.0);
            Ogrid[28] = TlPosition(-q1,  p1, 0.0);
            Ogrid[29] = TlPosition(-q1, -p1, 0.0);
            Ogrid[30] = TlPosition(q1, 0.0,  p1);
            Ogrid[31] = TlPosition(q1, 0.0, -p1);
            Ogrid[32] = TlPosition(-q1, 0.0,  p1);
            Ogrid[33] = TlPosition(-q1, 0.0, -p1);
            Ogrid[34] = TlPosition(0.0,  q1,  p1);
            Ogrid[35] = TlPosition(0.0,  q1, -p1);
            Ogrid[36] = TlPosition(0.0, -q1,  p1);
            Ogrid[37] = TlPosition(0.0, -q1, -p1);

            // Scale and Shift
            for (int i = 0; i < nOgrid; ++i) {
                Ogrid[i] *= r0;
                Ogrid[i] += core;
            }

        } else if (nOgrid == 86) {
            // ******Set the Lebedev grid of 86 points*****************
            // Set the Lebedev parameter

            const double A1 = 0.0115440115441;
            const double A3 = 0.0119439090859;
            const double B1 = 0.0111105557106;
            const double B2 = 0.0118765012945;
            const double C1 = 0.0118123037469;
            bm[0] = 0.852518311701;
            bm[1] = 0.189063552885;
            bl[0] = 0.369602846454;
            bl[1] = 0.694354006603;
            const double p1 = 0.927330657151;
            const double q1 = 0.374243039090;

            // Set the Lebedev weight
            for (int i = 0; i < 6; ++i) {
                w[i] *= A1;
            }
            for (int i = 6; i < 14; ++i) {
                w[i] *= A3;
            }
            for (int i = 14; i < 38; ++i) {
                w[i] *= B1;
            }
            for (int i = 38; i < 62; ++i) {
                w[i] *= B2;
            }
            for (int i = 62; i < 86; ++i) {
                w[i] *= C1;
            }

            // Set the grid coordinates
            Ogrid[ 0] = TlPosition(1.0,  0.0,  0.0);
            Ogrid[ 1] = TlPosition(0.0,  1.0,  0.0);
            Ogrid[ 2] = TlPosition(0.0,  0.0,  1.0);
            Ogrid[ 3] = TlPosition(-1.0,  0.0,  0.0);
            Ogrid[ 4] = TlPosition(0.0, -1.0,  0.0);
            Ogrid[ 5] = TlPosition(0.0,  0.0, -1.0);

            Ogrid[ 6] = TlPosition(SQ1_3,  SQ1_3,  SQ1_3);
            Ogrid[ 7] = TlPosition(-SQ1_3,  SQ1_3,  SQ1_3);
            Ogrid[ 8] = TlPosition(SQ1_3, -SQ1_3,  SQ1_3);
            Ogrid[ 9] = TlPosition(SQ1_3,  SQ1_3, -SQ1_3);
            Ogrid[10] = TlPosition(-SQ1_3, -SQ1_3,  SQ1_3);
            Ogrid[11] = TlPosition(-SQ1_3,  SQ1_3, -SQ1_3);
            Ogrid[12] = TlPosition(SQ1_3, -SQ1_3, -SQ1_3);
            Ogrid[13] = TlPosition(-SQ1_3, -SQ1_3, -SQ1_3);

            // set bik
            for (int i = 0; i < 2; ++i) {
                // bl[i] = sqrt( 1.0 - bm[i] * bm[i] ) / SQ2;
                Ogrid[ 14 + 24*i + 0 ] = TlPosition(bl[i],  bl[i],  bm[i]);
                Ogrid[ 14 + 24*i + 1 ] = TlPosition(-bl[i],  bl[i],  bm[i]);
                Ogrid[ 14 + 24*i + 2 ] = TlPosition(bl[i], -bl[i],  bm[i]);
                Ogrid[ 14 + 24*i + 3 ] = TlPosition(bl[i],  bl[i], -bm[i]);
                Ogrid[ 14 + 24*i + 4 ] = TlPosition(-bl[i], -bl[i],  bm[i]);
                Ogrid[ 14 + 24*i + 5 ] = TlPosition(-bl[i],  bl[i], -bm[i]);
                Ogrid[ 14 + 24*i + 6 ] = TlPosition(bl[i], -bl[i], -bm[i]);
                Ogrid[ 14 + 24*i + 7 ] = TlPosition(-bl[i], -bl[i], -bm[i]);
            }
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 8; ++j) {
                    TlPosition tmp = Ogrid[ 14 + 24*i + j ];
                    Ogrid[ 14 + 24*i +  8 + j ] = TlPosition(tmp.y(), tmp.z(), tmp.x());
                    Ogrid[ 14 + 24*i + 16 + j ] = TlPosition(tmp.z(), tmp.x(), tmp.y());
                }
            }

            // set ci1
            Ogrid[ 62 ] = TlPosition(p1,  q1, 0.0);
            Ogrid[ 63 ] = TlPosition(p1, -q1, 0.0);
            Ogrid[ 64 ] = TlPosition(-p1,  q1, 0.0);
            Ogrid[ 65 ] = TlPosition(-p1, -q1, 0.0);

            Ogrid[ 66 ] = TlPosition(p1, 0.0,  q1);
            Ogrid[ 67 ] = TlPosition(p1, 0.0, -q1);
            Ogrid[ 68 ] = TlPosition(-p1, 0.0,  q1);
            Ogrid[ 69 ] = TlPosition(-p1, 0.0, -q1);

            Ogrid[ 70 ] = TlPosition(0.0,  p1,  q1);
            Ogrid[ 71 ] = TlPosition(0.0,  p1, -q1);
            Ogrid[ 72 ] = TlPosition(0.0, -p1,  q1);
            Ogrid[ 73 ] = TlPosition(0.0, -p1, -q1);

            Ogrid[ 74 ] = TlPosition(q1,  p1, 0.0);
            Ogrid[ 75 ] = TlPosition(q1, -p1, 0.0);
            Ogrid[ 76 ] = TlPosition(-q1,  p1, 0.0);
            Ogrid[ 77 ] = TlPosition(-q1, -p1, 0.0);

            Ogrid[ 78 ] = TlPosition(q1, 0.0,  p1);
            Ogrid[ 79 ] = TlPosition(q1, 0.0, -p1);
            Ogrid[ 80 ] = TlPosition(-q1, 0.0,  p1);
            Ogrid[ 81 ] = TlPosition(-q1, 0.0, -p1);

            Ogrid[ 82 ] = TlPosition(0.0,  q1,  p1);
            Ogrid[ 83 ] = TlPosition(0.0,  q1, -p1);
            Ogrid[ 84 ] = TlPosition(0.0, -q1,  p1);
            Ogrid[ 85 ] = TlPosition(0.0, -q1, -p1);

            // Scale and Shift
            for (int i = 0; i < nOgrid; ++i) {
                Ogrid[i] *= r0;
                Ogrid[i] += core;
            }

        } else if (nOgrid == 194) {
            // ******Set the Lebedev grid of 194 points******************
            // Set the Lebedev parameter
            const double A1 = 0.00178234044724;
            const double A2 = 0.00571690594998;
            const double A3 = 0.00557338317884;
            const double B1 = 0.00551877146727;
            const double B2 = 0.00515823771181;
            const double B3 = 0.00560870408259;
            const double B4 = 0.00410677702817;
            const double C1 = 0.00505184606462;
            const double D1 = 0.00553024891623;
            bm[0]  = 0.777493219315;
            bm[1]  = 0.912509096867;
            bm[2]  = 0.314196994183;
            bm[3]  = 0.982972302707;
            bl[0]  = 0.444693317871;
            bl[1]  = 0.289246562758;
            bl[2]  = 0.671297344270;
            bl[3]  = 0.129933544765;

            const double p1 = 0.938319218138;
            const double q1 = 0.345770219761;
            const double dr = 0.836036015482;
            const double du = 0.159041710538;
            const double dw = 0.525118572443;

            // Set the Lebedev weight
            for (int i = 0; i < 6; ++i) {
                w[i] *= A1;
            }
            for (int i = 6; i < 18; ++i) {
                w[i] *= A2;
            }
            for (int i = 18; i < 26; ++i) {
                w[i] *= A3;
            }
            for (int i = 26; i < 50; ++i) {
                w[i] *= B1;
            }
            for (int i = 50; i < 74; ++i) {
                w[i] *= B2;
            }
            for (int i = 74; i < 98; ++i) {
                w[i] *= B3;
            }
            for (int i = 98; i < 122; ++i) {
                w[i] *= B4;
            }
            for (int i = 122; i < 146; ++i) {
                w[i] *= C1;
            }
            for (int i = 146; i < 194; ++i) {
                w[i] *= D1;
            }

            // Set the grid coordinates
            Ogrid[ 0] = TlPosition(1.0,  0.0,  0.0);
            Ogrid[ 1] = TlPosition(0.0,  1.0,  0.0);
            Ogrid[ 2] = TlPosition(0.0,  0.0,  1.0);
            Ogrid[ 3] = TlPosition(-1.0,  0.0,  0.0);
            Ogrid[ 4] = TlPosition(0.0, -1.0,  0.0);
            Ogrid[ 5] = TlPosition(0.0,  0.0, -1.0);

            Ogrid[ 6 ] = TlPosition(SQ1_2,  SQ1_2, 0.0);
            Ogrid[ 7 ] = TlPosition(-SQ1_2,  SQ1_2, 0.0);
            Ogrid[ 8 ] = TlPosition(SQ1_2, -SQ1_2, 0.0);
            Ogrid[ 9 ] = TlPosition(-SQ1_2, -SQ1_2, 0.0);

            for (int i = 0; i < 4; ++i) {
                const TlPosition tmp = Ogrid[ 6 + i ];
                Ogrid[ 6 + 4 + i ] = TlPosition(tmp.y(), tmp.z(), tmp.x());
                Ogrid[ 6 + 8 + i ] = TlPosition(tmp.z(), tmp.x(), tmp.y());
            }

            Ogrid[ 18 ] = TlPosition(SQ1_3,  SQ1_3,  SQ1_3);
            Ogrid[ 19 ] = TlPosition(-SQ1_3,  SQ1_3,  SQ1_3);
            Ogrid[ 20 ] = TlPosition(SQ1_3, -SQ1_3,  SQ1_3);
            Ogrid[ 21 ] = TlPosition(SQ1_3,  SQ1_3, -SQ1_3);
            Ogrid[ 22 ] = TlPosition(-SQ1_3, -SQ1_3,  SQ1_3);
            Ogrid[ 23 ] = TlPosition(-SQ1_3,  SQ1_3, -SQ1_3);
            Ogrid[ 24 ] = TlPosition(SQ1_3, -SQ1_3, -SQ1_3);
            Ogrid[ 25 ] = TlPosition(-SQ1_3, -SQ1_3, -SQ1_3);

            // set bik; no need to set cik
            for (int i = 0; i < 4; ++i) {
                // bl[i] = sqrt( 1.0 - bm[i] * bm[i] ) / SQ2;
                Ogrid[ 26 + 24*i + 0 ] = TlPosition(bl[i],  bl[i],  bm[i]);
                Ogrid[ 26 + 24*i + 1 ] = TlPosition(-bl[i],  bl[i],  bm[i]);
                Ogrid[ 26 + 24*i + 2 ] = TlPosition(bl[i], -bl[i],  bm[i]);
                Ogrid[ 26 + 24*i + 3 ] = TlPosition(bl[i],  bl[i], -bm[i]);
                Ogrid[ 26 + 24*i + 4 ] = TlPosition(-bl[i], -bl[i],  bm[i]);
                Ogrid[ 26 + 24*i + 5 ] = TlPosition(-bl[i],  bl[i], -bm[i]);
                Ogrid[ 26 + 24*i + 6 ] = TlPosition(bl[i], -bl[i], -bm[i]);
                Ogrid[ 26 + 24*i + 7 ] = TlPosition(-bl[i], -bl[i], -bm[i]);
            }

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 8; ++j) {
                    const TlPosition tmp = Ogrid[ 26 + 24*i + j ];
                    Ogrid[ 26 + 24*i +  8 + j ] = TlPosition(tmp.y(), tmp.z(), tmp.x());
                    Ogrid[ 26 + 24*i + 16 + j ] = TlPosition(tmp.z(), tmp.x(), tmp.y());
                }
            }

            // set ci1
            Ogrid[ 122 ] = TlPosition(p1,  q1, 0.0);
            Ogrid[ 123 ] = TlPosition(p1, -q1, 0.0);
            Ogrid[ 124 ] = TlPosition(-p1,  q1, 0.0);
            Ogrid[ 125 ] = TlPosition(-p1, -q1, 0.0);

            Ogrid[ 126 ] = TlPosition(p1, 0.0,  q1);
            Ogrid[ 127 ] = TlPosition(p1, 0.0, -q1);
            Ogrid[ 128 ] = TlPosition(-p1, 0.0,  q1);
            Ogrid[ 129 ] = TlPosition(-p1, 0.0, -q1);

            Ogrid[ 130 ] = TlPosition(0.0,  p1,  q1);
            Ogrid[ 131 ] = TlPosition(0.0,  p1, -q1);
            Ogrid[ 132 ] = TlPosition(0.0, -p1,  q1);
            Ogrid[ 133 ] = TlPosition(0.0, -p1, -q1);

            Ogrid[ 134 ] = TlPosition(q1,  p1, 0.0);
            Ogrid[ 135 ] = TlPosition(q1, -p1, 0.0);
            Ogrid[ 136 ] = TlPosition(-q1,  p1, 0.0);
            Ogrid[ 137 ] = TlPosition(-q1, -p1, 0.0);

            Ogrid[ 138 ] = TlPosition(q1, 0.0,  p1);
            Ogrid[ 139 ] = TlPosition(q1, 0.0, -p1);
            Ogrid[ 140 ] = TlPosition(-q1, 0.0,  p1);
            Ogrid[ 141 ] = TlPosition(-q1, 0.0, -p1);

            Ogrid[ 142 ] = TlPosition(0.0,  q1,  p1);
            Ogrid[ 143 ] = TlPosition(0.0,  q1, -p1);
            Ogrid[ 144 ] = TlPosition(0.0, -q1,  p1);
            Ogrid[ 145 ] = TlPosition(0.0, -q1, -p1);

            // set dik
            Ogrid[ 146 + 0 ] = TlPosition(dr,  du,  dw);
            Ogrid[ 146 + 1 ] = TlPosition(-dr,  du,  dw);
            Ogrid[ 146 + 2 ] = TlPosition(dr, -du,  dw);
            Ogrid[ 146 + 3 ] = TlPosition(dr,  du, -dw);
            Ogrid[ 146 + 4 ] = TlPosition(-dr, -du,  dw);
            Ogrid[ 146 + 5 ] = TlPosition(-dr,  du, -dw);
            Ogrid[ 146 + 6 ] = TlPosition(dr, -du, -dw);
            Ogrid[ 146 + 7 ] = TlPosition(-dr, -du, -dw);

            for (int i = 0; i < 8; ++i) {
                const TlPosition tmp = Ogrid[ 146 + i ];
                Ogrid[ 146 +  8 + i ] = TlPosition(tmp.x(), tmp.z(), tmp.y());
                Ogrid[ 146 + 16 + i ] = TlPosition(tmp.y(), tmp.x(), tmp.z());
                Ogrid[ 146 + 24 + i ] = TlPosition(tmp.y(), tmp.z(), tmp.x());
                Ogrid[ 146 + 32 + i ] = TlPosition(tmp.z(), tmp.y(), tmp.x());
                Ogrid[ 146 + 40 + i ] = TlPosition(tmp.z(), tmp.x(), tmp.y());
            }

            // Scale and Shift
            for (int i = 0; i < nOgrid; ++i) {
                Ogrid[i] *= r0;
                Ogrid[i] += core;
            }
        }
    } else if ((this->m_gridType == MEDIUM) ||
               (this->m_gridType == MEDIUM_FINE)) {
        const double A1 = 0.000599631368862;
        const double A2 = 0.00737299971862;
        const double A3 = 0.00721051536014;
        const double B1 = 0.00757439415905;
        const double B2 = 0.00675382948631;
        const double B3 = 0.00711635549312;
        const double D1 = 0.00699108735330;

        const double dr = 0.882270011260;
        const double du = 0.140355381171;
        const double dw = 0.449332832327;

        const double CORRP = 1.00;
        const double invCORRP = 1.0 / CORRP;

        bm[0] = 0.974888643677;
        bm[1] = 0.807089818360;
        bm[2] = 0.291298882210;

        // Set Lebedev weight
        // Note CORRP is a temporal correction.
        for (int i =  0; i <   6; ++i) {
            // w[i] *= A1 * 146.0;
            w[i] *= A1 * 146.0 * invCORRP;
        }
        for (int i =  6; i <  18; ++i) {
            // w[i] *= A2 * 146.0;
            w[i] *= A2 * 146.0 * invCORRP;
        }
        for (int i = 18; i <  26; ++i) {
            // w[i] *= A3 * 146.0;
            w[i] *= A3 * 146.0 * invCORRP;
        }
        for (int i = 26; i <  50; ++i) {
            // w[i] *= B1 * 146.0;
            w[i] *= B1 * 146.0 * invCORRP;
        }
        for (int i = 50; i <  74; ++i) {
            // w[i] *= B2 * 146.0;
            w[i] *= B2 * 146.0 * invCORRP;
        }
        for (int i = 74; i <  98; ++i) {
            // w[i] *= B3 * 146.0;
            w[i] *= B3 * 146.0 * invCORRP;
        }
        for (int i = 98; i < 146; ++i) {
            // w[i] *= D1 * 146.0;
            w[i] *= D1 * 146.0 * invCORRP;
        }

        // Generate Omega Grid
        // set aik
        Ogrid[ 0 ] = TlPosition(1.0,  0.0,  0.0);
        Ogrid[ 1 ] = TlPosition(0.0,  1.0,  0.0);
        Ogrid[ 2 ] = TlPosition(0.0,  0.0,  1.0);
        Ogrid[ 3 ] = TlPosition(-1.0,  0.0,  0.0);
        Ogrid[ 4 ] = TlPosition(0.0, -1.0,  0.0);
        Ogrid[ 5 ] = TlPosition(0.0,  0.0, -1.0);

        Ogrid[ 6 ] = TlPosition(SQ1_2,  SQ1_2, 0.0);
        Ogrid[ 7 ] = TlPosition(-SQ1_2,  SQ1_2, 0.0);
        Ogrid[ 8 ] = TlPosition(SQ1_2, -SQ1_2, 0.0);
        Ogrid[ 9 ] = TlPosition(-SQ1_2, -SQ1_2, 0.0);

        for (int i = 0; i < 4; ++i) {
            const TlPosition tmp = Ogrid[ 6 + i ];
            Ogrid[ 6 + 4 + i ] = TlPosition(tmp.y(), tmp.z(), tmp.x());
            Ogrid[ 6 + 8 + i ] = TlPosition(tmp.z(), tmp.x(), tmp.y());
        }

        Ogrid[ 18 ] = TlPosition(SQ1_3,  SQ1_3,  SQ1_3);
        Ogrid[ 19 ] = TlPosition(-SQ1_3,  SQ1_3,  SQ1_3);
        Ogrid[ 20 ] = TlPosition(SQ1_3, -SQ1_3,  SQ1_3);
        Ogrid[ 21 ] = TlPosition(SQ1_3,  SQ1_3, -SQ1_3);
        Ogrid[ 22 ] = TlPosition(-SQ1_3, -SQ1_3,  SQ1_3);
        Ogrid[ 23 ] = TlPosition(-SQ1_3,  SQ1_3, -SQ1_3);
        Ogrid[ 24 ] = TlPosition(SQ1_3, -SQ1_3, -SQ1_3);
        Ogrid[ 25 ] = TlPosition(-SQ1_3, -SQ1_3, -SQ1_3);

        // set bik; no need to set cik
        for (int i = 0; i < 3; ++i) {
            bl[i] = sqrt(1.0 - bm[i] * bm[i]) / SQ2;

            Ogrid[ 26 + 24*i + 0 ] = TlPosition(bl[i],  bl[i],  bm[i]);
            Ogrid[ 26 + 24*i + 1 ] = TlPosition(-bl[i],  bl[i],  bm[i]);
            Ogrid[ 26 + 24*i + 2 ] = TlPosition(bl[i], -bl[i],  bm[i]);
            Ogrid[ 26 + 24*i + 3 ] = TlPosition(bl[i],  bl[i], -bm[i]);
            Ogrid[ 26 + 24*i + 4 ] = TlPosition(-bl[i], -bl[i],  bm[i]);
            Ogrid[ 26 + 24*i + 5 ] = TlPosition(-bl[i],  bl[i], -bm[i]);
            Ogrid[ 26 + 24*i + 6 ] = TlPosition(bl[i], -bl[i], -bm[i]);
            Ogrid[ 26 + 24*i + 7 ] = TlPosition(-bl[i], -bl[i], -bm[i]);
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 8; ++j) {
                const TlPosition tmp = Ogrid[ 26 + 24*i + j ];
                Ogrid[ 26 + 24*i +  8 + j ] = TlPosition(tmp.y(), tmp.z(), tmp.x());
                Ogrid[ 26 + 24*i + 16 + j ] = TlPosition(tmp.z(), tmp.x(), tmp.y());
            }
        }

        // set dik
        Ogrid[ 98 + 0 ] = TlPosition(dr,  du,  dw);
        Ogrid[ 98 + 1 ] = TlPosition(-dr,  du,  dw);
        Ogrid[ 98 + 2 ] = TlPosition(dr, -du,  dw);
        Ogrid[ 98 + 3 ] = TlPosition(dr,  du, -dw);
        Ogrid[ 98 + 4 ] = TlPosition(-dr, -du,  dw);
        Ogrid[ 98 + 5 ] = TlPosition(-dr,  du, -dw);
        Ogrid[ 98 + 6 ] = TlPosition(dr, -du, -dw);
        Ogrid[ 98 + 7 ] = TlPosition(-dr, -du, -dw);

        for (int i = 0; i < 8; ++i) {
            const TlPosition tmp = Ogrid[ 98 + i ];
            Ogrid[ 98 +  8 + i ] = TlPosition(tmp.x(), tmp.z(), tmp.y());
            Ogrid[ 98 + 16 + i ] = TlPosition(tmp.y(), tmp.x(), tmp.z());
            Ogrid[ 98 + 24 + i ] = TlPosition(tmp.y(), tmp.z(), tmp.x());
            Ogrid[ 98 + 32 + i ] = TlPosition(tmp.z(), tmp.y(), tmp.x());
            Ogrid[ 98 + 40 + i ] = TlPosition(tmp.z(), tmp.x(), tmp.y());
        }

        // Scale and Shift
        for (int i = 0; i < nOgrid; ++i) {
            Ogrid[i] *= r0;
            Ogrid[i] += core;
        }
    }

}
