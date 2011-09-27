#include <ios>
#include <cmath>

#include "DfDmatrix.h"
#include "Fl_Out.h"
#include "Fl_Tbl_Orbital.h"
#include "Fl_Tbl_Orbital.h"
#include "TlUtils.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlStringTokenizer.h"
#include "TlFile.h"

/*********************************************************
MO_OVERLAP_ITER:
軌道の重なりの対応を使用しはじめるiteration回数を指定する。
**********************************************************/

DfDmatrix::DfDmatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;

    this->orbitalCorrespondenceMethod_ = OCM_NONE;
    const bool isOrbitalCorrespondence = pdfParam["model"]["orbital-correspondence"].getBoolean();
    const int startItr = pdfParam["model"]["orbital-correspondence-start"].getInt();
    if (isOrbitalCorrespondence == true) {
        const std::string method = TlUtils::toUpper(pdfParam["model"]["orbital-correspondence-method"].getStr());
        if ((method == "MO-OVERLAP") && (this->m_nIteration >= startItr)) {
            this->orbitalCorrespondenceMethod_ = OCM_OVERLAP;
        } else if ((method == "MO-PROJECTION") && (this->m_nIteration >= startItr)) {
            this->orbitalCorrespondenceMethod_ = OCM_PROJECTION;
        }
    }
}


DfDmatrix::~DfDmatrix()
{
}


void DfDmatrix::DfDmatrixMain()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main(RUN_RKS);
        break;

    case METHOD_UKS:
        this->main(RUN_UKS_ALPHA);
        this->main(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->main(RUN_ROKS);
        break;

    default:
        CnErr.abort();
        break;
    }
}


void DfDmatrix::main(const DfObject::RUN_TYPE runType)
{
    // occupation
    TlVector currOcc;
    switch (this->orbitalCorrespondenceMethod_) {
    case OCM_OVERLAP:
        this->logger(" orbital correspondence method: MO-overlap\n");
        currOcc = this->getOccupationUsingOverlap<TlMatrix>(runType);
        currOcc.save(this->getOccupationPath(runType));
        break;

    case OCM_PROJECTION:
        this->logger(" orbital correspondence method: MO-projection\n");
        currOcc = this->getOccupationUsingProjection<TlMatrix, TlSymmetricMatrix>(runType);
        currOcc.save(this->getOccupationPath(runType));
        break;

    default:
        this->logger(" orbital correspondence method: none\n");
        if (TlFile::isExist(this->getOccupationPath(runType)) == true) {
            currOcc.load(this->getOccupationPath(runType));
        } else {
            currOcc = this->createOccupation(runType);
            currOcc.save(this->getOccupationPath(runType));
        }
        break;
    }
    
    this->generateDensityMatrix<TlMatrix, TlSymmetricMatrix>(runType, currOcc);
}

TlVector DfDmatrix::getOccupation(const DfObject::RUN_TYPE runType)
{
    const std::string sFileName = this->getOccupationPath(runType);

    TlVector occ;
    occ.load(sFileName);
    assert(occ.getSize() == this->m_nNumOfMOs);

    return occ;
}

void DfDmatrix::checkOccupation(const TlVector& prevOcc, const TlVector& currOcc)
{
    const double xx = prevOcc.sum();
    const double yy = currOcc.sum();

    if (std::fabs(xx - yy) > 1.0e-10) {
        TlLogX& log = TlLogX::getInstance();
        this->logger(" #####   SUM pre_occ != SUM crr_occ #####\n");
        this->logger(TlUtils::format(" SUM pre_occ is %10.4lf,  SUM crr_occ is %10.4lf\n", xx, yy));
        this->logger("previous occupation\n");
        prevOcc.print(log);
        this->logger("current  occupation\n");
        currOcc.print(log);

        CnErr.abort("DfDmatrix", "", "", "SUM error for occ !!");
    }
}


void DfDmatrix::printOccupation(const TlVector& occ)
{
    TlLogX& log = TlLogX::getInstance();
    occ.print(log);
}


// print out Two Vectors' elements
void DfDmatrix::printTwoVectors(const TlVector& a, const TlVector& b,
                                const std::string& title, int pnumcol)
{
    this->logger(TlUtils::format("\n\n       %s\n\n", title.c_str()));

    if (a.getSize() != b.getSize()) {
        this->logger("DfDmatrix::printTwoVectors() : dimensions of vectors are not much\n");
        this->logger("                             : omit\n");
    }

    this->logger("       two vectors\n");

    const int number_of_emt = a.getSize();
    for (int ord = 0; ord < number_of_emt; ord += pnumcol) {
        this->logger("       ");
        for (int j = ord; ((j < ord + pnumcol) && (j < number_of_emt)); ++j) {
            this->logger(TlUtils::format("   %5d th", j+1));
        }
        this->logger("\n     ");

        for (int j = ord; ((j < ord + pnumcol) && (j < number_of_emt)); ++j) {
            this->logger("-----------");
        }
        this->logger("----\n       ");

        for (int j = ord; j < ord + pnumcol && j < number_of_emt; ++j) {
            const double aj = a[j];
            this->logger(TlUtils::format(" %6.0lf    ", aj + 1));
        }
        this->logger("\n\n       ");

        for (int j = ord; ((j < ord + pnumcol) && (j < number_of_emt)); ++j) {
            const double bj = b[j];
            this->logger(TlUtils::format(" %10.6lf", bj));
        }
        this->logger("\n\n");
    }
}


TlVector DfDmatrix::createOccupation(const DfObject::RUN_TYPE runType)
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    
    // construct guess occupations
    TlVector occ(this->m_nNumOfMOs);
    switch (runType) {
    case RUN_RKS: {
        std::vector<int> docLevel = this->getLevel(pdfParam["model"]["method/nsp/occlevel"].getStr());
        for (std::vector<int>::const_iterator p = docLevel.begin(); p != docLevel.end(); p++) {
            occ[*p -1] = 2.0;
        }
    }
    break;

    case RUN_UKS_ALPHA: {
        std::vector<int> aoocLevel = this->getLevel(pdfParam["model"]["method/sp/alpha-spin-occlevel"].getStr());
        for (std::vector<int>::const_iterator p = aoocLevel.begin(); p != aoocLevel.end(); p++) {
            occ[*p -1] = 1.0;
        }
    }
    break;

    case RUN_UKS_BETA: {
        std::vector<int> boocLevel = this->getLevel(pdfParam["model"]["method/sp/beta-spin-occlevel"].getStr());
        for (std::vector<int>::const_iterator p = boocLevel.begin(); p != boocLevel.end(); p++) {
            occ[*p -1] = 1.0;
        }
    }
    break;
    case RUN_ROKS: {
        std::vector<int> docLevel = this->getLevel(pdfParam["model"]["method/roks/closed-shell"].getStr());
        for (std::vector<int>::const_iterator p = docLevel.begin(); p != docLevel.end(); p++) {
            occ[*p -1] = 2.0;
        }

        std::vector<int> socLevel = this->getLevel(pdfParam["model"]["method/roks/open-shell"].getStr());
        for (std::vector<int>::const_iterator p = socLevel.begin(); p != socLevel.end(); p++) {
            // nsoc
            occ[*p -1] = 1.0;
        }
    }
    break;

    default:
        std::cerr << "program error. (DfDmatrix::createOccupation)" << std::endl;
        break;
    }

    return occ;
}

std::vector<int> DfDmatrix::getLevel(std::string sLevel)
{
    std::vector<int> answer;
    answer.clear();

    // "-" を " - " に置換
    TlUtils::replace(sLevel, "-", " - ");

    TlStringTokenizer token(sLevel);
    bool bRegionMode = false;
    int nPrevIndex = 0;
    while (token.hasMoreTokens()) {
        std::string tmp = token.nextToken();

        if (tmp == "nil") {
            continue;
        }

        if (tmp == "-") {
            if (nPrevIndex != 0) {
                bRegionMode = true;
                continue;
            } else {
                abort();
                //CnErr.abort("DfPreScf", "", "putdoclevel", "syntax error.");
            }
        }

        int nIndex = atoi(tmp.c_str());
        //std::cerr << "nIndex = " << nIndex << std::endl;
        if (nIndex > 0) {
            if (bRegionMode == true) {
                // 数字が xx - yy の形で入力
                for (int i = nPrevIndex +1; i <= nIndex; i++) {
                    answer.push_back(i);
                }
                bRegionMode = false;
                nPrevIndex = 0;
            } else {
                // 数字が単独で入力
                answer.push_back(nIndex);
                nPrevIndex = nIndex;
                continue;
            }
        } else {
            abort();
            //CnErr.abort("DfPreScf", "", "putdoclevel", "syntax error.");
        }
    }


    return answer;
}
