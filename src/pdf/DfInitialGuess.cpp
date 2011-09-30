#include <fstream>
#include "DfInitialGuess.h"
#include "DfInitialguess.h"
#include "DfInitialGuessHuckel.h"
#include "DfInitialGuessHarris.h"
#include "DfDmatrix.h"
#include "TlStringTokenizer.h"

DfInitialGuess::DfInitialGuess(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
}


DfInitialGuess::~DfInitialGuess()
{
}


void DfInitialGuess::exec()
{
    switch (this->initialGuessType_) {
    case GUESS_RHO:
        // go below
    case GUESS_FILE_RHO:
        this->createRho();
        this->createOccupation();
        break;

//         switch (this->m_nMethodType) {
//         case METHOD_RKS:
//             this->saveRho1(RUN_RKS);
//             break;

//         case METHOD_UKS:
//             this->saveRho1(RUN_UKS_ALPHA);
//             this->saveRho1(RUN_UKS_BETA);
//             break;

//         case METHOD_ROKS:
//             this->saveRho1(RUN_ROKS_CLOSE);
//             this->saveRho1(RUN_ROKS_OPEN);
//             break;
//         } 
//         this->createOccupation();
//         break;

    case GUESS_DENSITY:
        this->createOccupation();
        break;

    case GUESS_LCAO:
        this->createInitialGuessUsingLCAO();
        break;

    case GUESS_HUCKEL:
        this->createOccupation();
        this->createInitialGuessUsingHuckel();
        break;

    case GUESS_CORE:
        this->createOccupation();
        this->createInitialGuessUsingCore();
        break;

    case GUESS_HARRIS:
        this->createOccupation();
        this->createInitialGuessUsingHarris();
        break;

    default:
        this->logger("unknown initial guess parameter.\n");
        break;
    }
}

void DfInitialGuess::createRho()
{
    DfInitialguess diguess(this->pPdfParam_);
    diguess.dfGusMain();
}


void DfInitialGuess::saveRho1(const RUN_TYPE runType)
{
    this->loggerTime(" loading initial rho file.");

    TlVector rho;
    std::ifstream ifs;
    ifs.open("guess.rho", std::ios::in);
    if (!ifs) {
        this->logger(" could not open: guess.rho.\n");
    } else {
        std::size_t size;
        ifs >> size;
        rho.resize(size);

        double data = 0.0;
        for (std::size_t i = 0; i < size; ++i) {
            ifs >> data;
            rho[i] = data;
        }
    }
    ifs.close();
    
    this->loggerTime(" saving rho file.");
    rho.save(this->getRhoPath(runType, 1));
    this->loggerTime(" finished.");
}

void DfInitialGuess::createInitialGuessUsingHuckel()
{
    DfInitialGuessHuckel huckel(this->pPdfParam_);
}


void DfInitialGuess::createInitialGuessUsingCore()
{
    DfInitialGuessHuckel huckel(this->pPdfParam_);
}


void DfInitialGuess::createInitialGuessUsingHarris()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        {
            DfInitialGuessHarris harris(this->pPdfParam_);
            harris.main();
        }
        break;

    case METHOD_UKS:
        CnErr.abort("Sorry. harris method is not supported except RKS. stop.\n");
        break;

    case METHOD_ROKS:
        CnErr.abort("Sorry. harris method is not supported except RKS. stop.\n");
        break;

    default:
        CnErr.abort();
        break;
    }
}


void DfInitialGuess::createInitialGuessUsingLCAO()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->createInitialGuessUsingLCAO(RUN_RKS);
        break;

    case METHOD_UKS:
        this->createInitialGuessUsingLCAO(RUN_UKS_ALPHA);
        this->createInitialGuessUsingLCAO(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->createInitialGuessUsingLCAO(RUN_ROKS);
        break;

    default:
        abort();
        break;
    }
}


void DfInitialGuess::createInitialGuessUsingLCAO(const RUN_TYPE runType)
{
    // read guess lcao
    const TlMatrix LCAO = this->getLCAO<TlMatrix>(runType);
    this->saveC0(runType, LCAO);

    // read guess occupation
    const TlVector aOccupation = this->getOccupation(runType);
    this->saveOccupation(runType, aOccupation);

    {
        TlSerializeData tmpParam = *(this->pPdfParam_);
        tmpParam["orbital-correspondence"] = false;
        tmpParam["orbital-overlap-correspondence-method"] = "keep";
        tmpParam["control-iteration"] = 0;

        // 密度行列の作成
        DfDmatrix dfDmatrix(&tmpParam);
        dfDmatrix.DfDmatrixMain(); // RKS only?
    }
}


TlVector DfInitialGuess::getOccupation(const RUN_TYPE runType)
{
    TlVector occupation;
    const std::string sFile = std::string("./guess.occ.") + this->m_sRunTypeSuffix[runType];
    occupation.loadText(sFile.c_str());
//     if (this->m_nNumOfMOs < occupation.getSize()) {
//         this->logger("Occupation vector is shrinked.");
//         occupation.resize(this->m_nNumOfMOs);
//     }

    return occupation;
}


void DfInitialGuess::createOccupation()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->createOccupation(RUN_RKS);
        break;

    case METHOD_UKS:
        this->createOccupation(RUN_UKS_ALPHA);
        this->createOccupation(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->createOccupation(RUN_ROKS);
        break;

    default:
        CnErr.abort();
        break;
    }
}


void DfInitialGuess::createOccupation(const RUN_TYPE runType)
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);

    // construct guess occupations
    TlVector guess_occ(this->m_nNumOfMOs);
    switch (runType) {
    case RUN_RKS:
        {
            std::vector<int> docLevel = this->getLevel(pdfParam["method/nsp/occlevel"].getStr());
            for (std::vector<int>::const_iterator p = docLevel.begin(); p != docLevel.end(); p++) {
                guess_occ[*p -1] = 2.0;
            }
        }
        break;

    case RUN_UKS_ALPHA:
        {
            std::vector<int> aoocLevel = this->getLevel(pdfParam["method/sp/alpha-spin-occlevel"].getStr());
            for (std::vector<int>::const_iterator p = aoocLevel.begin(); p != aoocLevel.end(); p++) {
                guess_occ[*p -1] = 1.0;
            }
        }
        break;

    case RUN_UKS_BETA:
        {
            std::vector<int> boocLevel = this->getLevel(pdfParam["method/sp/beta-spin-occlevel"].getStr());
            for (std::vector<int>::const_iterator p = boocLevel.begin(); p != boocLevel.end(); p++) {
                guess_occ[*p -1] = 1.0;
            }
        }
        break;

    case RUN_ROKS:
        {
            std::vector<int> docLevel = this->getLevel(pdfParam["method/roks/closed-shell"].getStr());
            for (std::vector<int>::const_iterator p = docLevel.begin(); p != docLevel.end(); p++) {
                guess_occ[*p -1] = 2.0;
            }

            std::vector<int> socLevel = this->getLevel(pdfParam["method/roks/open-shell"].getStr());
            for (std::vector<int>::const_iterator p = socLevel.begin(); p != socLevel.end(); p++) {
                // nsoc
                guess_occ[*p -1] = 1.0;
            }
        }
        break;

    default:
        CnErr.abort();
        break;
    }

    // output occupation number to a files in fl_Work directory
    const std::string sOccFileName = this->getOccupationPath(runType);
    guess_occ.save(sOccFileName);
}


std::vector<int> DfInitialGuess::getLevel(std::string sLevel)
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


void DfInitialGuess::saveOccupation(const RUN_TYPE runType, const TlVector& rOccupation)
{
    const std::string sOccFileName = this->getOccupationPath(runType);
    rOccupation.save(sOccFileName);
}
