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

#include <fstream>
#include "CnError.h"
#include "DfInitialGuess.h"
#include "DfInitialGuessHuckel.h"
#include "DfInitialGuessHarris.h"
#include "DfDmatrix.h"
#include "DfPopulation.h"
#include "TlStringTokenizer.h"

DfInitialGuess::DfInitialGuess(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
    this->isNormalizeDensityMatrix_ = true;
    if ((*this->pPdfParam_)["guess/normalize_density_matrix"].getStr().empty() == false) {
        this->isNormalizeDensityMatrix_ = (*this->pPdfParam_)["guess/normalize_density_matrix"].getBoolean();
    }
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
        // this->createRho();
        // this->createOccupation();
        this->log_.critical("sorry. guess rho parameter is obsolete.");
        CnErr.abort();
        break;

    case GUESS_DENSITY:
        this->createInitialGuessUsingDensityMatrix();
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
        this->log_.warn("unknown initial guess parameter.");
        break;
    }
}

void DfInitialGuess::createInitialGuessUsingHuckel()
{
    DfInitialGuessHuckel huckel(this->pPdfParam_);
    huckel.createGuess();
}


void DfInitialGuess::createInitialGuessUsingCore()
{
    DfInitialGuessHuckel huckel(this->pPdfParam_);
    huckel.createGuess();
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
        this->createInitialGuessUsingLCAO(RUN_ROKS_CLOSED);
        this->createInitialGuessUsingLCAO(RUN_ROKS_OPEN);
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
        tmpParam["orbital-overlap-correspondence-method"] = "simple";
        tmpParam["num_of_iterations"] = 0;

        // 密度行列の作成
        DfDmatrix dfDmatrix(&tmpParam);
        dfDmatrix.DfDmatrixMain(); // RKS only?
    }
}


void DfInitialGuess::createInitialGuessUsingDensityMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->createInitialGuessUsingDensityMatrix(RUN_RKS);
        break;

    case METHOD_UKS:
        this->createInitialGuessUsingDensityMatrix(RUN_UKS_ALPHA);
        this->createInitialGuessUsingDensityMatrix(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->createInitialGuessUsingDensityMatrix(RUN_ROKS_CLOSED);
        this->createInitialGuessUsingDensityMatrix(RUN_ROKS_OPEN);
        break;

    default:
        abort();
        break;
    }

    this->createOccupation();
}


void DfInitialGuess::createInitialGuessUsingDensityMatrix(const RUN_TYPE runType)
{
    // read guess lcao
    TlSymmetricMatrix P = this->getInitialDensityMatrix<TlSymmetricMatrix>(runType);
    if (this->isNormalizeDensityMatrix_) {
        P = this->normalizeDensityMatrix<TlSymmetricMatrix, DfPopulation>(runType, P);
    }
    this->savePpqMatrix(runType, 0, P);

    // make occupation data
    this->createOccupation();
}


TlVector DfInitialGuess::getOccupation(const RUN_TYPE runType)
{
    TlVector occupation;
    const std::string binFile = TlUtils::format("./guess.occ.%s.vtr", this->m_sRunTypeSuffix[runType].c_str());
    const std::string txtFile = TlUtils::format("./guess.occ.%s.txt", this->m_sRunTypeSuffix[runType].c_str());

    if (TlFile::isExist(binFile)) {
        occupation.load(binFile);
    } else if (TlFile::isExist(txtFile)) {
        occupation.loadText(txtFile);
    } else {
        this->log_.warn(TlUtils::format("file not found.: %s", binFile.c_str()));
    }

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
        this->createOccupation(RUN_ROKS_CLOSED);
        this->createOccupation(RUN_ROKS_OPEN);
        break;

    default:
        CnErr.abort();
        break;
    }
}


TlVector DfInitialGuess::createOccupation(const RUN_TYPE runType)
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);

    // construct guess occupations
    const index_type numOfMOs = this->m_nNumOfMOs;
    TlVector guess_occ(numOfMOs);
    switch (runType) {
    case RUN_RKS:
        {
            const std::vector<int> docLevel = this->getLevel(pdfParam["method/rks/occlevel"].getStr());
            for (std::vector<int>::const_iterator p = docLevel.begin(); p != docLevel.end(); ++p) {
                const int level = *p -1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ[*p -1] = 2.0;
                }
            }
        }
        break;

    case RUN_UKS_ALPHA:
        {
            const std::vector<int> occLevel = this->getLevel(pdfParam["method/uks/alpha_occlevel"].getStr());
            for (std::vector<int>::const_iterator p = occLevel.begin(); p != occLevel.end(); ++p) {
                const int level = *p -1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ[level] = 1.0;
                }
            }
        }
        break;

    case RUN_UKS_BETA:
        {
            const std::vector<int> occLevel = this->getLevel(pdfParam["method/uks/beta_occlevel"].getStr());
            for (std::vector<int>::const_iterator p = occLevel.begin(); p != occLevel.end(); ++p) {
                const int level = *p -1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ[*p -1] = 1.0;
                }
            }
        }
        break;

    case RUN_ROKS_CLOSED:
        {
            const std::vector<int> occLevel_c = this->getLevel(pdfParam["method/roks/closed_occlevel"].getStr());
            for (std::vector<int>::const_iterator p = occLevel_c.begin(); p != occLevel_c.end(); ++p) {
                const int level = *p -1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ[*p -1] = 2.0;
                }
            }
        }
        break;

    case RUN_ROKS_OPEN:
        {
            const std::vector<int> occLevel_o = this->getLevel(pdfParam["method/roks/open_occlevel"].getStr());
            for (std::vector<int>::const_iterator p = occLevel_o.begin(); p != occLevel_o.end(); ++p) {
                const int level = *p -1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ[*p -1] = 1.0;
                }
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

    return guess_occ;
}

std::vector<int> DfInitialGuess::getLevel(const std::string& inputStr)
{
    std::vector<int> answer;

    // 構文解釈
    std::string numStr = "";
    std::vector<int> stack;
    const int len = inputStr.size();
    for (int i = 0; i < len; ++i) {
        const char c = inputStr[i];
        if (std::isdigit(c) != 0) {
            numStr.append(1, c);
        } else {
            if (numStr.size() > 0) {
                const int num = std::atoi(numStr.c_str());
                stack.push_back(num);
                numStr = "";
            }

            if (c == '-') {
                stack.push_back(-1);
            }
        }
    }
    if (numStr.empty() == false) {
        stack.push_back(std::atoi(numStr.c_str()));
    }

    // 翻訳
    const int stackSize = stack.size();
    for (int i = 0; i < stackSize; ++i) {
        const int v = stack[i];
        if (v > 0) {
            answer.push_back(v);
        } else if (v == -1) {
            const int i1 = i + 1;
            if ((i1 < stackSize) && (answer.size() > 0)) {
                const int end = stack[i1];
                const int start = answer.at(answer.size() -1);
                if (start > 0) { 
                    for (int v = start +1; v <= end; ++v) {
                        answer.push_back(v);
                    }
                }
                ++i;
            }
        }
    }

    return answer;
}


void DfInitialGuess::saveOccupation(const RUN_TYPE runType, const TlVector& rOccupation)
{
    const std::string sOccFileName = this->getOccupationPath(runType);
    rOccupation.save(sOccFileName);
}


void DfInitialGuess::makeDensityMatrix()
{
    TlSerializeData tmpParam = *(this->pPdfParam_);
    tmpParam["orbital-correspondence"] = false;
    tmpParam["orbital-overlap-correspondence-method"] = "simple";
    tmpParam["control-iteration"] = 0;

    DfDmatrix* pDfDmat = getDfDmatrixObject(&tmpParam);
    pDfDmat->DfDmatrixMain();
    delete pDfDmat;
    pDfDmat = NULL;
}

DfDmatrix* DfInitialGuess::getDfDmatrixObject(TlSerializeData* param)
{
    DfDmatrix* obj = new DfDmatrix(param);
    return obj;
}
