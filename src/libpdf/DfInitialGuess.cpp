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

#include "DfInitialGuess.h"

#include <fstream>

#include "CnError.h"
#include "DfDmatrix.h"
#include "DfInitialGuessHarris.h"
#include "DfInitialGuessHuckel.h"
#include "TlStringTokenizer.h"

#ifdef HAVE_EIGEN
#include "df_population_eigen.h"
#endif  // HAVE_EIGEN

#ifdef HAVE_LAPACK
#include "df_population_lapack.h"
#endif  // HAVE_LAPACK

DfInitialGuess::DfInitialGuess(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
    this->isNormalizeDensityMatrix_ = true;
    if ((*this->pPdfParam_)["guess/normalize_density_matrix"]
            .getStr()
            .empty() == false) {
        this->isNormalizeDensityMatrix_ =
            (*this->pPdfParam_)["guess/normalize_density_matrix"].getBoolean();
    }
}

DfInitialGuess::~DfInitialGuess() {}

void DfInitialGuess::exec() {
    unsigned int calcState = this->loadCalcState();

    if ((calcState & GUESS_DONE) == 0) {
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

        calcState |= GUESS_DONE;
        this->saveCalcState(calcState);
    }
}

unsigned int DfInitialGuess::loadCalcState() const {
    unsigned int cs = GUESS_UNPROCESSED;

    const bool isRestart = (*this->pPdfParam_)["restart"].getBoolean();
    if (isRestart) {
        this->logger("restart calculation is enabled.");
        cs = (*this->pPdfParam_)["control"]["guess_state"].getUInt();
    }

    return cs;
}

void DfInitialGuess::saveCalcState(unsigned int cs) {
    (*this->pPdfParam_)["control"]["guess_state"] = cs;

    // save PDF parameter
    const std::string pdfParamPath = (*this->pPdfParam_)["pdf_param_path"].getStr();
    TlMsgPack pdfParam_mpac(*this->pPdfParam_);
    pdfParam_mpac.save(pdfParamPath);
}

void DfInitialGuess::createInitialGuessUsingHuckel() {
    DfInitialGuessHuckel huckel(this->pPdfParam_);
    huckel.createGuess();
}

void DfInitialGuess::createInitialGuessUsingCore() {
    DfInitialGuessHuckel huckel(this->pPdfParam_);
    huckel.createGuess();
}

void DfInitialGuess::createInitialGuessUsingHarris() {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            DfInitialGuessHarris harris(this->pPdfParam_);
            harris.main();
        } break;

        case METHOD_UKS:
            CnErr.abort(
                "Sorry. harris method is not supported except RKS. stop.\n");
            break;

        case METHOD_ROKS:
            CnErr.abort(
                "Sorry. harris method is not supported except RKS. stop.\n");
            break;

        default:
            CnErr.abort();
            break;
    }
}

// ----------------------------------------------------------------------------
// LCAO
// ----------------------------------------------------------------------------
void DfInitialGuess::createInitialGuessUsingLCAO() {
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

void DfInitialGuess::createInitialGuessUsingLCAO(const RUN_TYPE runType) {
    this->createLcaoByFile(runType);

    this->createOccupationByFile(runType);

    this->makeDensityMatrix();
}

std::string DfInitialGuess::getLcaoPath_txt(const RUN_TYPE runType) {
    return TlUtils::format("./guess.lcao.%s.txt", DfObject::m_sRunTypeSuffix[runType].c_str());
}

std::string DfInitialGuess::getLcaoPath_bin(const RUN_TYPE runType) {
    return TlUtils::format("./guess.lcao.%s.mat", DfObject::m_sRunTypeSuffix[runType].c_str());
}

void DfInitialGuess::createLcaoByFile(const RUN_TYPE runType) {
    // MatrixType lcao;
    const std::string binFile = DfInitialGuess::getLcaoPath_bin(runType);
    const std::string txtFile = DfInitialGuess::getLcaoPath_txt(runType);

    if (TlFile::isExistFile(binFile)) {
        this->createLcaoByBinFile(runType);
    } else if (TlFile::isExistFile(txtFile)) {
        this->createLcaoByTxtFile(runType);
    } else {
        this->log_.warn(TlUtils::format("file not found.: %s", binFile.c_str()));
    }
}

void DfInitialGuess::createLcaoByBinFile(const RUN_TYPE runType) {
    const std::string binFile = DfInitialGuess::getLcaoPath_bin(runType);
    this->log_.info(TlUtils::format("check LCAO file(bin): %s", binFile.c_str()));

    TlMatrixObject::HeaderInfo headerInfo;
    const bool isLoadable = TlMatrixUtils::getHeaderInfo(binFile, &headerInfo);
    if (isLoadable == true) {
        std::string matrixType = "";
        switch (headerInfo.matrixType) {
            case TlMatrixObject::RLHD:
                matrixType = "symmetric";
                break;

            case TlMatrixObject::CSFD:
                matrixType = "normal (column-major)";
                break;

            case TlMatrixObject::RSFD:
                matrixType = "normal (row-major)";
                break;

            default:
                this->log_.critical(TlUtils::format("unknown matrix type: %s", binFile.c_str()));
        }
        this->log_.info(TlUtils::format("size: %d x %d", headerInfo.numOfRows, headerInfo.numOfCols));
    } else {
        CnErr.abort(TlUtils::format("cannot open file: %s", binFile.c_str()));
    }

    const std::string path = DfObject::getCMatrixPath(runType, 0);
    this->log_.info(TlUtils::format("copy the LCAO matrix: %s -> %s", binFile.c_str(), path.c_str()));
    TlFile::copy(binFile, path);
}

void DfInitialGuess::createLcaoByTxtFile(const RUN_TYPE runType) {
    const std::string txtFile = DfInitialGuess::getLcaoPath_txt(runType);

    this->log_.info(TlUtils::format("LCAO: loading: %s", txtFile.c_str()));
    std::ifstream fi;
    fi.open(txtFile.c_str(), std::ios::in);
    if ((fi.rdstate() & std::ifstream::failbit) != 0) {
        CnErr.abort(TlUtils::format("cannot open file %s.", txtFile.c_str()));
    }

    std::string dummy_line;
    fi >> dummy_line;

    int row_dimension, col_dimension;
    fi >> row_dimension >> col_dimension;
    if (row_dimension != this->m_nNumOfAOs) {
        CnErr.abort("DfInitialGuess", "", "prepare_occupation_and_or_mo",
                    "inputted guess lcao has illegal dimension");
    }

    const std::string path = DfObject::getCMatrixPath(runType, 0);
    TlDenseGeneralMatrix_mmap lcao(path);
    lcao.resize(row_dimension, col_dimension);

    const int maxRows = row_dimension;
    const int maxCols = col_dimension;
    for (int i = 0; i < maxRows; ++i) {
        for (int j = 0; j < maxCols; ++j) {
            double v;
            fi >> v;
            lcao.set(i, j, v);
        }
    }
}

// ----------------------------------------------------------------------------
// density matrix
// ----------------------------------------------------------------------------
void DfInitialGuess::createInitialGuessUsingDensityMatrix() {
    switch (this->linearAlgebraPackage_) {
        case LAP_VIENNACL:
            // use LAP_EIGEN
        case LAP_EIGEN: {
#ifdef HAVE_EIGEN
            {
                this->createInitialGuessUsingDensityMatrix_tmpl<TlDenseSymmetricMatrix_Eigen, DfPopulation_Eigen>();
            }
#else
            {
                CnErr.abort("linear algebra package mismatch.");
            }
#endif  // HAVE_EIGEN
        } break;

        case LAP_LAPACK: {
#ifdef HAVE_LAPACK
            {
                this->createInitialGuessUsingDensityMatrix_tmpl<TlDenseSymmetricMatrix_Lapack, DfPopulation_Lapack>();
            }
#else
            {
                CnErr.abort("linear algebra package mismatch.");
            }
#endif  // HAVE_LAPACK
        } break;

        default:
            CnErr.abort("linear algebra package mismatch.");
    }
}

// void DfInitialGuess::createInitialGuessUsingDensityMatrix(const RUN_TYPE runType) {
//     // read guess lcao
//     TlDenseSymmetricMatrix_Lapack P = this->getInitialDensityMatrix<TlDenseSymmetricMatrix_Lapack>(runType);
//     if (this->isNormalizeDensityMatrix_) {
//         P = this->normalizeDensityMatrix<TlDenseSymmetricMatrix_Lapack,
//                                          DfPopulation_Lapack>(runType, P);
//     }
//     this->savePpqMatrix(runType, 0, P);

//     // spin density matrix
//     const int iteration = 0;
//     if (runType == METHOD_RKS) {
//         this->saveSpinDensityMatrix(runType, iteration, 0.5 * P);
//     } else {
//         this->saveSpinDensityMatrix(runType, iteration, P);
//     }

//     // make occupation data
//     this->createOccupation();
// }

// ----------------------------------------------------------------------------
// occupation
// ----------------------------------------------------------------------------
void DfInitialGuess::createOccupationByFile(const RUN_TYPE runType) {
    TlDenseVector_Lapack occupation;
    const std::string binFile = TlUtils::format("./guess.occ.%s.vtr", this->m_sRunTypeSuffix[runType].c_str());
    const std::string txtFile = TlUtils::format("./guess.occ.%s.txt", this->m_sRunTypeSuffix[runType].c_str());

    if (TlFile::isExistFile(binFile)) {
        occupation.load(binFile);
    } else if (TlFile::isExistFile(txtFile)) {
        occupation.loadText(txtFile);
    } else {
        this->log_.warn(TlUtils::format("file not found.: %s", binFile.c_str()));
    }

    const std::string occFilePath = this->getOccupationPath(runType);
    occupation.save(occFilePath);
}

void DfInitialGuess::createOccupation() {
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

void DfInitialGuess::createOccupation(const RUN_TYPE runType) {
    const TlSerializeData& pdfParam = *(this->pPdfParam_);

    // construct guess occupations
    const index_type numOfMOs = this->m_nNumOfMOs;
    TlDenseVector_Lapack guess_occ(numOfMOs);
    switch (runType) {
        case RUN_RKS: {
            const std::vector<int> docLevel =
                this->getLevel(pdfParam["method/rks/occlevel"].getStr());
            for (std::vector<int>::const_iterator p = docLevel.begin();
                 p != docLevel.end(); ++p) {
                const int level = *p - 1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ.set(*p - 1, 2.0);
                }
            }
        } break;

        case RUN_UKS_ALPHA: {
            const std::vector<int> occLevel =
                this->getLevel(pdfParam["method/uks/alpha_occlevel"].getStr());
            for (std::vector<int>::const_iterator p = occLevel.begin();
                 p != occLevel.end(); ++p) {
                const int level = *p - 1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ.set(level, 1.0);
                }
            }
        } break;

        case RUN_UKS_BETA: {
            const std::vector<int> occLevel =
                this->getLevel(pdfParam["method/uks/beta_occlevel"].getStr());
            for (std::vector<int>::const_iterator p = occLevel.begin();
                 p != occLevel.end(); ++p) {
                const int level = *p - 1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ.set(*p - 1, 1.0);
                }
            }
        } break;

        case RUN_ROKS_CLOSED: {
            const std::vector<int> occLevel_c = this->getLevel(
                pdfParam["method/roks/closed_occlevel"].getStr());
            for (std::vector<int>::const_iterator p = occLevel_c.begin();
                 p != occLevel_c.end(); ++p) {
                const int level = *p - 1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ.set(*p - 1, 2.0);
                }
            }
        } break;

        case RUN_ROKS_OPEN: {
            const std::vector<int> occLevel_o =
                this->getLevel(pdfParam["method/roks/open_occlevel"].getStr());
            for (std::vector<int>::const_iterator p = occLevel_o.begin();
                 p != occLevel_o.end(); ++p) {
                const int level = *p - 1;
                if ((0 <= level) && (level < numOfMOs)) {
                    guess_occ.set(*p - 1, 1.0);
                }
            }
        } break;

        default:
            CnErr.abort();
            break;
    }

    // output occupation number to a files in fl_Work directory
    const std::string sOccFileName = this->getOccupationPath(runType);
    guess_occ.save(sOccFileName);
}

std::vector<int> DfInitialGuess::getLevel(const std::string& inputStr) {
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
                const int start = answer.at(answer.size() - 1);
                if (start > 0) {
                    for (int v = start + 1; v <= end; ++v) {
                        answer.push_back(v);
                    }
                }
                ++i;
            }
        }
    }

    return answer;
}

// ----------------------------------------------------------------------------
// make density matrix
// ----------------------------------------------------------------------------
void DfInitialGuess::makeDensityMatrix() {
    this->log_.info("make density matrix");

    TlSerializeData tmpParam = *(this->pPdfParam_);
    tmpParam["orbital-correspondence"] = false;
    tmpParam["orbital-overlap-correspondence-method"] = "simple";
    tmpParam["num_of_iterations"] = 0;

    DfDmatrix* pDfDmat = getDfDmatrixObject(&tmpParam);
    pDfDmat->run();

    delete pDfDmat;
    pDfDmat = NULL;

    this->log_.info("setup density matrix");
    this->copyDensityMatrix();
}

DfDmatrix* DfInitialGuess::getDfDmatrixObject(TlSerializeData* param) {
    DfDmatrix* obj = new DfDmatrix(param);
    return obj;
}

void DfInitialGuess::copyDensityMatrix() {
    switch (this->m_nMethodType) {
        case METHOD_RKS: {
            const std::string P0_path = DfObject::getPOutMatrixPath(RUN_RKS, 0);
            if (TlFile::isExistFile(P0_path)) {
                const std::string P1_path = DfObject::getPInMatrixPath(RUN_RKS, 1);
                this->log_.info("copy density matrix: " + P0_path + " -> " + P1_path + ".");
                TlFile::copy(P0_path, P1_path);
            } else {
                this->log_.warn("density matrix is not found: " + P0_path + ".");
            }
        } break;

        case METHOD_UKS: {
            {
                const std::string P0_path = DfObject::getPOutMatrixPath(RUN_UKS_ALPHA, 0);
                if (TlFile::isExistFile(P0_path)) {
                    const std::string P1_path = DfObject::getPInMatrixPath(RUN_UKS_ALPHA, 1);
                    TlFile::copy(P0_path, P1_path);
                }
            }
            {
                const std::string P0_path = DfObject::getPOutMatrixPath(RUN_UKS_BETA, 0);
                if (TlFile::isExistFile(P0_path)) {
                    const std::string P1_path = DfObject::getPInMatrixPath(RUN_UKS_BETA, 1);
                    TlFile::copy(P0_path, P1_path);
                }
            }
        } break;

        case METHOD_ROKS: {
            {
                const std::string P0_path = DfObject::getPOutMatrixPath(RUN_ROKS_CLOSED, 0);
                if (TlFile::isExistFile(P0_path)) {
                    const std::string P1_path = DfObject::getPInMatrixPath(RUN_ROKS_CLOSED, 1);
                    TlFile::copy(P0_path, P1_path);
                }
            }
            {
                const std::string P0_path = DfObject::getPOutMatrixPath(RUN_ROKS_OPEN, 0);
                if (TlFile::isExistFile(P0_path)) {
                    const std::string P1_path = DfObject::getPInMatrixPath(RUN_ROKS_OPEN, 1);
                    TlFile::copy(P0_path, P1_path);
                }
            }
        } break;

        default:
            CnErr.abort("program error.");
            break;
    }
}
