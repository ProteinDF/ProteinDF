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

#include <algorithm>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "CnError.h"
#include "CnFile_parallel.h"
#include "DfCD_Parallel.h"
#include "DfEriEngine.h"
#include "DfOverlapEngine.h"
#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"
#include "TlSystem.h"
#include "TlTime.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_mmap.h"
#include "tl_matrix_utils.h"

#ifdef HAVE_EIGEN
#include "tl_dense_symmetric_matrix_eigen.h"
#endif  // HAVE_EIGEN

#ifdef HAVE_LAPACK
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#endif  // HAVE_LAPACK

#ifdef HAVE_SCALAPACK
#include "tl_dense_general_matrix_scalapack.h"
#endif  // HAVE_SCALAPACK

#define TRANS_MEM_SIZE (1 * 1024 * 1024 * 1024)  // 1GB
// #define CD_DEBUG

DfCD_Parallel::DfCD_Parallel(TlSerializeData* pPdfParam)
    : DfCD(pPdfParam, false) {
    this->file_ = new CnFile_parallel();

    this->isDebugSaveL_ = false;
    if (!(*this->pPdfParam_)["debug/saveL"].getStr().empty()) {
        this->isDebugSaveL_ = (*this->pPdfParam_)["debug/saveL"].getBoolean();
    }
}

DfCD_Parallel::~DfCD_Parallel() {
    delete this->file_;
    this->file_ = NULL;
}

// ----------------------------------------------------------------------------
// calc L (for JK)
// ----------------------------------------------------------------------------
void DfCD_Parallel::calcCholeskyVectorsForJK() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc CholeskyVectors (parallel)");

    // for J & K
    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set"]);
    {
        // productive code
        this->log_.info("calc L(JK) start");
        this->createEngines<DfEriEngine>();

        switch (this->cdIntermediateFileFormat_) {
            case CD_INTERMEDIATE_FILE_FORMAT_MMAP: {
                this->log_.info("L_jk build on mmap");

                // prepare L
                const bool isTmp = (this->localTempPath_.empty() != true);
                std::string L_mat_path;
                TlDenseGeneralMatrix_mmap* pL = NULL;
                if (rComm.isMaster()) {
                    L_mat_path = DfObject::getLjkMatrixPath(this->localTempPath_, isTmp);

                    this->log_.info(TlUtils::format("L saved as %s", L_mat_path.c_str()));
                    if (TlFile::isExistFile(L_mat_path)) {
                        TlFile::remove(L_mat_path);
                    }
                    pL = new TlDenseGeneralMatrix_mmap(L_mat_path, 1, 1);
                }

                // calc CD
                this->calcCholeskyVectorsOnTheFlyS(orbInfo, this->getI2pqVtrPath(), this->epsilon_,
                                                   &DfCD::calcDiagonals, &DfCD::getSuperMatrixElements, pL);

                // move L
                if (rComm.isMaster()) {
                    delete pL;
                    pL = NULL;

                    if (isTmp) {
                        this->log_.info(TlUtils::format("L move to %s", DfObject::getLjkMatrixPath().c_str()));
                        TlFile::move(L_mat_path, DfObject::getLjkMatrixPath());
                    }
                }
            } break;

            case CD_INTERMEDIATE_FILE_FORMAT_ARRAY: {
                this->log_.info("L_jk build on arrays");

                TlDenseGeneralMatrix_arrays_RowOriented Ljk(rComm.getNumOfProcs(), 1, rComm.getNumOfProcs(), rComm.getRank());
                this->calcCholeskyVectorsOnTheFlyS_new(orbInfo, this->getI2pqVtrPath(), this->epsilon_, &DfCD::calcDiagonals,
                                                       &DfCD_Parallel::getSuperMatrixElements, &Ljk);

                if (this->optCdFile_) {
                    this->log_.info("optimize L matrix.");

                    this->saveL(Ljk, DfObject::getLjkMatrixPath());
                } else {
                    this->log_.info("skip optimize L matrix");
                    const int subunitID = rComm.getRank();
                    const std::string path = TlUtils::format("%s.part%d.mat", DfObject::getLjkMatrixPath().c_str(), subunitID);
                    Ljk.save(path);
                    this->log_.info("save partial L matrix");
                }

                // if (rComm.isMaster()) {
                //     if (!this->localTempPath_.empty()) {
                //         this->log_.info(TlUtils::format("L move to %s", DfObject::getLjkMatrixPath().c_str()));
                //         TlFile::move(L_mat_path, DfObject::getLjkMatrixPath());
                //     }
                // }
            } break;

            case CD_INTERMEDIATE_FILE_FORMAT_ARRAY_MMAP: {
                this->log_.info("L_jk build on arrays with mmap");

                const std::string L_basePath = DfObject::getLjkMatrixPath();

                // prepare L
                const bool isTmp = (this->localTempPath_.empty() != true);
                const std::string L_basePath_tmp = DfObject::getLjkMatrixPath(this->localTempPath_, isTmp);
                this->log_.info(TlUtils::format("L base name: %s", L_basePath_tmp.c_str()));

                std::string L_savedPath = "";
                {
                    TlDenseGeneralMatrix_arrays_mmap_RowOriented L(L_basePath_tmp, 1, 1,
                                                                   rComm.getNumOfProcs(), rComm.getRank());
                    // calc CD
                    this->calcCholeskyVectorsOnTheFlyS(orbInfo, this->getI2pqVtrPath(), this->epsilon_,
                                                       &DfCD::calcDiagonals, &DfCD_Parallel::getSuperMatrixElements, &L);
                    L_savedPath = L.getFilePath();
                }

                if (isTmp) {
                    const int subunitId = rComm.getRank();
                    const std::string new_L_savedPath = TlDenseGeneralMatrix_arrays_mmap_RowOriented::getFileName(L_basePath, subunitId);
                    this->log_.info(TlUtils::format("L move to %s", new_L_savedPath.c_str()));
                    TlFile::move(L_savedPath, new_L_savedPath);
                }

                // finalize
                rComm.barrier();
                this->log_.info("L_jk saved.");

                if (this->optCdFile_) {
                    this->log_.info("optimize L matrix.");
                    this->transpose2CSFD_mpi(L_basePath, DfObject::getLjkMatrixPath());
                } else {
                    this->log_.info("skip optimize L matrix");
                }
            } break;

            default: {
                this->log_.critical(TlUtils::format("program error: %s,%d", __FILE__, __LINE__));
            }
        }

        this->destroyEngines();
        this->log_.info("calc L(JK) end");
    }
}

bool DfCD_Parallel::transpose2CSFD_mpi(const std::string& rvmBasePath, const std::string& outputMatrixPath,
                                       const bool verbose, const bool showProgress) {
    bool answer = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int rank = rComm.getRank();

    // 初期データ読み込み
    TlMatrixObject::index_type numOfRows = 0;
    TlMatrixObject::index_type numOfCols = 0;
    int numOfSubunits = 0;
    int sizeOfChunk = 0;
    // bool isLoadable = false;
    {
        int subunitID = rank;
        const std::string inputPath0 = TlDenseMatrix_arrays_mmap_Object::getFileName(rvmBasePath, subunitID);

        TlMatrixObject::HeaderInfo headerInfo;
        const bool isLoadable = TlMatrixUtils::getHeaderInfo(inputPath0, &headerInfo);

        if (isLoadable != true) {
            std::cerr << TlUtils::format("can not open file: %s@%d", inputPath0.c_str(), rank) << std::endl;
            return false;
        }

        numOfRows = headerInfo.numOfVectors;
        numOfCols = headerInfo.sizeOfVector;
        numOfSubunits = headerInfo.numOfSubunits;
        sizeOfChunk = headerInfo.sizeOfChunk;

        if (verbose) {
            std::cerr << "rows: " << numOfRows << std::endl;
            std::cerr << "cols: " << numOfCols << std::endl;
            std::cerr << "units: " << numOfSubunits << std::endl;
            std::cerr << "chunk: " << sizeOfChunk << std::endl;
        }
    }

    // prepare CSFD file
    if (rComm.isMaster()) {
        if (TlFile::isExistFile(outputMatrixPath)) {
            if (verbose) {
                std::cerr << "file overwrite: " << outputMatrixPath << std::endl;
            }
            TlFile::remove(outputMatrixPath);
        }
        TlDenseGeneralMatrix_mmap outMat(outputMatrixPath, numOfRows, numOfCols);
        if (verbose) {
            std::cerr << "output matrix has been prepared by mmap." << std::endl;
        }
    }

    // convert
    // {
    //     const std::string L_path = TlDenseMatrix_arrays_mmap_Object::getFileName(rvmBasePath, rank);
    //     TlDenseGeneralMatrix_arrays_mmap_RowOriented L(L_path);

    //     const std::string tempCsfdMatPath = TlUtils::format("L_csfd.%d.mat", rank);
    //     L.convertMemoryLayout(tempCsfdMatPath);

    //     // @todo remove unused func `convert2csfd`
    //     // convert2csfd(rvmBasePath, rank, tempMatPath, verbose, showProgress);
    // }

    // @todo transfer tempCsfdMat to master note

    // 最終書き込み
    rComm.barrier();
    if (rComm.isMaster()) {
        TlDenseGeneralMatrix_mmap outMat(outputMatrixPath);

        for (int unit = 0; unit < numOfSubunits; ++unit) {
            this->log_.info(TlUtils::format("optimizing unit=%d", unit));
            const std::string inputPath = TlDenseMatrix_arrays_mmap_Object::getFileName(rvmBasePath, unit);
            TlDenseGeneralMatrix_arrays_mmap_RowOriented inMat(inputPath);

            const std::string tempCsfdMatPath = TlUtils::format("L_csfd.%d.mat", unit);
            inMat.convertMemoryLayout(tempCsfdMatPath, verbose, showProgress);
            inMat.set2csfd(&outMat, verbose, showProgress);

            TlFile::remove(tempCsfdMatPath);
        }
    }

    return answer;
}

// ----------------------------------------------------------------------------
// calc L (for K only)
// ----------------------------------------------------------------------------
void DfCD_Parallel::calcCholeskyVectorsForK() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc CholeskyVectors (K) (parallel)");

    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set"]);
    switch (this->fastCDK_mode_) {
        case FASTCDK_PRODUCTIVE_FULL: {
            this->log_.info("fast CDK routine(parallel; full).");
            this->createEngines<DfEriEngine>();

            TlDenseGeneralMatrix_arrays_RowOriented Lk(rComm.getNumOfProcs(), 1, rComm.getNumOfProcs(), rComm.getRank());
            this->calcCholeskyVectorsOnTheFlyS_new(orbInfo, this->getI2prVtrPath(), this->epsilon_K_, &DfCD::calcDiagonals_K_full,
                                                   &DfCD_Parallel::getSuperMatrixElements_K_full, &Lk);
            this->saveL(Lk, DfObject::getLkMatrixPath());

            this->destroyEngines();
            this->log_.info("");
        } break;

        case FASTCDK_PRODUCTIVE: {
            this->log_.info("fast CDK routine(parallel).");
            this->createEngines<DfEriEngine>();

            TlDenseGeneralMatrix_arrays_RowOriented Lk(rComm.getNumOfProcs(), 1, rComm.getNumOfProcs(), rComm.getRank());
            this->calcCholeskyVectorsOnTheFlyS_new(
                orbInfo, this->getI2prVtrPath(), this->epsilon_K_, &DfCD::calcDiagonals_K_half,
                &DfCD_Parallel::getSuperMatrixElements_K_half, &Lk);
            this->saveL(Lk, DfObject::getLkMatrixPath());

            this->destroyEngines();
            this->log_.info("");
        } break;

        case FASTCDK_NONE:
            this->log_.info("fast CDK routine is invalided.");
            break;

        default:
            this->log_.critical("program error");
            CnErr.abort();
            break;
    }
}

// ----------------------------------------------------------------------------
// calc L (for gridfree)
// ----------------------------------------------------------------------------
void DfCD_Parallel::calcCholeskyVectorsForGridFree() {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // for XC(gridfree)
    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set"]);

    if (this->isDedicatedBasisForGridFree_) {
        const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set_gridfree"]);
        // productive code
        const TlDenseGeneralMatrix_arrays_RowOriented Lxc =
            DfCD::calcCholeskyVectorsOnTheFly<DfOverlapEngine>(orbInfo_p, orbInfo_q, this->getI2pqVtrXCPath());
        this->saveL(Lxc, DfObject::getLxcMatrixPath());
    } else {
        // productive code
        this->log_.info("build Lxc matrix by on-the-fly method.");
        this->createEngines<DfOverlapEngine>();

        TlDenseGeneralMatrix_arrays_RowOriented Lxc(rComm.getNumOfProcs(), 1, rComm.getNumOfProcs(), rComm.getRank());
        this->calcCholeskyVectorsOnTheFlyS_new(orbInfo_p, this->getI2pqVtrXCPath(), this->epsilon_,
                                               &DfCD::calcDiagonals, &DfCD_Parallel::getSuperMatrixElements, &Lxc);

        this->saveL(Lxc, DfObject::getLxcMatrixPath());
        this->destroyEngines();
    }
}

DfTaskCtrl* DfCD_Parallel::getDfTaskCtrlObject() const {
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl_Parallel(this->pPdfParam_);
    return pDfTaskCtrl;
}

void DfCD_Parallel::finalize(TlDenseSymmetricMatrix_Lapack* pMat) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(pMat);
}

void DfCD_Parallel::finalize(TlSparseSymmetricMatrix* pMat) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    // rComm.gatherToMaster(*pMat);
    // rComm.broadcast(*pMat);
    rComm.allReduce_SUM(*pMat);
}

void DfCD_Parallel::finalize_I2PQ(PQ_PairArray* pI2PQ) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();

    // gather to master
    if (rComm.isMaster() == true) {
        const int slaves = numOfProcs - 1;
        int proc = 0;
        std::vector<index_type> shellArray;
        for (int i = 0; i < slaves; ++i) {
            rComm.receiveDataFromAnySource(shellArray, &proc);

            const std::size_t I2PQ_size = shellArray.size() / 2;
            PQ_PairArray i2pq_tmp(I2PQ_size);
            for (std::size_t i = 0; i < I2PQ_size; ++i) {
                i2pq_tmp[i] = Index2(shellArray[i * 2], shellArray[i * 2 + 1]);
            }

            pI2PQ->insert(pI2PQ->end(), i2pq_tmp.begin(), i2pq_tmp.end());
        }
        // std::sort(pI2PQ->begin(), pI2PQ->end(), PQ_Pair_less());
        std::sort(pI2PQ->begin(), pI2PQ->end());
    } else {
        const std::size_t I2PQ_size = pI2PQ->size();
        std::vector<index_type> shellArray(I2PQ_size * 2);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            shellArray[i * 2] = (*pI2PQ)[i].index1();
            shellArray[i * 2 + 1] = (*pI2PQ)[i].index2();
        }
        rComm.sendData(shellArray);
    }

    // broadcast
    {
        std::vector<index_type> tmp;
        if (rComm.isMaster() == true) {
            const std::size_t size = pI2PQ->size();
            tmp.resize(size * 2);
            for (std::size_t i = 0; i < size; ++i) {
                tmp[i * 2] = (*pI2PQ)[i].index1();
                tmp[i * 2 + 1] = (*pI2PQ)[i].index2();
            }
        }
        rComm.broadcast(tmp);
        if (rComm.isMaster() != true) {
            const std::size_t size = tmp.size() / 2;
            pI2PQ->resize(size);
            for (std::size_t i = 0; i < size; ++i) {
                (*pI2PQ)[i] = Index2(tmp[i * 2], tmp[i * 2 + 1]);
            }
        }
    }
}

void DfCD_Parallel::saveI2PQ(const PQ_PairArray& I2PQ, const std::string& filepath) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfCD::saveI2PQ(I2PQ, filepath);
    }
}

DfCD::PQ_PairArray DfCD_Parallel::getI2PQ(const std::string& filepath) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("distribute I2PQ table.");

    PQ_PairArray I2PQ;
    std::vector<index_type> shellArray;
    if (rComm.isMaster() == true) {
        I2PQ = DfCD::getI2PQ(filepath);

        const std::size_t I2PQ_size = I2PQ.size();
        shellArray.resize(I2PQ_size * 2);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            shellArray[i * 2] = I2PQ[i].index1();
            shellArray[i * 2 + 1] = I2PQ[i].index2();
        }
    }

    rComm.broadcast(shellArray);
    if (rComm.isMaster() != true) {
        const std::size_t I2PQ_size = shellArray.size() / 2;
        I2PQ.resize(I2PQ_size);
        for (std::size_t i = 0; i < I2PQ_size; ++i) {
            I2PQ[i] = Index2(shellArray[i * 2], shellArray[i * 2 + 1]);
        }
    }

    return I2PQ;
}

void DfCD_Parallel::divideCholeskyBasis(const index_type numOfCVs, index_type* pStart, index_type* pEnd) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const index_type range = (numOfCVs + numOfProcs - 1) / numOfProcs;
    *pStart = range * rank;
    *pEnd = std::min(range * (rank + 1), numOfCVs);
}

TlDenseSymmetricMatrix_Scalapack DfCD_Parallel::getCholeskyVector_distribute(const TlDenseVector_Lapack& L_col,
                                                                             const PQ_PairArray& I2PQ) {
    const index_type numOfItilde = L_col.getSize();
    TlDenseSymmetricMatrix_Scalapack answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(), I2PQ[i].index2(), L_col.get(i));
    }

    return answer;
}

TlDenseSymmetricMatrix_Scalapack DfCD_Parallel::getCholeskyVectorA_distribute(const TlOrbitalInfoObject& orbInfo_p,
                                                                              const TlOrbitalInfoObject& orbInfo_q,
                                                                              const TlDenseVector_Lapack& L_col,
                                                                              const PQ_PairArray& I2PQ) {
    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();

    const index_type numOfItilde = L_col.getSize();
    TlDenseGeneralMatrix_Scalapack answer(numOfOrbs_p, numOfOrbs_q);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(), I2PQ[i].index2(), L_col.get(i));
    }

    return answer;
}

// ----------------------------------------------------------------------------
// [integral] calc CD
// ----------------------------------------------------------------------------
void DfCD_Parallel::calcCholeskyVectorsOnTheFlyS_new(
    const TlOrbitalInfoObject& orbInfo, const std::string& I2PQ_path, const double threshold,
    CalcDiagonalsFunc calcDiagonalsFunc, GetSuperMatrixElementsFuncP getSuperMatrixElementsFunc,
    TlDenseGeneralMatrix_arrays_RowOriented* pL) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int myRank = rComm.getRank();

    this->log_.info("call on-the-fly Cholesky Decomposition routine (parallel; symmetric) (new)");
    assert(this->pEngines_ != NULL);

    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfAOs * (numOfAOs + 1) / 2;
    this->log_.info(TlUtils::format("number of orbitals: %d", numOfAOs));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));

    // CDAM
    PQ_PairArray I2PQ;
    std::vector<double> diagonals;  // 対角成分
    (this->*calcDiagonalsFunc)(orbInfo, &I2PQ, &diagonals);
    assert(diagonals.size() == I2PQ.size());

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);

    // debug
    // this->debug_I2PQ_ = I2PQ;

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", threshold));
    const index_type numOfPQtilde = I2PQ.size();
    // TlDenseGeneralMatrix_arrays_RowOriented L(numOfPQtilde, 1, rComm.getNumOfProcs(), myRank);

    double error = std::accumulate(diagonals.begin(), diagonals.end(), 0.0);
    std::vector<std::size_t> pivot(numOfPQtilde);
    for (index_type i = 0; i < numOfPQtilde; ++i) {
        pivot[i] = i;
    }

    std::vector<double> errors(numOfPQtilde);

    int progress = 0;
    index_type division = std::max<index_type>(numOfPQtilde * 0.01, 100);

    // ------------------------------------------------------------
    index_type reserveLCols = this->m_nNumOfAOs * 5;
    this->log_.info(TlUtils::format("resize L col: %d", reserveLCols));
    // pL->reserveColSize(division);
    pL->resize(numOfPQtilde, reserveLCols);
    // ------------------------------------------------------------

    index_type numOfCDVcts = 0;
    while ((error > threshold) && (numOfCDVcts < numOfPQtilde)) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", numOfCDVcts, numOfPQtilde, error));
#endif  // DEBUG_CD

        // progress
        if (numOfCDVcts >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e [RSS=% 8.3e MB]", numOfCDVcts, error, TlSystem::getMaxRSS()));
            ++progress;
        }

        // memory allocation
        if (numOfCDVcts >= reserveLCols) {
            reserveLCols += this->m_nNumOfAOs + this->m_nNumOfAOs;
            this->log_.info(TlUtils::format("resize L col: %d", reserveLCols));
            // pL->reserveColSize(division * progress);
            pL->resize(numOfPQtilde, reserveLCols);
        }

        // pivot
        {
            const std::size_t argmax = this->argmax_pivot(diagonals, pivot, numOfCDVcts);
            std::swap(pivot[numOfCDVcts], pivot[argmax]);
        }

        const index_type pivot_m = pivot[numOfCDVcts];
        error = diagonals[pivot_m];
        errors[numOfCDVcts] = error;
        if (error < threshold) {
            break;
        }

        // get super-matrix elements
        std::vector<double> G_pm_local;
        const index_type numOf_G_cols = numOfPQtilde - (numOfCDVcts + 1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type c = 0; c < numOf_G_cols; ++c) {
                const index_type pivot_i = pivot[(numOfCDVcts + 1) + c];  // from (m+1) to N
                G_col_list[c] = pivot_i;
            }
            (this->*getSuperMatrixElementsFunc)(orbInfo, pivot_m, G_col_list, I2PQ, &G_pm_local);
        }
        assert(static_cast<index_type>(G_pm_local.size()) == numOf_G_cols);

        // allReduce_SUM (G_pm)
        std::vector<double> G_pm(numOf_G_cols, 0.0);
        rComm.iAllReduce_SUM(&(G_pm_local[0]), &(G_pm[0]), numOf_G_cols);

        // CD calc
        // output:
        //   out_L_rows: row elements at the column numOfCDVcts(target) in L

        // TlDenseVector_Lapack L_pm(numOfCDVcts + 1);
        std::valarray<double> L_pm_x(0.0, numOfCDVcts + 1);
        {
            const int PE_in_charge = pL->getSubunitID(pivot_m);
            if (PE_in_charge == myRank) {
                const index_type copyCount_m = pL->getRowVector(pivot_m, &(L_pm_x[0]), numOfCDVcts + 1);
                assert(copyCount_m == numOfCDVcts + 1);
            }
            rComm.broadcast(L_pm_x, PE_in_charge);
        }

        std::valarray<double> out_L_rows_local(0.0, numOfPQtilde);
        const double l_m_pm = std::sqrt(diagonals[pivot_m]);
        const double inv_l_m_pm = 1.0 / l_m_pm;
        // std::vector<double> L_xm(numOf_G_cols, 0.0);
        // pL->set(pivot_m, numOfCDVcts, l_m_pm);
        if (rComm.isMaster()) {
            out_L_rows_local[pivot_m] = l_m_pm;
        }

        std::vector<double> update_diagonals_local(numOfPQtilde, 0.0);

        // wait (output_G_pm)
        rComm.wait(&(G_pm[0]));

#pragma omp parallel
        {
            std::valarray<double> L_pi_x(0.0, numOfCDVcts + 1);
#pragma omp for schedule(runtime)
            for (index_type i = 0; i < numOf_G_cols; ++i) {
                const index_type pivot_i = pivot[(numOfCDVcts + 1) + i];  // from (m+1) to N

                if (pL->getSubunitID(pivot_i) == myRank) {
                    const index_type copyCount_i = pL->getRowVector(pivot_i, &(L_pi_x[0]), numOfCDVcts + 1);
                    assert(copyCount_i == numOfCDVcts + 1);

                    // const double sum_ll = (L_pi.dotInPlace(L_pm)).sum();
                    // const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;
                    const double sum_ll = (L_pm_x * L_pi_x).sum();
                    const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

                    // #pragma omp atomic
                    //                     L_xm[i] += l_m_pi;
                    out_L_rows_local[pivot_i] = l_m_pi;

                    // #pragma omp atomic
                    update_diagonals_local[pivot_i] -= l_m_pi * l_m_pi;
                }
            }
        }
        // pL->setColVector(numOfCDVcts, out_L_rows);

        // allReduce_SUM (output_L_rows)
        // std::vector<double> output_L_xm(numOf_G_cols, 0.0);
        std::valarray<double> out_L_rows(0.0, numOfPQtilde);
        rComm.iAllReduce_SUM(&(out_L_rows_local[0]), &(out_L_rows[0]), numOfPQtilde);

        // allReduce_SUM (update_diagonals)
        std::vector<double> update_diagonals(numOfPQtilde, 0.0);
        rComm.iAllReduce_SUM(&(update_diagonals_local[0]), &(update_diagonals[0]), numOfPQtilde);

        rComm.wait(&(out_L_rows[0]));
        // for (index_type i = 0; i < numOf_G_cols; ++i) {
        //     const index_type pivot_i = pivot[(numOfCDVcts + 1) + i];  // from (m+1) to N
        //     pL->set(pivot_i, numOfCDVcts, output_L_xm[i]);
        // }
        pL->setColVector(numOfCDVcts, out_L_rows);

        rComm.wait(&(update_diagonals[0]));
        // diagonals += output_update_diagonals;
        std::transform(diagonals.begin(), diagonals.end(), update_diagonals.begin(), diagonals.begin(), std::plus<double>());

        // error = diagonals[pivot[numOfCDVcts]];

        ++numOfCDVcts;
    }

    {
        this->log_.info("save errors");
        if (rComm.isMaster()) {
            errors.resize(numOfCDVcts);
            TlDenseVector_Lapack e(errors);
            DfObject::saveLjkErrorsVector(e);
        }
    }

    this->log_.info(TlUtils::format("Cholesky Vectors: %d", numOfCDVcts));
    pL->resize(numOfPQtilde, numOfCDVcts);

    // return L;
}

void DfCD_Parallel::calcCholeskyVectorsOnTheFlyS(const TlOrbitalInfoObject& orbInfo, const std::string& I2PQ_path,
                                                 const double threshold, CalcDiagonalsFunc calcDiagonalsFunc,
                                                 GetSuperMatrixElementsFuncP getSuperMatrixElementsFunc,
                                                 TlDenseGeneralMatrix_arrays_mmap_RowOriented* pL) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int myRank = rComm.getRank();

    this->log_.info("call on-the-fly Cholesky Decomposition routine (parallel; symmetric)");
    this->log_.info("pass in: DfCD_Parallel::calcCholeskyVectorsOnTheFlyS()");
    assert(this->pEngines_ != NULL);

    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfAOs * (numOfAOs + 1) / 2;
    this->log_.info(TlUtils::format("number of orbitals: %d", numOfAOs));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));

    // CDAM
    PQ_PairArray I2PQ;
    std::vector<double> diagonals;  // 対角成分
    (this->*calcDiagonalsFunc)(orbInfo, &I2PQ, &diagonals);
    assert(diagonals.size() == I2PQ.size());

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);

    // debug
    // this->debug_I2PQ_ = I2PQ;

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", threshold));
    const index_type numOfPQtilde = I2PQ.size();

    double error = std::accumulate(diagonals.begin(), diagonals.end(), 0.0);
    std::vector<std::size_t> pivot(numOfPQtilde);
    for (index_type i = 0; i < numOfPQtilde; ++i) {
        pivot[i] = i;
    }

    std::vector<double> errors(numOfPQtilde);

    int progress = 0;
    const index_type division = std::max<index_type>(numOfPQtilde * 0.01, 100);

    // ------------------------------------------------------------
    index_type reserveLCols = this->m_nNumOfAOs * 5;
    this->log_.info(TlUtils::format("resize L col: %d", reserveLCols));
    pL->resize(numOfPQtilde, reserveLCols);
    // ------------------------------------------------------------

    index_type numOfCDVcts = 0;
    while ((error > threshold) && (numOfCDVcts < numOfPQtilde)) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", numOfCDVcts, numOfPQtilde, error));
#endif  // DEBUG_CD

        // progress
        if (numOfCDVcts >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e", numOfCDVcts, error));
            ++progress;
        }
        pL->resize(numOfPQtilde, numOfCDVcts + 1);

        // memory allocation
        if (numOfCDVcts >= reserveLCols) {
            reserveLCols = numOfCDVcts + this->m_nNumOfAOs;
            this->log_.info(TlUtils::format("resize L col: %d", reserveLCols));
            pL->resize(numOfPQtilde, reserveLCols);
        }

        // pivot
        {
            const std::size_t argmax = this->argmax_pivot(diagonals, pivot, numOfCDVcts);
            std::swap(pivot[numOfCDVcts], pivot[argmax]);
        }

        const index_type pivot_m = pivot[numOfCDVcts];
        error = diagonals[pivot_m];
        errors[numOfCDVcts] = error;
        if (error < threshold) {
            break;
        }

        // get supermatrix elements
        const index_type numOf_G_cols = numOfPQtilde - (numOfCDVcts + 1);
        std::vector<double> G_pm(numOf_G_cols, 0.0);
        {
            std::vector<double> G_pm_local(numOf_G_cols, 0.0);
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type c = 0; c < numOf_G_cols; ++c) {
                const index_type pivot_i = pivot[(numOfCDVcts + 1) + c];  // from (m+1) to N
                G_col_list[c] = pivot_i;
            }
            (this->*getSuperMatrixElementsFunc)(orbInfo, pivot_m, G_col_list, I2PQ, &G_pm_local);
            assert(static_cast<index_type>(G_pm_local.size()) == numOf_G_cols);
            rComm.iAllReduce_SUM(&(G_pm_local[0]), &(G_pm[0]), numOf_G_cols);
        }

        // CD calc
        // output:
        //   out_L_rows: row elements at the column numOfCDVcts(target) in L
        std::valarray<double> L_pm_x(0.0, numOfCDVcts + 1);
        {
            const int PE_in_charge = pL->getSubunitID(pivot_m);
            if (PE_in_charge == myRank) {
                const std::size_t copyCount_m = pL->getRowVector(pivot_m, &(L_pm_x[0]), numOfCDVcts + 1);
                // std::cerr << copyCount_m << ", " << numOfCDVcts << ", " << pL->getNumOfCols() << std::endl;
                assert(copyCount_m == static_cast<std::size_t>(numOfCDVcts + 1));
            }
            rComm.broadcast(L_pm_x, PE_in_charge);
        }

        std::valarray<double> out_L_rows(0.0, numOfPQtilde);
        const double l_m_pm = std::sqrt(diagonals[pivot_m]);
        const double inv_l_m_pm = 1.0 / l_m_pm;
        // L.set(pivot_m, numOfCDVcts, l_m_pm);
        if (rComm.isMaster()) {
            out_L_rows[pivot_m] = l_m_pm;
        }

        std::vector<double> update_diagonals(numOfPQtilde, 0.0);
        rComm.wait(&(G_pm[0]));
#pragma omp parallel
        {
            std::valarray<double> L_pi_x(0.0, numOfCDVcts + 1);
#pragma omp for schedule(runtime)
            for (index_type i = 0; i < numOf_G_cols; ++i) {
                const index_type pivot_i = pivot[(numOfCDVcts + 1) + i];  // from (m+1) to N

                if (pL->getSubunitID(pivot_i) == myRank) {
                    const std::size_t copyCount_i = pL->getRowVector(pivot_i, &(L_pi_x[0]), numOfCDVcts + 1);
                    assert(copyCount_i == static_cast<std::size_t>(numOfCDVcts + 1));

                    const double sum_ll = (L_pm_x * L_pi_x).sum();
                    const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

                    out_L_rows[pivot_i] += l_m_pi;
                    update_diagonals[pivot_i] -= l_m_pi * l_m_pi;
                }
            }
        }

        // allReduce_SUM (out_L_rows)
        std::valarray<double> output_L_xm(0.0, numOfPQtilde);
        rComm.iAllReduce_SUM(&(out_L_rows[0]), &(output_L_xm[0]), numOfPQtilde);

        // allReduce_SUM (update_diagonals)
        std::vector<double> output_update_diagonals(numOfPQtilde, 0.0);
        rComm.iAllReduce_SUM(&(update_diagonals[0]), &(output_update_diagonals[0]), numOfPQtilde);

        rComm.wait(&(output_L_xm[0]));
        pL->setColVector(numOfCDVcts, output_L_xm);

        rComm.wait(&(output_update_diagonals[0]));
        // diagonals += output_update_diagonals;
        std::transform(diagonals.begin(), diagonals.end(), output_update_diagonals.begin(), diagonals.begin(),
                       std::plus<double>());

        // error = diagonals[pivot[numOfCDVcts]];
        ++numOfCDVcts;
    }

    if (rComm.isMaster()) {
        this->log_.info("save errors");
        errors.resize(numOfCDVcts);
        TlDenseVector_Lapack e(errors);
        DfObject::saveLjkErrorsVector(e);
    }

    this->log_.info(TlUtils::format("Cholesky Vectors: %d", numOfCDVcts));
    pL->resize(numOfPQtilde, numOfCDVcts);
}

void DfCD_Parallel::calcCholeskyVectorsOnTheFlyS(const TlOrbitalInfoObject& orbInfo, const std::string& I2PQ_path,
                                                 const double threshold, CalcDiagonalsFunc calcDiagonalsFunc,
                                                 GetSuperMatrixElementsFunc getSuperMatrixElements,
                                                 TlDenseGeneralMatrix_mmap* pL) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    this->log_.info("call on-the-fly Cholesky Decomposition routine (symmetric, parallel)");
    assert(this->pEngines_ != NULL);

    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfAOs * (numOfAOs + 1) / 2;
    this->log_.info(TlUtils::format("number of orbitals: %d", numOfAOs));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));

    // CDAM
    assert(this->pEngines_ != NULL);
    PQ_PairArray I2PQ;
    std::vector<double> diagonals;  // 対角成分
    (this->*calcDiagonalsFunc)(orbInfo, &I2PQ, &diagonals);
    assert(diagonals.size() == I2PQ.size());
    // if (rComm.isMaster()) {
    //     TlDenseVector_Lapack v_diagonal = diagonals;
    //     v_diagonal.save("CDAM_diag.vtr");
    // }

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);

    // debug
    // this->debug_I2PQ_ = I2PQ;

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", threshold));
    const index_type numOfPQtilde = I2PQ.size();

    double error = std::accumulate(diagonals.begin(), diagonals.end(), 0.0);
    std::vector<std::size_t> pivot(numOfPQtilde);
    for (index_type i = 0; i < numOfPQtilde; ++i) {
        pivot[i] = i;
    }

    int progress = 0;
    const index_type division = std::max<index_type>(numOfPQtilde * 0.01, 100);
    index_type L_cols = this->m_nNumOfAOs * 5;
    this->log_.info(TlUtils::format("reservoed L col: %d", L_cols));
    if (rComm.isMaster()) {
        pL->resize(numOfPQtilde, L_cols);
    }

    // TlTime timeLoop;
    // TlTime timeERI;
    // TlTime timeCD;

    // timeLoop.start();
    {
        index_type numOfCDVcts = 0;
        while ((error > threshold) && (numOfCDVcts < numOfPQtilde)) {
#ifdef DEBUG_CD
            this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", numOfCDVcts, numOfPQtilde, error));
#endif  // DEBUG_CD

            // progress
            if (numOfCDVcts >= progress * division) {
                this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e", numOfCDVcts, error));
                ++progress;
            }
            // メモリの確保
            if (numOfCDVcts >= L_cols) {
                L_cols = numOfCDVcts + this->m_nNumOfAOs;
                if (rComm.isMaster()) {
                    this->log_.info(TlUtils::format("L col: %d", L_cols));
                    pL->resize(numOfPQtilde, L_cols);
                }
            }

            // pivot
            {
                const std::size_t argmax = this->argmax_pivot(diagonals, pivot, numOfCDVcts);
                std::swap(pivot[numOfCDVcts], pivot[argmax]);
            }

            const index_type pivot_m = pivot[numOfCDVcts];
            error = diagonals[pivot_m];
            if (error < threshold) {
                break;
            }

            const double l_m_pm = std::sqrt(diagonals[pivot_m]);
            const double inv_l_m_pm = 1.0 / l_m_pm;
            // if (rComm.isMaster()) {
            //     pL->set(pivot_m, numOfCDVcts, l_m_pm);
            // }

            // get supermatrix elements
            // timeERI.start();
            const index_type numOf_G_cols = numOfPQtilde - (numOfCDVcts + 1);
            std::vector<double> G_pm(numOf_G_cols);
            {
                std::vector<index_type> G_col_list(numOf_G_cols);
                for (index_type c = 0; c < numOf_G_cols; ++c) {
                    const index_type pivot_i = pivot[(numOfCDVcts + 1) + c];  // from (m+1) to N
                    G_col_list[c] = pivot_i;
                }
                (this->*getSuperMatrixElements)(orbInfo, pivot_m, G_col_list, I2PQ, &G_pm);
            }
            assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);
            rComm.allReduce_SUM(G_pm);
            // timeERI.stop();

            // CD calc
            // output:
            //   out_L_rows; row elements at the column numOfCDVcts(target) in L
            //   diagonals:
            // timeCD.start();
            if (rComm.isMaster()) {
                std::valarray<double> L_pm_x(0.0, numOfCDVcts + 1);
                const std::size_t copyCount_m = pL->getRowVector(pivot_m, &(L_pm_x[0]), numOfCDVcts + 1);
                assert(copyCount_m == static_cast<std::size_t>(numOfCDVcts + 1));

                std::valarray<double> out_L_rows(0.0, numOfPQtilde);
                out_L_rows[pivot_m] = l_m_pm;
#pragma omp parallel
                {
                    std::valarray<double> L_pi_x(0.0, numOfCDVcts + 1);
#pragma omp for schedule(runtime)
                    for (index_type i = 0; i < numOf_G_cols; ++i) {
                        const index_type pivot_i = pivot[(numOfCDVcts + 1) + i];  // from (m+1) to N

                        const std::size_t copyCount_i = pL->getRowVector(pivot_i, &(L_pi_x[0]), numOfCDVcts + 1);
                        assert(copyCount_i == static_cast<std::size_t>(numOfCDVcts + 1));

                        const double sum_ll = (L_pm_x * L_pi_x).sum();
                        const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

                        out_L_rows[pivot_i] = l_m_pi;
                        diagonals[pivot_i] -= l_m_pi * l_m_pi;
                    }
                }
                pL->setColVector(numOfCDVcts, out_L_rows);
            }
            // timeCD.stop();

            rComm.broadcast(diagonals);
            ++numOfCDVcts;
        }
        // timeLoop.stop();

        if (rComm.isMaster()) {
            pL->resize(numOfPQtilde, numOfCDVcts);
        }
        this->log_.info(TlUtils::format("Cholesky Vectors: %d", numOfCDVcts));

        // this->log_.info(
        //     TlUtils::format("  ERI:  %16.2f sec", timeERI.getElapseTime()));
        // this->log_.info(
        //     TlUtils::format("  CD:   %16.2f sec", timeCD.getElapseTime()));
        // this->log_.info(
        //     TlUtils::format("  LOOP: %16.2f sec", timeLoop.getElapseTime()));
    }
}

TlDenseGeneralMatrix_arrays_RowOriented DfCD_Parallel::calcCholeskyVectorsOnTheFlyA(
    const TlOrbitalInfoObject& orbInfo_p, const TlOrbitalInfoObject& orbInfo_q, const std::string& I2PQ_path) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info(
        "call on-the-fly Cholesky Decomposition routine (MPI parallel; "
        "asymmetric)");
    assert(this->pEngines_ != NULL);
    this->initializeCutoffStats(std::max(orbInfo_p.getMaxShellType(), orbInfo_q.getMaxShellType()));

    const index_type numOfOrbs_p = orbInfo_p.getNumOfOrbitals();
    const index_type numOfOrbs_q = orbInfo_q.getNumOfOrbitals();
    const std::size_t numOfPQs = numOfOrbs_p * numOfOrbs_q;
    this->log_.info(TlUtils::format("number of orbitals1: %d", numOfOrbs_p));
    this->log_.info(TlUtils::format("number of orbitals2: %d", numOfOrbs_q));
    this->log_.info(TlUtils::format("number of pair of orbitals: %ld", numOfPQs));
    PQ_PairArray I2PQ;
    TlSparseMatrix schwartzTable(numOfOrbs_p, numOfOrbs_q);
    std::vector<double> global_diagonals;  // 対角成分

    assert(this->pEngines_ != NULL);
    this->calcDiagonalsA(orbInfo_p, orbInfo_q, &I2PQ, &schwartzTable, &global_diagonals);

    this->log_.info(TlUtils::format("number of screened pairs of orbitals: %ld", I2PQ.size()));
    this->saveI2PQ(I2PQ, I2PQ_path);
    // this->ERI_cache_manager_.setMaxItems(I2PQ.size() * 2);

    // prepare variables
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));
    const double threshold = this->epsilon_;
    const index_type N = I2PQ.size();

    // prepare variables (parallel)
    TlDenseGeneralMatrix_arrays_RowOriented L(
        N, 1, rComm.getNumOfProcs(),
        rComm.getRank());  // 答えとなる行列Lは各PEに行毎に短冊状(行ベクトル)で分散して持たせる
    const index_type local_N = L.getNumOfLocalVectors();
    std::vector<double> L_pm(N);
    std::vector<int> global_pivot(N);   // ERI計算リストを作成するために必要
    std::vector<int> reverse_pivot(N);  // global_pivotの逆引き
    std::vector<int> local_pivot(local_N);
    std::vector<double> local_diagonals(local_N);
    const int myRank = rComm.getRank();

    double error = 0.0;
    index_type error_global_loc = 0;
    index_type error_local_loc = 0;
    {
        int local_i = 0;
        for (int global_i = 0; global_i < N; ++global_i) {
            global_pivot[global_i] = global_i;
            reverse_pivot[global_i] = global_i;
            if (L.getSubunitID(global_i) == myRank) {
                local_pivot[local_i] = global_i;
                local_diagonals[local_i] = global_diagonals[global_i];
                if (error < local_diagonals[local_i]) {
                    error_local_loc = local_i;
                    error_global_loc = local_pivot[local_i];
                    error = local_diagonals[local_i];
                }
                ++local_i;
            }
        }
        assert(local_i == local_N);
    }
    rComm.allReduce_MAXLOC(&error, &error_global_loc);

    index_type m = 0;
    index_type local_m = 0;
    if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
        std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
        ++local_m;
    }

    int progress = 0;
    index_type division = std::max<index_type>(N * 0.01, 100);
    while (error > threshold) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 16.10e", m, N, error));
#endif  // DEBUG_CD

        // progress
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e", m, error));
            ++progress;

            // メモリの確保
            L.reserveColSize(progress * division);
        }
        L.resize(N, m + 1);

        // pivot
        std::swap(global_pivot[m], global_pivot[error_global_loc]);
        reverse_pivot[global_pivot[m]] = m;
        reverse_pivot[global_pivot[error_global_loc]] = error_global_loc;

        const double l_m_pm = std::sqrt(error);
        const index_type pivot_m = global_pivot[m];
        L.set(pivot_m, m, l_m_pm);  // 通信発生せず。関係無いPEは値を捨てる。
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // calc
        std::vector<double> G_pm;
        // const index_type numOf_G_cols = N -(m+1);
        const index_type numOf_G_cols = N - (m + 1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type i = 0; i < numOf_G_cols; ++i) {
                const index_type pivot_i = global_pivot[m + 1 + i];  // from (m+1) to N
                G_col_list[i] = pivot_i;
            }
            G_pm = this->getSuperMatrixElementsA(orbInfo_p, orbInfo_q, pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(static_cast<index_type>(G_pm.size()) == numOf_G_cols);

        // CD calc
        {
            // 全PEに分配
            const int PEinCharge = L.getSubunitID(pivot_m);
            if (PEinCharge == rComm.getRank()) {
                // const index_type copySize = L.getVector(pivot_m, &(L_pm[0]),
                // m +1);
                L.getVector(pivot_m, &(L_pm[0]), m + 1);
                // assert(L_pm[0].size() == m +1);
            }
            rComm.broadcast(&(L_pm[0]), m + 1, PEinCharge);
        }

        error = 0.0;
#pragma omp parallel
        {
            std::vector<double> L_pi(m + 1);
            double my_error = 0.0;
            int my_error_global_loc = 0;
            int my_error_local_loc = 0;

#pragma omp for schedule(runtime)
            for (int i = local_m; i < local_N; ++i) {
                const int pivot_i = local_pivot[i];
                // const index_type copySize = L.getVector(pivot_i, &(L_pi[0]),
                // m +1);
                L.getVector(pivot_i, &(L_pi[0]), m + 1);
                // assert(L_pi == m +1);
                double sum_ll = 0.0;
                for (index_type j = 0; j < m; ++j) {
                    sum_ll += L_pm[j] * L_pi[j];
                }

                const int G_pm_index = reverse_pivot[pivot_i] - (m + 1);
                const double l_m_pi = (G_pm[G_pm_index] - sum_ll) * inv_l_m_pm;
                const double ll = l_m_pi * l_m_pi;

#pragma omp critical(DfCD_Parallel__calcCholeskyVectorsOnTheFlyA_updateL)
                {
                    L.set(pivot_i, m, l_m_pi);
                    global_diagonals[pivot_i] -= ll;
                }

                if (global_diagonals[pivot_i] > my_error) {
                    my_error = global_diagonals[pivot_i];
                    my_error_global_loc = reverse_pivot[pivot_i];  // == m +1 + i
                    my_error_local_loc = i;
                }
            }

#pragma omp critical(DfCD_Parallel__calcCholeskyVectorsOnTheFlyA_update_error)
            {
                if (error < my_error) {
                    error = my_error;
                    error_global_loc = my_error_global_loc;
                    error_local_loc = my_error_local_loc;
                }
            }
        }
        rComm.allReduce_MAXLOC(&error, &error_global_loc);
        global_diagonals[global_pivot[error_global_loc]] = error;

        ++m;
        if (error_global_loc == reverse_pivot[local_pivot[error_local_loc]]) {
            std::swap(local_pivot[local_m], local_pivot[error_local_loc]);
            ++local_m;
        }
    }
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->schwartzCutoffReport(std::max(orbInfo_p.getMaxShellType(), orbInfo_q.getMaxShellType()));

    return L;
}

void DfCD_Parallel::getSuperMatrixElements(const TlOrbitalInfoObject& orbInfo, const index_type G_row,
                                           const std::vector<index_type>& G_col_list, const PQ_PairArray& I2PQ,
                                           std::vector<double>* pSuperMatrixElements) {
    assert(pSuperMatrixElements != NULL);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    this->ERI_cache_.clear();

    const index_type sizeOfCols = G_col_list.size();
    const index_type block = (sizeOfCols + numOfProcs - 1) / numOfProcs;
    const index_type start = std::min(block * myRank, sizeOfCols);
    const index_type end = std::min(block * (myRank + 1), sizeOfCols);

    const std::vector<IndexPair4S> calcList = this->getCalcList(orbInfo, G_row, G_col_list, start, end, I2PQ);
    this->calcERIs(orbInfo, calcList);
    *pSuperMatrixElements = this->setERIs(orbInfo, G_row, G_col_list, start, end, I2PQ);

    assert(static_cast<index_type>(pSuperMatrixElements->size()) == sizeOfCols);

    // rComm.allReduce_SUM(answer);
}

void DfCD_Parallel::getSuperMatrixElements_K_full(const TlOrbitalInfoObject& orbInfo, const index_type G_row,
                                                  const std::vector<index_type>& G_col_list, const PQ_PairArray& I2PR,
                                                  std::vector<double>* pSuperMatrixElements) {
    assert(pSuperMatrixElements != NULL);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    this->ERI_cache_.clear();

    const index_type sizeOfCols = G_col_list.size();
    const index_type block = (sizeOfCols + numOfProcs - 1) / numOfProcs;
    const index_type start = std::min(block * myRank, sizeOfCols);
    const index_type end = std::min(block * (myRank + 1), sizeOfCols);

    const std::vector<IndexPair4S> calcList = DfCD::getCalcList_K_full(orbInfo, G_row, G_col_list, start, end, I2PR);
    DfCD::calcERIs_K(orbInfo, calcList);
    *pSuperMatrixElements = DfCD::setERIs_K_full(orbInfo, G_row, G_col_list, start, end, I2PR);

    assert(static_cast<index_type>(pSuperMatrixElements->size()) == sizeOfCols);

    // rComm.allReduce_SUM(answer);
}

void DfCD_Parallel::getSuperMatrixElements_K_half(const TlOrbitalInfoObject& orbInfo, const index_type G_row,
                                                  const std::vector<index_type>& G_col_list, const PQ_PairArray& I2PR,
                                                  std::vector<double>* pSuperMatrixElements) {
    assert(pSuperMatrixElements != NULL);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    this->ERI_cache_.clear();

    const index_type sizeOfCols = G_col_list.size();
    const index_type block = (sizeOfCols + numOfProcs - 1) / numOfProcs;
    const index_type start = std::min(block * myRank, sizeOfCols);
    const index_type end = std::min(block * (myRank + 1), sizeOfCols);

    const std::vector<IndexPair4S> calcList = DfCD::getCalcList_K_half(orbInfo, G_row, G_col_list, start, end, I2PR);
    DfCD::calcERIs_K(orbInfo, calcList);
    *pSuperMatrixElements = DfCD::setERIs_K_half(orbInfo, G_row, G_col_list, start, end, I2PR);

    assert(static_cast<index_type>(pSuperMatrixElements->size()) == sizeOfCols);

    // rComm.allReduce_SUM(answer);
}

void DfCD_Parallel::saveL(const TlDenseGeneralMatrix_arrays_RowOriented& L, const std::string& path) {
    L.save(path + ".rvm");

    switch (this->cdFileFormat_) {
        case CD_FILE_FORMAT_CSFD: {
            this->log_.info("save the Cholesky vectors. (mmap)");
            this->transLMatrix2mmap(L, path);
        } break;

        case CD_FILE_FORMAT_ABGD: {
            this->log_.info("transform the Cholesky vectors.");
            TlDenseGeneralMatrix_arrays_ColOriented colVecL = this->transLMatrix(L);

            this->log_.info("save the Cholesky vectors.");
            colVecL.save(path);
        } break;

        default: {
            this->log_.critical(TlUtils::format("program error: %s,%d", __FILE__, __LINE__));
        }
    }
}

TlDenseGeneralMatrix_Lapack DfCD_Parallel::mergeL(const TlDenseGeneralMatrix_arrays_RowOriented& L) {
    this->log_.info("merge L(row): start");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const index_type numOfRows = L.getNumOfRows();
    const index_type numOfCols = L.getNumOfCols();
    TlDenseGeneralMatrix_Lapack answer(numOfRows, numOfCols);

    const div_t turns = std::div(numOfRows, numOfProcs);
    const index_type localRows = turns.quot + ((rank < turns.rem) ? 1 : 0);
    for (index_type r = 0; r < localRows; ++r) {
        const index_type row = r * numOfProcs + rank;
        TlDenseVector_Lapack rowVec = L.getVector(row);
        for (index_type col = 0; col < numOfCols; ++col) {
            answer.set(row, col, rowVec.get(col));
        }
    }
    rComm.allReduce_SUM(&answer);

    this->log_.info("merge L: end");
    return answer;
}

TlDenseGeneralMatrix_Lapack DfCD_Parallel::mergeL(const TlDenseGeneralMatrix_arrays_ColOriented& L) {
    this->log_.info("merge L(col): start");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();

    const index_type numOfRows = L.getNumOfRows();
    const index_type numOfCols = L.getNumOfCols();
    TlDenseGeneralMatrix_Lapack answer(numOfRows, numOfCols);

    const div_t turns = std::div(numOfCols, numOfProcs);
    const index_type localCols = turns.quot + ((rank < turns.rem) ? 1 : 0);
    for (index_type c = 0; c < localCols; ++c) {
        const index_type col = c * numOfProcs + rank;
        TlDenseVector_Lapack colVec = L.getVector(col);
        for (index_type row = 0; row < numOfRows; ++row) {
            answer.set(row, col, colVec.get(row));
        }
    }
    rComm.allReduce_SUM(&answer);

    this->log_.info("merge L: end");
    return answer;
}

void DfCD_Parallel::getJ(TlDenseSymmetricMatrix_Lapack* pJ) {
    switch (this->cdFileFormat_) {
        case CD_FILE_FORMAT_CSFD: {
            this->getJ_mmap_DC(pJ);
        } break;

        case CD_FILE_FORMAT_ABGD: {
            this->getJ_cvm(pJ);
        } break;

        default: {
            this->log_.critical(TlUtils::format("program error: %s,%d", __FILE__, __LINE__));
        }
    }
}

void DfCD_Parallel::getJ_cvm(TlDenseSymmetricMatrix_Lapack* pJ) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc J by CD method (parallel).");

    // cholesky vector
    TlDenseGeneralMatrix_arrays_ColOriented L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLjkMatrixPath(), rComm.getRank());
    assert(L.getNumOfSubunits() == rComm.getNumOfProcs());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const TlDenseVector_Lapack vP = DfCD_Parallel::getScreenedDensityMatrix<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(
        DfCD_Parallel::getPMatrix<TlDenseSymmetricMatrix_Lapack>(), I2PQ);

    const index_type cvSize = L.getNumOfRows();
    assert(std::size_t(cvSize) == I2PQ.size());
    const index_type numOfCVs = L.getNumOfCols();

    TlDenseVector_Lapack vJ(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            const TlDenseVector_Lapack LI = L.getVector(I);

            TlDenseVector_Lapack tmpLI = LI;
            const double qi = tmpLI.dotInPlace(vP).sum();

            vJ += qi * LI;
        }
    }

    rComm.allReduce_SUM(&vJ);
    this->expandJMatrix(vJ, I2PQ, pJ);
}

void DfCD_Parallel::getJ_mmap_DC(TlDenseSymmetricMatrix_Lapack* pJ) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    this->log_.info("calc J by CD method (parallel; mmap).");

    // cholesky vector
    TlDenseGeneralMatrix_mmap L(DfObject::getLjkMatrixPath());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const TlDenseVector_Lapack vP = DfCD_Parallel::getScreenedDensityMatrix<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(
        DfCD_Parallel::getPMatrix<TlDenseSymmetricMatrix_Lapack>(), I2PQ);

    const index_type cvSize = L.getNumOfRows();
    assert(std::size_t(cvSize) == I2PQ.size());
    const index_type numOfCVs = L.getNumOfCols();

    const index_type numOfTasks = (numOfCVs + numOfProcs - 1) / numOfProcs;
    const index_type beginTask = numOfTasks * myRank;
    const index_type endTask = std::min(numOfTasks * (myRank + 1), numOfCVs);

    TlDenseVector_Lapack vJ(cvSize);
    for (index_type I = beginTask; I < endTask; ++I) {
        const TlDenseVector_Lapack LI = L.getColVector(I);

        TlDenseVector_Lapack tmpLI = LI;
        const double qi = tmpLI.dotInPlace(vP).sum();

        vJ += qi * LI;
    }

    rComm.allReduce_SUM(&vJ);
    this->expandJMatrix(vJ, I2PQ, pJ);
}

TlDenseSymmetricMatrix_Lapack DfCD_Parallel::getPMatrix(const RUN_TYPE runType, const int itr) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    TlDenseSymmetricMatrix_Lapack P;
    if (rComm.isMaster()) {
        P = DfCD::getPMatrix(runType, itr);
    }
    rComm.broadcast(&P);
    return P;
}

void DfCD_Parallel::getK_S_woCD(const RUN_TYPE runType, TlDenseSymmetricMatrix_Lapack* pK) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc K by CD method (parallel).");

    // cholesky vector
    TlDenseGeneralMatrix_arrays_ColOriented L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLjkMatrixPath(), rComm.getRank());
    assert(L.getNumOfSubunits() == rComm.getNumOfProcs());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();

    this->log_.info("load density matrix");
    TlDenseSymmetricMatrix_Lapack P = this->getPMatrix(runType, this->m_nIteration - 1);

    this->log_.info("calc by Cholesky Vectors");
    TlDenseVector_Lapack cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);

            const TlDenseSymmetricMatrix_Lapack l = this->getCholeskyVector<TlDenseSymmetricMatrix_Lapack>(cv, I2PQ);

            TlDenseGeneralMatrix_Lapack X = l * P;
            X *= l;

            *pK += X;
        }
    }
    rComm.allReduce_SUM(pK);

    *pK *= -1.0;
}

void DfCD_Parallel::getK_S_woCD_mmap_DC(const RUN_TYPE runType, TlDenseSymmetricMatrix_Lapack* pK) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    this->log_.info("calc K by CD method (parallel; mmap).");

    // cholesky vector
    TlDenseGeneralMatrix_mmap L(DfObject::getLjkMatrixPath());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();

    this->log_.info("load density matrix");
    TlDenseSymmetricMatrix_Lapack P = this->getPMatrix(runType, this->m_nIteration - 1);

    this->log_.info("calc by Cholesky Vectors");
    const index_type numOfTasks = (numOfCVs + numOfProcs - 1) / numOfProcs;
    const index_type beginTask = numOfTasks * myRank;
    const index_type endTask = std::min(numOfTasks * (myRank + 1), numOfCVs);

    TlDenseVector_Lapack cv(cvSize);
    for (index_type I = beginTask; I < endTask; ++I) {
        cv = L.getColVector(I);
        assert(cv.getSize() == cvSize);

        const TlDenseSymmetricMatrix_Lapack l = this->getCholeskyVector<TlDenseSymmetricMatrix_Lapack>(cv, I2PQ);

        TlDenseGeneralMatrix_Lapack X = l * P;
        X *= l;

        *pK += X;
    }
    rComm.allReduce_SUM(pK);

    *pK *= -1.0;
}

void DfCD_Parallel::getK_byLk(const RUN_TYPE runType, TlDenseSymmetricMatrix_Lapack* pK) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc K(fast) by CD method (parallel).");

    // cholesky vector
    TlDenseGeneralMatrix_arrays_ColOriented L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    // TlDenseGeneralMatrix_arrays_ColOriented L;
    L.load(DfObject::getLkMatrixPath(), rComm.getRank());
    assert(L.getNumOfSubunits() == rComm.getNumOfProcs());

    const PQ_PairArray I2PR = this->getI2PQ(this->getI2prVtrPath());
    const TlDenseVector_Lapack vP = DfCD_Parallel::getScreenedDensityMatrix<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(runType, I2PR);

    const index_type cvSize = L.getNumOfRows();
    assert(std::size_t(cvSize) == I2PR.size());
    const index_type numOfCVs = L.getNumOfCols();

    TlDenseVector_Lapack vK(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            const TlDenseVector_Lapack LI = L.getVector(I);

            TlDenseVector_Lapack tmpLI = LI;
            const double qi = tmpLI.dotInPlace(vP).sum();

            vK += qi * LI;
        }
    }
    vK *= -1.0;

    rComm.allReduce_SUM(&vK);
    this->expandKMatrix(vK, I2PR, pK);
}

void DfCD_Parallel::getJ_D(TlDenseSymmetricMatrix_Scalapack* pJ) {
    this->log_.info("calc J by CD method (parallel; distributed).");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    const TlDenseGeneralMatrix_mmap L(DfObject::getLjkMatrixPath());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const TlDenseVector_Lapack vP = this->getScreenedDensityMatrixD(I2PQ);
    // vP.save("vP.vtr");

    const index_type cvSize = L.getNumOfRows();
    assert(std::size_t(cvSize) == I2PQ.size());
    const index_type numOfCVs = L.getNumOfCols();

    const index_type numOfTasks = (numOfCVs + numOfProcs - 1) / numOfProcs;
    const index_type beginTask = numOfTasks * myRank;
    const index_type endTask = std::min(numOfTasks * (myRank + 1), numOfCVs);

    TlDenseVector_Lapack vJ(cvSize);
    for (index_type I = beginTask; I < endTask; ++I) {
        const TlDenseVector_Lapack LI = L.getColVector(I);

        TlDenseVector_Lapack tmpLI = LI;
        const double qi = tmpLI.dotInPlace(vP).sum();

        vJ += qi * LI;
    }

    rComm.allReduce_SUM(&vJ);
    // vJ.save("vJ.vtr");
    this->expandJMatrixD(vJ, I2PQ, pJ);
}

TlDenseVector_Lapack DfCD_Parallel::getScreenedDensityMatrixD(const PQ_PairArray& I2PQ) {
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const TlDenseSymmetricMatrix_Scalapack P = DfObject::getPInMatrix<TlDenseSymmetricMatrix_Scalapack>(RUN_RKS, this->m_nIteration);
    const std::size_t numOfI = I2PQ.size();
    TlDenseVector_Lapack answer(numOfI);

    for (std::size_t i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PQ[i];
        const index_type r = pair.index1();
        const index_type c = pair.index2();

        const double coef = (r != c) ? 2.0 : 1.0;
        // The elements of density matrix (P) are get by `getLocal()`.
        answer.set(i, coef * P.getLocal(r, c));
    }
    rComm.allReduce_SUM(&answer);

    return answer;
}

void DfCD_Parallel::expandJMatrixD(const TlDenseVector_Lapack& vJ, const PQ_PairArray& I2PQ,
                                   TlDenseSymmetricMatrix_Scalapack* pJ) {
    assert(pJ != NULL);
    const index_type numOfI = I2PQ.size();
    for (index_type i = 0; i < numOfI; ++i) {
        const Index2& pair = I2PQ[i];
        pJ->set(pair.index1(), pair.index2(), vJ.get(i));
    }
}

void DfCD_Parallel::getK_D(const RUN_TYPE runType) {
    TlDenseSymmetricMatrix_Scalapack K(this->m_nNumOfAOs);
    this->getK_S_woCD_D(runType, &K);

    this->file_->saveMatrix(DfObject::getHFxMatrixPath(runType, this->m_nIteration), K);
}

void DfCD_Parallel::getK_S_woCD_D(const RUN_TYPE runType, TlDenseSymmetricMatrix_Scalapack* pK) {
    this->log_.info("calc K by CD method (parallel; distributed).");
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const TlDenseGeneralMatrix_mmap L(DfObject::getLjkMatrixPath());

    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();

    this->log_.info("load density matrix");
    const TlDenseSymmetricMatrix_Scalapack P = 0.5 * DfObject::getPInMatrix<TlDenseSymmetricMatrix_Scalapack>(runType, this->m_nIteration);  // RKS

    TlTime time_all;
    TlTime time_bcast;
    TlTime time_translate;
    TlTime time_multimat1;
    TlTime time_multimat2;
    TlTime time_transpose;
    TlTime time_add;

    time_all.start();
    this->log_.info("calc by Cholesky Vectors");
    int progress = 0;
    const int division = numOfCVs * 0.1;
    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        // progress
        if (I >= progress * division) {
            const double rate = double(I) / double(numOfCVs) * 100.0;
            this->log_.info(TlUtils::format("K loop progress: %5.2f%%", rate));
            ++progress;
        }

        time_bcast.start();
        // const int PEinCharge = L.getSubunitID(I);
        // if (PEinCharge == rComm.getRank()) {
        //     cv = L.getVector(I);
        //     assert(cv.getSize() == cvSize);
        // }
        // rComm.broadcast(&(cv[0]), cvSize, PEinCharge);
        if (rComm.isMaster()) {
            cv = L.getColVector(I);
            assert(static_cast<TlMatrixObject::index_type>(cv.size()) == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, 0);
        time_bcast.stop();

        time_translate.start();
        const TlDenseSymmetricMatrix_Scalapack l = this->getCholeskyVector_distribute(cv, I2PQ);
        time_translate.stop();

        time_multimat1.start();
        TlDenseGeneralMatrix_Scalapack X = l * P;
        time_multimat1.stop();

        time_multimat2.start();
        X *= l;
        time_multimat2.stop();

        time_add.start();
        *pK += X;
        time_add.stop();
    }

    *pK *= -1.0;
    time_all.stop();

    // timing data
    this->log_.info(TlUtils::format("K all:       %10.1f sec.", time_all.getElapseTime()));
    this->log_.info(TlUtils::format("K bcast:     %10.1f sec.", time_bcast.getElapseTime()));
    this->log_.info(TlUtils::format("K translate: %10.1f sec.", time_translate.getElapseTime()));
    this->log_.info(TlUtils::format("K multi1:    %10.1f sec.", time_multimat1.getElapseTime()));
    this->log_.info(TlUtils::format("K multi2:    %10.1f sec.", time_multimat2.getElapseTime()));
    this->log_.info(TlUtils::format("K transpose: %10.1f sec.", time_transpose.getElapseTime()));
    this->log_.info(TlUtils::format("K add:       %10.1f sec.", time_add.getElapseTime()));
}

void DfCD_Parallel::getM(const TlDenseSymmetricMatrix_Lapack& P, TlDenseSymmetricMatrix_Lapack* pM) {
    this->log_.info("DfCD_Parallel::getM()");
    if (this->isDedicatedBasisForGridFree_) {
        this->getM_A(P, pM);
    } else {
        this->getM_S(P, pM);
    }
}

void DfCD_Parallel::getM_S(const TlDenseSymmetricMatrix_Lapack& P, TlDenseSymmetricMatrix_Lapack* pM) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (symmetric routine; parallel)");

    const TlOrbitalInfo orbInfo((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set"]);
    const index_type numOfAOs = orbInfo.getNumOfOrbitals();
    pM->resize(numOfAOs);

    // cholesky vector
    TlDenseGeneralMatrix_arrays_ColOriented L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLxcMatrixPath());
    // assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    // assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());
    this->log_.info(TlUtils::format("I2PQ size: %ld", I2PQ.size()));

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    this->log_.info(TlUtils::format("L size: %d x %d", cvSize, numOfCVs));

    TlDenseVector_Lapack cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            // const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);

            TlDenseSymmetricMatrix_Lapack LI = this->getCholeskyVector<TlDenseSymmetricMatrix_Lapack>(cv, I2PQ);
            assert(LI.getNumOfRows() == this->m_nNumOfAOs);
            assert(LI.getNumOfCols() == this->m_nNumOfAOs);

            TlDenseSymmetricMatrix_Lapack QI = LI;
            QI.dotInPlace(P);
            const double qi = QI.sum();

            *pM += qi * LI;
        }
    }

    rComm.allReduce_SUM(pM);
}

void DfCD_Parallel::getM_A(const TlDenseSymmetricMatrix_Lapack& P, TlDenseSymmetricMatrix_Lapack* pM) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (asymmetric routine; parallel)");

    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set_gridfree"]);
    const index_type dim_M = orbInfo_q.getNumOfOrbitals();
    pM->resize(dim_M);

    // cholesky vector
    TlDenseGeneralMatrix_arrays_ColOriented L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLxcMatrixPath());
    // assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    // assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrXCPath());
    this->log_.info(TlUtils::format("I2PQ size: %ld", I2PQ.size()));

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    this->log_.info(TlUtils::format("L size: %d x %d", cvSize, numOfCVs));

    TlDenseGeneralMatrix_Lapack C;
    P.pivotedCholeskyDecomposition(&C, this->epsilon_);

    TlDenseVector_Lapack cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            // const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(cv.getSize() == cvSize);

            const TlDenseGeneralMatrix_Lapack l = this->getCholeskyVectorA(orbInfo_p, orbInfo_q, cv, I2PQ);
            const TlDenseGeneralMatrix_Lapack lt = l.transpose();

            // TODO: ltを作らずにXtを計算できるはず。
            const TlDenseGeneralMatrix_Lapack X = lt * C;
            const TlDenseGeneralMatrix_Lapack Xt = X.transpose();

            TlDenseSymmetricMatrix_Lapack XX = X * Xt;
            *pM += XX;
        }
    }
    rComm.allReduce_SUM(pM);
}

void DfCD_Parallel::getM(const TlDenseSymmetricMatrix_Scalapack& P, TlDenseSymmetricMatrix_Scalapack* pM) {
    this->log_.critical("sorry, NO IMPLEMENTED!");
    this->log_.critical(
        "DfCD_Parallel::getM(const TlDenseSymmetricMatrix_Scalapack&, "
        "TlDenseSymmetricMatrix_Scalapack*)");
    abort();
    // if (this->isDedicatedBasisForGridFree_) {
    //     this->getM_A(P, pM);
    // } else {
    //     this->getM_S(P, pM);
    // }
}

void DfCD_Parallel::getM_S(const TlDenseSymmetricMatrix_Scalapack& P, TlDenseSymmetricMatrix_Scalapack* pM) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (symmetric routine; parallel; distributed)");

    // const TlDenseSymmetricMatrix_Scalapack P =
    //     DfObject::getPInMatrix<TlDenseSymmetricMatrix_Scalapack>(RUN_RKS,
    //     this->m_nIteration);

    // cholesky vector
    TlDenseGeneralMatrix_arrays_ColOriented L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLxcMatrixPath());
    // assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    // assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    this->log_.info(TlUtils::format("I2PQ size: %ld", I2PQ.size()));

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    this->log_.info(TlUtils::format("L size: %d x %d", cvSize, numOfCVs));

    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            // const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(static_cast<TlMatrixObject::index_type>(cv.size()) == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, PEinCharge);

        TlDenseSymmetricMatrix_Scalapack LI = this->getCholeskyVector_distribute(cv, I2PQ);
        assert(LI.getNumOfRows() == this->m_nNumOfAOs);
        assert(LI.getNumOfCols() == this->m_nNumOfAOs);

        TlDenseSymmetricMatrix_Scalapack QI = LI;
        QI.dotInPlace(P);
        const double qi = QI.sum();

        *pM += qi * LI;
    }
}

void DfCD_Parallel::getM_A(const TlDenseSymmetricMatrix_Scalapack& P, TlDenseSymmetricMatrix_Scalapack* pM) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    this->log_.info("calc M by CD method. (asymmetric routine; parallel; distributed)");

    const TlOrbitalInfo orbInfo_p((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set"]);
    const TlOrbitalInfo orbInfo_q((*this->pPdfParam_)["coordinates"], (*this->pPdfParam_)["basis_set_gridfree"]);
    const index_type dim_M = orbInfo_q.getNumOfOrbitals();
    pM->resize(dim_M);

    // cholesky vector
    TlDenseGeneralMatrix_arrays_ColOriented L(1, 1, rComm.getNumOfProcs(), rComm.getRank());
    L.load(DfObject::getLxcMatrixPath());
    // assert(L.getNumOfAllProcs() == rComm.getNumOfProcs());
    // assert(L.getRank() == rComm.getRank());
    const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());
    this->log_.info(TlUtils::format("I2PQ size: %ld", I2PQ.size()));

    const index_type cvSize = L.getNumOfRows();
    const index_type numOfCVs = L.getNumOfCols();
    this->log_.info(TlUtils::format("L size: %d x %d", cvSize, numOfCVs));

    this->log_.info("calc CD of density matrix");
    this->log_.info(TlUtils::format("epsilon = %8.3e", this->epsilon_));
    // TlDenseSymmetricMatrix_Scalapack P =
    //     0.5 *
    //     DfObject::getPInMatrix<TlDenseSymmetricMatrix_Scalapack>(RUN_RKS,
    //     this->m_nIteration); // RKS
    TlDenseGeneralMatrix_Scalapack C;
    P.pivotedCholeskyDecomposition(&C, this->epsilon_);

    this->log_.info("start loop");
    int progress = 0;
    const int division = numOfCVs * 0.1;
    std::vector<double> cv(cvSize);
    for (index_type I = 0; I < numOfCVs; ++I) {
        // progress
        if (I >= progress * division) {
            const double rate = double(I) / double(numOfCVs) * 100.0;
            this->log_.info(TlUtils::format("K loop progress: %5.2f%%", rate));
            ++progress;
        }

        const int PEinCharge = L.getSubunitID(I);
        if (PEinCharge == rComm.getRank()) {
            // const index_type copySize = L.getVector(I, &(cv[0]), cvSize);
            cv = L.getVector(I);
            assert(static_cast<TlMatrixObject::index_type>(cv.size()) == cvSize);
        }
        rComm.broadcast(&(cv[0]), cvSize, PEinCharge);

        TlDenseGeneralMatrix_Scalapack l = this->getCholeskyVectorA_distribute(orbInfo_p, orbInfo_q, cv, I2PQ);
        l.transposeInPlace();

        TlDenseGeneralMatrix_Scalapack X = l * C;
        TlDenseGeneralMatrix_Scalapack Xt = X;
        Xt.transposeInPlace();
        TlDenseSymmetricMatrix_Scalapack XX = X * Xt;
        *pM += XX;
    }
}

TlDenseGeneralMatrix_arrays_ColOriented DfCD_Parallel::transLMatrix(
    const TlDenseGeneralMatrix_arrays_RowOriented& rowVectorMatrix) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();
    const int tag = 12345;  // 適当

    const index_type numOfRows = rowVectorMatrix.getNumOfRows();
    const index_type numOfCols = rowVectorMatrix.getNumOfCols();

    TlDenseGeneralMatrix_arrays_ColOriented colVectorMatrix(numOfRows, numOfCols, numOfProcs, myRank);

    // prepare partial matrix
    std::vector<unsigned long> matrixElementsSizes(numOfProcs);
    std::vector<std::vector<TlMatrixObject::MatrixElement> > matrixElements(numOfProcs);
    {
        std::vector<TlSparseMatrix> partMat(numOfProcs, TlSparseMatrix(numOfRows, numOfCols));
        for (index_type row = 0; row < numOfRows; ++row) {
            if (rowVectorMatrix.getSubunitID(row) == myRank) {
                for (index_type col = 0; col < numOfCols; ++col) {
                    const double v = rowVectorMatrix.get(row, col);

                    const int unit = colVectorMatrix.getSubunitID(col);
                    partMat[unit].set(row, col,
                                      v);  // TODO: ベクトルのアクセスにするべき
                }
            }
        }

        for (int proc = 0; proc < numOfProcs; ++proc) {
            matrixElementsSizes[proc] = partMat[proc].getSize();
            matrixElements[proc] = partMat[proc].getMatrixElements();
        }
    }

    // send
    for (int proc = 0; proc < numOfProcs; ++proc) {
        if (proc != myRank) {
            rComm.iSendData(matrixElementsSizes[proc], proc, tag);

            if (matrixElementsSizes[proc] > 0) {
                rComm.iSendDataX(&(matrixElements[proc][0]), matrixElementsSizes[proc], proc, tag + 1);
            }
        }
    }

    // my partMat
    {
        std::vector<TlMatrixObject::MatrixElement>::const_iterator itEnd = matrixElements[myRank].end();
        for (std::vector<TlMatrixObject::MatrixElement>::const_iterator it = matrixElements[myRank].begin();
             it != itEnd; ++it) {
            const index_type r = it->row;
            const index_type c = it->col;
            colVectorMatrix.set(r, c, it->value);
        }
    }

    // recv
    {
        for (int i = 0; i < (numOfProcs - 1); ++i) {
            int proc = 0;
            unsigned long tmp_size = 0;
            rComm.receiveDataFromAnySource(tmp_size, &proc, tag);

            if (tmp_size > 0) {
                std::vector<TlMatrixObject::MatrixElement> tmp_elements(tmp_size);
                rComm.receiveDataX(&(tmp_elements[0]), tmp_size, proc, tag + 1);

                std::vector<TlMatrixObject::MatrixElement>::const_iterator itEnd = tmp_elements.end();
                for (std::vector<TlMatrixObject::MatrixElement>::const_iterator it = tmp_elements.begin(); it != itEnd;
                     ++it) {
                    colVectorMatrix.set(it->row, it->col, it->value);
                }
            }
        }
    }

    // wait
    for (int proc = 0; proc < numOfProcs; ++proc) {
        if (proc != myRank) {
            rComm.wait(&(matrixElementsSizes[proc]));
            rComm.wait(&(matrixElements[proc][0]));
        }
    }

    return colVectorMatrix;
}

void DfCD_Parallel::transLMatrix2mmap(const TlDenseGeneralMatrix_arrays_RowOriented& rowVectorMatrix,
                                      const std::string& savePath) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int myRank = rComm.getRank();

    const index_type numOfRows = rowVectorMatrix.getNumOfRows();
    const index_type numOfCols = rowVectorMatrix.getNumOfCols();
    const int sizeOfChunk = rowVectorMatrix.getSizeOfChunk();

    // const int root = 0;
    // const int endMessageTag = -1;

    if (rComm.isMaster()) {
        // std::cout << TlUtils::format("mmap: %d, %d (chunk=%d)", numOfRows,
        //                              numOfCols, sizeOfChunk)
        //           << std::endl;
        TlDenseGeneralMatrix_mmap output(savePath, numOfRows, numOfCols);

        // host data
        {
            std::vector<double> chunkBuf(numOfCols * sizeOfChunk);
            std::vector<double> transBuf(numOfCols * sizeOfChunk);
            for (index_type row = 0; row < numOfRows; row += sizeOfChunk) {
                if (rowVectorMatrix.getSubunitID(row) == 0) {
                    // std::cout << TlUtils::format("[0] row=%d", row) <<
                    // std::endl;
                    rowVectorMatrix.getChunk(row, &(chunkBuf[0]), numOfCols * sizeOfChunk);

                    // change memory layout
                    const index_type readRowChunks = std::min(sizeOfChunk, numOfRows - row);
                    TlUtils::changeMemoryLayout(&(chunkBuf[0]), readRowChunks, numOfCols, &(transBuf[0]));

                    // write to matrix
#if defined HAVE_EIGEN
                    TlDenseGeneralMatrix_Eigen tmpMat(readRowChunks, numOfCols, &(transBuf[0]));
#elif defined HAVE_LAPACK
                    TlDenseGeneralMatrix_Lapack tmpMat(readRowChunks, numOfCols, &(transBuf[0]));
#else
#error "not implemented matrix type"
#endif  //
                    output.block(row, 0, tmpMat);
                }
            }
        }

        // remote data
        if (numOfProcs > 1) {
            std::vector<double> recvBuf(numOfCols * sizeOfChunk);
            std::vector<int> recvCounts(numOfProcs, 0);

            // std::cout << "begin loop:" << std::endl;
            int src, tag;
            bool endOfLoop = false;
            while (endOfLoop == false) {
                // recv
                rComm.iReceiveDataFromAnySourceAnyTag(&(recvBuf[0]), numOfCols * sizeOfChunk);
                rComm.wait(&(recvBuf[0]), &src, &tag);

                // write to matrix
                const index_type row = sizeOfChunk * src + recvCounts[src] * sizeOfChunk * numOfProcs;
                // std::cout << TlUtils::format("[0] <- [%d] row=%d tag=%d",
                // src,
                //                              row, tag)
                //           << std::endl;
                assert(row == tag);
                const index_type readRowChunks = std::min(sizeOfChunk, numOfRows - row);
                {
#if defined HAVE_EIGEN
                    TlDenseGeneralMatrix_Eigen tmpMat(readRowChunks, numOfCols, &(recvBuf[0]));
#elif defined HAVE_LAPACK
                    TlDenseGeneralMatrix_Lapack tmpMat(readRowChunks, numOfCols, &(recvBuf[0]));
#else
#error "not implemented matrix type"
#endif  //
                    output.block(row, 0, tmpMat);
                }

                ++recvCounts[src];

                endOfLoop = true;
                for (int proc = 1; proc < numOfProcs; ++proc) {
                    // How much data was received?
                    const index_type recvIndex = sizeOfChunk * proc + (recvCounts[proc] - 1) * sizeOfChunk * numOfProcs;

                    if ((recvIndex + sizeOfChunk * numOfProcs) < numOfRows) {
                        endOfLoop = false;
                        // std::cout
                        //     << TlUtils::format("[0] end range @%d %d < %d
                        //     [NG]",
                        //                        proc, recvIndex, numOfRows)
                        //     << std::endl;
                        break;
                    } else {
                        this->log_.debug(TlUtils::format("end range @%d %d > %d [OK]", proc, recvIndex, numOfRows));
                    }
                }
            }
            // std::cout << "end loop:" << std::endl;
        }
    } else {
        std::vector<double> chunkBuf(numOfCols * sizeOfChunk);
        std::vector<double> sendBuf(numOfCols * sizeOfChunk);
        for (index_type row = 0; row < numOfRows; row += sizeOfChunk) {
            // if (myRank == 1) {
            //     std::cout << TlUtils::format("row = %d, subunitid = %d", row,
            //                                  rowVectorMatrix.getSubunitID(row))
            //               << std::endl;
            // }

            // get chunk
            if (rowVectorMatrix.getSubunitID(row) == myRank) {
                const int copied = rowVectorMatrix.getChunk(row, &(chunkBuf[0]), numOfCols * sizeOfChunk);
                assert(copied <= numOfCols * sizeOfChunk);

                // change memory layout
                const index_type readChunks = std::min(sizeOfChunk, numOfRows - row);
                std::fill(sendBuf.begin(), sendBuf.end(), 0.0);
                TlUtils::changeMemoryLayout(&(chunkBuf[0]), readChunks, numOfCols, &(sendBuf[0]));

                // send
                rComm.iSendDataX((double*)(&(sendBuf[0])), numOfCols * sizeOfChunk, 0, row);
                rComm.wait(&(sendBuf[0]));
                // std::cout << TlUtils::format("[%d] send tag=%d", myRank, row)
                //           << std::endl;
            }
        }

        // std::cout << TlUtils::format("[%d] send finish", myRank) <<
        // std::endl;
    }

    assert(rComm.checkNonBlockingCommunications());
}
