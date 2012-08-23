#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <set>
#include "DfCD.h"
#include "DfEriEngine.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlTime.h"
#include "TlRowVectorMatrix2.h"
#include "TlUtils.h"
#include "TlSystem.h"

DfCD::DfCD(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam), pEriEngines_(NULL),
      orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_sets"])
{
    this->numOfPQs_ = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

    this->cutoffThreshold_ = 1.0E-10;
    if ((*pPdfParam)["cut-value"].getStr().empty() != true) {
        this->cutoffThreshold_ = (*pPdfParam)["cut-value"].getDouble();
    }    
    this->cutoffEpsilon3_ = this->cutoffThreshold_ * 0.01;
    if ((*pPdfParam)["cutoff_epsilon3"].getStr().empty() != true) {
        this->cutoffEpsilon3_ = (*pPdfParam)["cutoff_epsilon3"].getDouble();
    }    

    this->CDAM_tau_ = 1.0E-4;
    if ((*pPdfParam)["CDAM_tau"].getStr().empty() != true) {
        this->CDAM_tau_ = (*pPdfParam)["CDAM_tau"].getDouble();
    }    

    this->epsilon_ = 1.0E-10;
    if ((*pPdfParam)["CD_epsilon"].getStr().empty() != true) {
        this->epsilon_ = (*pPdfParam)["CD_epsilon"].getDouble();
    }    

    this->isStoreERIs_ = true;
    if ((*pPdfParam)["CD_store_ERIs"].getStr().empty() != true) {
        this->isStoreERIs_ = (*pPdfParam)["CD_store_ERIs"].getBoolean();
    }    
}

DfCD::~DfCD()
{
}

void DfCD::createEngines()
{
    assert(this->pEriEngines_ == NULL);
    
#ifdef _OPENMP
    {
        const int numOfThreads = omp_get_max_threads();
        this->pEriEngines_ = new DfEriEngine[numOfThreads];
    }
#else
    this->pEriEngines_ = new DfEriEngine[1];
#endif // _OPENMP
}


void DfCD::destroyEngines()
{
    if (this->pEriEngines_ != NULL) {
        delete[] this->pEriEngines_;
    }
    this->pEriEngines_ = NULL;
}

void DfCD::calcCholeskyVectors()
{
    this->calcCholeskyVectors_onTheFly();
}

void DfCD::saveI2PQ(const PQ_PairArray& I2PQ) 
{
    std::string filepath = this->getI2pqVtrPath();
    std::ofstream ofs;
    ofs.open(filepath.c_str(), std::ofstream::out | std::ofstream::binary);

    const std::size_t size = I2PQ.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(std::size_t));
    for (std::size_t i = 0; i < size; ++i) {
        const index_type index1 = I2PQ[i].index1();
        const index_type index2 = I2PQ[i].index2();
        ofs.write(reinterpret_cast<const char*>(&index1), sizeof(index_type));
        ofs.write(reinterpret_cast<const char*>(&index2), sizeof(index_type));
    }

    ofs.close();
}

DfCD::PQ_PairArray DfCD::getI2PQ()
{
    std::string filepath = this->getI2pqVtrPath();
    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
    if (ifs.fail()) {
        abort();
    }

    std::size_t size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    PQ_PairArray answer(size);
    index_type shellIndex1 = 0;
    index_type shellIndex2 = 0;
    for (std::size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&shellIndex1), sizeof(index_type));
        ifs.read(reinterpret_cast<char*>(&shellIndex2), sizeof(index_type));
        answer[i] = IndexPair2(shellIndex1, shellIndex2);
    }

    ifs.close();
    return answer;
}


void DfCD::saveL(const TlMatrix& L)
{
    DfObject::saveLMatrix(L);
}


TlMatrix DfCD::getL()
{
    TlMatrix L = DfObject::getLMatrix<TlMatrix>();
    return L;
}


DfTaskCtrl* DfCD::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    // pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    // pDfTaskCtrl->setCutoffEpsilon_density(0.0);
    // pDfTaskCtrl->setCutoffEpsilon_distribution(this->CDAM_tau_);

    return pDfTaskCtrl;
}

void DfCD::finalize(TlSymmetricMatrix* pMat)
{
    // do nothing
}

void DfCD::finalize(TlSparseSymmetricMatrix *pMat) 
{
    // do nothing
}

void DfCD::finalize_I2PQ(PQ_PairArray *pI2PQ)
{
    std::sort(pI2PQ->begin(), pI2PQ->end());
}


TlSymmetricMatrix DfCD::getCholeskyVector(const TlVector& L_col,
                                          const PQ_PairArray& I2PQ)
{
    const index_type numOfItilde = L_col.getSize();
    TlSymmetricMatrix answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(),
                   I2PQ[i].index2(),
                   L_col[i]);
    }

    return answer;
}

void DfCD::getJ(TlSymmetricMatrix* pJ)
{
    const TlSymmetricMatrix P = this->getPMatrix();

    // cholesky vector
    TlMatrix L = this->getL();
    const index_type numOfCBs = L.getNumOfCols();

    const PQ_PairArray I2PQ = this->getI2PQ();
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlSymmetricMatrix LI = this->getCholeskyVector(L.getColVector(I), I2PQ);
        assert(LI.getNumOfRows() == this->m_nNumOfAOs);
        assert(LI.getNumOfCols() == this->m_nNumOfAOs);
        
        TlMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();

        *pJ += qi*LI;
    }

    this->finalize(pJ);
}


void DfCD::divideCholeskyBasis(const index_type numOfCBs,
                               index_type *pStart, index_type *pEnd)
{
    *pStart = 0;
    *pEnd = numOfCBs;
}


void DfCD::getK(const RUN_TYPE runType,
                TlSymmetricMatrix *pK)
{
    TlMatrix L = this->getL();
    const index_type numOfCBs = L.getNumOfCols();
    
    TlSymmetricMatrix P = 0.5 * this->getPMatrix(); // RKS
    const TlMatrix C = P.choleskyFactorization2(this->epsilon_);
    
    const PQ_PairArray I2PQ = this->getI2PQ();
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlSymmetricMatrix l = this->getCholeskyVector(L.getColVector(I), I2PQ);
    
        TlMatrix X = l * C;
        TlMatrix Xt = X;
        Xt.transpose();
        
        TlSymmetricMatrix XX = X * Xt;
        *pK += XX;
    }
    
    *pK *= -1.0;
    this->finalize(pK);
}


TlSymmetricMatrix DfCD::getPMatrix()
{
    TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    return P;
}


void DfCD::calcCholeskyVectors_onTheFly()
{
    // timing data
    TlTime CD_all_time;
    TlTime CD_diagonals_time;
    TlTime CD_resizeL_time;
    TlTime CD_pivot_time;
    TlTime CD_ERI_time;
    // TlTime CD_Lpm_time;
    TlTime CD_calc_time;
    // TlTime CD_d_time;
    TlTime CD_save_time;

    CD_all_time.start();

    this->createEngines();
    this->initializeCutoffStats();

    this->log_.info(TlUtils::format("# of PQ dimension: %d", int(this->numOfPQs_)));
    TlSparseSymmetricMatrix schwartzTable(this->m_nNumOfAOs);
    PQ_PairArray I2PQ;
    TlVector d; // 対角成分

    CD_diagonals_time.start();
    this->calcDiagonals(&schwartzTable, &I2PQ, &d);
    CD_diagonals_time.stop();

    this->log_.info(TlUtils::format("# of I~ dimension: %d", int(I2PQ.size())));
    this->saveI2PQ(I2PQ);
    // this->ERI_cache_manager_.setMaxItems(I2PQ.size() * 2);

    const index_type N = I2PQ.size();
    double error = d.getMaxAbsoluteElement();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    // prepare variables
    const bool isUsingMemManager = this->isEnableMmap_;
    TlRowVectorMatrix2 L(N, 1, 1, 0, isUsingMemManager);
    const double threshold = this->epsilon_;
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));

    int progress = 0;
    index_type division =  index_type(N * 0.01);
    L.reserve_cols(division);
    index_type m = 0;
    while (error > threshold) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 8.3e", m, N, error));
#endif //DEBUG_CD

        CD_resizeL_time.start();
        // progress 
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e, local mem:%8.1f MB",
                                            m, error, TlSystem::getMaxRSS()));
            ++progress;

            // メモリの確保
            L.reserve_cols(division * progress);
        }
        L.resize(N, m+1);
        CD_resizeL_time.stop();

        // pivot
        CD_pivot_time.start();
        {
            std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                           pivot.end());
            const index_type i = it - pivot.begin();
            std::swap(pivot[m], pivot[i]);
        }
        CD_pivot_time.stop();

        error = d[pivot[m]];
        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(pivot[m], m, l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // ERI
        CD_ERI_time.start();
        const index_type pivot_m = pivot[m];
        std::vector<double> G_pm;
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type c = 0; c < numOf_G_cols; ++c) {
                const index_type pivot_i = pivot[m+1 +c]; // from (m+1) to N
                G_col_list[c] = pivot_i;
            }
            G_pm = this->getSuperMatrixElements(pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(G_pm.size() == numOf_G_cols);
        CD_ERI_time.stop();

        // CD calc
        CD_calc_time.start();
        const TlVector L_pm = L.getRowVector(pivot_m);
        std::vector<double> L_xm(numOf_G_cols);
#pragma omp parallel for schedule(runtime)
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            TlVector L_pi = L.getRowVector(pivot_i);
            const double sum_ll = (L_pi.dot(L_pm)).sum();
            const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

            //L.set(pivot_i, m, l_m_pi);
#pragma omp atomic
            L_xm[i] += l_m_pi; // for OpenMP
            
#pragma omp atomic
            d[pivot_i] -= l_m_pi * l_m_pi;
        }
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            L.set(pivot_i, m, L_xm[i]);
        }
        CD_calc_time.stop();

        ++m;
    }
    L.resize(N, m);
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->destroyEngines();
    this->schwartzCutoffReport();

    CD_save_time.start();
    this->saveL(L.getTlMatrix());
    CD_save_time.stop();

    CD_all_time.stop();

    // timing data
    this->log_.info(TlUtils::format("CD all:       %12.2f sec.", CD_all_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD diagonals: %12.2f sec.", CD_diagonals_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD resize L:  %12.2f sec.", CD_resizeL_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD pivot:     %12.2f sec.", CD_pivot_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD ERI:       %12.2f sec.", CD_ERI_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD calc:      %12.2f sec.", CD_calc_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD save:      %12.2f sec.", CD_save_time.getElapseTime()));
}


void DfCD::calcDiagonals(TlSparseSymmetricMatrix *pSchwartzTable,
                         PQ_PairArray *pI2PQ,
                         TlVector *pDiagonals)
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(numOfAOs == this->orbitalInfo_.getNumOfOrbitals());

    const double tau = this->CDAM_tau_;
    this->log_.info(TlUtils::format("CDAM tau: %e", tau));
    this->log_.info(TlUtils::format("primitive GTO quartet threshold: %e", this->cutoffEpsilon3_));

    // initialize
    pI2PQ->clear();
    pI2PQ->reserve(this->numOfPQs_);
    pSchwartzTable->clear();
    pSchwartzTable->resize(numOfAOs);
    TlSparseSymmetricMatrix diagonalMat(numOfAOs);

    // task
    this->log_.info("diagonal calculation: start");
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(this->orbitalInfo_,
                                          true,
                                          this->grainSize_,
                                          &taskList, true);
    while (hasTask == true) {
        this->calcDiagonals_kernel(taskList,
                                   pSchwartzTable,
                                   &diagonalMat, pI2PQ);
        hasTask = pDfTaskCtrl->getQueue2(this->orbitalInfo_,
                                         true,
                                         this->grainSize_,
                                         &taskList);
    }
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    
    // finalize
    this->log_.info("diagonal calculation: finalize");
    this->finalize_I2PQ(pI2PQ);
    this->finalize(&diagonalMat);
    this->finalize(pSchwartzTable);

    // set diagonals
    const index_type numOfI = pI2PQ->size();
    pDiagonals->resize(numOfI);
    for (index_type i = 0; i < numOfI; ++i) {
        const index_type row = (*pI2PQ)[i].index1();
        const index_type col = (*pI2PQ)[i].index2();
        const double value = diagonalMat.get(row, col);
        (*pDiagonals)[i] = value;
    }
}


void DfCD::calcDiagonals_kernel(const std::vector<DfTaskCtrl::Task2>& taskList,
                                TlSparseSymmetricMatrix *pSchwartzTable,
                                TlSparseSymmetricMatrix *pDiagonalMat,
                                PQ_PairArray *pI2PQ)
{
    const index_type numOfAOs = this->orbitalInfo_.getNumOfOrbitals();
    pDiagonalMat->resize(numOfAOs);

    const double tau = this->CDAM_tau_;
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        PQ_PairArray local_I2PQ;
        TlSparseSymmetricMatrix local_diagMat(numOfAOs);
        TlSparseSymmetricMatrix local_schwartzTable(numOfAOs);
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = this->orbitalInfo_.getShellType(shellIndexP);
            const int shellTypeQ = this->orbitalInfo_.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;

            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::CGTO_Pair PQ = 
                this->pEriEngines_[threadID].getCGTO_pair(this->orbitalInfo_,
                                                          shellIndexP,
                                                          shellIndexQ,
                                                          pairwisePGTO_cutoffThreshold);
            this->pEriEngines_[threadID].calc(queryPQ, queryPQ, PQ, PQ);
                
            const int maxStepsPQ = maxStepsP * maxStepsQ;
            double maxValue = 0.0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                for (int q = 0; q < maxStepsQ; ++q) {
                    const index_type indexQ = shellIndexQ + q;
                    
                    if ((shellIndexP != shellIndexQ) || (indexP >= indexQ)) {
                        const int pq_index = p * maxStepsQ + q;
                        const int pqpq_index = pq_index * maxStepsPQ + pq_index;
                        
                        const double value = this->pEriEngines_[threadID].WORK[pqpq_index];

                        // for schwartz
                        maxValue = std::max(maxValue, std::fabs(value));
                        
                        // for I~ to pq table
                        if (std::fabs(value) > tau) {
                            if (value > 0) {
                                local_diagMat.set(indexP, indexQ, value);
                                local_I2PQ.push_back(IndexPair2(indexP, indexQ));
                            } else {
                                this->log_.warn(TlUtils::format("pqpq_value: (%d %d)=% e is spoiled.",
                                                                indexP, indexQ, value));
                            }
                        }
                    }
                }
            }
            local_schwartzTable.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
        }

        // add up
#ifdef _OPENMP
        {
            const int numOfThreads = omp_get_num_threads();
            for (int i = 0; i < numOfThreads; ++i) {
                if (threadID == i) {
                    pI2PQ->insert(pI2PQ->end(),
                                  local_I2PQ.begin(), local_I2PQ.end());
                    pDiagonalMat->merge(local_diagMat);
                    pSchwartzTable->merge(local_schwartzTable);
                }
#pragma omp barrier                
            }
        }
#else
        {
            *pI2PQ = local_I2PQ;
            *pDiagonalMat = local_diagMat;
            *pSchwartzTable = local_schwartzTable;
        }
#endif // _OPENMP
    }
}


bool DfCD::isAliveBySchwartzCutoff(const index_type shellIndexP,
                                   const index_type shellIndexQ,
                                   const index_type shellIndexR,
                                   const index_type shellIndexS,
                                   const int shellQuartetType,
                                   const TlSparseSymmetricMatrix& schwarzTable,
                                   const double threshold)
{
    bool answer = false;

    const double sqrt_pqpq = schwarzTable.get(shellIndexP, shellIndexQ);
    const double sqrt_rsrs = schwarzTable.get(shellIndexR, shellIndexS);

    if ((sqrt_pqpq * sqrt_rsrs) >= threshold) {
        answer = true;

#pragma omp atomic
        ++(this->cutoffAlive_schwartz_[shellQuartetType]);
    }

#pragma omp atomic
    ++(this->cutoffAll_schwartz_[shellQuartetType]);

    return answer;
}


void DfCD::initializeCutoffStats()
{
    // clear cutoff stats
    const int maxShellType = this->orbitalInfo_.getMaxShellType();
    const int numOfShellPairType = maxShellType* maxShellType;
    const int numOfShellQuartetType = numOfShellPairType * numOfShellPairType;
    this->cutoffAll_schwartz_.clear();
    this->cutoffAlive_schwartz_.clear();
    this->cutoffAll_schwartz_.resize(numOfShellQuartetType, 0);
    this->cutoffAlive_schwartz_.resize(numOfShellQuartetType, 0);
}


void DfCD::schwartzCutoffReport()
{
    static const char typeStr4[][5] = {
        "SSSS", "SSSP", "SSSD", "SSPS", "SSPP", "SSPD", "SSDS", "SSDP", "SSDD",
        "SPSS", "SPSP", "SPSD", "SPPS", "SPPP", "SPPD", "SPDS", "SPDP", "SPDD",
        "SDSS", "SDSP", "SDSD", "SDPS", "SDPP", "SDPD", "SDDS", "SDDP", "SDDD",
        "PSSS", "PSSP", "PSSD", "PSPS", "PSPP", "PSPD", "PSDS", "PSDP", "PSDD",
        "PPSS", "PPSP", "PPSD", "PPPS", "PPPP", "PPPD", "PPDS", "PPDP", "PPDD",
        "PDSS", "PDSP", "PDSD", "PDPS", "PDPP", "PDPD", "PDDS", "PDDP", "PDDD",
        "DSSS", "DSSP", "DSSD", "DSPS", "DSPP", "DSPD", "DSDS", "DSDP", "DSDD",
        "DPSS", "DPSP", "DPSD", "DPPS", "DPPP", "DPPD", "DPDS", "DPDP", "DPDD",
        "DDSS", "DDSP", "DDSD", "DDPS", "DDPP", "DDPD", "DDDS", "DDDP", "DDDD",
    };
    const int maxShellType = this->orbitalInfo_.getMaxShellType();

    // cutoff report for schwarz
    bool hasCutoffSchwarz = false;
    for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
        for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
            const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
            for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                    const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                    if (this->cutoffAll_schwartz_[shellTypeABCD] != 0) {
                        hasCutoffSchwarz = true;
                        break;
                    }
                }
            }
        }
    }
    if (hasCutoffSchwarz == true) {
        this->log_.info("schwarz cutoff report");
        this->log_.info(TlUtils::format("threshold: %e", this->CDAM_tau_));
        this->log_.info("type: alive / all (ratio)");
        for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
            for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
                const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
                for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                    const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                    for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                        const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                        
                        if (this->cutoffAll_schwartz_[shellTypeABCD] > 0) {
                            const double ratio = (double)this->cutoffAlive_schwartz_[shellTypeABCD]
                                / (double)this->cutoffAll_schwartz_[shellTypeABCD]
                                * 100.0;
                            this->log_.info(TlUtils::format(" %4s: %12ld / %12ld (%6.2f%%)",
                                                            typeStr4[shellTypeABCD],
                                                            this->cutoffAlive_schwartz_[shellTypeABCD],
                                                            this->cutoffAll_schwartz_[shellTypeABCD],
                                                            ratio));
                        }
                    }
                }
            }
        }
    }
}


std::vector<double>
DfCD::getSuperMatrixElements(const index_type G_row,
                             const std::vector<index_type>& G_col_list,
                             const PQ_PairArray& I2PQ,
                             const TlSparseSymmetricMatrix& schwartzTable)
{
    this->ERI_cache_.clear();

    const std::vector<IndexPair4> calcList = this->getCalcList(G_row, G_col_list, I2PQ);
    this->calcERIs(calcList, schwartzTable);
    const std::vector<double> answer = this->setERIs(G_row, G_col_list, I2PQ);

    return answer;
}


std::vector<DfCD::IndexPair4> 
DfCD::getCalcList(const index_type G_row,
                  const std::vector<index_type>& G_col_list,
                  const PQ_PairArray& I2PQ)
{
    std::set<IndexPair4> calcSet;

    const index_type indexP = I2PQ[G_row].index1();
    const index_type indexQ = I2PQ[G_row].index2();
    const index_type shellIndexP = this->orbitalInfo_.getShellIndex(indexP);
    const index_type shellIndexQ = this->orbitalInfo_.getShellIndex(indexQ);

    const index_type numOf_G_cols = G_col_list.size();
#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < numOf_G_cols; ++i) {
        const index_type G_col = G_col_list[i];

        const index_type indexR = I2PQ[G_col].index1();
        const index_type indexS = I2PQ[G_col].index2();
        const index_type shellIndexR = this->orbitalInfo_.getShellIndex(indexR);
        const index_type shellIndexS = this->orbitalInfo_.getShellIndex(indexS);

        IndexPair4 indexPair4(shellIndexP, shellIndexQ, shellIndexR, shellIndexS);
#pragma omp critical(DfCD__getCalcList) 
        {
            calcSet.insert(indexPair4);
        }
    }

    std::vector<IndexPair4> calcList(calcSet.size());
    std::copy(calcSet.begin(), calcSet.end(), calcList.begin());

    return calcList;
}


void DfCD::calcERIs(const std::vector<IndexPair4>& calcList,
                    const TlSparseSymmetricMatrix& schwartzTable) 
{
    const int maxShellType = this->orbitalInfo_.getMaxShellType();
    const double threshold = this->CDAM_tau_;
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

    const int numOfList = calcList.size();
#pragma omp parallel
    {
        int threadID = 0;
        ERI_CacheType local_cache;

#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < numOfList; ++i) {
            const index_type shellIndexP = calcList[i].index1();
            const index_type shellIndexQ = calcList[i].index2();
            const index_type shellIndexR = calcList[i].index3();
            const index_type shellIndexS = calcList[i].index4();
            
            const int shellTypeP = this->orbitalInfo_.getShellType(shellIndexP);
            const int shellTypeQ = this->orbitalInfo_.getShellType(shellIndexQ);
            const int shellTypeR = this->orbitalInfo_.getShellType(shellIndexR);
            const int shellTypeS = this->orbitalInfo_.getShellType(shellIndexS);
            
            const int shellQuartetType =
                ((shellTypeP * maxShellType + shellTypeQ) * maxShellType + shellTypeP) * maxShellType + shellTypeQ;
            const bool isAlive = this->isAliveBySchwartzCutoff(shellIndexP, shellIndexQ,
                                                               shellIndexR, shellIndexS,
                                                               shellQuartetType,
                                                               schwartzTable,
                                                               threshold);
            if (isAlive == true) {
                const int maxStepsP = 2 * shellTypeP + 1;
                const int maxStepsQ = 2 * shellTypeQ + 1;
                const int maxStepsR = 2 * shellTypeR + 1;
                const int maxStepsS = 2 * shellTypeS + 1;
                
                const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(this->orbitalInfo_,
                                                                                            shellIndexP,
                                                                                            shellIndexQ,
                                                                                            pairwisePGTO_cutoffThreshold);
                const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(this->orbitalInfo_,
                                                                                            shellIndexR,
                                                                                            shellIndexS,
                                                                                            pairwisePGTO_cutoffThreshold);
                const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
                const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
                
                this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
                
                const int steps = maxStepsP * maxStepsQ * maxStepsR * maxStepsS;
                std::vector<double> buf(steps);
                std::copy(this->pEriEngines_[threadID].WORK, this->pEriEngines_[threadID].WORK + steps,
                          buf.begin());

                local_cache[calcList[i]] = buf;
            }
        }

        // merge cache
#pragma omp critical(DfCD__calcERIs)
        {
            this->ERI_cache_.insert(local_cache.begin(), local_cache.end());
        }
    }
}


std::vector<double>
DfCD::setERIs(const index_type G_row,
              const std::vector<index_type> G_col_list,
              const PQ_PairArray& I2PQ)
{
    const index_type indexP_orig = I2PQ[G_row].index1();
    const index_type indexQ_orig = I2PQ[G_row].index2();
    const index_type shellIndexP_orig = this->orbitalInfo_.getShellIndex(indexP_orig);
    const index_type shellIndexQ_orig = this->orbitalInfo_.getShellIndex(indexQ_orig);

    const index_type numOf_G_cols = G_col_list.size();
    std::vector<double> answer(numOf_G_cols);
#pragma omp parallel for
    for (index_type i = 0; i < numOf_G_cols; ++i) {
        index_type indexP = indexP_orig;
        index_type indexQ = indexQ_orig;
        index_type shellIndexP = shellIndexP_orig;
        index_type shellIndexQ = shellIndexQ_orig;

        const index_type G_col = G_col_list[i];

        index_type indexR = I2PQ[G_col].index1();
        index_type indexS = I2PQ[G_col].index2();
        index_type shellIndexR = this->orbitalInfo_.getShellIndex(indexR);
        index_type shellIndexS = this->orbitalInfo_.getShellIndex(indexS);

        // swap
        if (shellIndexP > shellIndexQ) {
            std::swap(shellIndexP, shellIndexQ);
            std::swap(indexP, indexQ);
        }
        if (shellIndexR > shellIndexS) {
            std::swap(shellIndexR, shellIndexS);
            std::swap(indexR, indexS);
        }
        {
            IndexPair2 pq(shellIndexP, shellIndexQ);
            IndexPair2 rs(shellIndexR, shellIndexS);
            if (pq.index() > rs.index()) {
                std::swap(shellIndexP, shellIndexR);
                std::swap(indexP, indexR);
                std::swap(shellIndexQ, shellIndexS);
                std::swap(indexQ, indexS);
            }
        }

        std::vector<double> values;
#pragma omp critical(DfCD__setERIs)
        {
            values = this->ERI_cache_[IndexPair4(shellIndexP, shellIndexQ,
                                                 shellIndexR, shellIndexS)];
        }

        if (values.empty() != true) {
            const int basisTypeP = indexP - shellIndexP;
            const int basisTypeQ = indexQ - shellIndexQ;
            const int basisTypeR = indexR - shellIndexR;
            const int basisTypeS = indexS - shellIndexS;
            
            const int shellTypeP = this->orbitalInfo_.getShellType(shellIndexP);
            const int shellTypeQ = this->orbitalInfo_.getShellType(shellIndexQ);
            const int shellTypeR = this->orbitalInfo_.getShellType(shellIndexR);
            const int shellTypeS = this->orbitalInfo_.getShellType(shellIndexS);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            const int maxStepsS = 2 * shellTypeS + 1;
            
            const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
            assert(values.size() > index);
            
#pragma omp atomic
            answer[i] += values.at(index);
        }
    }

    return answer;
}

