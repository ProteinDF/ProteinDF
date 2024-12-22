#ifndef DF_CDK_MATRIX_PARALLEL_H
#define DF_CDK_MATRIX_PARALLEL_H

#include "TlCommunicate.h"
#include "df_cdk_matrix.h"

class DfCdkMatrix_Parallel : public DfCdkMatrix {
public:
    DfCdkMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfCdkMatrix_Parallel();

public:
    virtual void getK();

protected:
    template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector, typename TempSymmetricMatrix>
    void getK_byLjk_useDenseMatrix(const RUN_TYPE runType);

    // ---------------------------------------------------------------------------
    // task control
    // ---------------------------------------------------------------------------
protected:
    DfTaskCtrl* getDfTaskCtrlObject() const;

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
protected:
    // virtual void finalize(TlDenseSymmetricMatrixObject* pK);
};

// template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector, typename TempSymmetricMatrix>
// void DfCdkMatrix_Parallel::getK_byLjk_useDenseMatrix_master(const RUN_TYPE runType) {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     assert(rComm.isMaster() == true);

//     this->log_.info("calc K by CD method");
//     this->log_.info("transformed by dens matrix");

//     const Ljk_MatrixType L(DfObject::getLjkMatrixPath());
//     this->log_.info(TlUtils::format("L(K): %d x %d", L.getNumOfRows(), L.getNumOfCols()));
//     index_type calcNumOfCVs = L.getNumOfCols();

//     TempSymmetricMatrix P;
//     if ((this->isUpdateMethod_ == true) && (this->m_nIteration > 1)) {
//         this->log_.info("update method is applied.");
//         P = DfObject::getDiffDensityMatrix<SymmetricMatrix>(runType, this->m_nIteration);

//         double threshold = 1.0E-10;
//         const double diffP = (*this->pPdfParam_)["DfDiffDensityMatrix"]["max_abs"][this->m_nIteration - 1].getDouble();
//         this->log_.info(TlUtils::format("dP max: %e", diffP));
//         if (diffP < 1.0) {
//             threshold /= diffP;
//             threshold = std::sqrt(threshold);
//         }
//         this->log_.info(TlUtils::format("threshold: %e", threshold));

//         std::vector<double> errors;
//         {
//             TlDenseVector_Lapack errs_vtr = DfObject::loadLjkErrorsVector<TlDenseVector_Lapack>();
//             errors = errs_vtr;
//         }
//         std::vector<double>::iterator it = std::find_if(errors.begin(), errors.end(),
//                                                         std::bind(std::less<double>(), std::placeholders::_1, threshold));
//         calcNumOfCVs = static_cast<index_type>(std::distance(errors.begin(), it));

//     } else {
//         P = this->getPInMatrix<SymmetricMatrix>(runType, this->m_nIteration);
//     }
//     if (runType == RUN_RKS) {
//         P *= 0.5;
//     }
//     rComm.broadcast(P);

//     const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());

//     int numOfThreads = 1;
// #ifdef _OPENMP
//     numOfThreads = omp_get_num_procs();
// #endif  // _OPENMP
//     const int taskSize = 100 * numOfThreads;
//     this->log_.info(TlUtils::format("task size: %d", taskSize));

//     DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
//     std::vector<std::size_t> tasks;
//     bool hasTasks = pTaskCtrl->getQueue(calcNumOfCVs, taskSize, &tasks, true);

//     TempSymmetricMatrix K(this->m_nNumOfAOs);
//     while (hasTasks == true) {
//         const int numOfTasks = tasks.size();
//         // #pragma omp parallel for schedule(runtime)
//         for (int i = 0; i < numOfTasks; ++i) {
//             const TempSymmetricMatrix l = this->getCholeskyVector<SymmetricMatrix, Vector>(L.getColVector(tasks[i]), I2PQ);
//             // this->log_.info(TlUtils::format("l: %s", TlUtils::getTypeName(l).c_str()));
//             assert(l.getNumOfRows() == this->m_nNumOfAOs);

//             // #pragma omp critical
//             { K += l * P * l; }
//         }
//         hasTasks = pTaskCtrl->getQueue(calcNumOfCVs, taskSize, &tasks);
//     }

//     K *= -1.0;

//     SymmetricMatrix Kout = K;
//     if ((this->isUpdateMethod_ == true) && (this->m_nIteration > 1)) {
//         SymmetricMatrix prevK = this->getHFxMatrix<SymmetricMatrix>(runType, this->m_nIteration - 1);
//         Kout += prevK;
//     }

//     this->log_.info("finalize");
//     // this->finalize(&K);

//     delete pTaskCtrl;
//     pTaskCtrl = NULL;

//     this->saveHFxMatrix(runType, this->m_nIteration, Kout);
// }

// template <typename Ljk_MatrixType, typename SymmetricMatrix, typename Vector, typename TempSymmetricMatrix>
// void DfCdkMatrix_Parallel::getK_byLjk_useDenseMatrix_worker(const RUN_TYPE runType) {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     assert(rComm.isMaster() != true);

//     TempSymmetricMatrix P;
//     rComm.broadcast(P);

//     const PQ_PairArray I2PQ = this->getI2PQ(this->getI2pqVtrPath());

//     int numOfThreads = 1;
// #ifdef _OPENMP
//     numOfThreads = omp_get_num_procs();
// #endif  // _OPENMP
//     const int taskSize = 100 * numOfThreads;
//     this->log_.info(TlUtils::format("task size: %d", taskSize));

//     DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();
//     std::vector<std::size_t> tasks;
//     bool hasTasks = pTaskCtrl->getQueue(calcNumOfCVs, taskSize, &tasks, true);

//     TempSymmetricMatrix K(this->m_nNumOfAOs);
//     while (hasTasks == true) {
//         const int numOfTasks = tasks.size();
//         // #pragma omp parallel for schedule(runtime)
//         for (int i = 0; i < numOfTasks; ++i) {
//             const TempSymmetricMatrix l = this->getCholeskyVector<SymmetricMatrix, Vector>(L.getColVector(tasks[i]), I2PQ);
//             // this->log_.info(TlUtils::format("l: %s", TlUtils::getTypeName(l).c_str()));
//             assert(l.getNumOfRows() == this->m_nNumOfAOs);

//             // #pragma omp critical
//             { K += l * P * l; }
//         }
//         hasTasks = pTaskCtrl->getQueue(calcNumOfCVs, taskSize, &tasks);
//     }

//     K *= -1.0;

//     SymmetricMatrix Kout = K;
//     if ((this->isUpdateMethod_ == true) && (this->m_nIteration > 1)) {
//         SymmetricMatrix prevK = this->getHFxMatrix<SymmetricMatrix>(runType, this->m_nIteration - 1);
//         Kout += prevK;
//     }

//     this->log_.info("finalize");
//     // this->finalize(&K);

//     delete pTaskCtrl;
//     pTaskCtrl = NULL;

//     this->saveHFxMatrix(runType, this->m_nIteration, Kout);
// }


#endif  // DF_CDK_MATRIX_PARALLEL_H
