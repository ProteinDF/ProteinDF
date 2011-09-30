#include "DfDiffDensityMatrix.h"
#include "TlFile.h"

DfDiffDensityMatrix::DfDiffDensityMatrix(TlSerializeData* pPdfData)
  : DfObject(pPdfData)
{
    this->isSaveDiffMatrix_ = (TlUtils::toUpper((*pPdfData)["save_diff_density_matrix"].getStr()) == "YES");
}


DfDiffDensityMatrix::~DfDiffDensityMatrix()
{
}


void DfDiffDensityMatrix::exec()
{
    // check memory
    const std::size_t needMem = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) * sizeof(double);
    if ((this->isWorkOnDisk_ == true) ||
        (this->procMaxMemSize_ < needMem)) {
        this->logger(" The differencial density matrix is build on disk.\n");
        TlMatrix::useMemManager(true);
    } else {
        this->logger(" The differencial density matrix is build on memory.\n");
        TlMatrix::useMemManager(false);
    }

    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->calc(RUN_RKS, this->m_nIteration);
        break;
        
    case METHOD_UKS:
        this->calc(RUN_UKS_ALPHA, this->m_nIteration);
        this->calc(RUN_UKS_BETA, this->m_nIteration);
        break;
        
    case METHOD_ROKS:
        this->calc_ROKS<TlSymmetricMatrix>();
        abort();
        break;
        
    default:
        std::cerr << "program error. " << __FILE__ << ":" << __LINE__ << std::endl;
        abort();
        break;
    }
}


void DfDiffDensityMatrix::calc(const DfObject::RUN_TYPE runType, const int iteration)
{
    TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration -1);
    if (TlFile::isExist(this->getPpqMatrixPath(runType, iteration -2)) == true) {
        P -= (this->getPpqMatrix<TlSymmetricMatrix>(runType, iteration -2));
    }

    this->saveDiffDensityMatrix(runType, iteration, P);
}


