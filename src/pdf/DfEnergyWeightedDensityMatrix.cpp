#include "DfEnergyWeightedDensityMatrix.h"

DfEnergyWeightedDensityMatrix::DfEnergyWeightedDensityMatrix()
    : DfObject() {
}


DfEnergyWeightedDensityMatrix::~DfEnergyWeightedDensityMatrix()
{
}


void DfEnergyWeightedDensityMatrix::calc(RUN_TYPE runType)
{
    const int itr = this->numOfIterations_;
    const int numOfMOs = this->getNumOfMOs_;

    TlVector eps;
    eps.load(this->getOccupationPath(runType));
    assert(occ.getSize() == numOfMOs);
    {
        TlVector eig;
        eig.load(this->getEigvalPath(runType, itr));
        assert(eig.getSize() == numOfMOs);

        eps.dot(eig);
    }

    TlMatrix C = this->getCMatrix(runType, itr);
    TlMatrix Ct = C;
    Ct.transpose();
    
    TlMatrix W = Ct * eps * C;
    return W;
}



