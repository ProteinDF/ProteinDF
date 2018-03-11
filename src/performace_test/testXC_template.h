#ifndef TEST_XC_TEMPLATE_H
#define TEST_XC_TEMPLATE_H

#include <string>
#include "DfFunctional_B3LYP.h"
#include "TlSerializeData.h"
#include "TlUtils.h"

template <class SymmetricMatrixType, class DfCalcGridClass>
void testXC(const TlSerializeData& inParam,
            const std::string& densityMatrixPath,
            const std::string& KMatrixPath, bool verbose = false) {
  TlSerializeData param = inParam;

  TlOrbitalInfo orbitalInfo(param["coordinates"], param["basis_set"]);
  const int numOfAOs = orbitalInfo.getNumOfOrbitals();
  if (verbose == true) {
    std::cerr << TlUtils::format("number of AOs: %d", numOfAOs) << std::endl;
  }

  SymmetricMatrixType P;
  if (verbose == true) {
    std::cerr << TlUtils::format("density matrix: %s",
                                 densityMatrixPath.c_str())
              << std::endl;
  }
  P.load(densityMatrixPath);
  if (verbose == true) {
    std::cerr << TlUtils::format("dimensions of density matrix: %d",
                                 P.getNumOfRows())
              << std::endl;
  }
  assert(P.getNumOfRows() == numOfAOs);

  // calc
  DfCalcGridClass dfCalcGrid(&param);
  DfFunctional_B3LYP b3lyp;
  SymmetricMatrixType K(numOfAOs);
  dfCalcGrid.calcXCIntegForFockAndEnergy(P, SymmetricMatrixType(), &b3lyp, &K,
                                         NULL);

  // save
  K.save(KMatrixPath);
}

#endif  // TEST_XC_TEMPLATE_H
