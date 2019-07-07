#ifndef TL_SYSTEM_EIGEN_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <string>

#ifdef HAVE_EIGEN
#include <Eigen/Core>

class TlSystem_Eigen {
   public:
    static std::string getSimdInstructionSetsInUse() {
        return Eigen::SimdInstructionSetsInUse();
    }
};
#endif //HAVE_EIGEN

#endif  // TL_SYSTEM_EIGEN_H
