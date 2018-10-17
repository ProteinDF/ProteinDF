#include <iostream>
#include "tl_mint_object.h"
#include "Fl_Geometry.h"
#include "TlOrbitalInfo.h"
#include "TlUtils.h"

TlMintObject::TlMintObject(const Fl_Geometry& geom, const TlOrbitalInfo& orbInfo) : geom_(geom), orbInfo_(orbInfo) {

}
TlMintObject::~TlMintObject() {}

bool TlMintObject::checkAO(int aoIndex) const {
    bool answer = false;
    int type = this->orbInfo_.getBasisType(aoIndex);
    if ((type == 0) || (type == 1) || (type == 4) || (type == 9)) {
      answer = true;
    } else {
      std::cerr << TlUtils::format(
                       "The AO index, %d, is not the beginning of shells.",
                       aoIndex)
                << std::endl;
    }

    return answer;
  }


