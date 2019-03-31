#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include "TlGetopt.h"
#include "TlSerializeData.h"
#include "TlMsgPack.h"
#include "TlUtils.h"

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hb:p:");

    std::string pdfParamPath = "pdfparam.mpac";
    if (! opt["p"].empty()) {
        pdfParamPath = opt["p"];
    }

    std::string orbInput = "";
    if (! opt["b"].empty()) {
        const std::string opt_b = opt["b"];
        orbInput = opt_b;
    }

  // パラメータファイルの読み込み
  TlSerializeData param;
  {
    TlMsgPack mpac;
    mpac.load(pdfParamPath);
    param = mpac.getSerializeData();
  }


    std::vector<int> orb_list = TlUtils::vector_notation(orbInput);
    std::vector<int>::const_iterator itEnd = orb_list.end();
    for (std::vector<int>::const_iterator it = orb_list.begin(); it != itEnd; ++it) {
        std::cout << *it << std::endl;
    }

    return EXIT_SUCCESS;
}
