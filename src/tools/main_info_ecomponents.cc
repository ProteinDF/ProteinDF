#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlUtils.h"

#include "df_total_energy_eigen.h"

int main(int argc, char* argv[]) {
  // parse args
  TlGetopt opt(argc, argv, "hi:p:v");
  const bool verbose = (opt["v"] == "defined");

  std::string pdfParamPath = "pdfparam.mpac";
  if (!opt["p"].empty()) {
    pdfParamPath = opt["p"];
  }

  // パラメータファイルの読み込み
  TlSerializeData param;
  {
    TlMsgPack mpac;
    mpac.load(pdfParamPath);
    param = mpac.getSerializeData();
  }

  //
  int iteration = -1;
  {
    if (!opt["i"].empty()) {
      iteration = std::atoi(opt["i"].c_str());
    } else {
      DfObject dfObj(&param);
      iteration = dfObj.iteration();
    }
  }
  if (verbose) {
    std::cerr << TlUtils::format("iteration: %d", iteration) << std::endl;
  }

  // parse args
  const int numOfArgs = opt.getCount();
  std::vector<std::vector<int> > group;
  for (int i = 1; i < numOfArgs; ++i) {
    const std::string input = opt[i];
    const std::vector<int> inputArray =
        TlUtils::nonreduntant_vector(TlUtils::vector_notation(input));
    if (verbose) {
      std::cerr << TlUtils::format("group[%d]: %s", i, input.c_str())
                << std::endl;
      std::cerr << TlUtils::vector2str(inputArray) << std::endl;
    }
    group.push_back(inputArray);
  }

  DfTotalEnergy_Eigen dfTotalEnergy(&param);
  dfTotalEnergy.calc(iteration);
  dfTotalEnergy.output();

  return EXIT_SUCCESS;
}
