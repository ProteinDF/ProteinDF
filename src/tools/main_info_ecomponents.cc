#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"
#include "TlUtils.h"

#include "df_total_energy_eigen.h"

std::vector<int> atom2AO(const TlOrbitalInfo& orbInfo,
                         const std::vector<int>& atoms) {
  std::vector<int> AOs;
  const int numOfAOs = orbInfo.getNumOfOrbitals();
  for (int ao = 0; ao < numOfAOs; ++ao) {
    const int atomIndex = orbInfo.getAtomIndex(ao);
    
    std::vector<int>::const_iterator itEnd = atoms.end();
    for (std::vector<int>::const_iterator it = atoms.begin(); it != itEnd; ++it) {
      if (*it == atomIndex) {
        AOs.push_back(ao);
        break;
      }
    }
  }
  return AOs;
}

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
  std::vector<std::vector<int> > groups_atom;
  for (int i = 1; i < numOfArgs; ++i) {  // opt[0] is program name
    const std::string input = opt[i];
    std::cout << TlUtils::format("input %d: %s", i, input.c_str()) << std::endl;
    const std::vector<int> inputArray =
        TlUtils::nonreduntant_vector(TlUtils::vector_notation(input));
    groups_atom.push_back(inputArray);
  }
  if (groups_atom.size() == 0) {
    DfObject dfObj(&param);
    int numOfAtoms = dfObj.getNumOfAtoms();
    const std::string allAtoms = TlUtils::format("0-%d", numOfAtoms - 1);
    const std::vector<int> allAtomsArray =
        TlUtils::nonreduntant_vector(TlUtils::vector_notation(allAtoms));
    groups_atom.push_back(allAtomsArray);
  }
  if (verbose) {
    const int numOfGroups = groups_atom.size();
    for (int i = 0; i < numOfGroups; ++i) {
      std::cerr << TlUtils::format("group[%d] (atom): %s", i,
                                   TlUtils::vector2str(groups_atom[i]).c_str())
                << std::endl;
    }
  }

  // TlOrbitalInfo
  // const TlOrbitalInfo orbInfo(param["coordinates"], param["basis_set"]);

  // atom -> AO
  // std::vector<std::vector<int> > groups_AO;
  // std::vector<std::vector<int> >::const_iterator itEnd = groups_atom.end();
  // for (std::vector<std::vector<int> >::const_iterator it = groups_atom.begin(); it != itEnd; ++it) {
  //   std::vector<int> AOs = atom2AO(orbInfo, *it);
  //   groups_AO.push_back(AOs);
  // }
  // if (verbose) {
  //   const int numOfGroups = groups_AO.size();
  //   for (int i = 0; i < numOfGroups; ++i) {
  //     std::cerr << TlUtils::format("group[%d] (AO): %s", i,
  //                                  TlUtils::vector2str(groups_AO[i]).c_str())
  //               << std::endl;
  //   }
  // }

  //
  DfTotalEnergy_Eigen dfTotalEnergy(&param);
  dfTotalEnergy.calc(iteration);
  dfTotalEnergy.output();

  //
  double IE_total = 0.0;
  const int numOfGroups = groups_atom.size();
  for (int i = 0; i < numOfGroups; ++i) {
    {
      std::cout << TlUtils::format(">>>> group: %d <-> group: %d", i, i)
                << std::endl;
      double IE = dfTotalEnergy.get_IE(groups_atom[i]);
      IE_total += IE;
    }
    for (int j = i +1; j < numOfGroups; ++j) {
      std::cout << TlUtils::format(">>>> group: %d <-> group: %d", i, j)
                << std::endl;
      // std::cout << TlUtils::format("group: %d [%s]", i,
      //                              TlUtils::vector2str(groups_atom[i]).c_str())
      //           << std::endl;
      // std::cout << TlUtils::format("group: %d [%s]", j,
      //                              TlUtils::vector2str(groups_atom[j]).c_str())
      //           << std::endl;
      double IE = dfTotalEnergy.get_IE(groups_atom[i], groups_atom[j]);
      IE_total += IE;
    }
  }
  std::cout << TlUtils::format("IE (total) = %lf", IE_total) << std::endl;

  return EXIT_SUCCESS;
}
