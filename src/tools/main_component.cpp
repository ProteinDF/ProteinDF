// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <iostream>
#include <vector>

#include "DfObject.h"
#include "TlFile.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"
#include "TlSymmetricMatrix.h"
#include "TlUtils.h"

#define AU_BOHR 1.889762

struct LCAO {
 public:
  int index;
  double value;
};

struct LCAO_cmp {
 public:
  bool operator()(const LCAO& a, const LCAO& b) {
    return (std::fabs(a.value) > std::fabs(b.value));
  }
};

enum EXEC_MODE { UNDEFINED, MO_MODE, ATOM_MODE, SAVE_MODE };

void showHelp(const std::string& name) {
  std::cout << TlUtils::format("%s <OPTIONS> -M MO_LEVEL NUM_OF_ITEMS",
                               name.c_str())
            << std::endl;
  std::cout << " or " << std::endl;
  std::cout << TlUtils::format("%s <OPTIONS> -A ATOM_INDEX", name.c_str())
            << std::endl;
  std::cout << std::endl;
  std::cout << "OPTIONS:" << std::endl;
  std::cout
      << "  -d DIRECTORY    set directory to read (default: current directory)"
      << std::endl;
  std::cout << "  -c path         specify LCAO matrix (default: guessed)"
            << std::endl;
  std::cout << "  -s path         specify overlap matrix (default: guessed)"
            << std::endl;
  std::cout << "  -h:             show help(this)" << std::endl;
  std::cout << "  -v:             enable verbose mode" << std::endl;
  std::cout << "  <MO mode>" << std::endl;
  std::cout << "  -M MO_LEVEL     show AO components for the MO level"
            << std::endl;
  std::cout << "  -a:             output all basis functions" << std::endl;
  std::cout << "  -n NUM          output number of items" << std::endl;
  std::cout << "  <atom mode>" << std::endl;
  std::cout << "  -A ATOM_INDEX   output MO levels in which AOs are major"
            << std::endl;
  std::cout << "  -t THRESHOLD    threshold of w_atom" << std::endl;
  std::cout << "  <save mode>" << std::endl;
  std::cout << "  -S path         save intermediate files" << std::endl;
}

/// MO モード
/// topに-1が指定された場合は、全軌道が出力される
int exec_MO_mode(const TlMatrix& SC, const TlMatrix& C,
                 const TlOrbitalInfo& readOrbInfo, int level, int top,
                 bool isVerbose = false) {
  // MO mode -------------------------------------------------------------
  const int numOfAOs = C.getNumOfRows();
  const int numOfMOs = C.getNumOfCols();

  if (level < numOfMOs) {
    std::cerr << TlUtils::format(
                     "input level is out of range. input=%d, max=%d", level,
                     C.getNumOfCols())
              << std::endl;
    return EXIT_FAILURE;
  }

  if (isVerbose == true) {
    std::cerr << "MO level: " << level << std::endl;
  }

  // MO
  const TlVector MO = C.getColVector(level);
  std::vector<LCAO> AOs(numOfAOs);
  double sum = 0.0;
  for (int i = 0; i < numOfAOs; ++i) {
    AOs[i].index = i;
    const double w = std::fabs(MO[i] * SC.get(i, level));
    sum += w;
    AOs[i].value = w;
  }

  // sort
  std::sort(AOs.begin(), AOs.end(), LCAO_cmp());

  // display
  const char shellTypes[] = {'s', 'p', 'd'};
  const char* angularMomentumTypes[] = {"s",   "px",  "py",     "pz", "dxy",
                                        "dyz", "dzx", "dxx-yy", "dzz"};

  std::cout << TlUtils::format("MO level: %d", level) << std::endl;

  if (top == -1) {
    top = numOfAOs;
  }

  double accumulate = 0.0;
  for (int rank = 0; rank < top; ++rank) {
    const int AO = AOs[rank].index;
    const double w = AOs[rank].value;
    const TlPosition pos = readOrbInfo.getPosition(AO);

    const double rate = w / sum * 100.0;
    accumulate += rate;
    std::cout << TlUtils::format(
                     "No.%2d %5d %2s (% 8.3f, % 8.3f, % 8.3f) %c(%6s) coef=% "
                     "8.3f w=% 8.3f(%6.2f%%; %6.2f%%)",
                     rank + 1, readOrbInfo.getAtomIndex(AO),
                     readOrbInfo.getAtomName(AO).c_str(), pos.x() / AU_BOHR,
                     pos.y() / AU_BOHR, pos.z() / AU_BOHR,
                     shellTypes[readOrbInfo.getShellType(AO)],
                     angularMomentumTypes[readOrbInfo.getBasisType(AO)],
                     MO.get(AO), w, rate, accumulate)
              << std::endl;
  }

  return EXIT_SUCCESS;
}

int exec_atom_mode(const TlMatrix& SC, const TlMatrix& C,
                   const TlOrbitalInfo& readOrbInfo, int inputAtomIndex,
                   const double threshold, bool isVerbose = false) {
  const int numOfAOs = C.getNumOfRows();
  const int numOfMOs = C.getNumOfCols();

  if (isVerbose) {
    std::cerr << TlUtils::format("threshold: %8.2e", threshold) << std::endl;
  }

  std::vector<int> pickupMO;
  for (int MO = 0; MO < numOfMOs; ++MO) {
    const TlVector MO_vec = C.getColVector(MO);

    if (isVerbose) {
      std::cerr << TlUtils::format("check %5dth MO ...", MO) << std::endl;
    }
    double sum = 0.0;
    double w_atom = 0.0;
    for (int i = 0; i < numOfAOs; ++i) {
      const double w = std::fabs(MO_vec[i] * SC.get(i, MO));
      sum += w;

      const int atomIndex = readOrbInfo.getAtomIndex(i);
      if (atomIndex == inputAtomIndex) {
        w_atom += w;
      }
    }

    const double ratio = w_atom / sum;
    if (ratio > threshold) {
      std::cout << TlUtils::format("MO=%5d w_atom=% 8.3f", MO, w_atom)
                << std::endl;
    }
  }

  return EXIT_SUCCESS;
}

int exec_save_mode(const TlMatrix& SC, const TlMatrix& C,
                   const std::string& savePath, bool isVerbose = false) {
  if (isVerbose) {
    std::cerr << "calc matrix..." << std::endl;
  }

  TlMatrix CSC = SC;
  CSC.dot(C);

  if (isVerbose) {
    std::cerr << TlUtils::format("save matrix: %s", savePath.c_str())
              << std::endl;
  }
  CSC.save(savePath);

  return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "A:M:S:d:c:s:hvan:t:");

  // parameters - common
  if (opt["h"] == "defined") {
    showHelp(opt[0]);
    return EXIT_FAILURE;
  }

  const bool isVerbose = (opt["v"] == "defined");

  std::string readDir = "./";
  if (opt["d"].empty() == false) {
    readDir = opt["d"];
  }

  std::string LCAO_path = "";
  if (opt["c"].empty() == false) {
    LCAO_path = opt["c"];
  }

  std::string Spq_path = "";
  if (opt["s"].empty() == false) {
    Spq_path = opt["s"];
  }

  // parameters - exec mode
  int MO_level = 0;
  int inputAtomIndex = 0;
  std::string savePath = "";
  EXEC_MODE execMode = UNDEFINED;
  if (opt["M"].empty() == false) {
    execMode = MO_MODE;
    MO_level = std::atoi(opt["M"].c_str());
  } else if (opt["A"].empty() == false) {
    execMode = ATOM_MODE;
    inputAtomIndex = std::atoi(opt["A"].c_str());
  } else if (opt["S"].empty() == false) {
    execMode = SAVE_MODE;
    savePath = opt["S"];
  }

  // parameters - MO mode
  int top = 10;
  if (opt["n"].empty() == false) {
    top = std::atoi(opt["n"].c_str());
  }
  if (opt["a"] == "defined") {
    top = -1;
  }

  // parameters - atom mode
  double threshold = 0.5;
  if (opt["t"].empty() == false) {
    threshold = std::atof(opt["t"].c_str());
  }

  // for reading object
  const std::string readParamPath =
      TlUtils::format("%s/pdfparam.mpac", readDir.c_str());
  TlMsgPack readMsgPack;
  readMsgPack.load(readParamPath);
  TlSerializeData pdfparam = readMsgPack.getSerializeData();
  const TlOrbitalInfo readOrbInfo(pdfparam["coordinates"],
                                  pdfparam["basis_set"]);
  const int lastIteration = pdfparam["num_of_iterations"].getInt();
  DfObject dfObj(&pdfparam);

  // load S
  if (Spq_path.empty()) {
    Spq_path = dfObj.getSpqMatrixPath();
  }
  TlSymmetricMatrix S;
  if (isVerbose == true) {
    std::cout << TlUtils::format("read Spq: %s", Spq_path.c_str()) << std::endl;
  }
  S.load(Spq_path);

  // load C
  if (LCAO_path.empty()) {
    LCAO_path = dfObj.getCMatrixPath(DfObject::RUN_RKS, lastIteration);
  }
  TlMatrix C;
  if (isVerbose == true) {
    std::cout << TlUtils::format("read LCAO: %s", LCAO_path.c_str())
              << std::endl;
  }
  C.load(LCAO_path);
  if (isVerbose == true) {
    std::cout << TlUtils::format("row = %d, col = %d", C.getNumOfRows(),
                                 C.getNumOfCols())
              << std::endl;
  }

  const TlMatrix SC = S * C;

  int answer = 0;
  switch (execMode) {
    case MO_MODE:
      answer = exec_MO_mode(SC, C, readOrbInfo, MO_level, top, isVerbose);
      break;
    case ATOM_MODE:
      answer = exec_atom_mode(SC, C, readOrbInfo, inputAtomIndex, threshold,
                              isVerbose);
      break;
    case SAVE_MODE:
      answer = exec_save_mode(SC, C, savePath, isVerbose);
      break;
    default:
      std::cerr << "please specify mode: " << std::endl;
      showHelp(opt[0]);
      break;
  }

  return answer;
}
