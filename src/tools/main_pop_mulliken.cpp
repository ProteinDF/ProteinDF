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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "DfPopulation.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"

void usage();

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "hi:vs:");

  std::string paramPath = "pdfparam.mpac";
  if (opt["p"].empty() != true) {
    paramPath = opt["p"];
  }
  const std::string savePath = opt["s"];

  const bool isVerbose = (opt["v"] == "defined");
  const bool showHelp = (opt["h"] == "defined");

  if (showHelp == true) {
    usage();
    return EXIT_SUCCESS;
  }

  TlMsgPack mpac;
  mpac.load(paramPath);
  TlSerializeData param = mpac.getSerializeData();

  int iteration = 0;
  if (opt["i"].empty() != true) {
    iteration = std::atoi(opt["i"].c_str());
  } else {
    iteration = param["num_of_iterations"].getInt();
  }

  if (isVerbose == true) {
    std::cerr << "iteration = " << iteration << std::endl;
  }

  DfPopulation dfPop(&param);

  if (savePath.empty() != true) {
    const TlMatrix mtx = dfPop.getAtomPopData(iteration);
    if (isVerbose == true) {
      std::cerr << "save Mulliken Data as " << savePath << std::endl;
    }
    mtx.save(savePath);
  } else {
    dfPop.getReport(iteration, std::cout);
  }

  return EXIT_SUCCESS;
}

void usage() {
  std::cerr << "usage: pdf-pop-mulliken [options] " << std::endl;
  std::cerr << std::endl;
  std::cerr << "calculation Mulliken Population." << std::endl;
  std::cerr << "  -p param:\n";
  std::cerr << "  -i num:    set SCF iteration to get Mulliken Population\n";
  std::cerr << "  -s path:   save atom population data as pdf matrix file\n";
  std::cerr << "  -h:        show help(this)\n";
  std::cerr << "  -v:        verbose\n";
}
