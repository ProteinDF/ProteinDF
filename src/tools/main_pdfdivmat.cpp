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

#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlDistributeVector.h"
#include "TlGetopt.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

void showHelp() {
  std::cout << "pdfdivmat [options] row col path" << std::endl;
  std::cout << "divide matrix" << std::endl;
  std::cout << " OPTIONS:" << std::endl;
  std::cout << "  -h:      show help" << std::endl;
  std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[]) {
  TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

  TlGetopt opt(argc, argv, "b:hv");
  int scalapackBlockSize = 64;
  if (!opt["b"].empty()) {
    scalapackBlockSize = std::atoi(opt["b"].c_str());
  }
  if (opt["h"] == "defined") {
    if (rComm.isMaster()) {
      showHelp();
    }
    rComm.finalize();
    return EXIT_SUCCESS;
  }
  const bool isVerbose = (opt["v"] == "defined");

  if (opt.getCount() < 2) {
    if (rComm.isMaster()) {
      showHelp();
    }
    rComm.finalize();
    return EXIT_FAILURE;
  }

  const std::string loadFilePath = opt[1];
  const std::string saveFilePath =
      TlUtils::format("%s.%d", loadFilePath.c_str(), rComm.getRank());

  TlDistributeMatrix::setSystemBlockSize(scalapackBlockSize);
  TlDistributeVector::setSystemBlockSize(scalapackBlockSize);

  if (TlSymmetricMatrix::isLoadable(loadFilePath) == true) {
    TlDistributeSymmetricMatrix M;
    TlDistributeMatrix::setUsingPartialIO(false);
    if (isVerbose) {
      if (rComm.isMaster()) {
        std::cerr << TlUtils::format("load symmetric matrix: %s.",
                                     loadFilePath.c_str())
                  << std::endl;
      }
    }
    M.load(loadFilePath);

    TlDistributeMatrix::setUsingPartialIO(true);
    if (isVerbose) {
      if (rComm.isMaster()) {
        std::cerr << TlUtils::format("save distribute matrix on 0: %s.",
                                     saveFilePath.c_str())
                  << std::endl;
      }
    }
    M.save(saveFilePath);
  } else if (TlMatrix::isLoadable(loadFilePath) == true) {
    TlDistributeMatrix M;
    TlDistributeMatrix::setUsingPartialIO(false);
    if (isVerbose) {
      if (rComm.isMaster()) {
        std::cerr << TlUtils::format("load general matrix: %s.",
                                     loadFilePath.c_str())
                  << std::endl;
      }
    }
    M.load(loadFilePath);

    TlDistributeMatrix::setUsingPartialIO(true);
    if (isVerbose) {
      if (rComm.isMaster()) {
        std::cerr << TlUtils::format("save distribute matrix on 0: %s.",
                                     saveFilePath.c_str())
                  << std::endl;
      }
    }
    M.save(saveFilePath);
  } else {
    std::cerr << "unknown file type: " << loadFilePath << std::endl;
    rComm.finalize();
    return EXIT_FAILURE;
  }

  rComm.finalize();
  return EXIT_SUCCESS;
}
