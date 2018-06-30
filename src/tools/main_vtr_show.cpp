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

#include "TlGetopt.h"
#include "tl_dense_vector_blas.h"
#include "tl_vector_utils.h"

void showHelp(const std::string& progname) {
  std::cout << TlUtils::format("Usage: %s [options]... FILE", progname.c_str())
            << std::endl;
  std::cout << "display ProteinDF vector file" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << " -g           show guess mode" << std::endl;
  std::cout << " -l           show list mode" << std::endl;
  std::cout << " -h           show help message (this)." << std::endl;
  std::cout << " -v           show message verbosely." << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "ghlv");

  const bool bGuessMode = (opt["g"] == "defined");
  const bool bListMode = (opt["l"] == "defined");
  const bool bVerbose = (opt["v"] == "defined");

  if (opt["h"] == "defined") {
    showHelp(opt[0]);
  }

  std::string sPath = opt[1];
  if (bVerbose) {
    std::cerr << "loading... " << sPath << std::endl;
  }

  TlVector_BLAS v;
  bool bIsLoadable = TlVectorUtils::isLoadable(sPath);
  if (bIsLoadable == false) {
    std::cerr << "could not open: " << sPath << std::endl;
    return EXIT_FAILURE;
  }
  v.load(sPath);

  if (bGuessMode == true) {
    v.outputText(std::cout);
  } else if (bListMode == true) {
    const int nSize = v.getSize();
    for (int i = 0; i < nSize; ++i) {
      std::cout << v[i] << std::endl;
    }
  } else {
    v.print(std::cout);
  }

  return EXIT_SUCCESS;
}
