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

#include <iostream>

#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_vector_lapack.h"

void help(const std::string& prgname) {
  std::cout << TlUtils::format("Usage: %s [options]... FILE1 FILE2",
                               prgname.c_str())
            << std::endl;
  std::cout << "compare ProteinDF vector files" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << " -h           show help message (this)." << std::endl;
  std::cout << " -v           show message verbosely." << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "hv");

  bool bVerbose = (opt["v"] == "defined");
  if (opt.getCount() < 2) {
    help(opt[0]);
    std::exit(1);
  }

  const std::string sPath1 = opt[1];
  const std::string sPath2 = opt[2];
  if (bVerbose) {
    std::cerr << "loading... " << sPath1 << std::endl;
    std::cerr << "loading... " << sPath2 << std::endl;
  }

  int nErrorCode = 0;

  TlDenseVector_Lapack v1, v2;
  if (v1.load(sPath1) == false) {
    std::cerr << "could not open: " << sPath1 << std::endl;
    return 1;
  }
  if (v2.load(sPath2) == false) {
    std::cerr << "could not open: " << sPath2 << std::endl;
    return 1;
  }

  v1 -= v2;
  std::cout << v1 << std::endl;

  return nErrorCode;
}
