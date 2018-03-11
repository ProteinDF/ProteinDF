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
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

void showHelp(const std::string& name) {
  std::cout << TlUtils::format("%s [options] input_path output_path",
                               name.c_str())
            << std::endl;
  std::cout << " OPTIONS:" << std::endl;
  std::cout << "  -t row: top row" << std::endl;
  std::cout << "  -b row: bottom row" << std::endl;
  std::cout << "  -l col: left col" << std::endl;
  std::cout << "  -r col: right col" << std::endl;
  std::cout << "  -h:      show help" << std::endl;
  std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "t:b:l:r:hv");

  if (opt["h"] == "defined") {
    showHelp(opt[0]);
    return EXIT_SUCCESS;
  }

  const bool isVerbose = (opt["v"] == "defined");
  if (opt.getCount() <= 2) {
    showHelp(opt[0]);
    return EXIT_FAILURE;
  }
  std::string inputMatrixPath = opt[1];
  std::string outputMatrixPath = opt[2];

  if (isVerbose == true) {
    std::cerr << "load matrix: " << inputMatrixPath << std::endl;
  }

  TlMatrix A;
  if (TlMatrix::isLoadable(inputMatrixPath)) {
    A.load(inputMatrixPath);
  } else if (TlSymmetricMatrix::isLoadable(inputMatrixPath)) {
    TlSymmetricMatrix tmp;
    tmp.load(inputMatrixPath);
    A = tmp;
  } else {
    std::cerr << TlUtils::format("cannot load: %s", inputMatrixPath.c_str())
              << std::endl;
    return EXIT_FAILURE;
  }

  TlMatrix::index_type top = 0;
  if (!opt["t"].empty()) {
    top = std::atoi(opt["t"].c_str());
    if ((top < 0) || (A.getNumOfRows() <= top)) {
      std::cerr << TlUtils::format("top: out of range(0 - %d): %d",
                                   A.getNumOfRows() - 1, top)
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  TlMatrix::index_type bottom = A.getNumOfRows();
  if (!opt["t"].empty()) {
    bottom = std::atoi(opt["b"].c_str());
    if ((bottom < 0) || (A.getNumOfRows() <= bottom)) {
      std::cerr << TlUtils::format("bottom: out of range(0 - %d): %d",
                                   A.getNumOfRows() - 1, bottom)
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (top >= bottom) {
    std::cerr << TlUtils::format("no left row lines: top: %d bottom: %d", top,
                                 bottom)
              << std::endl;
    return EXIT_FAILURE;
  }

  TlMatrix::index_type left = 0;
  if (!opt["l"].empty()) {
    left = std::atoi(opt["l"].c_str());
    if ((left < 0) || (A.getNumOfCols() <= left)) {
      std::cerr << TlUtils::format("left: out of range(0 - %d): %d",
                                   A.getNumOfCols() - 1, left)
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  TlMatrix::index_type right = A.getNumOfCols();
  if (!opt["t"].empty()) {
    right = std::atoi(opt["r"].c_str());
    if ((right < 0) || (A.getNumOfCols() <= right)) {
      std::cerr << TlUtils::format("right: out of range(0 - %d): %d",
                                   A.getNumOfCols() - 1, right)
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (left >= right) {
    std::cerr << TlUtils::format("no left col lines: left: %d right: %d", left,
                                 right)
              << std::endl;
    return EXIT_FAILURE;
  }

  TlMatrix::index_type numOfRows = bottom - top + 1;
  TlMatrix::index_type numOfCols = right - left + 1;
  TlMatrix answer(numOfRows, numOfCols);
  for (TlMatrix::index_type r = 0; r < numOfRows; ++r) {
    for (TlMatrix::index_type c = 0; c < numOfCols; ++c) {
      double v = A.get(top + r, left + c);
      answer.set(r, c, v);
    }
  }

  if (isVerbose == true) {
    std::cerr << "save matrix: " << outputMatrixPath << std::endl;
  }
  answer.save(outputMatrixPath);

  return EXIT_SUCCESS;
}
