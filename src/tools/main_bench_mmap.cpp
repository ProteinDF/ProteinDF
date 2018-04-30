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

#include <stdlib.h>
#include <time.h>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

#include "TlGetopt.h"
#include "TlMemManager.h"

double bench_rand(std::size_t count, TlMatrix* pMtrx) {
  const int maxRow = pMtrx->getNumOfRows();
  const int maxCol = pMtrx->getNumOfCols();

  srand((unsigned int)time(NULL));

  clock_t startTime = clock();
  for (std::size_t i = 0; i < count; ++i) {
    const int x1 = rand() % maxRow;
    const int y1 = rand() % maxCol;
    double v = pMtrx->get(x1, y1) + 1.0;

    const int x2 = rand() % maxRow;
    const int y2 = rand() % maxCol;
    pMtrx->set(x2, y2, v);
  }
  clock_t endTime = clock();

  return double(endTime - startTime) / double(CLOCKS_PER_SEC);
}

double bench_seq(std::size_t count, TlMatrix* pMtrx) {
  const int maxRow = pMtrx->getNumOfRows();
  const int maxCol = pMtrx->getNumOfCols();

  srand((unsigned int)time(NULL));

  std::size_t c = 0;
  clock_t startTime = clock();
  for (int x1 = 0; x1 < maxRow; ++x1) {
    for (int y1 = 0; y1 < maxCol; ++y1) {
      double v = pMtrx->get(x1, y1) + 1.0;
      pMtrx->set(x1, y1, v);
      ++c;
      if (c >= count) {
        break;
      }
    }
  }
  clock_t endTime = clock();

  return double(endTime - startTime) / double(CLOCKS_PER_SEC);
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "d:f:mrs");

  std::size_t count = 10 * 1024 * 1024;

  int dim = 1000;
  if (!opt["d"].empty()) {
    dim = std::atoi(opt["d"].c_str());
  }

  std::string mmap_file = "/tmp/mmap_bench.map";
  if (!opt["f"].empty()) {
    mmap_file = opt["f"];
  }

  bool isSymmetric = false;
  if (!opt["s"].empty()) {
    isSymmetric = true;
  }

  bool isMMap = false;
  if (!opt["m"].empty()) {
    isMMap = true;
  }

  bool isRandom = false;
  if (!opt["r"].empty()) {
    isRandom = true;
  }

  // show information
  std::cout << TlUtils::format("(%6d, %6d) ", dim, dim);
  if (isSymmetric == true) {
    std::cout << TlUtils::format(
        "[sym] mem_size=%8.2fMB ",
        double(dim * (dim + 1) * sizeof(double)) / (2.0 * 1024.0 * 1024.0));
  } else {
    std::cout << TlUtils::format(
        "[---] mem_size=%8.2fMB ",
        double(dim * dim * sizeof(double)) / (1024.0 * 1024.0));
  }
  if (isMMap == true) {
    std::cout << "[mmap  ] ";
  } else {
    std::cout << "[memory] ";
  }
  if (isRandom == true) {
    std::cout << "[random R/W    ] " << std::endl;
  } else {
    std::cout << "[sequential R/W] " << std::endl;
  }

  // bench
  double elapse = 0.0;
  if (isMMap == true) {
    std::size_t needMemSize =
        dim * dim * sizeof(double) + 100;  // 100は念のため
    TlMemManager::setParam(needMemSize, mmap_file);
    TlMatrix::useMemManager(true);
  }

  if (isSymmetric == true) {
    TlSymmetricMatrix m(dim);
    if (isRandom == true) {
      elapse = bench_rand(count, &m);
    } else {
      elapse = bench_seq(count, &m);
    }
  } else {
    TlMatrix m(dim, dim);
    if (isRandom == true) {
      elapse = bench_rand(count, &m);
    } else {
      elapse = bench_seq(count, &m);
    }
  }

  // results
  std::cout << "access time = " << elapse << " [s]" << std::endl;
  std::size_t dataSize = count * sizeof(double) * 2;  // *2 = get, set
  double ratio = dataSize / elapse;
  std::cout << "ratio = " << (ratio / (1024.0 * 1024.0)) << " MB/sec"
            << std::endl;

  return EXIT_SUCCESS;
}
