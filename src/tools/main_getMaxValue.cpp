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
#include "TlMatrix.h"

enum FileType { FT_NONE, FT_SYMMAT, FT_MAT };

struct ElementInfo {
  ElementInfo(TlMatrixObject::index_type r, TlMatrixObject::index_type c,
              double v)
      : row(r), col(c), value(v) {}

  TlMatrixObject::index_type row;
  TlMatrixObject::index_type col;
  double value;
};

struct CompElementInfo {
  bool operator()(const ElementInfo& rhs1, const ElementInfo& rhs2) const {
    return (rhs1.value > rhs2.value);  // 大きい順
  }
};

struct CompElementInfo_ABS {
  bool operator()(const ElementInfo& rhs1, const ElementInfo& rhs2) const {
    return (rhs1.value > rhs2.value);
  }
};

void showHelp();
FileType getFileType(const std::string& filePath);

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "ac:hs:v");

  int errorCode = 0;
  bool verbose = (opt["v"] == "defined");
  int count = 10;
  if (opt["c"].size() != 0) {
    count = std::atoi(opt["c"].c_str());
  }
  bool abs_mode = false;
  if (opt["a"] == "defined") {
    abs_mode = true;
  }

  const std::string filePath = opt[1];
  if (verbose) {
    std::cerr << "loading... " << filePath << std::endl;
  }

  std::ifstream ifs;
  ifs.open(filePath.c_str());
  if (ifs.fail()) {
    std::cerr << "could not open file. " << filePath << std::endl;
    errorCode = 1;
    return errorCode;
  }

  const FileType fileType = getFileType(filePath);
  TlMatrix* pMat = NULL;
  switch (fileType) {
    case FT_SYMMAT:
      pMat = new TlDenseSymmetricMatrix_BLAS_Old();
      break;

    case FT_MAT:
      pMat = new TlMatrix();
      break;

    default:
      errorCode = 2;
      break;
  }

  if (pMat != NULL) {
    if (pMat->load(filePath) == true) {
      const TlMatrixObject::index_type numOfRows = pMat->getNumOfRows();
      const TlMatrixObject::index_type numOfCols = pMat->getNumOfCols();

      std::vector<ElementInfo> answer;
      answer.reserve(count + 1);

      if (fileType == FT_SYMMAT) {
        for (TlMatrixObject::index_type row = 0; row < numOfRows; ++row) {
          for (TlMatrixObject::index_type col = 0; col <= row; ++col) {
            const ElementInfo ei(row, col, pMat->get(row, col));
            std::vector<ElementInfo>::iterator it;
            if (abs_mode == true) {
              it = std::lower_bound(answer.begin(), answer.end(), ei,
                                    CompElementInfo());
            } else {
              it = std::lower_bound(answer.begin(), answer.end(), ei,
                                    CompElementInfo_ABS());
            }
            answer.insert(it, ei);
            if (answer.size() > (std::size_t)count) {
              answer.pop_back();
            }
          }
        }
      } else {
        for (TlMatrixObject::index_type row = 0; row < numOfRows; ++row) {
          for (TlMatrixObject::index_type col = 0; col < numOfCols; ++col) {
            const ElementInfo ei(row, col, pMat->get(row, col));
            std::vector<ElementInfo>::iterator it;
            if (abs_mode == true) {
              it = std::lower_bound(answer.begin(), answer.end(), ei,
                                    CompElementInfo());
            } else {
              it = std::lower_bound(answer.begin(), answer.end(), ei,
                                    CompElementInfo_ABS());
            }
            answer.insert(it, ei);
            if (answer.size() > (std::size_t)count) {
              answer.pop_back();
            }
          }
        }
      }

      // output
      std::cout << "row, col, value" << std::endl;
      for (std::vector<ElementInfo>::const_iterator it = answer.begin();
           it != answer.end(); ++it) {
        std::cout << it->row << ", " << it->col << ", " << it->value
                  << std::endl;
      }
    }
  }

  return errorCode;
}

void showHelp() {
  std::cout << "getMaxValue [options] <file1> <file2>" << std::endl;
}

FileType getFileType(const std::string& filePath) {
  FileType fileType = FT_NONE;
  if (TlDenseSymmetricMatrix_BLAS_Old::isLoadable(filePath) == true) {
    fileType = FT_SYMMAT;
  } else if (TlMatrix::isLoadable(filePath) == true) {
    fileType = FT_MAT;
  }

  return fileType;
}
