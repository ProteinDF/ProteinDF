#include <cstdlib>
#include <iostream>
#include <string>

#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlTime.h"
#include "TlUtils.h"

#define MATRIX_FILENAME_FOR_INVERSE "benchForInverse.matrix"
#define MATRIX_FILENAME_FOR_DIAGONAL "benchForDiagonal.matrix"

void msg(const std::string& str);
void inverse();
void diagonal();

int main(int argc, char* argv[]) {
  TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

  inverse();

  diagonal();

  rComm.finalize();
  return EXIT_SUCCESS;
}

void msg(const std::string& msg) {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  if (rComm.isMaster() == true) {
    std::cout << msg << std::endl;
  }
}

void inverse() {
  msg("start benckmark 'inverse'");
  TlTime tlTime;
  TlDistributeSymmetricMatrix matrix;
  matrix.load(MATRIX_FILENAME_FOR_INVERSE);
  msg(TlUtils::format("elapse time for loading: %f [sec]",
                      tlTime.getElapseTime()));

  tlTime.start();
  matrix.inverse();

  msg(TlUtils::format("elapse time for inverse: %f [sec]",
                      tlTime.getElapseTime()));
  msg("end of benchmark 'inverse'\n");
}

void diagonal() {
  msg("start benckmark 'diagonal'");
  TlTime tlTime;
  TlDistributeSymmetricMatrix matrix;
  matrix.load(MATRIX_FILENAME_FOR_DIAGONAL);
  msg(TlUtils::format("elapse time for loading: %f [sec]",
                      tlTime.getElapseTime()));

  tlTime.start();
  const int size = matrix.getNumOfRows();
  TlVector eigVal(size);
  TlDistributeMatrix eigVec(size, size);
  matrix.diagonal(&eigVal, &eigVec);

  msg(TlUtils::format("elapse time for diagonal: %f [sec]",
                      tlTime.getElapseTime()));
  msg("end of benchmark 'diagonal'\n");
}
