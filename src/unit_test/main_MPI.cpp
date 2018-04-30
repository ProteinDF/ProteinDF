#include "TlCommunicate.h"
#include "gtest/gtest.h"

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);

  // MPI_Init(&argc, &argv);
  TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

  int result = RUN_ALL_TESTS();

  // MPI_Finalize();
  rComm.finalize();

  return result;
}
