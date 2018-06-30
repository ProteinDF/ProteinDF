#include "TlCommunicate.h"
#include "gtest/gtest.h"

TEST(TlCommunicate, sendrecv_int) {
  TlCommunicate& rComm = TlCommunicate::getInstance();
  const int proc = rComm.getNumOfProc();
  const int rank = rComm.getRank();

  for (int i = 1; i < proc; i++) {
    if (rank == i) {
      rComm.sendData(100 + i);
    } else if (rComm.isMaster()) {
      int rData = 0;
      rComm.receiveData(rData, i);

      EXPECT_EQ(100 + i, rData);
    }
  }

  rComm.barrier();
}
