#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <string>

#include "TlCommunicate.h"
#include "TlTime.h"
#include "TlUtils.h"

int main(int argc, char* argv[]) {
  TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);
  const int numOfProcs = rComm.getNumOfProc();
  const int tag = 999;

  if (rComm.isMaster() == true) {
    std::vector<int> msg(numOfProcs, 0);
    for (int i = 1; i < numOfProcs; ++i) {
      rComm.iReceiveData(msg[i], i, tag);
    }

    bool breakLoop = false;
    std::vector<bool> finished(numOfProcs, false);
    while (breakLoop != true) {
      for (int i = 1; i < numOfProcs; ++i) {
        if (finished[i] != true) {
          if (rComm.test(msg[i]) == true) {
            rComm.wait(msg[i]);
            std::cout << TlUtils::format(
                             "[OK] rank0: receive message from rank%d", i)
                      << std::endl;
            finished[i] = true;
          }
        } else {
          std::cout << TlUtils::format(
                           "[--] rank0: waiting message from rank%d", i)
                    << std::endl;
        }
      }

      breakLoop = true;
      for (int i = 1; i < numOfProcs; ++i) {
        if (finished[i] != true) {
          breakLoop = false;
          break;
        }
      }
    }
  } else {
    TlTime::sleep(4000);  // wait 4 sec.

    int msg = rComm.getRank();
    rComm.sendData(msg, 0, tag);
  }

  rComm.finalize();
  return EXIT_SUCCESS;
}
