#include <iostream>

#include "TlCommunicate.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlUtils.h"

void showResultMessage(const std::string& sFunction, bool bIsPassed);
void showResultMessageAll(const std::string& sFunction, bool bIsPassed);
void testCommunicate_msgpack();

int main(int argc, char* argv[]) {
    // initialize
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

    // ===================================================================
    testCommunicate_msgpack();
    // ===================================================================

    // finalize
    rComm.finalize();
    return EXIT_SUCCESS;
}

void showResultMessage(const std::string& sFunction, bool bIsPassed) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        std::cout << "TEST: " << sFunction << "() ";
        std::cout << ((bIsPassed == true) ? "." : "F") << std::endl;
    }
    rComm.barrier();
}

void showResultMessageAll(const std::string& sFunction, bool bIsPassed) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    int nProc = rComm.getNumOfProc();
    int nRank = rComm.getRank();

    rComm.barrier();
    if (rComm.isMaster() == true) {
        std::cout << "TEST: " << sFunction << "() ";
    }
    for (int i = 0; i < nProc; ++i) {
        if (i == nRank) {
            std::cout << ((bIsPassed == true) ? "." : "F");
            std::cout.flush();
        }
        rComm.barrier();
    }
    if (rComm.isMaster() == true) {
        std::cout << std::endl;
    }

    rComm.barrier();
}

void testCommunicate_msgpack() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    TlMsgPack msgpack;
    std::string buf;
    if (rComm.isMaster() == true) {
        msgpack.load("sample.mpac");
        TlSerializeData so = msgpack.getSerializeData();
        // std::cerr << so.str() << std::endl;
        buf = msgpack.dump();
        // std::cerr << "buf size=" << buf.size() << std::endl;
    }

    rComm.broadcast(buf);

    // std::cerr << TlUtils::format("[%d] buf size=%d", rank, buf.size()) <<
    // std::endl;

    for (int i = 0; i < proc; ++i) {
        if (rank == i) {
            msgpack.pack(buf);
            const TlSerializeData so = msgpack.getSerializeData();
            std::cerr << TlUtils::format("[%d] %s", rank, so.str().c_str())
                      << std::endl;
        }
        rComm.barrier();
    }

    showResultMessage("testCommunicate_msgpack", bIsPassed);
}
