#include <cstdlib>
#include <iostream>
#include <string>

#include <iostream>
#include "ProteinDF.h"

#include "TlDenseSymmetricMatrix_BLAS_Old.h"
#include "tl_sparse_matrix.h"
#include "tl_sparse_symmetric_matrix.h"

#include "TlCommunicate.h"
#include "TlFile.h"
#include "TlMemManager.h"

#define ASSERT_EQ(value, expect)                                          \
    do {                                                                  \
        if (std::fabs(value - expect) > 1.0E-5) {                         \
            std::cerr                                                     \
                << TlUtils::format(                                       \
                       "ASSERTION ERROR at %s:%d: value=% f, expect=% f", \
                       __FILE__, __LINE__, value, expect)                 \
                << std::endl;                                             \
        }                                                                 \
    } while (0);

#define HOSTNAME_LEN 256

void showResultMessage(const std::string& sFunction, bool bIsPassed);
void showResultMessageAll(const std::string& sFunction, bool bIsPassed);

void checkAllProcess();
void testTlLog_Parallel();
void testCommunicate_sendInt();
void testCommunicate_isend_int();
void testCommunicate_sendLong();
void testCommunicate_sendDouble();
void testCommunicate_sendVectorInt();
void testCommunicate_sendVectorLong();
void testCommunicate_sendVectorDouble();
void testCommunicate_sendString();
void testCommunicate_send_TlSparseMatrix();
void testCommunicate_isend_TlSparseMatrix();
void testCommunicate_recvFronAnySource_TlSparseMatrix();

void testGatherToMaster_TlSparseMatrix();

void testBroadcast_bool();
void testBroadcast_int();
void testBroadcast_string();
void testBroadcast_vectorString();
void testBroadcast_valarray();
void testBroadcast_TlVector();
void testBroadcast_TlSparseMatrix();
void testBroadcast_TlMmapSymmetricMatrix();

void testAllreduceSum_double();
void testAllreduceSum_vectorDouble();
void testAllreduceSum_TlVector();
void testAllreduceSum_matrix_symmetric();
void testAllreduceSum_TlMmapSymmetricMatrix();

void test_nonblocking();

int main(int argc, char* argv[]) {
    // initialize
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

    TlMemManager::setParam(1024 * 1024 * 1024,  // 1 GB
                           "/tmp/parallelTest");

    //
    {
        char* pHostName = new char[HOSTNAME_LEN];
        (void)gethostname(pHostName, HOSTNAME_LEN);
        const std::string hostname(pHostName);
        delete[] pHostName;
        pHostName = NULL;

        const int numOfProcs = rComm.getNumOfProcs();
        for (int i = 0; i < numOfProcs; ++i) {
            if (i == rComm.getRank()) {
                std::cerr << TlUtils::format("PID %d on %s as rank %d",
                                             getpid(), hostname.c_str(),
                                             rComm.getRank())
                          << std::endl;
            }
            rComm.barrier();
        }
        rComm.barrier();

        const int waitingTime = 4000;
        rComm.barrier();
        if (rComm.isMaster() == true) {
            std::cerr << TlUtils::format("waiting %d msec.", waitingTime)
                      << std::endl;
        }
        TlTime::sleep(waitingTime);
        rComm.barrier();
    }

    // ===================================================================
    checkAllProcess();

    testCommunicate_sendInt();
    testCommunicate_isend_int();
    testCommunicate_sendLong();
    testCommunicate_sendDouble();
    testCommunicate_sendVectorInt();
    testCommunicate_sendVectorLong();
    testCommunicate_sendVectorDouble();
    testCommunicate_sendString();

    testCommunicate_send_TlSparseMatrix();
    testCommunicate_isend_TlSparseMatrix();
    testCommunicate_recvFronAnySource_TlSparseMatrix();
    testGatherToMaster_TlSparseMatrix();

    testBroadcast_bool();
    testBroadcast_int();
    testBroadcast_string();
    testBroadcast_vectorString();
    testBroadcast_valarray();
    testBroadcast_TlVector();
    testBroadcast_TlSparseMatrix();
    // testBroadcast_TlMmapSymmetricMatrix();

    testAllreduceSum_double();
    testAllreduceSum_vectorDouble();
    testAllreduceSum_TlVector();
    testAllreduceSum_matrix_symmetric();
    // testAllreduceSum_TlMmapSymmetricMatrix();
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

void checkAllProcess() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    // const int rank = rComm.getRank();

    for (int i = 0; i < proc; ++i) {
        rComm.barrier();
    }

    showResultMessage("checkAllProcess", bIsPassed);
}

// void testTlLog_Parallel(){
//   TlLog_Parallel info(TlLog::INFO);

//   info << "this message is written by only master process." << std::endl;
// }

void testCommunicate_sendInt() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            rComm.sendData(100 + i);
        } else if (rComm.isMaster()) {
            int rData = 0;
            rComm.receiveData(rData, i);

            if (rData != (100 + i)) {
                std::cout << "recv(int) = " << rData
                          << ", expect = " << (100 + i) << std::endl;
                bIsPassed = false;
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendInt", bIsPassed);
}

void testCommunicate_isend_int() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    // const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();
    const int tag = 999;

    if (rank == 1) {
        int msg = 123;
        rComm.iSendData(msg, 0, tag);
        rComm.wait(&msg);
    } else if (rank == 0) {
        int msg = 0;
        rComm.receiveData(msg, 1, tag);

        ASSERT_EQ(msg, 123);
    }

    showResultMessage("testCommunicate_isend_int", bIsPassed);
}

void testCommunicate_sendLong() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            long sData = 1000 + i;
            rComm.sendData(sData);
        } else if (rComm.isMaster()) {
            long rData = 0;
            rComm.receiveData(rData, i);

            if (rData != (1000 + i)) {
                std::cout << "recv(long) = " << rData
                          << ", expect = " << (1000 + i) << std::endl;
                bIsPassed = false;
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendLong", bIsPassed);
}

void testCommunicate_sendString() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            std::string sData = "Hello, world!";
            for (int j = 0; j < i; j++) {
                sData += " hoge!";
            }
            rComm.sendData(sData);
        } else if (rComm.isMaster()) {
            std::string rData;
            rComm.receiveData(rData, i);

            std::string sExpectData = "Hello, world!";
            for (int j = 0; j < i; j++) {
                sExpectData += " hoge!";
            }

            if (rData != sExpectData) {
                std::cout << "recv(string) = " << rData
                          << ", expect = " << sExpectData << std::endl;
                bIsPassed = false;
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendString", bIsPassed);
}

void testCommunicate_sendDouble() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            rComm.sendData(double(100 + i));
        } else if (rComm.isMaster()) {
            double rData = 0.0;
            rComm.receiveData(rData, i);

            if (std::fabs(rData - (100 + i)) > 1.0E-16) {
                std::cout << "recv(double) = " << rData
                          << ", expect = " << (100 + i) << std::endl;
                bIsPassed = false;
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendDouble", bIsPassed);
}

void testCommunicate_sendVectorInt() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            std::vector<int> v(proc);
            for (int j = 0; j < proc; ++j) {
                v[j] = 100 + j;
            }
            rComm.sendData(v);
        } else if (rComm.isMaster()) {
            std::vector<int> rData;
            rComm.receiveData(rData, i);

            for (int j = 0; j < proc; ++j) {
                if (rData[j] != (100 + j)) {
                    std::cout << "recv(vector int) = " << rData[j]
                              << ", expect = " << (100 + j) << std::endl;
                    bIsPassed = false;
                }
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendVectorInt", bIsPassed);
}

void testCommunicate_sendVectorLong() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            std::vector<long> v(proc);
            for (int j = 0; j < proc; ++j) {
                v[j] = 100 + j;
            }
            rComm.sendData(v);
        } else if (rComm.isMaster()) {
            std::vector<long> rData;
            rComm.receiveData(rData, i);

            for (int j = 0; j < proc; ++j) {
                if (rData[j] != (100 + j)) {
                    std::cout << "recv(vector long) = " << rData[j]
                              << ", expect = " << (100 + j) << std::endl;
                    bIsPassed = false;
                }
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendVectorLong", bIsPassed);
}

void testCommunicate_sendVectorDouble() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            std::vector<double> v(proc);
            for (int j = 0; j < proc; ++j) {
                v[j] = double(100 + j);
            }
            rComm.sendData(v);
        } else if (rComm.isMaster()) {
            std::vector<double> rData;
            rComm.receiveData(rData, i);

            for (int j = 0; j < proc; ++j) {
                if (std::fabs(rData[j] - (100 + j)) > 1.0E-16) {
                    std::cout << "recv(vector double) = " << rData[j]
                              << ", expect = " << (100 + j) << std::endl;
                    bIsPassed = false;
                }
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendVectorDouble", bIsPassed);
}

void testCommunicate_send_TlSparseMatrix() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            TlSparseMatrix m(100, 200);
            m.set(rank, rank * 2, 1.23 * rank);
            rComm.sendData(m);
        } else if (rComm.isMaster()) {
            TlSparseMatrix m;
            rComm.receiveData(m, i);

            ASSERT_EQ(m.getNumOfRows(), 100);
            ASSERT_EQ(m.getNumOfCols(), 200);

            ASSERT_EQ(m.get(i, i * 2), 1.23 * i);
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_send_TlSparseMatrix", bIsPassed);
}

void testCommunicate_isend_TlSparseMatrix() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    // const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();
    const int tag = 999;

    if (rank == 1) {
        TlSparseMatrix m(10, 10);
        m.set(1, 1, 1.0);
        m.set(2, 2, 2.0);
        m.set(3, 3, 3.0);

        const std::size_t numOfSize = m.getSize();
        const std::vector<TlMatrixElement> elements = m.getMatrixElements();

        rComm.iSendData(numOfSize, 0, tag);
        rComm.iSendDataX(&(elements[0]), numOfSize, 0, tag);

        rComm.wait(&numOfSize);
        rComm.wait(&(elements[0]));
    } else if (rank == 0) {
        TlSparseMatrix m(10, 10);
        std::size_t numOfSize = 0;
        std::vector<TlMatrixElement> elements;

        rComm.receiveData(numOfSize, 1, tag);
        elements.resize(numOfSize);
        rComm.receiveDataX(&(elements[0]), numOfSize, 1, tag);

        m.add(elements);
        ASSERT_EQ(m.getNumOfRows(), 10);
        ASSERT_EQ(m.getNumOfCols(), 10);
        ASSERT_EQ(m.get(1, 1), 1.0);
        ASSERT_EQ(m.get(2, 2), 2.0);
        ASSERT_EQ(m.get(3, 3), 3.0);
    }

    showResultMessage("testCommunicate_isend_TlSparseMatrix", bIsPassed);
}

void testCommunicate_recvFronAnySource_TlSparseMatrix() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    const int tag = 12345;
    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            TlSparseMatrix m(100, 200);
            m.set(rank, rank * 2, 1.23 * rank);
            rComm.sendData(m, 0, tag);
        } else if (rComm.isMaster()) {
            TlSparseMatrix m;
            int src = 0;
            rComm.receiveDataFromAnySource(&m, &src, tag);

            ASSERT_EQ(m.getNumOfRows(), 100);
            ASSERT_EQ(m.getNumOfCols(), 200);
            ASSERT_EQ(m.get(src, src * 2), 1.23 * src);
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendTlSparseMatrix", bIsPassed);
}

void testGatherToMaster_TlSparseMatrix() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    // 準備
    TlSparseSymmetricMatrix sm(proc);
    sm(rank, rank) = rank + 1;

    // test
    rComm.allReduce_SUM(sm);

    bool bIsPassed = true;
    if (rComm.isMaster() == true) {
        for (int i = 0; i < proc; ++i) {
            // 対角項
            const double check = std::fabs(sm(i, i) - (i + 1.0));
            if (check > 1.0E-16) {
                bIsPassed = false;
            }

            // 非対角項
            for (int j = 0; j < i; ++j) {
                if (std::fabs(sm(i, j) - 0.0) > 1.0E-16) {
                    bIsPassed = false;
                }
            }
            for (int j = i + 1; j < proc; ++j) {
                if (std::fabs(sm(i, j) - 0.0) > 1.0E-16) {
                    bIsPassed = false;
                }
            }
        }
    }

    showResultMessage("testGatherToMaster_TlSparseMatrix", bIsPassed);
}

void testBroadcast_bool() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    bool b = false;
    if (rComm.isMaster() == true) {
        b = true;
    }

    // do test
    rComm.broadcast(b);

    if (b != true) {
        bIsPassed = false;
    }

    showResultMessageAll("testBroadcast_bool", bIsPassed);
}

void testBroadcast_int() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    int num = 0;
    if (rComm.isMaster() == true) {
        num = 317;
    }

    rComm.broadcast(num);

    if (num != 317) {
        bIsPassed = false;
    }

    showResultMessageAll("testBroadcast_int", bIsPassed);
}

void testBroadcast_string() {
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    std::string str = "slave";
    if (rComm.isMaster() == true) {
        str = "master";
    }

    rComm.broadcast(str);

    if (str != "master") {
        bIsPassed = false;
    }

    showResultMessageAll("testBroadcast_string", bIsPassed);
}

void testBroadcast_vectorString() {
    bool bIsPassed = true;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    std::vector<std::string> aStr;
    if (rComm.isMaster() == true) {
        aStr.resize(3);
        aStr[0] = "hoge";
        aStr[1] = "bar";
        aStr[2] = "foo";
    } else {
        aStr.resize(1);
        aStr[0] = "nothing";
    }

    rComm.broadcast(aStr);

    const int nMax = aStr.size();
    if (nMax != 3) {
        bIsPassed = false;
    } else {
        if ((aStr[0] != "hoge") || (aStr[1] != "bar") || (aStr[2] != "foo")) {
            bIsPassed = false;
        }
    }

    showResultMessageAll("testBroadcast_vectorString", bIsPassed);
}

void testBroadcast_valarray() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    std::valarray<double> v;
    if (rComm.isMaster() == true) {
        v.resize(10, 0.0);

        for (int i = 0; i < 10; ++i) {
            v[i] = 1.0 * i;
        }
    }

    rComm.broadcast(v);

    if (v.size() == 10) {
        bool tmp = true;
        for (int i = 0; i < 10; ++i) {
            if (std::fabs(v[i] - i) > 1.0E-16) {
                tmp = false;
            }
        }

        if (tmp == true) {
            bIsPassed = true;
        }
    }

    showResultMessageAll("testBroadcast_valarray", bIsPassed);
}

void testBroadcast_TlVector() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    TlVector v;
    if (rComm.isMaster() == true) {
        v.resize(10);
        for (int i = 0; i < 10; ++i) {
            v[i] = 1.0 * i;
        }
    }

    rComm.broadcast(v);

    if (v.getSize() == 10) {
        bool tmp = true;
        for (int i = 0; i < 10; ++i) {
            if (std::fabs(v[i] - i) > 1.0E-16) {
                tmp = false;
            }
        }

        if (tmp == true) {
            bIsPassed = true;
        }
    }

    showResultMessageAll("testBroadcast_TlVector", bIsPassed);
}

void testBroadcast_TlSparseMatrix() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    TlSparseMatrix sm(9, 9);
    if (rComm.isMaster() == true) {
        sm(0, 0) = 1.0;
        sm(3, 1) = 2.0;
        sm(5, 2) = 4.0;
    }

    rComm.broadcast(sm);

    if ((std::fabs(sm(0, 0) - 1.0) < 1.0E-16) &&
        (std::fabs(sm(3, 1) - 2.0) < 1.0E-16) &&
        (std::fabs(sm(5, 2) - 4.0) < 1.0E-16)) {
        bIsPassed = true;
    }

    //   for (int i = 0; i < rComm.getNumOfProc(); ++i) {
    //     if (i == rComm.getRank()) {
    //       sm.print(std::cout);
    //     }
    //     rComm.barrier();
    //   }

    showResultMessageAll("testBroadcast_TlSparseMatrix", bIsPassed);
}

void testBroadcast_TlMmapSymmetricMatrix() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    TlMatrix::useMemManager(true);
    TlDenseSymmetricMatrix_BLAS_Old sm(9);
    if (rComm.isMaster() == true) {
        sm(0, 0) = 1.0;
        sm(3, 1) = 2.0;
        sm(5, 2) = 4.0;
    }
    rComm.broadcast(sm);

    if ((std::fabs(sm(0, 0) - 1.0) < 1.0E-16) &&
        (std::fabs(sm(3, 1) - 2.0) < 1.0E-16) &&
        (std::fabs(sm(5, 2) - 4.0) < 1.0E-16)) {
        bIsPassed = true;
    }

    for (int i = 0; i < rComm.getNumOfProc(); ++i) {
        if (i == rComm.getRank()) {
            sm.print(std::cout);
        }
        rComm.barrier();
    }
    TlMatrix::useMemManager(false);

    showResultMessageAll("testBroadcast_TlMmapSymmetricMatrix", bIsPassed);
}

void testAllreduceSum_double() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    double value = double(rank);

    rComm.allReduce_SUM(value);

    double exact = 0.0;
    for (int i = 0; i < proc; ++i) {
        exact += i;
    }

    // bool check = true;
    if (std::fabs(value - exact) > 1.0E-12) {
        std::cout << TlUtils::format(
                         "testAllreduceSum_double: value=%f, exact=%f", value,
                         exact)
                  << std::endl;
        bIsPassed = false;
    } else {
        bIsPassed = true;
    }

    showResultMessageAll("testAllreduceSum_double", bIsPassed);
}

void testAllreduceSum_vectorDouble() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    int nSize = 5;
    std::vector<double> s_data(nSize);
    for (int i = 0; i < nSize; ++i) {
        s_data[i] = 1.0 * (rank + 1);
    }

    rComm.allReduce_SUM(s_data);

    bool check = true;
    for (int i = 0; i < 5; ++i) {
        if (std::fabs(s_data[i] - ((1 + proc) * proc / 2)) > 1.0E-16) {
            check = false;
        }
    }
    if (check == true) {
        bIsPassed = true;
    }

    showResultMessageAll("testAllreduceSum_vectorDouble", bIsPassed);
}

void testAllreduceSum_TlVector() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    int nSize = 100;
    TlVector s_data(nSize);
    for (int i = 0; i < nSize; ++i) {
        s_data[i] = 1.0 * (rank + 1);
    }

    rComm.allReduce_SUM(s_data);

    bool check = true;
    for (int i = 0; i < 5; ++i) {
        if (std::fabs(s_data[i] - ((1 + proc) * proc / 2)) > 1.0E-16) {
            check = false;
        }
    }
    if (check == true) {
        bIsPassed = true;
    }

    showResultMessageAll("testAllreduceSum_TlVector", bIsPassed);
}

void testAllreduceSum_matrix_symmetric() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    int nSize = 5;
    TlDenseSymmetricMatrix_BLAS_Old s_data(nSize);

    for (int x = 0; x < nSize; ++x) {
        for (int y = x; y < nSize; ++y) {
            s_data(x, y) = rank + 1.0;
        }
    }

    rComm.allReduce_SUM(s_data);

    bool check = true;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (std::fabs(s_data(i, j) - ((1 + proc) * proc / 2)) > 1.0E-16) {
                check = false;
            }
        }
    }
    if (check == true) {
        bIsPassed = true;
    }

    showResultMessageAll("testAllreduceSum_matrix_symmetric", bIsPassed);
}

void testAllreduceSum_TlMmapSymmetricMatrix() {
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    int nSize = 1000;
    TlMatrix::useMemManager(true);
    TlDenseSymmetricMatrix_BLAS_Old s_data(nSize);

    for (int x = 0; x < nSize; ++x) {
        for (int y = x; y < nSize; ++y) {
            s_data(x, y) = rank + 1.0;
        }
    }

    rComm.allReduce_SUM(s_data);

    bool check = true;
    for (int i = 0; i < nSize; ++i) {
        for (int j = 0; j < nSize; ++j) {
            if (std::fabs(s_data(i, j) - ((1 + proc) * proc / 2)) > 1.0E-16) {
                check = false;
            }
        }
    }
    if (check == true) {
        bIsPassed = true;
    }

    TlMatrix::useMemManager(false);
    showResultMessageAll("testAllreduceSum_TlMmapSymmetricMatrix", bIsPassed);
}

void test_nonblocking() {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    static const int VECTOR_SIZE = 10;

    if (rComm.getRank() == 1) {
        int foo = 12;
        std::vector<int> v(VECTOR_SIZE);
        for (int i = 0; i < VECTOR_SIZE; ++i) {
            v[i] = i;
        }

        rComm.iSendData(foo, 0);
        rComm.iSendData(v, 0);

        rComm.wait(foo);
        rComm.wait(v);

    } else if (rComm.getRank() == 0) {
        int foo = 0;
        std::vector<int> v;

        rComm.iReceiveData(foo, 1);
        rComm.iReceiveData(v, 1);

        rComm.wait(foo);
        rComm.wait(v);
    }
}
