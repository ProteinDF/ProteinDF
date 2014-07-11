#include <cstdlib>
#include <iostream>
#include <string>

#include "ProteinDF.h"

#include "TlSymmetricMatrix.h"
#include "TlSparseMatrix.h"
#include "TlSparseSymmetricMatrix.h"

#include "TlCommunicate.h"
#include "TlMemManager.h"
#include "TlFile.h"

void showResultMessage(const std::string& sFunction, bool bIsPassed);
void showResultMessageAll(const std::string& sFunction, bool bIsPassed);

void checkAllProcess();
void testTlLog_Parallel();
void testCommunicate_sendInt();
void testCommunicate_sendLong();
void testCommunicate_sendDouble();
void testCommunicate_sendVectorInt();
void testCommunicate_sendVectorLong();
void testCommunicate_sendVectorDouble();
void testCommunicate_sendString();

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

int main(int argc, char *argv[])
{
    // initialize
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

    TlMemManager::setParam(1024 * 1024 * 1024, // 1 GB
                           "/tmp/parallelTest");

    // ===================================================================
    checkAllProcess();

    testCommunicate_sendInt();
    testCommunicate_sendLong();
    testCommunicate_sendDouble();
    testCommunicate_sendVectorInt();
    testCommunicate_sendVectorLong();
    testCommunicate_sendVectorDouble();
    testCommunicate_sendString();

    //testGatherToMaster_TlSparseMatrix();

    testBroadcast_bool();
    testBroadcast_int();
    testBroadcast_string();
    testBroadcast_vectorString();
    testBroadcast_valarray();
    testBroadcast_TlVector();
    testBroadcast_TlSparseMatrix();
    //testBroadcast_TlMmapSymmetricMatrix();

    testAllreduceSum_double();
    testAllreduceSum_vectorDouble();
    testAllreduceSum_TlVector();
    testAllreduceSum_matrix_symmetric();
    //testAllreduceSum_TlMmapSymmetricMatrix();
    // ===================================================================

    // finalize
    rComm.finalize();
    return EXIT_SUCCESS;
}

void showResultMessage(const std::string& sFunction, bool bIsPassed)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        std::cout << "TEST: " << sFunction << "() ";
        std::cout << ((bIsPassed == true) ? "." : "F") << std::endl;
    }
    rComm.barrier();
}

void showResultMessageAll(const std::string& sFunction, bool bIsPassed)
{
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

void checkAllProcess()
{
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    //const int rank = rComm.getRank();

    for (int i = 0; i < proc; ++i) {
        rComm.barrier();
    }

    showResultMessage("checkAllProcess", bIsPassed);
}

// void testTlLog_Parallel(){
//   TlLog_Parallel info(TlLog::INFO);

//   info << "this message is written by only master process." << std::endl;
// }

void testCommunicate_sendInt()
{
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            rComm.sendData(100 +i);
        } else if (rComm.isMaster()) {
            int rData = 0;
            rComm.receiveData(rData, i);

            if (rData != (100 +i)) {
                std::cout << "recv(int) = " << rData << ", expect = " << (100 +i) << std::endl;
                bIsPassed = false;
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendInt", bIsPassed);
}

void testCommunicate_sendLong()
{
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

            if (rData != (1000 +i)) {
                std::cout << "recv(long) = " << rData << ", expect = " << (1000 +i) << std::endl;
                bIsPassed = false;
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendLong", bIsPassed);
}

void testCommunicate_sendString()
{
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
                std::cout << "recv(string) = " << rData << ", expect = " << sExpectData << std::endl;
                bIsPassed = false;
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendString", bIsPassed);
}

void testCommunicate_sendDouble()
{
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            rComm.sendData(double(100 +i));
        } else if (rComm.isMaster()) {
            double rData = 0.0;
            rComm.receiveData(rData, i);

            if (std::fabs(rData - (100 +i)) > 1.0E-16) {
                std::cout << "recv(double) = " << rData << ", expect = " << (100 +i) << std::endl;
                bIsPassed = false;
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendDouble", bIsPassed);
}

void testCommunicate_sendVectorInt()
{
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            std::vector<int> v(proc);
            for (int j = 0; j < proc; ++j) {
                v[j] = 100 +j;
            }
            rComm.sendData(v);
        } else if (rComm.isMaster()) {
            std::vector<int> rData;
            rComm.receiveData(rData, i);

            for (int j = 0; j < proc; ++j) {
                if (rData[j] != (100 +j)) {
                    std::cout << "recv(vector int) = " << rData[j] << ", expect = " << (100 +j) << std::endl;
                    bIsPassed = false;
                }
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendVectorInt", bIsPassed);
}

void testCommunicate_sendVectorLong()
{
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            std::vector<long> v(proc);
            for (int j = 0; j < proc; ++j) {
                v[j] = 100 +j;
            }
            rComm.sendData(v);
        } else if (rComm.isMaster()) {
            std::vector<long> rData;
            rComm.receiveData(rData, i);

            for (int j = 0; j < proc; ++j) {
                if (rData[j] != (100 +j)) {
                    std::cout << "recv(vector long) = " << rData[j] << ", expect = " << (100 +j) << std::endl;
                    bIsPassed = false;
                }
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendVectorLong", bIsPassed);
}

void testCommunicate_sendVectorDouble()
{
    bool bIsPassed = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    for (int i = 1; i < proc; i++) {
        if (rank == i) {
            std::vector<double> v(proc);
            for (int j = 0; j < proc; ++j) {
                v[j] = double(100 +j);
            }
            rComm.sendData(v);
        } else if (rComm.isMaster()) {
            std::vector<double> rData;
            rComm.receiveData(rData, i);

            for (int j = 0; j < proc; ++j) {
                if (std::fabs(rData[j] - (100 +j)) > 1.0E-16) {
                    std::cout << "recv(vector double) = " << rData[j] << ", expect = " << (100 +j) << std::endl;
                    bIsPassed = false;
                }
            }
        }

        rComm.barrier();
    }

    showResultMessage("testCommunicate_sendVectorDouble", bIsPassed);
}

void testGatherToMaster_TlSparseMatrix()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    // $B=`Hw(B
    TlSparseSymmetricMatrix sm(proc);
    sm(rank, rank) = rank +1;

    // test
    rComm.gatherToMaster(sm);

    bool bIsPassed = true;
    if (rComm.isMaster() == true) {
        for (int i = 0; i < proc; ++i) {
            // $BBP3Q9`(B
            const double check = std::fabs(sm(i, i) - (i +1.0));
            if (check > 1.0E-16) {
                bIsPassed = false;
            }

            // $BHsBP3Q9`(B
            for (int j = 0; j < i; ++j) {
                if (std::fabs(sm(i, j) - 0.0) > 1.0E-16) {
                    bIsPassed = false;
                }
            }
            for (int j = i+1; j < proc; ++j) {
                if (std::fabs(sm(i, j) - 0.0) > 1.0E-16) {
                    bIsPassed = false;
                }
            }
        }
    }

    showResultMessage("testGatherToMaster_TlSparseMatrix", bIsPassed);
}

void testBroadcast_bool()
{
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

void testBroadcast_int()
{
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

void testBroadcast_string()
{
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

void testBroadcast_vectorString()
{
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
        if ((aStr[0] != "hoge") ||
                (aStr[1] != "bar") ||
                (aStr[2] != "foo")) {
            bIsPassed = false;
        }
    }

    showResultMessageAll("testBroadcast_vectorString", bIsPassed);
}

void testBroadcast_valarray()
{
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

void testBroadcast_TlVector()
{
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

void testBroadcast_TlSparseMatrix()
{
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


void testBroadcast_TlMmapSymmetricMatrix()
{
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // prepare
    TlMatrix::useMemManager(true);
    TlSymmetricMatrix sm(9);
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


void testAllreduceSum_double()
{
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
    
    bool check = true;
    if (std::fabs(value - exact) > 1.0E-12) {
        std::cout << TlUtils::format("testAllreduceSum_double: value=%f, exact=%f",
                                     value, exact)
                  << std::endl;
        bIsPassed = false;
    } else {
        bIsPassed = true;
    }

    showResultMessageAll("testAllreduceSum_double", bIsPassed);
}



void testAllreduceSum_vectorDouble()
{
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    int nSize = 5;
    std::vector<double> s_data(nSize);
    for (int i = 0; i < nSize; ++i) {
        s_data[i] = 1.0 * (rank +1);
    }

    rComm.allReduce_SUM(s_data);

    bool check = true;
    for (int i = 0; i < 5; ++i) {
        if (std::fabs(s_data[i] - ((1+proc)*proc/2)) > 1.0E-16) {
            check = false;
        }
    }
    if (check == true) {
        bIsPassed = true;
    }

    showResultMessageAll("testAllreduceSum_vectorDouble", bIsPassed);
}

void testAllreduceSum_TlVector()
{
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    int nSize = 100;
    TlVector s_data(nSize);
    for (int i = 0; i < nSize; ++i) {
        s_data[i] = 1.0 * (rank +1);
    }

    rComm.allReduce_SUM(s_data);

    bool check = true;
    for (int i = 0; i < 5; ++i) {
        if (std::fabs(s_data[i] - ((1+proc)*proc/2)) > 1.0E-16) {
            check = false;
        }
    }
    if (check == true) {
        bIsPassed = true;
    }

    showResultMessageAll("testAllreduceSum_TlVector", bIsPassed);
}


void testAllreduceSum_matrix_symmetric()
{
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    int nSize = 5;
    TlSymmetricMatrix s_data(nSize);

    for (int x = 0; x < nSize; ++x) {
        for (int y = x; y < nSize; ++y) {
            s_data(x, y) = rank +1.0;
        }
    }

    rComm.allReduce_SUM(s_data);

    bool check = true;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (std::fabs(s_data(i, j) - ((1+proc)*proc/2)) > 1.0E-16) {
                check = false;
            }
        }
    }
    if (check == true) {
        bIsPassed = true;
    }

    showResultMessageAll("testAllreduceSum_matrix_symmetric", bIsPassed);
}


void testAllreduceSum_TlMmapSymmetricMatrix()
{
    bool bIsPassed = false;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getNumOfProc();
    const int rank = rComm.getRank();

    int nSize = 1000;
    TlMatrix::useMemManager(true);
    TlSymmetricMatrix s_data(nSize);

    for (int x = 0; x < nSize; ++x) {
        for (int y = x; y < nSize; ++y) {
            s_data(x, y) = rank +1.0;
        }
    }

    rComm.allReduce_SUM(s_data);

    bool check = true;
    for (int i = 0; i < nSize; ++i) {
        for (int j = 0; j < nSize; ++j) {
            if (std::fabs(s_data(i, j) - ((1+proc)*proc/2)) > 1.0E-16) {
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


void test_nonblocking()
{
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
