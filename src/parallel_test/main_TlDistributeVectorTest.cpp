#include <cstdlib>
#include <iostream>
#include <string>

#include "TlTime.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlDistributeVector.h"
#include "TlVector.h"
#include "TlCommunicate.h"

void showResultMessage(const std::string& sFunction, bool bIsPassed);
void showResultMessageAll(const std::string& sFunction, bool bIsPassed);

void testConstructer();
void testCopyConstructer();
void testSet();
void testOperatorPlusEqual();
void testOperatorMultiEqual();
void testSave();
void testLoad();
void testMatrixVectorOperation();
void testVectorVectorOperation();

bool check(double actual, double expect, double range,
           const std::string& fileStr, int line)
{
    if (std::fabs(actual - expect) > range) {
        std::string str = TlUtils::format("[FAIL] %s:%d actual=%f, expect=%f",
                                          fileStr.c_str(), line,
                                          actual, expect);
        std::cerr << str << std::endl;
        return false;
    } else {
        return true;
    }
}

int main(int argc, char *argv[])
{
    // initialize
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);
    TlDistributeVector::setSystemBlockSize(10);
    
    // ===================================================================
    testConstructer();
    testCopyConstructer();
    testSet();
    testOperatorPlusEqual();
    testOperatorMultiEqual();
    testSave();
    testLoad();
    testMatrixVectorOperation();
    testVectorVectorOperation();
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


void testConstructer()
{
    TlDistributeVector v(100);

    std::size_t size = v.getSize();

    bool bIsPassed = (size == 100) ? true : false;

    showResultMessageAll("testConstructer", bIsPassed);
}


void testCopyConstructer()
{
    TlDistributeVector v(100);
    v[10] =  10.0;
    v[17] =  56.0;
    v[21] = - 2.5;
    v[90] =  12.0;

    const TlDistributeVector w(v);

    bool bIsPassed = true;
    if (w.getSize() != 100) {
        bIsPassed = false;
    }
    if (std::fabs(w.get(10) -  10.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(w.get(17) -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(w.get(21) - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(w.get(90) -   12.0) > 1.0E-16) {
        bIsPassed = false;
    }

    showResultMessageAll("testCopyConstructer", bIsPassed);
}


void testSet()
{
    const int size = 100;
    TlDistributeVector v(size);

    for (int i = 0; i < size; ++i) {
        v.set(i, rand());
    }

    v[10] = 10.0;
    v[17] = 56.0;
    v[21] =- 2.5;
    v[90] = 12.0;

    bool bIsPassed = true;
    if (v.getSize() != size) {
        bIsPassed = false;
    }
    if (std::fabs(v[10] -  10.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(v[17] -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(v[21] - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(v[90] -  12.0) > 1.0E-16) {
        bIsPassed = false;
    }

    showResultMessageAll("testSet", bIsPassed);
}


void testOperatorPlusEqual()
{
    TlDistributeVector v(100);
    v[10] = 10.0;
    v[17] = 56.0;
    v[21] =- 2.5;
    v[90] = 12.0;

    TlDistributeVector w(100);
    w[10] =   5.0;
    w[20] = -12.0;
    w[30] =  51.0;

    v += w;

    bool bIsPassed = true;
    if (v.getSize() != 100) {
        bIsPassed = false;
    }
    if (std::fabs(v[10] -  15.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(v[17] -  56.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(v[20] - (-12.0)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(v[21] - (-2.5)) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(v[30] -  51.0) > 1.0E-16) {
        bIsPassed = false;
    }
    if (std::fabs(v[90] -  12.0) > 1.0E-16) {
        bIsPassed = false;
    }

    showResultMessageAll("testOperatorPlusEqual", bIsPassed);
}


void testOperatorMultiEqual()
{
    TlDistributeVector v(100);
    v[10] = 10.0;
    v[17] = 56.0;
    v[21] =- 2.5;
    v[90] = 12.0;

    v *= 3.0;
    
    bool bIsPassed = true;
    if (v.getSize() != 100) {
        bIsPassed = false;
    }

    bIsPassed &= check(v[10], 10.0*3.0, 1.0E-16, __FILE__, __LINE__);
    bIsPassed &= check(v[17], 56.0*3.0, 1.0E-16, __FILE__, __LINE__);
    bIsPassed &= check(v[21], -2.5*3.0, 1.0E-16, __FILE__, __LINE__);
    bIsPassed &= check(v[90], 12.0*3.0, 1.0E-16, __FILE__, __LINE__);

    showResultMessageAll("testOperatorMultiEqual", bIsPassed);
}


void testSave()
{
    bool bIsPassed = true;

    const int size = 1000;
    TlDistributeVector v(size);
    TlVector w(size);

    int count = 0;
    for (int i = 0; i < size; ++i) {
        double d = double(count);
        v[i] = d;
        w[i] = d;
        ++count;
    }

    v.save("test.vector");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        w.save("test.TlVector.vector");
    }
    
    showResultMessageAll("testSave", bIsPassed);
    rComm.barrier();
}


void testLoad()
{
    bool bIsPassed = true;

    TlDistributeVector v;
    v.load("test.vector");
    TlVector w;
    w.load("test.TlVector.vector");

    if (v.getSize() != w.getSize()) {
        std::cout << "size = " << v.getSize() << std::endl;
        bIsPassed = false;
    }

    const int size = v.getSize();
    for (int i = 0; i < size; ++i) {
        bIsPassed &= check(v[i], w[i], 1.0E-16, __FILE__, __LINE__);
        if (bIsPassed == false) {
            break;
        }
    }

    showResultMessageAll("testLoad", bIsPassed);
}


void testMatrixVectorOperation()
{
    bool bIsPassed = true;

    TlSymmetricMatrix S;
    S.load("Sab.mtx");
    TlVector R;
    R.load("rho.vct");

    TlVector A = S * R;

    TlDistributeSymmetricMatrix s;
    s.load("Sab.mtx");
    TlDistributeVector r;
    r.load("rho.vct");

    TlDistributeVector a = s * r;

    const int size = A.getSize();
    for (int i = 0; i < size; ++i) {
        bIsPassed &= check(a[i], A[i], 1.0E-10, __FILE__, __LINE__);
        if (bIsPassed == false) {
            break;
        }
    }

    showResultMessageAll("MatrixVectorOperation", bIsPassed);
}


void testVectorVectorOperation()
{
    bool bIsPassed = true;

    TlVector R;
    R.load("rho.vct");

    double A = R * R;

    TlDistributeVector r;
    r.load("rho.vct");

    double a = r * r;

    bIsPassed &= check(a, A, 1.0E-10, __FILE__, __LINE__);

    showResultMessageAll("VectorVectorOperation", bIsPassed);
}

