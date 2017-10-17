#include "gtest/gtest.h"
#include "TlFileSymmetricMatrix.h"
#include "TlFile.h"
#include "TlSymmetricMatrix.h"

// const double EPS = 1.0E-10; // std::numeric_limits<double>::epsilon();
static const std::string mat_path = "temp.filesymmat.mat";

static void cleanup()
{
    if (TlFile::isExistFile(mat_path)) {
        TlFile::remove(mat_path);
    }
}


TEST(TlFileSymmetricMatrix, constructer)
{
    {
        cleanup();
        TlFileSymmetricMatrix fsm(mat_path, 100);
    }

    {
        TlSymmetricMatrix a;
        a.load(mat_path);
        EXPECT_EQ(100, a.getNumOfRows());
        EXPECT_EQ(100, a.getNumOfCols());
        for (int r = 0; r < 100; ++r) {
            for (int c = 0; c <= r; ++c) {
                EXPECT_DOUBLE_EQ(0.0, a.get(r, c));
            }
        }
    }
}


TEST(TlFileSymmetricMatrix, constructer_by_existed)
{
    TlSymmetricMatrix sm(100);
    {
        sm( 3, 17) =  51.0;
        sm( 0, 28) = - 1.0;
        sm.save(mat_path);
    }
    
    TlFileSymmetricMatrix fsm(mat_path);
}


TEST(TlFileSymmetricMatrix, get)
{
    TlSymmetricMatrix sm(100);
    sm( 3, 17) =  51.0;
    sm( 0, 28) = - 1.0;
    sm.save(mat_path);
    
    TlFileSymmetricMatrix fsm(mat_path);
    
    EXPECT_DOUBLE_EQ( 0.0, fsm.get( 0,  0));
    EXPECT_DOUBLE_EQ( 0.0, fsm.get( 0, 99));
    EXPECT_DOUBLE_EQ( 0.0, fsm.get(99,  0));
    EXPECT_DOUBLE_EQ( 0.0, fsm.get(99, 99));
    EXPECT_DOUBLE_EQ(51.0, fsm.get( 3, 17));
    EXPECT_DOUBLE_EQ(-1.0, fsm.get( 0, 28));
}


TEST(TlFileSymmetricMatrix, set)
{
    cleanup();
    {
        TlFileSymmetricMatrix fsm(mat_path, 100);
        
        fsm.set( 3, 17, 51.0);
        fsm.set( 0, 28, -1.0);
    }
    
    TlSymmetricMatrix sm;
    sm.load(mat_path);
    
    EXPECT_EQ(100, sm.getNumOfRows());
    EXPECT_EQ(100, sm.getNumOfCols());
    EXPECT_DOUBLE_EQ( 0.0, sm( 0,  0));
    EXPECT_DOUBLE_EQ( 0.0, sm( 0, 99));
    EXPECT_DOUBLE_EQ( 0.0, sm(99,  0));
    EXPECT_DOUBLE_EQ( 0.0, sm(99, 99));
    EXPECT_DOUBLE_EQ(51.0, sm( 3, 17));
    EXPECT_DOUBLE_EQ(-1.0, sm( 0, 28));    
}


TEST(TlFileSymmetricMatrix, add)
{
    {
        TlSymmetricMatrix sm(100);
        sm( 3, 17) =  51.0;
        sm( 0, 28) = - 1.0;
        sm.save(mat_path);
    }
    
    {
        TlFileSymmetricMatrix fsm(mat_path);
        fsm.add(3, 17, 12.3);
    }
    
    TlSymmetricMatrix sm;
    sm.load(mat_path);
    EXPECT_EQ(100, sm.getNumOfRows());
    EXPECT_EQ(100, sm.getNumOfCols());
    EXPECT_DOUBLE_EQ( 0.0, sm( 0,  0));
    EXPECT_DOUBLE_EQ( 0.0, sm( 0, 99));
    EXPECT_DOUBLE_EQ( 0.0, sm(99,  0));
    EXPECT_DOUBLE_EQ( 0.0, sm(99, 99));
    EXPECT_DOUBLE_EQ(63.3, sm( 3, 17));
    EXPECT_DOUBLE_EQ(-1.0, sm( 0, 28));    
}


TEST(TlFileSymmetricMatrix, resize)
{
    cleanup();

    // prepare
    TlFileSymmetricMatrix m(mat_path, 100);
    m.set(30, 40, 3.14);
    EXPECT_EQ(100, m.getNumOfRows());
    EXPECT_EQ(100, m.getNumOfCols());

    // resize
    m.resize(150);
    EXPECT_EQ(150, m.getNumOfRows());
    EXPECT_EQ(150, m.getNumOfCols());
    EXPECT_DOUBLE_EQ(3.14, m.get(30, 40));
    EXPECT_DOUBLE_EQ(3.14, m.get(40, 30));
}


