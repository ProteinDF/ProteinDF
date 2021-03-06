#include <cppunit/TextOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "TlDenseSymmetricMatrix_BLAS_OldTest.h"
#include "TlFileMatrixTest.h"
#include "TlFileSymmetricMatrixTest.h"
#include "TlMatrixTest.h"
#include "TlMemManagerTest.h"
#include "TlMsgPackTest.h"
#include "TlPartialSymmetricMatrixTest.h"
#include "TlSerializeDataTest.h"
#include "TlSparseMatrixTest.h"
#include "TlSparseSymmetricMatrixTest.h"
#include "TlSparseVectorTest.h"
#include "TlStlUtilsTest.h"
#include "TlStringTokenizerTest.h"
#include "TlUtilsTest.h"
#include "TlVectorTest.h"

#include "TlDenseGeneralMatrix_arrays_ColOrientedTest.h"
#include "TlDenseGeneralMatrix_arrays_RowOrientedTest.h"

//#include "Fl_UserinputXTest.h"

#include "DfFunctional_B3LYPTest.h"
#include "DfFunctional_Becke88Test.h"
#include "DfFunctional_LYPTest.h"
#include "DfFunctional_PW91XTest.h"
#include "DfFunctional_SlaterTest.h"
#include "DfFunctional_VWNTest.h"

int main(int argc, char* argv[]) {
    // TestRunnerを生成してrun()を実行する
    CppUnit::TextUi::TestRunner runner;

    // ここにテストを追加していく
    // runner.addTest(CppUnit::TestFactoryRegistry::getRegistry().makeTest());

    // Tl class
    runner.addTest(TlUtilsTest::suite());
    runner.addTest(TlStringTokenizerTest::suite());
    runner.addTest(TlMemManagerTest::suite());
    runner.addTest(TlMsgPackTest::suite());
    runner.addTest(TlSerializeDataTest::suite());

    runner.addTest(TlVectorTest::suite());
    runner.addTest(TlMatrixTest::suite());
    runner.addTest(TlDenseSymmetricMatrix_BLAS_OldTest::suite());
    runner.addTest(TlSparseMatrixTest::suite());
    runner.addTest(TlSparseSymmetricMatrixTest::suite());
    runner.addTest(TlSparseVectorTest::suite());
    runner.addTest(TlFileMatrixTest::suite());
    runner.addTest(TlFileSymmetricMatrixTest::suite());
    runner.addTest(TlPartialSymmetricMatrixTest::suite());
    runner.addTest(TlStlUtilsTest::suite());

    runner.addTest(TlDenseGeneralMatrix_arrays_RowOrientedTest::suite());
    runner.addTest(TlDenseGeneralMatrix_arrays_ColOrientedTest::suite());

    //
    // runner.addTest(Fl_UserinputXTest::suite());

    // Df class
    runner.addTest(DfFunctional_SlaterTest::suite());
    runner.addTest(DfFunctional_Becke88Test::suite());
    runner.addTest(DfFunctional_VWNTest::suite());
    runner.addTest(DfFunctional_LYPTest::suite());
    runner.addTest(DfFunctional_B3LYPTest::suite());
    runner.addTest(DfFunctional_PW91XTest::suite());

    CppUnit::Outputter* outputter =
        new CppUnit::TextOutputter(&runner.result(), std::cout);
    runner.setOutputter(outputter);

    return runner.run() ? 0 : 1;
}
