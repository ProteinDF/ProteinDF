#include <vector>

#include "TlCommunicate.h"
#include "TlMath.h"
#include "scalapack.h"
#include "tl_scalapack_context.h"

TlScalapackContext* TlScalapackContext::m_pTlScalapackContextInstance = NULL;
int TlScalapackContext::m_nContext = 0;
int TlScalapackContext::m_nProc = 0;
int TlScalapackContext::m_nRank = 0;
int TlScalapackContext::m_nProcGridRow = 0;
int TlScalapackContext::m_nProcGridCol = 0;
int TlScalapackContext::systemBlockSize_ = 64;

TlScalapackContext::TlScalapackContext() {
  TlCommunicate& rComm = TlCommunicate::getInstance();

  // initialize
  Cblacs_pinfo(&(this->m_nRank), &(this->m_nProc));
  if (this->m_nProc < 1) {
    if (this->m_nRank == 0) {
      this->m_nProc = rComm.getNumOfProc();
    }
    Cblacs_setup(&(this->m_nRank), &(this->m_nProc));
  }
  Cblacs_get(-1, 0, &(this->m_nContext));

  // // process grid
  TlScalapackContext::m_nProc = rComm.getNumOfProc();
  {
    std::vector<int> f = TlMath::factor(TlScalapackContext::m_nProc);
    int r = 1;
    int c = 1;
    int f_size = f.size();
    for (int i = 0; i < f_size; ++i) {
      if ((i % 2) == 0) {
        c *= f[i];
      } else {
        r *= f[i];
      }
    }
    TlScalapackContext::m_nProcGridRow = r;
    TlScalapackContext::m_nProcGridCol = c;
  }

  // "Row-major" is default
  Cblacs_gridinit(&(TlScalapackContext::m_nContext), "Row-major",
                  TlScalapackContext::m_nProcGridRow,
                  TlScalapackContext::m_nProcGridCol);
}

TlScalapackContext::~TlScalapackContext() {
  Cblacs_gridexit(TlScalapackContext::m_nContext);
  TlScalapackContext::m_nContext = 0;
  Cblacs_exit(1);
}

void TlScalapackContext::getData(int& rContext, int& rProc, int& rRank,
                                 int& rProcGridRow, int& rProcGridCol) {
  if (TlScalapackContext::m_pTlScalapackContextInstance == NULL) {
    TlScalapackContext::m_pTlScalapackContextInstance =
        new TlScalapackContext();
  }
  assert(TlScalapackContext::m_pTlScalapackContextInstance != NULL);

  rContext = TlScalapackContext::m_nContext;
  rProc = TlScalapackContext::m_nProc;
  rRank = TlScalapackContext::m_nRank;
  rProcGridRow = TlScalapackContext::m_nProcGridRow;
  rProcGridCol = TlScalapackContext::m_nProcGridCol;
}

int TlScalapackContext::getBlockSize() {
  return TlScalapackContext::systemBlockSize_;
}

void TlScalapackContext::setBlockSize(int blockSize) {
  TlScalapackContext::systemBlockSize_ = blockSize;
}

void TlScalapackContext::finalize() {
  if (TlScalapackContext::m_pTlScalapackContextInstance != NULL) {
    delete TlScalapackContext::m_pTlScalapackContextInstance;
    TlScalapackContext::m_pTlScalapackContextInstance = NULL;
  }
}
