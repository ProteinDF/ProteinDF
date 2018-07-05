#ifndef TL_SCALAPACK_CONTEXT_H
#define TL_SCALAPACK_CONTEXT_H

class TlScalapackContext {
 public:
  static void getData(int& rContext, int& rProc, int& rRank, int& rProcGridRow,
                      int& rProcGridCol);

  static int getBlockSize();
  static void setBlockSize(int blockSize);

  static void finalize();

 private:
  TlScalapackContext();
  ~TlScalapackContext();

 private:
  static TlScalapackContext* m_pTlScalapackContextInstance;
  static int m_nContext;
  static int m_nProc;
  static int m_nRank;
  static int m_nProcGridRow;
  static int m_nProcGridCol;

  static int systemBlockSize_;
};

#endif  // TL_SCALAPACK_CONTEXT_H
