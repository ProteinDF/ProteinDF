#ifndef TL_DENSE_GENERAL_MATRIX_IMPL_SCALAPACK_H
#define TL_DENSE_GENERAL_MATRIX_IMPL_SCALAPACK_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

#include <bitset>
#include <list>
#include <vector>

#include "tl_dense_matrix_impl_object.h"
#include "tl_dense_vector_object.h"
#include "tl_sparse_matrix.h"

class TlDenseGeneralMatrixObject;
class TlDenseSymmetricMatrix_ImplScalapack;
class TlDenseVector_ImplLapack;
class TlDenseVector_ImplScalapack;
class TlDenseMatrix_IO_object;

#define MAX_SESSION_ID (100)

class TlDenseGeneralMatrix_ImplScalapack : public TlDenseMatrix_ImplObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_ImplScalapack(
      TlMatrixObject::index_type row = 1, TlMatrixObject::index_type col = 1);
  TlDenseGeneralMatrix_ImplScalapack(
      const TlDenseGeneralMatrix_ImplScalapack& rhs);
  // TlDenseGeneralMatrix_ImplScalapack(const
  // TlDenseSymmetricMatrix_ImplScalapack& rhs);
  // TlDenseGeneralMatrix_ImplScalapack(const TlDistributedVector& rhs,
  //                                  const int row, const int col);
  // TlDenseGeneralMatrix_ImplScalapack(
  // const TlDenseGeneralMatrix_arrays_RowOriented& rhs);
  virtual ~TlDenseGeneralMatrix_ImplScalapack();

  // ---------------------------------------------------------------------------
  // static
  // ---------------------------------------------------------------------------
  // static void setSystemBlockSize(int blockSize);
  // static void setUsingPartialIO(bool isUsePIO);

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlMatrixObject::index_type getNumOfRows() const;
  virtual TlMatrixObject::index_type getNumOfCols() const;
  virtual TlMatrixObject::size_type getNumOfElements() const;

  virtual void resize(TlMatrixObject::index_type row,
                      TlMatrixObject::index_type col);

  virtual double get(TlMatrixObject::index_type row,
                     TlMatrixObject::index_type col) const;
  virtual void set(TlMatrixObject::index_type row,
                   TlMatrixObject::index_type col, double value);
  virtual void add(TlMatrixObject::index_type row,
                   TlMatrixObject::index_type col, double value);

  /// ローカルから該当する要素があれば値を返す
  ///
  /// 範囲外であれば0.0を返す
  /// すべてのプロセスが同時に呼ぶ必要はない
  virtual double getLocal(const TlMatrixObject::index_type row,
                          const TlMatrixObject::index_type col) const;

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  // need to call by all processes.
  TlDenseGeneralMatrix_ImplScalapack& operator=(
      const TlDenseGeneralMatrix_ImplScalapack& rhs);

  TlDenseGeneralMatrix_ImplScalapack& operator+=(
      const TlDenseGeneralMatrix_ImplScalapack& rhs);
  TlDenseGeneralMatrix_ImplScalapack& operator-=(
      const TlDenseGeneralMatrix_ImplScalapack& rhs);
  TlDenseGeneralMatrix_ImplScalapack& operator*=(double coef);
  TlDenseGeneralMatrix_ImplScalapack& operator/=(double coef);
  TlDenseGeneralMatrix_ImplScalapack& operator*=(
      const TlDenseGeneralMatrix_ImplScalapack& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  virtual std::vector<double> diagonals() const;

  virtual double sum() const;
  virtual double trace() const;
  virtual double getRMS() const;
  virtual double getMaxAbsoluteElement(
      TlMatrixObject::index_type* pOutRow = NULL,
      TlMatrixObject::index_type* pOutCol = NULL) const;

  TlDenseGeneralMatrix_ImplScalapack transpose() const;
  virtual void transposeInPlace();

  TlDenseGeneralMatrix_ImplScalapack dot(const TlDenseGeneralMatrix_ImplScalapack& rhs) const;
  const TlDenseGeneralMatrix_ImplScalapack& dotInPlace(
      const TlDenseGeneralMatrix_ImplScalapack& rhs);

  TlDenseGeneralMatrix_ImplScalapack inverse() const;

  // TlDenseGeneralMatrix_ImplScalapack getLeastSquaresSolution(
  //     const TlDenseGeneralMatrix_ImplScalapack& B) const;

  bool getSparseMatrix(TlSparseMatrix* pMatrix, bool isFinalize = false) const;

  /// 各ノードが与えた疎行列を大域行列に加算する。
  ///
  /// 全ノードがこの関数を呼び出す必要がある。
  void mergeSparseMatrix(const TlSparseMatrix& M);

  std::vector<TlMatrixObject::index_type> getRowIndexTable() const;
  std::vector<TlMatrixObject::index_type> getColIndexTable() const;
  void getLocalMatrix(TlDenseGeneralMatrixObject* pOutMatrix) const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  virtual bool load(const std::string& filePath);
  virtual bool save(const std::string& filePath) const;

  // #ifdef HAVE_HDF5
  //   virtual bool loadHdf5(const std::string& filepath, const std::string&
  //   h5path);
  //   virtual bool saveHdf5(const std::string& filepath,
  //                         const std::string& h5path) const;
  // #endif  // HAVE_HDF5

  static void setUsingPartialIO(bool isUsePIO);

 public:
//   void dump(TlDenseVector_ImplScalapack* v) const;
//   void restore(const TlDenseVector_ImplScalapack& v);

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------
 protected:
  virtual void initialize();
  virtual TlMatrixObject::size_type getNumOfMyElements() const;

  // ブロックサイズを返す
  int getBlockSize() const;

  // int getProcIdForIndex(const TlMatrixObject::index_type globalRow,
  //                       const TlMatrixObject::index_type globalCol) const;

  // ローカル行列のインデックスを返す。
  // 保持していない場合は-1を返す。
  TlMatrixObject::size_type getIndex0(
      TlMatrixObject::index_type globalRow,
      TlMatrixObject::index_type globalCol) const;

  // 担当rank, ローカル行列のインデックスを返す。
  void getGlobalRowCol2LocalRowCol(
      const TlMatrixObject::index_type globalRow,
      const TlMatrixObject::index_type globalCol, int* pRank,
      TlMatrixObject::index_type* pLocalRow = NULL,
      TlMatrixObject::index_type* pLocalCol = NULL) const;

  // get global matrix index from local matrix index
  // void getLocalRowCol2GlobalRowCol(
  //     const TlMatrixObject::index_type localRow,
  //     const TlMatrixObject::index_type localCol,
  //     TlMatrixObject::index_type* pGlobalRow,
  //     TlMatrixObject::index_type* pGlobalCol) const;

  TlMatrixObject::index_type getLocal2Global_row(
      const TlMatrixObject::index_type localIndex) const;
  TlMatrixObject::index_type getLocal2Global_col(
      const TlMatrixObject::index_type localIndex) const;
  TlMatrixObject::index_type getLocal2Global(
      const TlMatrixObject::index_type localIndex, const int procID,
      const int numOfProcs) const;

  // get local matrix array index
  TlMatrixObject::index_type getLocalIndex(
      const TlMatrixObject::index_type localRow,
      const TlMatrixObject::index_type localCol) const;

  //
  /// ローカル行列における要素の絶対値の最大値を返す
  ///
  /// @param[out] outRow 該当要素の行
  /// @param[out] outCol 該当要素の列
  /// @return 要素の絶対値の最大値
  virtual double getLocalMaxAbsoluteElement(
      TlMatrixObject::index_type* pOutRow = NULL,
      TlMatrixObject::index_type* pOutCol = NULL) const;

  void getPartialMatrix_ServerTasks(const bool isFinalize) const;
  void getSparseMatrix_registerTask(TlSparseMatrix* pMatrix) const;
  int getPartialMatrix_getSessionId(const int proc) const;
  bool getPartialMatrix_ClientTasks(TlMatrixObject* pMatrix) const;

  void mergeMatrix_common(
      const std::vector<std::vector<TlDenseVectorObject::index_type> >&
          indexArrays,
      const std::vector<std::vector<double> >& values);

  virtual bool load(std::fstream& filePath);

  // #ifdef HAVE_HDF5
  //  protected:
  //   bool saveHdf5(const std::string& filepath, const std::string& h5path,
  //                 const int saveMatType) const;
  //   virtual size_type getArrayIndex(const index_type row,
  //                                   const index_type col) const;
  //   void saveElements(
  //       TlHdf5Utils* pH5, const std::string& path,
  //       const std::vector<TlMatrixObject::MatrixElement>& elements) const;
  // #endif  // HAVE_HDF5

  // std::vector<TlMatrixObject::MatrixElement> getMatrixElementsInLocal() const;
  virtual std::vector<TlMatrixObject::MatrixElement> getMatrixElementsInLocal2() const;

  void saveElements(
      TlDenseMatrix_IO_object* pFileMatrix,
      const std::vector<TlMatrixObject::MatrixElement>& elements) const;

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  static const std::size_t FILE_BUFFER_SIZE;
  static bool isUsingPartialIO;

 protected:
  int context_;
  int pDESC_[9];
  TlMatrixObject::index_type rows_;    // 大域行列の行数
  TlMatrixObject::index_type cols_;    // 大域行列の列数
  TlMatrixObject::index_type myRows_;  // ローカル行列の行数
  TlMatrixObject::index_type myCols_;  // ローカル行列の列数

  int rank_;         // プロセスのランク
  int proc_;         // 総プロセス数
  int procGridRow_;  // プロセスグリッドの行数
  int procGridCol_;  // プロセスグリッドの列数
  int m_nMyProcRow;  // プロセスグリッドにおける自分の行数
  int m_nMyProcCol;  // プロセスグリッドにおける自分の列数
  int blockSize_;    // ブロックサイズ

  // index table
  // usage: m_RowIndexTable[local_index] = global_index
  std::vector<TlMatrixObject::index_type> m_RowIndexTable;
  std::vector<TlMatrixObject::index_type> m_ColIndexTable;

  double* pData_;

 protected:
  // static int systemBlockSize_;

 protected:
  /// MPI通信タグ
  enum {
    // for copy constructor
    TAG_CONSTRUCTOR_SIZE = 10001,
    TAG_CONSTRUCTOR_INDEX = 10002,
    TAG_CONSTRUCTOR_VALUE = 10003,
    // for LOAD
    TAG_LOAD_SIZE = 10011,
    TAG_LOAD_ROWCOLS = 10012,
    TAG_LOAD_VALUES = 10013,
    TAG_LOAD_END = 10014,
    // for save
    TAG_SAVE_HANDSHAKE = 10021,
    TAG_SAVE_HANDSHAKE_OK = 10022,
    TAG_SAVE_DATA_ROWS = 10023,
    TAG_SAVE_DATA_COLS = 10024,
    TAG_SAVE_DATA_ROWINDEX = 10025,
    TAG_SAVE_DATA_COLINDEX = 10026,
    TAG_SAVE_DATA = 10027,
    // for getPartialMatrix / getSparseMatrix
    TAG_GPM_SESSION_ID = 10100,
    TAG_GPM_NUM_OF_COMPONENTS = 10200,
    TAG_GPM_COMPONENTS = 10300,
    TAG_GPM_ELEMENT_VALUES = 10400,
    // for MergeMatrix
    TAG_MERGE_MATRIX_SESSION_ID = 10500,
    TAG_MERGE_MATRIX_NUM_OF_CONTENTS = 10600,
    TAG_MERGE_MATRIX_INDECES = 10700,
    TAG_MERGE_MATRIX_VALUES = 10800
  };

  enum {
    TAG_TLDISTRIBUTE_MATRIX__MERGE__NUM_OF_INDECES = 2001,
    TAG_TLDISTRIBUTE_MATRIX__MERGE__INDECES = 2002,
    TAG_TLDISTRIBUTE_MATRIX__MERGE__VALUES = 2003
  };

 protected:
  // GPM
  // -----------------------------------------------------------------------
  bool isDebugOut_GPM_;

  struct GPM_ServerTask {
   public:
    GPM_ServerTask()
        : state(0), sessionId(0), requestProc(0), numOfComponents(0) {}

   public:
    unsigned int state;
    int sessionId;
    int requestProc;
    std::size_t numOfComponents;
    std::vector<TlMatrixObject::index_type> components;
    std::vector<double> elementValues;
  };

  struct GPM_ClientTask {
    std::vector<unsigned int> state;
    std::vector<int> sessionIds;

    /// プロセス毎へリクエストする要素番号情報
    /// プロセス番号を添字として、行番号、列番号の順に1次元配列に格納する
    std::vector<std::vector<TlMatrixObject::index_type> > components;

    /// プロセス毎へリクエストする要素番号情報の総数
    std::vector<std::size_t> numOfComponents;

    /// 受信用行列データ
    std::vector<std::vector<double> > elementValues;

    bool isFinished;
  };

  typedef std::list<GPM_ServerTask> GpmServerTasksType;
  mutable GpmServerTasksType gpmServerTasks_;

  typedef std::map<TlMatrixObject*, GPM_ClientTask> GpmClientTasks;
  mutable GpmClientTasks gpmClientTasks_;

  mutable std::vector<std::bitset<MAX_SESSION_ID> > sessionTable_;

  /// プロセス毎の通信整理
  mutable std::vector<TlMatrixObject*> trafficControl_;

  mutable bool isWaitingRequestHandShake_;
  mutable int sessionId_;

  // リクエスト受信-行列要素送信用
  enum GpmServerState {
    GPM_SERVER_WAIT = 0,
    GPM_SERVER_RECV_NUM_OF_COMPONENTS = 1,
    GPM_SERVER_WAIT_NUM_OF_COMPONENTS = 2,
    GPM_SERVER_RECV_COMPONENTS = 4,
    GPM_SERVER_WAIT_COMPONENTS = 8,
    GPM_SERVER_SEND_ELEMENT_VALUES = 16,
    GPM_SERVER_WAIT_ELEMENT_VALUES = 32
  };

  enum {
    GPM_CLIENT_COUNT_COMPONENTS = 1,
    GPM_CLIENT_SEND_SESSION_ID = 2,
    GPM_CLIENT_WAIT_SESSION_ID = 4,
    GPM_CLIENT_SEND_NUM_OF_COMPONENTS = 8,
    GPM_CLIENT_WAIT_NUM_OF_COMPONENTS = 16,
    GPM_CLIENT_SEND_COMPONENTS = 32,
    GPM_CLIENT_WAIT_COMPONENTS = 64,
    GPM_CLIENT_RECV_ELEMENT_VALUES = 128,
    GPM_CLIENT_WAIT_ELEMENT_VALUES = 256
  };

  // ---------------------------------------------------------------------------
  // friends
  // ---------------------------------------------------------------------------
  friend TlDenseGeneralMatrix_ImplScalapack operator*(
      const TlDenseGeneralMatrix_ImplScalapack& X,
      const TlDenseGeneralMatrix_ImplScalapack& Y);

  friend TlDenseVector_ImplScalapack operator*(
      const TlDenseGeneralMatrix_ImplScalapack& A,
      const TlDenseVector_ImplScalapack& X);
  friend TlDenseVector_ImplScalapack operator*(
      const TlDenseVector_ImplScalapack& X,
      const TlDenseGeneralMatrix_ImplScalapack& A);

  friend bool diagonalByScaLapack_QR(
      const TlDenseSymmetricMatrix_ImplScalapack& inMatrix,
      TlDenseVector_ImplLapack* outEigVal,
      TlDenseGeneralMatrix_ImplScalapack* outEigVec);
  friend bool diagonalByScaLapack_DC(
      const TlDenseSymmetricMatrix_ImplScalapack& inMatrix,
      TlDenseVector_ImplLapack* outEigVal,
      TlDenseGeneralMatrix_ImplScalapack* outEigVec);

  // ---------------------------------------------------------------------------
};

#endif  // TL_DENSE_GENERAL_MATRIX_IMPL_SCALAPACK_H
