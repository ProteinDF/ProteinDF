#ifndef TL_DENSE_SCALAPACK_OBJECT_H
#define TL_DENSE_SCALAPACK_OBJECT_H

#include "tl_matrix_object.h"

class TlDenseScalapackObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlDenseScalapackObject(TlMatrixObject::index_type row = 1,
                                    TlMatrixObject::index_type col = 1);
    TlDenseScalapackObject(const TlDenseScalapackObject& rhs);
    virtual ~TlDenseScalapackObject();

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
    virtual void mul(TlMatrixObject::index_type row,
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
    TlDenseScalapackObject& operator=(const TlDenseScalapackObject& rhs);

    TlDenseScalapackObject& operator+=(const TlDenseScalapackObject& rhs);
    TlDenseScalapackObject& operator-=(const TlDenseScalapackObject& rhs);
    TlDenseScalapackObject& operator*=(double coef);
    TlDenseScalapackObject& operator/=(double coef);

    // ---------------------------------------------------------------------------
    // I/O
    // ---------------------------------------------------------------------------
   public:
    // ---------------------------------------------------------------------------
    // protected
    // ---------------------------------------------------------------------------
   protected:
    virtual void initialize();

    TlMatrixObject::index_type getNumOfMyRows() const;
    TlMatrixObject::index_type getNumOfMyCols() const;
    virtual TlMatrixObject::size_type getNumOfMyElements() const;

    // 担当rank, ローカル行列のインデックスを返す。
    void getGlobalRowCol2LocalRowCol(
        const TlMatrixObject::index_type globalRow,
        const TlMatrixObject::index_type globalCol, int* pRank,
        TlMatrixObject::index_type* pLocalRow = NULL,
        TlMatrixObject::index_type* pLocalCol = NULL) const;

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

    virtual std::vector<TlMatrixObject::MatrixElement>
    getMatrixElementsInLocal() const;

    // ---------------------------------------------------------------------------
    // variables
    // ---------------------------------------------------------------------------
   protected:
    static const std::size_t FILE_BUFFER_SIZE;

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
    // std::vector<TlMatrixObject::index_type> m_RowIndexTable;
    // std::vector<TlMatrixObject::index_type> m_ColIndexTable;

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

    // ---------------------------------------------------------------------------
};

#endif  // TL_DENSE_SCALAPACK_OBJECT_H
