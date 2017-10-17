// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef TLDISTRIBUTEMATRIX_H
#define TLDISTRIBUTEMATRIX_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <vector>
#include <list>
#include <bitset>
#include "TlMatrixObject.h"
#include "TlVector.h"
#include "TlSparseMatrix.h"
#include "TlRowVectorMatrix.h"
#include "TlUtils.h"
#include "TlLogging.h"
#include "TlCommunicate.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif // HAVE_HDF5

#define MAX_SESSION_ID (100)

class TlDistributeSymmetricMatrix;
class TlDistributeVector;

class TlScalapackContext {
public:
    static void getData(int& rContext, int& rProc, int& rRank,
                        int& rProcGridRow, int& rProcGridCol);
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
};


/// 分散クラス
class TlDistributeMatrix : public TlMatrixObject {
protected:
    typedef std::valarray<double> DataType;

public:
    struct LocalMatrixHeader {
    public:
        LocalMatrixHeader(int inType =0,
                          index_type inGlobalRow =0,
                          index_type inGlobalCol =0,
                          index_type inMyRows =0,
                          index_type inMyCols =0,
                          int inRank =0, int inNumOfProcs =0,
                          int inProcGridRow =0, int inProcGridCol =0,
                          int inMyProcRow =0, int inMyProcCol =0,
                          int inBlockSize =0)
            :type(inType), globalRow(inGlobalRow), globalCol(inGlobalCol), myRows(inMyRows), myCols(inMyCols),
             rank(inRank), numOfProcs(inNumOfProcs), procGridRow(inProcGridRow), procGridCol(inProcGridCol),
             myProcRow(inMyProcRow), myProcCol(inMyProcCol), blockSize(inBlockSize) {
        }

        void load(std::ifstream* pIs);
        void save(std::ofstream* pOs);
        
    public:
        int type;
        index_type globalRow;
        index_type globalCol;
        index_type myRows;
        index_type myCols;
        int rank;
        int numOfProcs;
        int procGridRow;
        int procGridCol;
        int myProcRow;
        int myProcCol;
        int blockSize;
    };
    
public:
    /// 行数・列数を指定してオブジェクトを作成する
    explicit TlDistributeMatrix(int row =1, int col =1);

    /// コピーコンストラクタ
    TlDistributeMatrix(const TlDistributeMatrix& rhs);

    /// 対称行列からの変換コンストラクタ
    TlDistributeMatrix(const TlDistributeSymmetricMatrix& rhs);

    /// ベクトルからの変換コンストラクタ
    TlDistributeMatrix(const TlDistributeVector& rhs,
                       const int row, const int col);

    TlDistributeMatrix(const TlRowVectorMatrix& rhs);

    virtual ~TlDistributeMatrix();

    static void setSystemBlockSize(int blockSize);
    static void setUsingPartialIO(bool isUsePIO);

public:
    /// 行数を返す
    ///
    /// @return 行数
    index_type getNumOfRows() const;


    /// 列数を返す
    ///
    /// @return 列数
    index_type getNumOfCols() const;


    /// サイズを変更する
    void resize(index_type row, index_type col);

protected:
    virtual size_type getNumOfElements() const;
    
    /// ブロックサイズを返す
    int getBlockSize() const;


public:
    /// 行列要素をベクトルにして返す
    ///
    /// 要素全体を返すのではなく、このプロセスが保持しているデータを返すことに注意
    virtual TlDistributeVector getVector() const;

    /// 指定した行の要素から構成されるベクトルを返す
    ///
    /// @param[in] nRow 指定する行
    virtual TlVector getRowVector(index_type row) const;

    /// 指定した列の要素から構成されるベクトルを返す
    ///
    /// @param[in] nCol 指定する列
    virtual TlVector getColumnVector(index_type col) const;

    /// 対角要素の和を返す
    ///
    /// @return 対角要素の和
    double trace() const;

    /// 要素の絶対値の最大値を返す
    ///
    /// @param[out] outRow 該当要素の行
    /// @param[out] outCol 該当要素の列
    /// @return 要素の絶対値の最大値
    virtual double getMaxAbsoluteElement(index_type* pOutRow =NULL, index_type* pOutCol =NULL) const;

protected:
    /// ローカル行列における要素の絶対値の最大値を返す
    ///
    /// @param[out] outRow 該当要素の行
    /// @param[out] outCol 該当要素の列
    /// @return 要素の絶対値の最大値
    virtual double getLocalMaxAbsoluteElement(index_type* pOutRow =NULL, index_type* pOutCol =NULL) const;

public:
    /// calc RMS
    virtual double getRMS() const;

    /// 全行列を各プロセスに均等分割された疎行列を返す
    // TlSparseMatrix getPartialMatrix(double threshold = 1.0E-16) const;

    bool getSparseMatrixX(TlSparseMatrix* pMatrix, bool isFinalize = false) const;

    /// 指定された要素を持つ疎行列を返す
    // void getPartialMatrix(TlSparseMatrix& ioMatrix) const;

    /// 各ノードが与えた疎行列を大域行列に加算する。
    /// 
    /// 全ノードがこの関数を呼び出す必要がある。
    void mergeSparseMatrix(const TlSparseMatrix& M);


    /// dot積を返す
    const TlDistributeMatrix& dot(const TlDistributeMatrix& X);

    /// 全ての要素の合計を返す
    virtual double sum() const;

    /// 指定された要素(グローバル)がどのプロセスが所有しているかを返す
    virtual int getProcIdForIndex(index_type globalRow, index_type globalCol) const;
    
    std::vector<index_type> getRowIndexTable() const;
    std::vector<index_type> getColIndexTable() const;

    TlMatrix getLocalMatrix() const;
    
    // need to call by all processes.
    TlDistributeMatrix& operator=(const TlDistributeMatrix& rhs);
    TlDistributeMatrix& operator=(const TlDistributeSymmetricMatrix& rhs);

    TlDistributeMatrix& operator+=(const TlDistributeMatrix& rhs);
    TlDistributeMatrix& operator-=(const TlDistributeMatrix& rhs);
    TlDistributeMatrix& operator*=(const TlDistributeMatrix& rhs);
    TlDistributeMatrix& operator*=(const TlDistributeSymmetricMatrix& rhs);
    TlDistributeMatrix& operator*=(double dCoef);
    TlDistributeMatrix& operator/=(double dCoef);
    
public:
    /// 指定された入力ストリームが読み込み可能かどうかを返す
    ///
    /// @param[in,out] ifs 入力ファイルストリーム
    /// @retval true 読み取り可能
    /// @retval false 読み取り不可能
    static bool isLoadable(std::ifstream& ifs);
    static bool isLoadable(const std::string& rFilePath);

    virtual bool load(const std::string& sFilePath);
    virtual bool load(std::ifstream& ifs);
    virtual bool save(const std::string& sFilePath) const;

#ifdef HAVE_HDF5
public:
    virtual bool saveHdf5(const std::string& filepath, const std::string& h5path) const;
    virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);

protected:
    bool saveHdf5(const std::string& filepath, const std::string& h5path, const int saveMatType) const;
    virtual size_type getArrayIndex(const index_type row, const index_type col) const;
    void saveElements(TlHdf5Utils* pH5, const std::string& path,
                      const std::vector<TlMatrixObject::MatrixElement>& elements) const;
#endif // HAVE_HDF5

public:
    virtual bool saveText(const std::string& sFilePath) const;
    virtual bool saveText(std::ofstream& ofs) const;

protected:
    virtual std::vector<TlMatrixObject::MatrixElement> getMatrixElementsInLocal() const;
    virtual void saveElements(TlFileMatrix* pFileMatrix, const std::vector<TlMatrixObject::MatrixElement>& elements) const;
    
public:
    /// 使用しているメモリのサイズを返す
    virtual std::size_t getMemSize() const;

public:
    // need to call by all processes.
    virtual double get(index_type row, index_type col) const;

    // need to call by all processes.
    virtual void set(index_type row, index_type col, double dValue);

    // need to call by all processes.
    virtual double operator()(index_type row, index_type col) const;
    virtual double& operator()(index_type row, index_type col);

    virtual void add(index_type row, index_type col, double value);

    virtual void addByList(const index_type* pIndexPairs,
                           const double* pValues,
                           const std::size_t size);

public:
    /// ローカルから該当する要素があれば値を返す
    ///
    /// 範囲外であれば0.0を返す
    /// すべてのプロセスが同時に呼ぶ必要はない
    double getLocal(index_type row, index_type col) const;

public:
    virtual bool inverse();

    virtual TlVector getDiagonalElements() const;
    
    /// 転置行列にする
    ///
    /// @return 転置行列となったこのオブジェクト
    virtual const TlDistributeMatrix& transpose();

public:
    // need to call by all processes.
    template <typename T>
    void print(T& out) const;

    // need to call by all processes.
    //virtual void update() const;

protected:
    virtual void initialize();

    
    /// 必要な行列要素の数を返す。
    /// initialize()などから呼び出される。
    virtual index_type getNumOfMyElements() const;


    /// ローカル行列のインデックスを返す。
    /// 保持していない場合は-1を返す。
    virtual size_type getIndex(index_type globalRow, index_type globalCol) const;


    /// 各ノードに配列保持された行列要素をマージする
    ///
    /// @param[in] indexArrays プロセス番号ごとに行列インデックス(行、列)を順に格納する。(i1, j1, i2, j2, ...)
    /// @param[in] values 行列要素の値を格納する。
    void mergeMatrix_common(const std::vector<std::vector<index_type> >& indexArrays,
                            const std::vector<std::vector<double> >& values);
    
protected:
    virtual bool loadLocal(const std::string& filePath);
    virtual bool saveLocal(const std::string& filePath) const;

protected:
    TlLogging& log_;

    int m_nContext;
    int m_pDESC[9];

    index_type m_nRows; // 大域行列の行数
    index_type m_nCols; // 大域行列の列数

    index_type m_nMyRows; // ローカル行列の行数
    index_type m_nMyCols; // ローカル行列の列数

    int m_nRank; // プロセスのランク
    int m_nProc; // 総プロセス数
    int m_nProcGridRow; // プロセスグリッドの行数
    int m_nProcGridCol; // プロセスグリッドの列数
    int m_nMyProcRow; // プロセスグリッドにおける自分の行数
    int m_nMyProcCol; // プロセスグリッドにおける自分の列数
    int m_nBlockSize; // ブロックサイズ

    double* pData_;

    std::vector<index_type> m_RowIndexTable; // m_RowIndexTable[local_index] = global_index
    std::vector<index_type> m_ColIndexTable;
    std::vector<int> processMatrix_;
    
    double m_dTempVar;

protected:
    static int systemBlockSize_;
    static const std::size_t FILE_BUFFER_SIZE;
    static const size_type MAX_LOOP;
    static bool isUsingPartialIO; // 分割保存形式を使う(true, defaultはfalse)
    
protected:
    /// MPI通信タグ
    enum {
        // for copy constructor
        TAG_CONSTRUCTOR_SIZE   = 10001,
        TAG_CONSTRUCTOR_INDEX  = 10002,
        TAG_CONSTRUCTOR_VALUE  = 10003,
        // for LOAD
        TAG_LOAD_SIZE          = 10011,
        TAG_LOAD_ROWCOLS       = 10012,
        TAG_LOAD_VALUES        = 10013,
        TAG_LOAD_END           = 10014,
        // for save
        TAG_SAVE_HANDSHAKE     = 10021,
        TAG_SAVE_HANDSHAKE_OK  = 10022,
        TAG_SAVE_DATA_ROWS     = 10023,
        TAG_SAVE_DATA_COLS     = 10024,
        TAG_SAVE_DATA_ROWINDEX = 10025,
        TAG_SAVE_DATA_COLINDEX = 10026,
        TAG_SAVE_DATA          = 10027,
        // for getPartialMatrix / getSparseMatrix
        TAG_GPM_SESSION_ID        = 10100,
        TAG_GPM_NUM_OF_COMPONENTS = 10200,
        TAG_GPM_COMPONENTS        = 10300,
        TAG_GPM_ELEMENT_VALUES    = 10400,
        // for MergeMatrix
        TAG_MERGE_MATRIX_SESSION_ID      = 10500,
        TAG_MERGE_MATRIX_NUM_OF_CONTENTS = 10600,
        TAG_MERGE_MATRIX_INDECES         = 10700,
        TAG_MERGE_MATRIX_VALUES          = 10800
    };
    
protected:
    bool isDebugOut_GPM_;
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
    struct GPM_ClientTask {
        std::vector<unsigned int> state;
        std::vector<int> sessionIds;
        
        /// プロセス毎へリクエストする要素番号情報
        /// プロセス番号を添字として、行番号、列番号の順に1次元配列に格納する
        std::vector<std::vector<index_type> > components; 

        /// プロセス毎へリクエストする要素番号情報の総数
        std::vector<std::size_t> numOfComponents;

        /// 受信用行列データ
        std::vector<std::vector<double> > elementValues;

        bool isFinished;
    };
    typedef std::map<TlMatrixObject*, GPM_ClientTask> GpmClientTasks;
    mutable GpmClientTasks gpmClientTasks_;
    
    mutable std::vector<std::bitset<MAX_SESSION_ID> > sessionTable_;

    
    /// プロセス毎の通信整理
    mutable std::vector<TlMatrixObject*> trafficControl_;

    bool getPartialMatrix_ClientTasks(TlMatrixObject* pMatrix) const;
    void getSparseMatrixX_registerTask(TlSparseMatrix* pMatrix) const;
    int getPartialMatrix_getSessionId(const int proc) const;
    void getPartialMatrix_ServerTasks(bool isFinalize) const;
    

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
    struct GPM_ServerTask {
    public:
        GPM_ServerTask()
            : state(0), sessionId(0), requestProc(0), numOfComponents(0) {
        }
    public:
        unsigned int state;
        int sessionId;
        int requestProc;
        std::size_t numOfComponents;
        std::vector<index_type> components;
        std::vector<double> elementValues;
    };
    typedef std::list<GPM_ServerTask> GpmServerTasksType;
    mutable GpmServerTasksType gpmServerTasks_;

    mutable bool isWaitingRequestHandShake_;
    mutable int sessionId_;

    // ===========
    enum {
        TAG_TLDISTRIBUTE_MATRIX__MERGE__NUM_OF_INDECES = 2001,
        TAG_TLDISTRIBUTE_MATRIX__MERGE__INDECES = 2002,
        TAG_TLDISTRIBUTE_MATRIX__MERGE__VALUES = 2003
    };

    // merge matrix ============================================================
public:
    virtual void mergeSparseMatrixAsync(const TlSparseMatrix* pMatrix,
                                        bool isFinalize = false);

protected:
    void mergeMatrixAsync_send(const std::vector<std::vector<index_type> >& indexArrays,
                               const std::vector<std::vector<double> >& values);
    void mergeMatrixAsync_recv(bool isFinalize);

    bool mm_isWaitingSessionId_;
    int mm_waitingSessionId_;
    
    struct MergeMatrixRecvTask {
    public:
        MergeMatrixRecvTask(int src = -1, int id = -1)
            : srcProc(src), sessionId(id), state(0), numOfContents(0) {
        }
        
    public:
        int srcProc;
        int sessionId;
        unsigned int state;
        std::size_t numOfContents;
        std::vector<index_type> indeces;
        std::vector<double> values;
    };

    typedef std::list<MergeMatrixRecvTask> MergeMatrixRecvTasks;
    MergeMatrixRecvTasks mergeMatrixRecvTasks_;

    enum {
        MM_RECV_NUM_OF_CONTENTS = 1,
        MM_WAIT_NUM_OF_CONTENTS = 2,
        MM_RECV_INDECES = 4,
        MM_WAIT_INDECES = 8,
        MM_RECV_VALUES = 16,
        MM_WAIT_VALUES = 32,
        MM_FINISHED = 64
    };
    
    // friend ------------------------------------------------------------------
    friend TlDistributeMatrix operator*(const TlDistributeMatrix& X, const TlDistributeMatrix& Y);
    friend TlDistributeVector operator*(const TlDistributeMatrix& A, const TlDistributeVector& X);
    //friend TlDistributeVector operator*(const TlDistributeVector& X, const TlDistributeMatrix& A);
    friend TlDistributeMatrix operator*(double X, const TlDistributeMatrix& Y);
    friend TlDistributeMatrix operator*(const TlDistributeMatrix& X, double Y);

    friend TlDistributeMatrix operator+(const TlDistributeMatrix& X, const TlDistributeMatrix& Y);
    friend TlDistributeMatrix operator-(const TlDistributeMatrix& X, const TlDistributeMatrix& Y);

    friend TlDistributeMatrix multiplicationByScaLapack(const TlDistributeSymmetricMatrix& X,
                                                        const TlDistributeMatrix& Y);
    friend TlDistributeMatrix multiplicationByScaLapack(const TlDistributeMatrix& X,
                                                        const TlDistributeSymmetricMatrix& Y);

    friend bool diagonalByScaLapack_QR(const TlDistributeSymmetricMatrix& inMatrix,
                                       TlVector* outEigVal, TlDistributeMatrix* outEigVec);
    friend bool diagonalByScaLapack_DC(const TlDistributeSymmetricMatrix& inMatrix,
                                       TlVector* outEigVal, TlDistributeMatrix* outEigVec);

    friend bool inverseByScaLapack(TlDistributeMatrix& X);
    friend bool inverseByScaLapack(TlDistributeSymmetricMatrix& X);
};

// =============================================================================
// inline double TlDistributeMatrix::getDirect(const index_type row, const index_type col) const
// {
//     const std::size_t index = row + col * this->m_nMyRows;
//     assert(index < this->getNumOfMyElements());

//     return this->pData_[index];
// }


template <typename T>
void TlDistributeMatrix::print(T& out) const
{
    //this->update();

    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int nNumOfRows = this->getNumOfRows();
    const int nNumOfCols = this->getNumOfCols();

    for (int ord = 0; ord < nNumOfCols; ord += 10) {
        if (rComm.isMaster() == true) {
            out << "       ";
            for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
                out << TlUtils::format("   %5d th", j+1);
            }
            out << "\n ----";

            for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
                out << "-----------";
            }
            out << "----\n";
        }

        for (int i = 0; i < nNumOfRows; ++i) {
            if (rComm.isMaster() == true) {
                out << TlUtils::format(" %5d  ", i+1);
            }

            for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
                std::string s = TlUtils::format(" %10.6lf", (*this)(i, j));

                if (rComm.isMaster() == true) {
                    out << s;
                }
            }

            if (rComm.isMaster() == true) {
                out << "\n";
            }
        }

        if (rComm.isMaster() == true) {
            out << "\n\n";
        }
    }

    if (rComm.isMaster() == true) {
        out.flush();
    }
}

#endif // TLDISTRIBUTEMATRIX_H
