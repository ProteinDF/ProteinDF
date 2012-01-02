#ifndef TLCOMMUNICATE_H
#define TLCOMMUNICATE_H

#include <mpi.h>
#include <stdint.h>
#include <vector>
#include <valarray>
#include <map>
#include <string>

#include "TlPartialMatrix.h"
#include "TlPartialSymmetricMatrix.h"
#include "TlTime.h"
#include "TlLogging.h"

//#define ANY_TAG_OK

class TlVector;
class TlMatrix;
class TlSymmetricMatrix;
class TlSparseMatrix;
class TlSparseSymmetricMatrix;
class TlFileMatrix;
class TlFileSymmetricMatrix;
class TlMmapMatrix;
class TlMmapSymmetricMatrix;
class TlParameter;
class TlSerializeData;

// シングルトンパターンで実装
class TlCommunicate {
public:
    static TlCommunicate& getInstance(int argc, char* argv[]);
    static TlCommunicate& getInstance();

private:
    // for Singleton
    TlCommunicate();
    TlCommunicate(const TlCommunicate& rhs);
    ~TlCommunicate();
    const TlCommunicate& operator=(const TlCommunicate& rhs);

public:
    /// 作業用メモリサイズを設定する
    ///
    /// @param workMemSize [in] 作業用メモリサイズ(byte単位)
    void setWorkMemSize(std::size_t workMemSize);

    /// 作業用メモリサイズを返す
    ///
    /// @return 作業用メモリサイズ(byte単位)
    std::size_t getWorkMemSize() const;

    int getNumOfProcs() const;

    // absolute
    int getNumOfProc() const {
        return this->getNumOfProcs();
    }
    
    int getRank() const;

    std::string getReport() const;
    
    /**
     *  rank が 0 だったらtrueを返す
     */
    bool isMaster() const;
    bool isSlave() const;

    /**
     *  同期を取る
     */
    int barrier(bool isDebugOut = false) const;

    // データ通信 ========================================================
    // 集団通信
    int reduce_SUM(unsigned int* pData, std::size_t size, int root = 0);
    int reduce_SUM(unsigned long* pData, std::size_t size, int root = 0);
    
    int allReduce_SUM(int& rData);
    int allReduce_SUM(unsigned int& rData);
    int allReduce_SUM(long& rData);
    int allReduce_SUM(unsigned long& rData);
    int allReduce_SUM(double& rData);
    int allReduce_SUM(std::vector<int>& rData);
    int allReduce_SUM(std::vector<unsigned int>& rData);
    int allReduce_SUM(std::vector<long>& rData);
    int allReduce_SUM(std::vector<unsigned long>& rData);
    int allReduce_SUM(std::vector<double>& rData);
    int allReduce_SUM(std::valarray<double>& rData);
    int allReduce_SUM(int* pData, std::size_t length);
    int allReduce_SUM(double* pData, std::size_t length);
    int allReduce_SUM(TlVector& rVector);
    int allReduce_SUM(TlMatrix& rMatrix);
    int allReduce_SUM(TlSymmetricMatrix& rMatrix);
    int allReduce_SUM(TlMmapMatrix& rMatrix);
    int allReduce_SUM(TlMmapSymmetricMatrix& rMatrix);
    int allReduce_AND(bool& rData);
    int allReduce_MAX(int& rData);
    int allReduce_MIN(int& rData);

    int gatherToMaster(TlSparseMatrix& rMatrix);

    /** master に保存されたデータをslaveにコピーします
     *
     *  @param[in/out] rData 送信(master)・受信(slave)オブジェクトの参照
     */
    int broadcast(bool& rData);
    int broadcast(int& data, int nRoot = 0);
    int broadcast(unsigned int& data, int nRoot = 0);
    int broadcast(long& data, int root = 0);
    int broadcast(unsigned long& data, int root = 0);
    int broadcast(double& data, int root = 0);
    int broadcast(std::string& rData);
    int broadcast(std::vector<int>& data, int root = 0);
    int broadcast(std::vector<long>& data, int root = 0);
    int broadcast(std::vector<unsigned long>& rData, int nRoot = 0);
    int broadcast(std::vector<double>& rData, int nRoot = 0);
    int broadcast(std::valarray<double>& rData, int nRoot = 0);
    int broadcast(std::vector<std::string>& rData);
    int broadcast(TlVector& rData, int root = 0);
    int broadcast(TlMatrix& rData, int root = 0);
    int broadcast(TlSymmetricMatrix& rData, int root = 0);
    int broadcast(TlSparseMatrix& rData, int nRoot = 0);
    int broadcast(TlMmapMatrix& data, int root = 0);
    int broadcast(TlMmapSymmetricMatrix& data, int root = 0);
    int broadcast(TlParameter& rParam);
    int broadcast(TlSerializeData& data);

    int broadcast(double* p, const std::size_t size, const int root);

    // 1:1通信
    int sendData(bool data, int destination = 0, int tag = 0);
    int sendData(int nData, int nDestination = 0, int nTag = 0);
    int sendData(unsigned int nData, int nDestination = 0, int nTag = 0);
    int sendData(long nData, int nDestination = 0, int nTag = 0);
    int sendData(unsigned long nData, int nDestination = 0, int nTag = 0);
    int sendData(double nData, int nDestination = 0, int nTag = 0);
    int sendData(const std::vector<int>& data, int nDestination = 0, int nTag = 0);
    int sendData(const std::vector<unsigned int>& data, int nDestination = 0, int nTag = 0);
    int sendData(const std::vector<long>& data, int nDestination = 0, int nTag = 0);
    int sendData(const std::vector<unsigned long>& data, int nDestination = 0, int nTag = 0);
    int sendData(const std::vector<double>& data, int nDestination = 0, int nTag = 0);
    int sendData(const std::valarray<double>& data, int nDestination = 0, int nTag = 0);
    int sendData(const std::string& data, int nDestination = 0, int nTag = 0);
    int sendData(const TlVector& data, int destination = 0, int tag = 0);
    int sendData(const TlMatrix& data, int destination = 0, int tag = 0);
    int sendData(const TlSymmetricMatrix& rData, int nDestination = 0, int nTag = 0);
    int sendData(const TlSparseSymmetricMatrix& rData, int nDestination = 0, int nTag = 0);
    int sendData(const TlPartialSymmetricMatrix& rData, int nDestination = 0, int nTag = 0);

    int receiveData(bool& data, int src, int tag = 0);
    int receiveData(int& rData, int nSrc, int nTag = 0);
    int receiveData(unsigned int& rData, int nSrc, int nTag = 0);
    int receiveData(long& rData, int nSrc, int nTag = 0);
    int receiveData(unsigned long& rData, int nSrc, int nTag = 0);
    int receiveData(double& rData, int nSrc, int nTag = 0);
    int receiveData(std::vector<int>& rData, int nSrc, int nTag = 0);
    int receiveData(std::vector<unsigned int>& rData, int nSrc, int nTag = 0);
    int receiveData(std::vector<long>& rData, int nSrc, int nTag = 0);
    int receiveData(std::vector<unsigned long>& rData, int nSrc, int nTag = 0);
    int receiveData(std::vector<double>& rData, int nSrc, int nTag = 0);
    int receiveData(std::valarray<double>& rData, int nSrc, int nTag = 0);
    int receiveData(std::string& pData, int nSrc, int nTag = 0);
    int receiveData(TlVector& data, int src, int tag = 0);
    int receiveData(TlMatrix& data, int src, int tag = 0);
    int receiveData(TlSymmetricMatrix& rData, int nSrc, int nTag = 0);
    int receiveData(TlSparseSymmetricMatrix& rData, int nSrc, int nTag = 0);
    int receiveData(TlPartialSymmetricMatrix& rData, int nSrc, int nTag = 0);

    int receiveDataFromAnySource(int& data, int* pSrc, const int tag);

    int receiveDataFromAnySource(int& data, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(unsigned int& data, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(long& data, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(unsigned long& data, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(std::vector<int>& rData, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(std::vector<unsigned int>& rData, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(std::vector<long>& rData, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(std::vector<unsigned long>& rData, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(std::vector<double>& rData, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(TlSparseSymmetricMatrix& rData, int* pSrc, int* pTag = NULL);
    int receiveDataFromAnySource(TlPartialSymmetricMatrix& rData, int* pSrc, int* pTag = NULL);

    // タグ付き
    int receiveDataFromAnySource(TlMatrix& rData, int* pSrc, int tag = 0);
    
    // 1:1通信 非ブロッキング
    int iSendData(const int& data, int destination = 0, int tag = 0);
    int iSendData(const unsigned int& data, int destination = 0, int tag = 0);
    int iSendData(const long& data, int destination = 0, int tag = 0);
    int iSendData(const unsigned long& data, int destination = 0, int tag = 0);
    int iSendData(const std::vector<int>& data, int destination = 0, int tag = 0);
    int iSendData(const std::vector<double>& data, int destination = 0, int tag = 0);

    int iReceiveData(int& data, int src, int tag = 0);
    int iReceiveData(unsigned int& data, int src, int tag = 0);
    int iReceiveData(long& data, int src, int tag = 0);
    int iReceiveData(unsigned long& data, int src, int tag = 0);
    int iReceiveData(std::vector<int>& data, int src, int tag = 0);
    int iReceiveData(std::vector<double>& data, int src, int tag = 0);

    int iReceiveDataFromAnySource(int& data, int tag = 0);
    int iReceiveDataFromAnySource(unsigned int& data, int tag = 0);
    int iReceiveDataFromAnySource(long& data, int tag = 0);
    int iReceiveDataFromAnySource(unsigned long& data, int tag = 0);
    int iReceiveDataFromAnySource(std::vector<int>& data, int tag = 0);
    int iReceiveDataFromAnySource(std::vector<double>& data, int tag = 0);

    int sendDataX(const int* pData, const std::size_t size,
                  const int dest, const int tag = 0);
    int sendDataX(const unsigned int* pData, const std::size_t size,
                  const int dest, const int tag = 0);
    int sendDataX(const double* pData, const std::size_t size,
                  const int dest, const int tag = 0);
    
    int receiveDataX(int* pData, const std::size_t size,
                     const int src, const int tag = 0);
    int receiveDataX(unsigned int* pData, const std::size_t size,
                     const int src, const int tag = 0);
    int receiveDataX(double* pData, const std::size_t size,
                     const int src, const int tag = 0);
    
    int receiveDataFromAnySourceX(int* pData, std::size_t size, int* pSrc, int tag = 0);

    int iSendDataX(const int* pData, const std::size_t size,
                   const int dest, const int tag = 0);
    int iSendDataX(const double* pData, const std::size_t size,
                   const int dest, const int tag = 0);
    
    int iReceiveDataX(int* pData, const std::size_t size,
                      const int src, const int tag = 0);
    int iReceiveDataX(double* pData, const std::size_t size,
                      const int src, const int tag = 0);

    int iReceiveDataFromAnySourceX(int* pData, std::size_t size, int tag = 0);

    /// 非同期通信をキャンセルする
    /// 内部でwait()関数を呼んでいるので、再度wait()/test()を呼ぶ必要はない。
    template<typename T>
    bool cancel(const T& data) {
        return this->cancel((void*)&data);
    }

    template<typename T>
    bool cancel(T* pData) {
        return this->cancel((void*)pData);
    }
    
    /// 転送が完了したかどうかをテストする
    /// iReceiveDataFromAnySourceを使う場合は、
    /// wait()関数ではなくtest()関数に送信元rankが指定されるので注意。
    template<typename T>
    bool test(const T& data, int* pSrc = NULL) {
        return this->test((void*)&data, pSrc);
    }

    template<typename T>
    bool test(T* pData, int* pSrc = NULL) {
        return this->test((void*)pData, pSrc);
    }
    
    template<typename T>
    bool test(const std::vector<T>& data, int* pSrc = NULL) {
        return this->test((void*)&(data[0]), pSrc);
    }

    /// wait
    template<typename T>
    int wait(const T& data, int* pSrc = NULL) {
        return this->wait((void*)&data, pSrc);
    }

    template<typename T>
    int wait(T* pData, int* pSrc = NULL) {
        return this->wait((void*)pData, pSrc);
    }

    template<typename T>
    int wait(const std::vector<T>& data, int* pSrc = NULL) {
        return this->wait((void*)&(data[0]), pSrc);
    }

    // 集計
    int allReduce_SUM(const TlFileMatrix& fromLocalMatrix,
                      const std::string& toMatrixFilePath);
    int allReduce_SUM(const TlFileSymmetricMatrix& fromLocalMatrix,
                      const std::string& toMatrixFilePath);

    // 終了処理
    // プログラム終了時は必ず呼ぶこと
    int finalize();

    /// プロセスのアボート
    int abort(const int errorCode);
    
    ///
//     int getNumOfWaiting() const {
//         return this->nonBlockingCommParamTable_.size();
//     }

    /// 非同期通信が残っていないかチェックする
    bool checkNonBlockingCommunications() const;

    
public:
    double getTime();

public:
    // for debug
    unsigned long getBarrierCount() {
        return this->counter_barrier_;
    };

protected:
    template<typename T>
    int reduce(T* pData, const MPI_Datatype mpiType,
               const std::size_t start, const std::size_t end,
               const MPI_Op mpiOp, const int root);

    template<typename T>
    int allReduce(T* pData, const MPI_Datatype mpiType,
                  const std::size_t start, const std::size_t end,
                  const MPI_Op mpiOp);
    
    template<typename T>
    int allReduce_SUM(std::vector<T>& data, const MPI_Datatype mpiType,
                      const std::size_t start, const std::size_t end);

    int allReduce_SUM(std::valarray<double>& data,
                      const std::size_t start, const std::size_t end);

//     int allReduce_SUM(double* pData,
//                       const std::size_t start, const std::size_t end);

    template<typename T>
    int sendData(const std::vector<T>& data, const MPI_Datatype mpiType,
                 const std::size_t start, const std::size_t end,
                 int destination, int tag);

    int sendData(const std::valarray<double>& data,
                 const std::size_t start, const std::size_t end,
                 int destination, int tag);

    template<typename T>
    int sendDataX(const T* pData, const MPI_Datatype mpiType,
                  const std::size_t start, const std::size_t end,
                  const int dest, const int tag);


    template<typename T>
    int receiveData(std::vector<T>& data, const MPI_Datatype mpiType,
                    const std::size_t start, const std::size_t end,
                    const int src, const int tag);

    int receiveData(std::valarray<double>& data,
                    const std::size_t start, const std::size_t end,
                    const int src, const int tag);

    template<typename T>
    int receiveDataX(T* pData, const MPI_Datatype mpiType,
                     const std::size_t start, const std::size_t end,
                     const int src, const int tag);

    template<typename T>
    int receiveDataFromAnySource(T& data, const MPI_Datatype mpiType, int* pSrc, const int tag);

    template<typename T>
    int receiveDataFromAnySource(T& data, const MPI_Datatype mpiType, int* pSrc, int* pTag);

    template<typename T>
    int receiveDataFromAnySource(std::vector<T>& data, const MPI_Datatype mpiType,
                                 int* pSrc, int*pTag);

    template<typename T>
    int receiveDataFromAnySourceX(T* pData, const MPI_Datatype mpiType,
                                  const std::size_t start, const std::size_t end,
                                  int* pSrc, const int tag);

    
    template<typename T>
    int iSendData(const T& data, const MPI_Datatype mpiType,
                  const int destination, const int tag);

    template<typename T>
    int iSendData(const std::vector<T>& data, const MPI_Datatype mpiType,
                  const std::size_t start, const std::size_t end,
                  const int destination, const int tag);

    template<typename T>
    int iReceiveData(T& data, const MPI_Datatype mpiType, const int src, const int tag);

    template<typename T>
    int iSendDataX(const T* pData, const MPI_Datatype mpiType,
                   const std::size_t start, const std::size_t end,
                   const int dest, const int tag);

    template<typename T>
    int iReceiveData(std::vector<T>& data, const MPI_Datatype mpiType,
                     const std::size_t start, const std::size_t end,
                     const int src, const int tag);

    template<typename T>
    int iReceiveDataX(T* pData, const MPI_Datatype mpiType,
                      const std::size_t start, const std::size_t end,
                      const int src, const int tag);

    template<typename T>
    int iReceiveDataFromAnySource(T& data,
                                  const MPI_Datatype mpiType,
                                  const int tag);
   
    template<typename T>
    int iReceiveDataFromAnySource(std::vector<T>& data,
                                  const MPI_Datatype mpiType,
                                  const int tag);

    template<typename T>
    int iReceiveDataFromAnySourceX(T* pData, const MPI_Datatype mpiType,
                                   const std::size_t start, const std::size_t end,
                                   const int tag);

    int broadcast(std::valarray<double>& data,
                  const std::size_t start, const std::size_t end,
                  int root);

private:
    int initialize(int argc, char* argv[]);

    bool test(void* pData, int* pSrc = NULL);
    int wait(void* pData, int* pSrc = NULL);
    bool cancel(void* pData);

private:
    template<typename T>
    int broadcast(T& data, const MPI_Datatype mpiType, const int root);

    template<typename T>
    int broadcast(T* p, const MPI_Datatype mpiType,
                  const std::size_t start, const std::size_t end, const int root);
    
    template<typename T>
    int broadcast(std::vector<T>& data, const MPI_Datatype mpiType,
                  const std::size_t start, const std::size_t end,
                  int root);

   
private:
    static TlCommunicate* m_pTlCommunicateInstance;

    /// 総プロセッサ数
    int m_nProc;

    /// 自分のランク数
    int m_nRank;

    std::string m_sProcName;

    /// 作業用メモリサイズ
    std::size_t workMemSize_;

    //static int m_nBarrierCount;
    mutable unsigned long counter_barrier_;
    mutable unsigned long counter_test_;
    mutable unsigned long counter_wait_;
    mutable unsigned long counter_allreduce_;
    
    /// 時間計測用変数(MPI_Barrier)
    mutable TlTime time_barrier_;
    mutable TlTime time_test_;
    mutable TlTime time_wait_;
    mutable TlTime time_allreduce_;

    TlLogging& log_;
    
private:
    // for non-blocking communication
    struct NonBlockingCommParam {
    public:
        NonBlockingCommParam()
            : requests(std::vector<uintptr_t>()), requestStates(),
              tag(0), property(0), source(-1) {
        }
        NonBlockingCommParam(const std::vector<uintptr_t>& reqs,
                             const int in_tag,
                             const unsigned int in_property)
            : requests(reqs), requestStates(reqs.size(), 0),
              tag(in_tag), property(in_property), source(-1) {
        }
        
    public:
        enum Property {
            RECV = 0,
            SEND = 1,
            COMPLETE = 2
        };

    public:
        std::vector<uintptr_t> requests;
        std::vector<unsigned int> requestStates;
        int tag;
        unsigned int property;
        int source;
    };
    typedef std::map<uintptr_t, NonBlockingCommParam> NonBlockingCommParamTableType;
    static NonBlockingCommParamTableType nonBlockingCommParamTable_;

private:
    /// 非同期通信テーブルの衝突をチェックする
    bool checkNonBlockingTableCollision(uintptr_t key,
                                        const NonBlockingCommParam& param,
                                        int line) const;
};


// memo
/*
int MPI_Init(int *argc, char **argv)
    MPIの実行環境の初期化を行う。
    argc      コマンド行の引数の数
    argv      コマンド行の引数

int MPI_Comm_size(MPI_Comm comm, int *size)
    通信を行うグループのサイズを決める。
    comm      通信を行うグループの指定
    size      グループ内のタスクの数を受け取る

int MPI_Comm_rank(MPI_Comm comm,int *rank)
    通信を行うグループのプロセスにタスク番号を与える。
    comm      通信を行うグループの指定。
    rank      commの中でのタスク番号を受け取る(0,1,2,..)

int MPI_Barrier(MPI_Comm comm)
    バリア同期を取る。
    comm      同期を取るグループの指定

int MP_Finalize(void)
    MPIの実行環境を終了する。

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest,
            int tag, MPI_Comm comm)
    送信関数
    buf       送信データバッファ
    count     送信データの個数
    datatype  データタイプ
    dest      メッセージの送信先を指定
    tag       メッセージタグ
    comm      通信を行うグループの指定

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source,
            int tag, MPI_Comm comm, MPI_Status *status)
    受信関数
    buf       受信データバッファ
    count     受信データの個数
    datatype  データタイプ
    source    メッセージの送信元のタスク番号を指定
              MPI_ANY_SOURCEで任意の送信元を指定
    tag       メッセージタグ
              MPI_ANY_TAGで任意のタグを指定
    comm      通信を行うグループの指定
    status    構造体MPI_Statusで受信状況を返す
              送信元、タグ、メッセージの大きさなど

Ｃプログラムでのデータタイプの引数とそれに対応するＣでのデータの型は、以下の通り。
    MPI Datatype       C Datatype
    MPI_CHAR           char
    MPI_SHORT          short
    MPI_INT            int
    MPI_LONG           long
    MPI_UNSIGNED_CHAR  unsigned char
    MPI_UNSIGNED_SHORT unsinged short
    MPI_UNSIGNED       unsinged int
    MPI_UNSIGNED_LONG  unsinged long
    MPI_FLOAT          float
    MPI_DOUBLE         double
    MPI_LONG_DOUBLE    long double
    MPI_BYTE           対応する型はありません
    MPI_PACKED         対応する型はありません

int MPI_Allreduce(void *sendbuf,void *recvbuf,int count,
                    MPI_Datatype datatype,MPI_Op op,MPI_Comm comm);

  sendbuf   受信バッファの初期アドレス(IN)
  recvbuf   送信バッファの初期アドレス(OUT)
  count     要素数(IN)
  datatype  受信バッファ要素のデータタイプのハンドル(IN)
  op        操作
  comm      コミュニケータ(IN)

  全プロセッサー間で合計、最大値、最小値などを求め、全プロセッサーにその
結果を伝える関数。op には、MPI_SUM、MPI_MIN、MPI_MAX などを指定する。

*/

#endif // TLCOMMUNICATE_H
