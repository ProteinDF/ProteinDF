#include <algorithm>
#include <cassert>
#include <functional>

#include "TlCommunicate.h"
#include "TlUtils.h"
#include "scalapack.h"
#include "tl_dense_vector_impl_lapack.h"
#include "tl_dense_vector_impl_scalapack.h"
#include "tl_dense_vector_lapack.h"
#include "tl_scalapack_context.h"

#include "TlTime.h"

// const std::size_t TlDenseVector_ImplScalapack::FILE_BUFFER_SIZE =
//     100 * 1024 * 1024;  // 100 MB

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVector_ImplScalapack::TlDenseVector_ImplScalapack(
    const TlDenseVectorObject::index_type nSize)
    : TlDenseScalapackObject(nSize, 1) {}

TlDenseVector_ImplScalapack::TlDenseVector_ImplScalapack(
    const TlDenseVector_ImplScalapack& rhs)
    : TlDenseScalapackObject(rhs) {}

TlDenseVector_ImplScalapack::TlDenseVector_ImplScalapack(
    const TlDenseVector_ImplLapack& rhs)
    : TlDenseScalapackObject(rhs.getSize()) {
    TlDenseVectorObject::size_type size = rhs.getSize();
    for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
        this->set(i, rhs.get(i));
    }
}

TlDenseVector_ImplScalapack::~TlDenseVector_ImplScalapack() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlDenseVectorObject::size_type TlDenseVector_ImplScalapack::getSize() const {
    return TlDenseScalapackObject::getNumOfRows();
}

void TlDenseVector_ImplScalapack::resize(
    const TlDenseVectorObject::index_type size) {
    assert(size > 0);
    TlDenseScalapackObject::resize(size, 1);
}

double TlDenseVector_ImplScalapack::get(
    const TlDenseVectorObject::size_type index) const {
    return TlDenseScalapackObject::get(index, 0);
}

void TlDenseVector_ImplScalapack::set(
    const TlDenseVectorObject::size_type index, const double value) {
    TlDenseScalapackObject::set(index, 0, value);
}

void TlDenseVector_ImplScalapack::add(
    const TlDenseVectorObject::size_type index, const double value) {
    TlDenseScalapackObject::add(index, 0, value);
}

void TlDenseVector_ImplScalapack::mul(
    const TlDenseVectorObject::size_type index, const double value) {
    TlDenseScalapackObject::mul(index, 0, value);
}

std::vector<double> TlDenseVector_ImplScalapack::getVector() const {
    std::vector<double> ans(this->getSize(), 0.0);
    const TlMatrixObject::index_type numOfGlobalCols = this->getNumOfCols();

    const TlMatrixObject::index_type maxLocalRow = this->getNumOfMyRows();
    const TlMatrixObject::index_type maxLocalCol = this->getNumOfMyCols();
    for (TlMatrixObject::index_type localRow = 0; localRow < maxLocalRow;
         ++localRow) {
        const TlMatrixObject::index_type globalRow =
            this->getLocal2Global_row(localRow);

        for (TlMatrixObject::index_type localCol = 0; localCol < maxLocalCol;
             ++localCol) {
            const TlMatrixObject::index_type globalCol =
                this->getLocal2Global_col(localCol);
            TlMatrixObject::index_type localIndex =
                this->getLocalIndex(localRow, localCol);

            const TlMatrixObject::index_type globalIndex =
                globalRow * numOfGlobalCols + globalCol;
            ans[globalIndex] = this->pData_[localIndex];
        }
    }

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(ans);

    return ans;
}

// ---------------------------------------------------------------------------
// operators
// ---------------------------------------------------------------------------
TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator=(
    const TlDenseVector_ImplScalapack& rhs) {
    if (&rhs != this) {
        TlDenseScalapackObject::operator=(rhs);
    }

    return *this;
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator+=(
    const TlDenseVector_ImplScalapack& rhs) {
    TlDenseScalapackObject::operator+=(rhs);

    return (*this);
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator-=(
    const TlDenseVector_ImplScalapack& rhs) {
    TlDenseScalapackObject::operator-=(rhs);

    return (*this);
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator*=(
    const double coef) {
    TlDenseScalapackObject::operator*=(coef);

    return (*this);
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::operator/=(
    const double coef) {
    return this->operator*=(1.0 / coef);
}

double TlDenseVector_ImplScalapack::operator*(
    const TlDenseVector_ImplScalapack& rhs) const {
    return this->dot(rhs);
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
void TlDenseVector_ImplScalapack::sortByGreater() {
    // TODO: implement
}

double TlDenseVector_ImplScalapack::dot(
    const TlDenseVector_ImplScalapack& rhs) const {
    // this const_cast is requiered for PGI compiler
    // "error: expression must be an lvalue or a function designator"
    TlDenseVector_ImplScalapack& Xtmp =
        const_cast<TlDenseVector_ImplScalapack&>(*this);
    TlDenseVector_ImplScalapack& Ytmp =
        const_cast<TlDenseVector_ImplScalapack&>(rhs);

    assert(this->getSize() == rhs.getSize());

    const int N = this->getSize();
    double DOT = 0.0;
    const int IX = 1;
    const int JX = 1;
    const int INCX = 1;
    const int IY = 1;
    const int JY = 1;
    const int INCY = 1;

    pddot_(&N, &DOT, Xtmp.pData_, &IX, &JX, Xtmp.pDESC_, &INCX, Ytmp.pData_,
           &IY, &JY, Ytmp.pDESC_, &INCY);

    // INCX=1の場合、列のプロセスにのみ値が返る。
    // したがって、全プロセスに値を送信する必要がある。
    {
        TlCommunicate& rComm = TlCommunicate::getInstance();
        rComm.broadcast(DOT);
    }

    return DOT;
}

TlDenseVector_ImplScalapack& TlDenseVector_ImplScalapack::dotInPlace(
    const TlDenseVector_ImplScalapack& rhs) {
    assert(this->getSize() == rhs.getSize());

    const TlDenseVectorObject::index_type localSize =
        this->getNumOfMyRows() * this->getNumOfMyCols();
    std::transform(this->pData_, this->pData_ + localSize, rhs.pData_,
                   this->pData_, std::multiplies<double>());

    return *this;
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// I/O
// ----------------------------------------------------------------------------
bool TlDenseVector_ImplScalapack::load(const std::string& sFilePath) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    bool answer = false;

    std::ifstream fs;
    if (rComm.isMaster() == true) {
        fs.open(sFilePath.c_str(), std::ios::in | std::ios::binary);
    }

    bool bIsFail = false;
    if (rComm.isMaster() == true) {
        bIsFail = fs.fail();
    }
    rComm.broadcast(bIsFail);

    if (bIsFail) {
        if (rComm.isMaster() == true) {
            this->log_.critical(
                TlUtils::format("[error] TlDistributedVector::load(): could "
                                "not open file. : %s",
                                sFilePath.c_str()));
        }
        abort();
    }

    answer = this->load(fs);

    if (answer != true) {
        this->log_.critical(
            TlUtils::format("Not supported file type: %s @%s.%d",
                            sFilePath.c_str(), __FILE__, __LINE__));
        std::abort();
    }

    if (rComm.isMaster() == true) {
        fs.close();
    }

    return answer;
}

bool TlDenseVector_ImplScalapack::load(std::ifstream& ifs) {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int myRank = rComm.getRank();
    assert(rComm.checkNonBlockingCommunications());

    const int numOfProcs = rComm.getNumOfProc();

    // read header
    bool answer = true;
    int globalSize = 0;
    if (rComm.isMaster() == true) {
        ifs.read((char*)&(globalSize), sizeof(int));
    }
    rComm.broadcast(globalSize);
    this->resize(globalSize);
    // std::cerr << TlUtils::format("[%d] resize: %d", rComm.getRank(),
    // globalSize) << std::endl;

    // read contents
    if (rComm.isMaster() == true) {
        static const std::size_t bufferCount =
            FILE_BUFFER_SIZE / sizeof(double);
        std::vector<double> buf(bufferCount, 0.0);

        TlVectorObject::index_type count = 0;
        const TlVectorObject::index_type maxCount = globalSize;
        bool isFinished = false;

        std::vector<TlMatrixObject::index_type> sizeLists(numOfProcs, 0);
        std::vector<std::vector<TlVectorObject::VectorElement> > elementsBuf(
            numOfProcs);
        std::vector<bool> isSendData(numOfProcs, false);

        while (isFinished == false) {
            // buffer分を一度に読み込み
            ifs.read((char*)(&buf[0]), sizeof(double) * bufferCount);

            // 各プロセスのバッファに振り分ける
            std::vector<std::vector<TlVectorObject::VectorElement> > elements(
                numOfProcs);
            for (std::size_t i = 0; i < bufferCount; ++i) {
                int rank = 0;
                this->getGlobalRowCol2LocalRowCol(count, 0, &rank);
                elements[rank].push_back(
                    TlVectorObject::VectorElement(count, buf[i]));

                ++count;
                if (count >= maxCount) {
                    isFinished = true;
                    break;
                }
            }

            // データを送信
            for (int proc = 1; proc < numOfProcs;
                 ++proc) {  // proc == 0 は送信しない
                if (isSendData[proc] == true) {
                    rComm.wait(sizeLists[proc]);
                    if (sizeLists[proc] > 0) {
                        rComm.wait(&(elementsBuf[proc][0]));
                    }
                    isSendData[proc] = false;
                }

                sizeLists[proc] = elements[proc].size();
                elementsBuf[proc] = elements[proc];

                rComm.iSendData(sizeLists[proc], proc, TAG_LOAD_SIZE);
                // std::cerr << TlUtils::format("[0]->[%d]: TAG_LOAD_SIZE=%d",
                // proc,
                //                              sizeLists[proc])
                //           << std::endl;
                if (sizeLists[proc] > 0) {
                    rComm.iSendDataX(&(elementsBuf[proc][0]), sizeLists[proc],
                                     proc, TAG_LOAD_VALUES);
                }
                isSendData[proc] = true;
            }

            // proc=0分データの書き込み
            {
                const TlVectorObject::index_type sizeList = elements[0].size();
                for (TlVectorObject::index_type i = 0; i < sizeList; ++i) {
                    const TlVectorObject::index_type globalIndex =
                        elements[0][i].index;

                    int rank = 0;
                    TlMatrixObject::index_type myRow;
                    TlMatrixObject::index_type myCol;
                    this->getGlobalRowCol2LocalRowCol(globalIndex, 0, &rank,
                                                      &myRow, &myCol);
                    assert(rank == myRank);
                    TlVectorObject::size_type index =
                        this->getLocalIndex(myRow, myCol);
                    // std::cerr << TlUtils::format("[0] vec[%d]=% 8.3f",
                    // globalIndex,
                    //                              elements[0][i].value)
                    //           << std::endl;
                    this->pData_[index] = elements[0][i].value;
                }
            }
        }  // end while

        for (int proc = 1; proc < numOfProcs;
             ++proc) {  // proc == 0 は送信しない
            if (isSendData[proc] == true) {
                rComm.wait(sizeLists[proc]);
                if (sizeLists[proc] > 0) {
                    rComm.wait(&(elementsBuf[proc][0]));
                }
                isSendData[proc] = false;
            }
        }
        assert(rComm.checkNonBlockingCommunications());

        // 終了メッセージを全ノードに送る
        // sizeList=-1 が送られると終了
        std::vector<int> endMsg(numOfProcs, -1);
        for (int proc = 1; proc < numOfProcs;
             ++proc) {  // proc == 0 は送信しない
            rComm.iSendData(endMsg[proc], proc, TAG_LOAD_SIZE);
        }
        for (int proc = 1; proc < numOfProcs;
             ++proc) {  // proc == 0 は送信しない
            rComm.wait(endMsg[proc]);
        }
    } else {
        // slave
        const int root = 0;
        TlMatrixObject::index_type sizeList = 0;
        std::vector<TlVectorObject::VectorElement> elements;
        int endMsg = 0;

        rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
        bool isLoopBreak = false;
        while (isLoopBreak == false) {
            if (rComm.test(sizeList) == true) {
                rComm.wait(sizeList);
                // std::cerr << TlUtils::format("[%d] TAG_LOAD_SIZE=%d",
                // rComm.getRank(),
                //                              sizeList)
                //           << std::endl;
                if (sizeList >= 0) {
                    if (sizeList > 0) {
                        elements.resize(sizeList);
                        rComm.receiveDataX(&(elements[0]), sizeList, root,
                                           TAG_LOAD_VALUES);

                        for (TlMatrixObject::index_type i = 0; i < sizeList;
                             ++i) {
                            const TlMatrixObject::index_type globalIndex =
                                elements[i].index;

                            int rank = 0;
                            TlMatrixObject::index_type myRow, myCol;
                            this->getGlobalRowCol2LocalRowCol(
                                globalIndex, 0, &rank, &myRow, &myCol);
                            assert(rank == myRank);
                            TlVectorObject::size_type index =
                                this->getLocalIndex(myRow, myCol);
                            // std::cerr << TlUtils::format("[%d]
                            // vec[%d]=% 8.3f", rComm.getRank(),
                            //                              globalIndex,
                            //                              elements[i].value)
                            //           << std::endl;
                            this->pData_[index] = elements[i].value;
                        }
                    }
                    rComm.iReceiveData(sizeList, root, TAG_LOAD_SIZE);
                } else {
                    assert(sizeList < 0);
                    isLoopBreak = true;
                    // std::cerr << TlUtils::format("[%d] done.",
                    // rComm.getRank())
                    //           << std::endl;
                }
            }
        }
    }

    assert(rComm.checkNonBlockingCommunications());
    return answer;
}

bool TlDenseVector_ImplScalapack::save(const std::string& sFilePath) const {
    bool answer = true;

    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications());

    if (rComm.isMaster() == true) {
        // master
        const TlVectorObject::index_type globalSize = this->getSize();
        TlDenseVector_Lapack v(globalSize);

        // store local matrix
        {
            const std::vector<TlVectorObject::VectorElement> buf =
                this->getVectorElementsInLocal();
            this->saveElements(&v, buf);
        }

        // recive submatrix & write
        // std::cerr << TlUtils::format("0: master start") << std::endl;
        const int numOfSlaves = rComm.getNumOfProc() - 1;
        for (int i = 0; i < numOfSlaves; ++i) {
            int src = 0;
            TlVectorObject::size_type bufSize = 0;
            rComm.iReceiveDataFromAnySource(bufSize, TAG_SAVE_HANDSHAKE);
            rComm.wait(&bufSize, &src);

            if (bufSize > 0) {
                std::vector<TlVectorObject::VectorElement> buf(bufSize);
                rComm.receiveDataX(&(buf[0]), bufSize, src, TAG_SAVE_DATA);
                this->saveElements(&v, buf);
            }
            // std::cerr << TlUtils::format("0: recv data src=%d", src) <<
            // std::endl;
        }

        answer = v.save(sFilePath);
    } else {
        // slave: send submatrix
        const int root = 0;
        const std::vector<TlVectorObject::VectorElement> buf =
            this->getVectorElementsInLocal();
        const TlVectorObject::size_type bufSize = buf.size();
        rComm.sendData(bufSize, root, TAG_SAVE_HANDSHAKE);
        if (bufSize > 0) {
            rComm.iSendDataX(&(buf[0]), buf.size(), root, TAG_SAVE_DATA);
            rComm.wait(&(buf[0]));
        }
    }

    rComm.broadcast(answer);
    assert(rComm.checkNonBlockingCommunications());
    return answer;
}

// ---------------------------------------------------------------------------
// protected
// ---------------------------------------------------------------------------
std::vector<TlVectorObject::VectorElement>
TlDenseVector_ImplScalapack::getVectorElementsInLocal() const {
    const TlVectorObject::index_type numOfRows = this->getNumOfRows();
    const TlVectorObject::index_type numOfCols = this->getNumOfCols();
    const TlVectorObject::index_type maxLocalRows = this->getNumOfMyRows();
    const TlVectorObject::index_type maxLocalCols = this->getNumOfMyCols();

    const TlVectorObject::size_type numOfMyElements =
        this->getNumOfMyElements();
    std::vector<TlVectorObject::VectorElement> answer;
    answer.reserve(numOfMyElements);

    for (TlVectorObject::index_type localRow = 0; localRow < maxLocalRows;
         ++localRow) {
        const TlVectorObject::index_type globalRow =
            this->getLocal2Global_row(localRow);
        if ((globalRow < 0) || (numOfRows <= globalRow)) {
            continue;
        }

        for (TlVectorObject::index_type localCol = 0; localCol < maxLocalCols;
             ++localCol) {
            const TlVectorObject::index_type globalCol =
                this->getLocal2Global_col(localCol);
            if ((globalCol < 0) || (numOfCols <= globalCol)) {
                continue;
            }

            const TlMatrixObject::index_type localIndex =
                this->getLocalIndex(localRow, localCol);
            answer.push_back(
                TlVectorObject::VectorElement(globalRow, pData_[localIndex]));
        }
    }

    // {
    //   TlCommunicate& rComm = TlCommunicate::getInstance();
    //   std::cerr << TlUtils::format("[%d] getVectorElementsInLocal(): %d x %d
    //   (%d x %d) [%d/%d]x[%d/%d] -> %d",
    //    rComm.getRank(),
    //   this->getNumOfRows(), this->getNumOfCols(), this->getNumOfMyRows(),
    //   this->getNumOfMyCols(), this->m_nMyProcRow, this->procGridRow_,
    //   this->m_nMyProcCol, this->procGridCol_, answer.size()) << std::endl;
    // }

    return answer;
}

void TlDenseVector_ImplScalapack::saveElements(
    TlDenseVector_Lapack* pVector,
    const std::vector<TlVectorObject::VectorElement>& elements) const {
    std::vector<TlVectorObject::VectorElement>::const_iterator itEnd =
        elements.end();
    for (std::vector<TlVectorObject::VectorElement>::const_iterator it =
             elements.begin();
         it != itEnd; ++it) {
        assert(0 <= it->index);
        assert(it->index < this->getSize());

        pVector->set(it->index, it->value);
    }
}
