#ifndef TLMATRIXUTILS_H
#define TLMATRIXUTILS_H

#include <fstream>
#include <string>

#include "tl_matrix_object.h"

class TlMatrixUtils {
public:
    typedef std::char_traits<char>::off_type FileSize;

public:
    static bool isLoadable(const std::string& filePath, const TlMatrixObject::MatrixType matrixType);

public:
    /// ヘッダ情報を読み取る
    /// header情報を読み取れた場合はtrueを返す
    /// 読み取り位置はヘッダの終了を指し示す。
    static bool getHeaderInfo(const std::string& filePath, TlMatrixObject::HeaderInfo* pHeaderInfo = NULL);

    /// ヘッダ情報を読み取る
    /// header情報を読み取れた場合はtrueを返す
    /// 読み取り位置はヘッダの終了を指し示す。
    static bool getHeaderInfo(std::ifstream& ifs, TlMatrixObject::HeaderInfo* pHeaderInfo = NULL);
    static bool getHeaderInfo(std::fstream& fs, TlMatrixObject::HeaderInfo* pHeaderInfo = NULL);

public:
    // static TlMatrixObject::size_type loadMatrix(
    //     const std::string& filepath,
    //     const TlMatrixObject::MatrixType inMatrixType, double* pBuf,
    //     const TlMatrixObject::size_type numOfItems, const FileSize startPos =
    //     0);

    static bool saveMatrix(const std::string& filepath, const TlMatrixObject::MatrixType matrixType,
                           const TlMatrixObject::index_type rows, const TlMatrixObject::index_type cols,
                           const double* pBuf, const std::size_t numOfItems);

public:
    static void CSFD2RSFD(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                          const double* pBufIn, double* pBufOut);
    static void RSFD2CSFD(const TlMatrixObject::index_type row, const TlMatrixObject::index_type col,
                          const double* pBufIn, double* pBufOut);

public:
    /// リストからまとめて要素を加算する
    ///
    /// @param[in] pIndexPairs インデックスのリスト(row1, col1, row2, col2, ...)
    /// @param[in] pValues 値のリスト
    /// @param[in] size 要素数
    template <class MatrixClass>
    static void addByList(const TlMatrixObject::index_type* pIndexPairs, const double* pValues, const std::size_t size,
                          MatrixClass* pMatrix);

protected:
    /// header情報を読み取れた場合はtrueを返す
    template <typename stream>
    static bool getHeaderSize_template(stream& s, TlMatrixObject::HeaderInfo* pHeaderInfo = NULL);

    /// header情報を読み取れた場合はtrueを返す
    template <typename StreamType, typename MatrixType, typename IndexType>
    static bool getHeaderSizeByType_template(StreamType& s, TlMatrixObject::HeaderInfo* pHeaderInfo = NULL);

    static std::size_t estimateFileSize(const TlMatrixObject::HeaderInfo& headerInfo);
};

template <class MatrixClass>
void TlMatrixUtils::addByList(const TlMatrixObject::index_type* pIndexPairs, const double* pValues,
                              const std::size_t size, MatrixClass* pMatrix) {
    for (std::size_t i = 0; i < size; ++i) {
        const TlMatrixObject::index_type globalRow = pIndexPairs[i * 2];
        const TlMatrixObject::index_type globalCol = pIndexPairs[i * 2 + 1];
        const double value = pValues[i];

        pMatrix->add(globalRow, globalCol, value);
    }
}

#endif  // TLMATRIXUTILS_H
