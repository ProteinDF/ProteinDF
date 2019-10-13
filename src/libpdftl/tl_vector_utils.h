#ifndef TL_VECTOR_UTILS_H
#define TL_VECTOR_UTILS_H

#include <fstream>
#include <string>
#include "tl_vector_abstract.h"

class TlVectorUtils {
   public:
    typedef std::char_traits<char>::off_type FileSize;

   public:
    static bool isLoadable(const std::string& filePath);

    /// ヘッダ情報を読み取る
    /// header情報を読み取れた場合はヘッダサイズを返す
    /// 読み取れなかった場合は0を返す
    static FileSize getHeaderInfo(const std::string& filePath,
                                  TlVectorAbstract::index_type* pSize = NULL);

   public:
    static bool load(const std::string& filepath, double* pBuf,
                     const TlVectorAbstract::size_type numOfItems,
                     const TlVectorAbstract::size_type startPos = 0);

    static bool save(const std::string& filepath,
                     const TlVectorAbstract::index_type length,
                     const double* pBuf, const std::size_t numOfItems);

   protected:
    /// header情報を読み取れた場合はヘッダサイズを返す
    /// 読み取れなかった場合は0を返す
    template <typename stream>
    static FileSize getHeaderSize_templ1(
        stream& s, TlVectorAbstract::index_type* pSize = NULL);

    /// header情報を読み取れた場合はヘッダサイズを返す
    /// 読み取れなかった場合は0を返す
    template <typename StreamType, typename VectorType>
    static FileSize getHeaderSize_templ2(StreamType& s,
                                         TlVectorAbstract::index_type* pSize);
};

#endif  // TL_VECTOR_UTILS_H
