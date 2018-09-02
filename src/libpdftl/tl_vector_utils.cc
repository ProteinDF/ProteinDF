#include "tl_vector_utils.h"
#include <iostream>
#include "TlUtils.h"

// -----------------------------------------------------------------------------
// public
// -----------------------------------------------------------------------------
bool TlVectorUtils::isLoadable(const std::string& filePath) {
  const TlVectorUtils::FileSize headerSize =
      TlVectorUtils::getHeaderInfo(filePath);

  bool answer = false;
  if (headerSize > 0) {
    answer = true;
  }
  return answer;
}

TlVectorUtils::FileSize TlVectorUtils::getHeaderInfo(
    const std::string& filepath, TlVectorAbstract::index_type* pSize) {
  FileSize headerSize = 0;

  std::ios_base::sync_with_stdio(false);
  std::fstream fs;
  fs.open(filepath.c_str(), std::ios::in | std::ios::binary);
  if (!fs.fail()) {
    headerSize = TlVectorUtils::getHeaderSize_templ1<std::fstream>(fs, pSize);
  } else {
    std::cerr << TlUtils::format("cannot open vector file: %s @%s:%d",
                                 filepath.c_str(), __FILE__, __LINE__)
              << std::endl;
  }
  fs.close();

  return headerSize;
}

bool TlVectorUtils::load(const std::string& filepath, double* pBuf,
                         const TlVectorAbstract::size_type numOfItems,
                         const TlVectorAbstract::size_type startPos) {
  std::ios_base::sync_with_stdio(false);
  std::fstream fs;
  fs.open(filepath.c_str(), std::ios::in | std::ios::binary);
  if (!fs.fail()) {
    const FileSize headerSize =
        TlVectorUtils::getHeaderSize_templ1<std::fstream>(fs);

    fs.seekg(headerSize + startPos);
    fs.read(reinterpret_cast<char*>(pBuf), sizeof(double) * numOfItems);
  } else {
    std::cerr << TlUtils::format("cannnot open vector file: %s @%s:%d",
                                 filepath.c_str(), __FILE__, __LINE__)
              << std::endl;
  }
  fs.close();

  return true;
}

bool TlVectorUtils::save(const std::string& filepath,
                         const TlVectorAbstract::index_type length,
                         const double* pBuf, const std::size_t numOfItems) {
  std::ios_base::sync_with_stdio(false);

  std::fstream fs(filepath.c_str(), std::ios::out | std::ios::binary);
  fs.write(reinterpret_cast<const char*>(&length),
           sizeof(TlVectorAbstract::index_type));
  fs.write(reinterpret_cast<const char*>(&pBuf), sizeof(double) * numOfItems);
  fs.close();

  return true;
}

// -----------------------------------------------------------------------------
// protected
// -----------------------------------------------------------------------------
template <typename StreamType>
TlVectorUtils::FileSize TlVectorUtils::getHeaderSize_templ1(
    StreamType& s, TlVectorAbstract::index_type* pSize) {
  TlVectorAbstract::index_type size = 0;

  FileSize headerSize = 0;
  // char case:
  {
    headerSize =
        TlVectorUtils::getHeaderSize_templ2<StreamType, char>(s, &size);
  }

  // int case:
  if (headerSize == 0) {
    headerSize = TlVectorUtils::getHeaderSize_templ2<StreamType, int>(s, &size);
  }

  // long case:
  if (headerSize == 0) {
    headerSize =
        TlVectorUtils::getHeaderSize_templ2<StreamType, long>(s, &size);
  }

  if (headerSize > 0) {
    if (pSize != NULL) {
      *pSize = size;
    }
  }

  return headerSize;
}

template <typename StreamType, typename VectorIndexType>
TlVectorUtils::FileSize TlVectorUtils::getHeaderSize_templ2(
    StreamType& s, TlVectorAbstract::index_type* pSize) {
  VectorIndexType size = 0;

  // get file size
  FileSize fileSize = 0;
  {
    s.seekg(0, std::ios_base::end);
    fileSize = s.tellg();
  }

  s.seekg(0, std::ios_base::beg);
  s.read((char*)&(size), sizeof(VectorIndexType));

  const FileSize headerSize = sizeof(VectorIndexType);
  const FileSize estimatedFileSize = headerSize + sizeof(double) * size;

  FileSize answer = 0;
  if (estimatedFileSize == fileSize) {
    answer = headerSize;
  }
  if (pSize != NULL) {
    *pSize = static_cast<TlVectorAbstract::index_type>(size);
  }

  return answer;
}
