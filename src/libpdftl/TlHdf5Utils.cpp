#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "TlFile.h"
#include "TlHdf5Utils.h"
#include "TlUtils.h"

const int TlHdf5Utils::initStrLength_ = 256;
const char TlHdf5Utils::delim_ = '/';

TlHdf5Utils::TlHdf5Utils(const std::string& path) : file_() {
  if (TlFile::isExistFile(path)) {
    this->file_ = H5::H5File(path.c_str(), H5F_ACC_RDWR | H5F_ACC_DEBUG);
    // this->file_ = H5::H5File(path, H5F_ACC_TRUNC | H5F_ACC_DEBUG);
  } else {
    this->file_ = H5::H5File(path.c_str(), H5F_ACC_EXCL | H5F_ACC_DEBUG);
  }
}

TlHdf5Utils::~TlHdf5Utils() { this->file_.close(); }

void TlHdf5Utils::flush() const { this->file_.flush(H5F_SCOPE_LOCAL); }

// =====================================================================
// create dataset
// =====================================================================
void TlHdf5Utils::createDataSet_int(const std::string& path,
                                    const std::size_t size) {
  this->createDataSet_common(path, H5::PredType::NATIVE_INT, size);
}

void TlHdf5Utils::createDataSet_double(const std::string& path,
                                       const std::size_t size) {
  this->createDataSet_common(path, H5::PredType::NATIVE_DOUBLE, size);
}

void TlHdf5Utils::createDataSet_common(const std::string& path,
                                       const H5::PredType& predType,
                                       const std::size_t size) {
  assert(size > 0);
  this->prepareToWrite(path);

  H5::DataType memType(predType);
  const int rank = 1;
  const hsize_t dims = size;
  H5::DataSpace memSpace(rank, &dims);

  H5::DataSet dataSet = this->file_.createDataSet(path.c_str(), memType, memSpace);
}

// =====================================================================
// write
// =====================================================================
template <typename T>
void TlHdf5Utils::write_common(const std::string& path, const T value,
                               const H5::PredType& predType) {
  this->prepareToWrite(path);

  H5::DataType memType(predType);
  const int rank = 1;
  const hsize_t dims = 1;
  H5::DataSpace memSpace(rank, &dims);

  H5::DataSet dataSet = this->file_.createDataSet(path.c_str(), memType, memSpace);
  dataSet.write(&value, memType);

  this->flush();
}

void TlHdf5Utils::prepareToWrite(const std::string& path) {
  this->createGroup(this->getDirName(path));
  if (this->hasDataSet(path)) {
    this->removeObject(path);
  }
}

void TlHdf5Utils::write(const std::string& path, const char value) {
  this->write_common(path, value, H5::PredType::NATIVE_CHAR);
}

void TlHdf5Utils::write(const std::string& path, const unsigned char value) {
  this->write_common(path, value, H5::PredType::NATIVE_UCHAR);
}

void TlHdf5Utils::write(const std::string& path, const int value) {
  this->write_common(path, value, H5::PredType::NATIVE_INT);
}

void TlHdf5Utils::write(const std::string& path, const unsigned int value) {
  this->write_common(path, value, H5::PredType::NATIVE_UINT);
}

void TlHdf5Utils::write(const std::string& path, const long value) {
  this->write_common(path, value, H5::PredType::NATIVE_LONG);
}

void TlHdf5Utils::write(const std::string& path, const unsigned long value) {
  this->write_common(path, value, H5::PredType::NATIVE_ULONG);
}

void TlHdf5Utils::write(const std::string& path, const float value) {
  this->write_common(path, value, H5::PredType::NATIVE_FLOAT);
}

void TlHdf5Utils::write(const std::string& path, const double value) {
  this->write_common(path, value, H5::PredType::NATIVE_DOUBLE);
}

void TlHdf5Utils::write(const std::string& path, const std::string& values) {
  this->prepareToWrite(path);

  // H5::StrType dataType(H5::PredType::C_S1, H5T_VARIABLE);
  H5::StrType dataType(H5::PredType::C_S1, values.size());

  const int rank = 1;
  const hsize_t dims = 1;
  H5::DataSpace dataSpace(rank, &dims);
  H5::DataSet dataSet = this->file_.createDataSet(path.c_str(), dataType, dataSpace);
  dataSet.write(values.c_str(), dataType);
}

template <typename T>
void TlHdf5Utils::write_common(const std::string& path, const T* pValues,
                               const H5::PredType& predType,
                               const std::size_t size) {
  this->prepareToWrite(path);

  H5::IntType memType(predType);
  const int rank = 1;
  const hsize_t dims = size;
  H5::DataSpace memSpace(rank, &dims);

  H5::DataSet dataSet = this->file_.createDataSet(path.c_str(), memType, memSpace);
  dataSet.write(pValues, memType);

  this->flush();
}

void TlHdf5Utils::write(const std::string& name, const int* pValues,
                        std::size_t size) {
  this->write_common<int>(name, pValues, H5::PredType::NATIVE_INT, size);
}

void TlHdf5Utils::write(const std::string& name, const unsigned int* pValues,
                        std::size_t size) {
  this->write_common<unsigned int>(name, pValues, H5::PredType::NATIVE_ULONG,
                                   size);
}

void TlHdf5Utils::write(const std::string& name, const long* pValues,
                        std::size_t size) {
  this->write_common<long>(name, pValues, H5::PredType::NATIVE_LONG, size);
}

void TlHdf5Utils::write(const std::string& name, const unsigned long* pValues,
                        std::size_t size) {
  this->write_common<unsigned long>(name, pValues, H5::PredType::NATIVE_ULONG,
                                    size);
}

void TlHdf5Utils::write(const std::string& name, const float* pValues,
                        std::size_t size) {
  this->write_common<float>(name, pValues, H5::PredType::NATIVE_FLOAT, size);
}

void TlHdf5Utils::write(const std::string& name, const double* pValues,
                        std::size_t size) {
  this->write_common<double>(name, pValues, H5::PredType::NATIVE_DOUBLE, size);
}

template <typename T>
void TlHdf5Utils::write_common(const std::string& path,
                               const std::vector<T>& values,
                               const H5::PredType& predType) {
  this->write_common<T>(path, values.data(), predType, values.size());
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<char>& values) {
  this->write_common<char>(path, values, H5::PredType::NATIVE_CHAR);
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<unsigned char>& values) {
  this->write_common<unsigned char>(path, values, H5::PredType::NATIVE_UCHAR);
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<int>& values) {
  this->write_common<int>(path, values, H5::PredType::NATIVE_INT);
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<unsigned int>& values) {
  this->write_common<unsigned int>(path, values, H5::PredType::NATIVE_UINT);
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<long>& values) {
  this->write_common<long>(path, values, H5::PredType::NATIVE_LONG);
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<unsigned long>& values) {
  this->write_common<unsigned long>(path, values, H5::PredType::NATIVE_ULONG);
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<float>& values) {
  this->write_common<float>(path, values, H5::PredType::NATIVE_FLOAT);
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<double>& values) {
  this->write_common<double>(path, values, H5::PredType::NATIVE_DOUBLE);
}

void TlHdf5Utils::write(const std::string& path,
                        const std::vector<std::string>& values) {
  this->prepareToWrite(path);

  H5::StrType dataType(H5::PredType::C_S1, H5T_VARIABLE);
  const int rank = 1;
  const hsize_t dims = values.size();
  H5::DataSpace dataSpace(rank, &dims);
  H5::DataSet dataSet = this->file_.createDataSet(path.c_str(), dataType, dataSpace);

  std::vector<const char*> pBuf(dims);
  for (hsize_t i = 0; i < dims; ++i) {
    pBuf[i] = values[i].c_str();
  }

  dataSet.write(pBuf.data(), dataType);
}

// =====================================================================
// get
// =====================================================================
template <typename T>
void TlHdf5Utils::get(const std::string& path, const H5::PredType& predType,
                      T* pOut) const {
  assert(pOut != NULL);
  const H5::DataSet dataSet = this->file_.openDataSet(path.c_str());

  H5::DataType memType(predType);
  const int rank = 1;
  const hsize_t dims[] = {1};
  H5::DataSpace dataSpace(rank, dims);

  dataSet.read(pOut, memType, dataSpace);
}

void TlHdf5Utils::get(const std::string& path, char* pOut) const {
  this->get<char>(path, H5::PredType::NATIVE_CHAR, pOut);
}

void TlHdf5Utils::get(const std::string& path, unsigned char* pOut) const {
  this->get<unsigned char>(path, H5::PredType::NATIVE_UCHAR, pOut);
}

void TlHdf5Utils::get(const std::string& path, int* pOut) const {
  this->get<int>(path, H5::PredType::NATIVE_INT, pOut);
}

void TlHdf5Utils::get(const std::string& path, unsigned int* pOut) const {
  this->get<unsigned int>(path, H5::PredType::NATIVE_UINT, pOut);
}

void TlHdf5Utils::get(const std::string& path, long* pOut) const {
  this->get<long>(path, H5::PredType::NATIVE_LONG, pOut);
}

void TlHdf5Utils::get(const std::string& path, unsigned long* pOut) const {
  this->get<unsigned long>(path, H5::PredType::NATIVE_ULONG, pOut);
}

void TlHdf5Utils::get(const std::string& path, float* pOut) const {
  this->get<float>(path, H5::PredType::NATIVE_FLOAT, pOut);
}

void TlHdf5Utils::get(const std::string& path, double* pOut) const {
  this->get<double>(path, H5::PredType::NATIVE_DOUBLE, pOut);
}

void TlHdf5Utils::get(const std::string& path, std::string* pOut) const {
  assert(pOut != NULL);
  const H5::DataSet dataSet = this->file_.openDataSet(path.c_str());

  H5::DataSpace dataSpace = dataSet.getSpace();
  const int rank = dataSpace.getSimpleExtentNdims();
  assert(rank > 0);

  // std::vector<hsize_t> dims(rank);
  // const int nDims = dataSpace.getSimpleExtentDims(dims.data(), NULL);
  // assert(rank == nDims);
  // hsize_t bufferSize = 1;
  // for (int dim = 0; dim < nDims; ++dim) {
  //   bufferSize *= dims[dim];
  // }

  H5::DataSpace memSpace(dataSpace);
  hsize_t count =  1; // bufferSize;
  hsize_t offset = 0;
  memSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
  dataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);

  const std::size_t dataSize = dataSet.getInMemDataSize();

  // H5::StrType memType(H5::PredType::C_S1, H5T_VARIABLE);
  H5::StrType memType(H5::PredType::C_S1, dataSize);

  char* pBuf = new char[dataSize +1];

  dataSet.read(pBuf, memType, memSpace);
  *pOut = std::string(pBuf);

  delete[] pBuf;
  pBuf = NULL;
}


template <typename T>
void TlHdf5Utils::get_common(const std::string& path,
                             const H5::PredType& predType, T* pOut,
                             const std::size_t size) const {
  const H5::DataSet dataSet = this->file_.openDataSet(path.c_str());

  H5::DataSpace dataSpace = dataSet.getSpace();
  const int rank = dataSpace.getSimpleExtentNdims();
  assert(rank > 0);

  std::vector<hsize_t> dims(rank);
  const int nDims = dataSpace.getSimpleExtentDims(dims.data(), NULL);
  assert(rank == nDims);

  hsize_t bufferSize = 1;
  for (int dim = 0; dim < nDims; ++dim) {
    bufferSize *= dims[dim];
  }

  H5::DataSpace memSpace(dataSpace);
  hsize_t count = size;
  hsize_t offset = 0;
  memSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
  dataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);

  H5::DataType memType(predType);
  dataSet.read(pOut, memType, memSpace, dataSpace);
}

void TlHdf5Utils::get(const std::string& path, double* pOut,
                      const std::size_t size) const {
  this->get_common(path, H5::PredType::NATIVE_DOUBLE, pOut, size);
}

template <typename T>
void TlHdf5Utils::get(const std::string& path, const H5::PredType& predType,
                      std::vector<T>* pOut) const {
  const H5::DataSet dataSet = this->file_.openDataSet(path.c_str());

  H5::DataSpace dataSpace = dataSet.getSpace();
  const int rank = dataSpace.getSimpleExtentNdims();
  assert(rank > 0);

  std::vector<hsize_t> dims(rank);
  const int nDims = dataSpace.getSimpleExtentDims(dims.data(), NULL);
  assert(rank == nDims);

  hsize_t bufferSize = 1;
  for (int dim = 0; dim < nDims; ++dim) {
    bufferSize *= dims[dim];
  }

  pOut->resize(bufferSize);
  H5::DataType memType(predType);
  dataSet.read(&((*pOut)[0]), memType);
}

void TlHdf5Utils::get(const std::string& path, std::vector<char>* pOut) const {
  this->get<char>(path, H5::PredType::NATIVE_CHAR, pOut);
}

void TlHdf5Utils::get(const std::string& path,
                      std::vector<unsigned char>* pOut) const {
  this->get<unsigned char>(path, H5::PredType::NATIVE_UCHAR, pOut);
}

void TlHdf5Utils::get(const std::string& path, std::vector<int>* pOut) const {
  this->get<int>(path, H5::PredType::NATIVE_INT, pOut);
}

void TlHdf5Utils::get(const std::string& path,
                      std::vector<unsigned int>* pOut) const {
  this->get<unsigned int>(path, H5::PredType::NATIVE_UINT, pOut);
}

void TlHdf5Utils::get(const std::string& path, std::vector<long>* pOut) const {
  this->get<long>(path, H5::PredType::NATIVE_LONG, pOut);
}

void TlHdf5Utils::get(const std::string& path,
                      std::vector<unsigned long>* pOut) const {
  this->get<unsigned long>(path, H5::PredType::NATIVE_ULONG, pOut);
}

void TlHdf5Utils::get(const std::string& path, std::vector<float>* pOut) const {
  this->get<float>(path, H5::PredType::NATIVE_FLOAT, pOut);
}

void TlHdf5Utils::get(const std::string& path,
                      std::vector<double>* pOut) const {
  this->get<double>(path, H5::PredType::NATIVE_DOUBLE, pOut);
}

void TlHdf5Utils::get(const std::string& path,
                      std::vector<std::string>* pOut) const {
  const H5::DataSet dataSet = this->file_.openDataSet(path.c_str());

  H5::DataSpace dataSpace = dataSet.getSpace();
  const int rank = dataSpace.getSimpleExtentNdims();
  assert(rank == 1);

  std::vector<hsize_t> dims(rank);
  const int nDims = dataSpace.getSimpleExtentDims(dims.data(), NULL);
  assert(rank == nDims);

  hsize_t bufferSize = 1;
  for (int dim = 0; dim < nDims; ++dim) {
    bufferSize *= dims[dim];
  }

  std::vector<char*> pStrs(bufferSize);
  pOut->resize(bufferSize);
  const H5::StrType memType(H5::PredType::C_S1, H5T_VARIABLE);
  dataSet.read(pStrs.data(), memType);
  for (hsize_t i = 0; i < bufferSize; ++i) {
    (*pOut)[i] = std::string(pStrs[i]);
  }
}

// ---------------------------------------------------------------------------
template <typename T>
void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const H5::PredType& predType,
                                      const std::vector<T>& data) {
  assert(coord.size() == data.size());

  const H5::DataSet dataSet = this->file_.openDataSet(path.c_str());

  H5::DataSpace dataSpace = dataSet.getSpace();

  const std::size_t coordSize = coord.size();
  dataSpace.selectElements(H5S_SELECT_SET, coordSize, &(coord[0]));

  H5::DataType memType(predType);

  const int memRank = 1;
  std::vector<hsize_t> memDims(memRank);
  memDims[0] = coordSize;
  H5::DataSpace memSpace(memRank, memDims.data());

  dataSet.write(&(data[0]), memType, memSpace, dataSpace);
  this->flush();
}

void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const std::vector<char>& data) {
  this->setSelectedElements<>(path, coord, H5::PredType::NATIVE_CHAR, data);
}

void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const std::vector<unsigned char>& data) {
  this->setSelectedElements<>(path, coord, H5::PredType::NATIVE_UCHAR, data);
}

void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const std::vector<int>& data) {
  this->setSelectedElements<>(path, coord, H5::PredType::NATIVE_INT, data);
}

void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const std::vector<unsigned int>& data) {
  this->setSelectedElements<>(path, coord, H5::PredType::NATIVE_UINT, data);
}

void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const std::vector<long>& data) {
  this->setSelectedElements<>(path, coord, H5::PredType::NATIVE_LONG, data);
}

void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const std::vector<unsigned long>& data) {
  this->setSelectedElements<>(path, coord, H5::PredType::NATIVE_ULONG, data);
}

void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const std::vector<float>& data) {
  this->setSelectedElements<>(path, coord, H5::PredType::NATIVE_FLOAT, data);
}

void TlHdf5Utils::setSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const std::vector<double>& data) {
  this->setSelectedElements<>(path, coord, H5::PredType::NATIVE_DOUBLE, data);
}

// ------------------------------------------------------------------
template <typename T>
void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      const H5::PredType& predType,
                                      std::vector<T>* pOut) {
  const H5::DataSet dataSet = this->file_.openDataSet(path.c_str());

  H5::DataSpace dataSpace = dataSet.getSpace();

  const std::size_t coordSize = coord.size();
  dataSpace.selectElements(H5S_SELECT_SET, coordSize, &(coord[0]));

  H5::DataType memType(predType);

  const int memRank = 1;
  std::vector<hsize_t> memDims(memRank);
  memDims[0] = coordSize;
  H5::DataSpace memSpace(memRank, memDims.data());

  pOut->resize(coordSize);

  dataSet.read(&((*pOut)[0]), memType, memSpace, dataSpace);
}

void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      std::vector<char>* pOut) {
  this->getSelectedElements(path, coord, H5::PredType::NATIVE_CHAR, pOut);
}

void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      std::vector<unsigned char>* pOut) {
  this->getSelectedElements(path, coord, H5::PredType::NATIVE_UCHAR, pOut);
}

void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      std::vector<int>* pOut) {
  this->getSelectedElements(path, coord, H5::PredType::NATIVE_INT, pOut);
}

void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      std::vector<unsigned int>* pOut) {
  this->getSelectedElements(path, coord, H5::PredType::NATIVE_UINT, pOut);
}

void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      std::vector<long>* pOut) {
  this->getSelectedElements(path, coord, H5::PredType::NATIVE_LONG, pOut);
}

void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      std::vector<unsigned long>* pOut) {
  this->getSelectedElements(path, coord, H5::PredType::NATIVE_ULONG, pOut);
}

void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      std::vector<float>* pOut) {
  this->getSelectedElements(path, coord, H5::PredType::NATIVE_FLOAT, pOut);
}

void TlHdf5Utils::getSelectedElements(const std::string& path,
                                      const std::vector<hsize_t>& coord,
                                      std::vector<double>* pOut) {
  this->getSelectedElements(path, coord, H5::PredType::NATIVE_DOUBLE, pOut);
}

// ---------------------------------------------------------------------------
// attribute
// ---------------------------------------------------------------------------
template <typename T>
void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const H5::PredType& predType, const T value) {
  H5::DataType memType(predType);
  const int rank = 1;
  const hsize_t dims[] = {1};
  H5::DataSpace memSpace(rank, dims);

  H5::DataSet dataSet = this->file_.openDataSet(path.c_str());
  H5::Attribute attr = dataSet.createAttribute(attrName.c_str(), memType, memSpace);
  attr.write(memType, &value);
}

void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const char value) {
  this->setAttr(path, attrName, H5::PredType::NATIVE_CHAR, value);
}

void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const unsigned char value) {
  this->setAttr(path, attrName, H5::PredType::NATIVE_UCHAR, value);
}

void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const int value) {
  this->setAttr(path, attrName, H5::PredType::NATIVE_INT, value);
}

void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const unsigned int value) {
  this->setAttr(path, attrName, H5::PredType::NATIVE_UINT, value);
}

void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const long value) {
  this->setAttr(path, attrName, H5::PredType::NATIVE_LONG, value);
}

void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const unsigned long value) {
  this->setAttr(path, attrName, H5::PredType::NATIVE_ULONG, value);
}

void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const float value) {
  this->setAttr(path, attrName, H5::PredType::NATIVE_FLOAT, value);
}

void TlHdf5Utils::setAttr(const std::string& path, const std::string& attrName,
                          const double value) {
  this->setAttr(path, attrName, H5::PredType::NATIVE_DOUBLE, value);
}

template <typename T>
void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          const H5::PredType& predType, T* pValue) {
  assert(pValue != NULL);

  H5::DataSet dataSet = this->file_.openDataSet(path.c_str());
  H5::Attribute attr = dataSet.openAttribute(attrName.c_str());

  H5::DataType memType(predType);
  // const int rank = 1;
  // const hsize_t dims[] = {1};
  // H5::DataSpace memSpace(rank, dims);

  attr.read(memType, pValue);
}

void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          char* pValue) {
  this->getAttr(path, attrName, H5::PredType::NATIVE_CHAR, pValue);
}

void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          unsigned char* pValue) {
  this->getAttr(path, attrName, H5::PredType::NATIVE_UCHAR, pValue);
}

void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          int* pValue) {
  this->getAttr(path, attrName, H5::PredType::NATIVE_INT, pValue);
}

void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          unsigned int* pValue) {
  this->getAttr(path, attrName, H5::PredType::NATIVE_UINT, pValue);
}

void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          long* pValue) {
  this->getAttr(path, attrName, H5::PredType::NATIVE_LONG, pValue);
}

void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          unsigned long* pValue) {
  this->getAttr(path, attrName, H5::PredType::NATIVE_ULONG, pValue);
}

void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          float* pValue) {
  this->getAttr(path, attrName, H5::PredType::NATIVE_FLOAT, pValue);
}

void TlHdf5Utils::getAttr(const std::string& path, const std::string& attrName,
                          double* pValue) {
  this->getAttr(path, attrName, H5::PredType::NATIVE_DOUBLE, pValue);
}

// ---------------------------------------------------------------------------
void TlHdf5Utils::createGroup(const std::string& path) {
  if (path.length() > 0) {
    const std::string dirName = this->getDirName(path);
    if (dirName.length() > 0) {
      if (this->hasGroup(dirName) != true) {
        this->createGroup(dirName);
      }
    }

    if (this->hasGroup(path) != true) {
      H5::Group group = this->file_.createGroup(path.c_str());
    }
  }
}

// ---------------------------------------------------------------------------
bool TlHdf5Utils::hasGroup(const std::string& path) {
  return this->hasObject(&(this->file_), path, H5O_TYPE_GROUP);
}

// ---------------------------------------------------------------------------
bool TlHdf5Utils::hasDataSet(const std::string& path) {
  return this->hasObject(&(this->file_), path, H5O_TYPE_DATASET);
}

bool TlHdf5Utils::hasGroup(const H5::Group* pParent,
                           const std::string& name) {
  return this->hasObject(pParent, name, H5O_TYPE_GROUP);
}

bool TlHdf5Utils::hasDataSet(const H5::Group* pParent,
                             const std::string& name) {
  return this->hasObject(pParent, name, H5O_TYPE_DATASET);
}

bool TlHdf5Utils::hasObject(const H5::Group* pParent,
                            const std::string& name,
                            const H5O_type_t requestType) {
  bool answer = false;

  std::string rest = "";
  std::string cname = this->getChildName(name, &rest);
  // std::cerr << "hasObject():" << name << "> " << cname << ", " << rest <<
  // std::endl;

  if (rest.length() > 0) {
    if (this->hasGroup(pParent, cname)) {
      H5::Group grp = this->getGroup(pParent, cname);
      answer = this->hasObject(&grp, rest, requestType);
    } else {
      // NOT FOUND!
      // std::cerr << "path not found: " << name << "(" << cname << ")" <<
      // std::endl;
    }
  } else {
    answer = this->hasChildObject(pParent, cname, requestType);
  }

  return answer;
}

bool TlHdf5Utils::hasChildObject(const H5::Group* pParent,
                                 const std::string& name,
                                 const H5O_type_t requestType) {
  bool answer = false;
  const hsize_t numOfObjs = pParent->getNumObjs();
  for (hsize_t i = 0; i < numOfObjs; ++i) {
    // const std::string objName = pParent->getObjnameByIdx(i);
    const std::string objName = this->getObjnameByIdx(pParent, i, TlHdf5Utils::initStrLength_);
    if (objName == name) {
      const H5O_type_t objType = pParent->childObjType(objName.c_str());
      if (objType == requestType) {
        answer = true;
      }
      break;
    }
  }

  return answer;
}

std::string TlHdf5Utils::getObjnameByIdx(const H5::Group* pParent, const hsize_t index, const std::size_t size) {
  std::string answer;
  char* pBuf = new char[size +1];

  const std::size_t copied = pParent->getObjnameByIdx(index, pBuf, size);
  if (copied < size) {
    answer = std::string(pBuf);
    delete[] pBuf;
    pBuf = NULL;
  } else {
    delete[] pBuf;
    pBuf = NULL;
    answer = this->getObjnameByIdx(pParent, index, size * 2);
  }

  return answer;
}

// ---------------------------------------------------------------------------
H5::Group TlHdf5Utils::getGroup(const H5::Group* pParent,
                                const std::string& name) {
  assert(this->hasGroup(pParent, name));

  return pParent->openGroup(name.c_str());
}

// ---------------------------------------------------------------------------
void TlHdf5Utils::removeObject(const std::string& name) {
  this->file_.unlink(name.c_str());
}

// ---------------------------------------------------------------------------
// パスにおいて、直下の名前を返す
//
// "/aaa/bbb/ccc" -> "aaa" (rest="bbb/ccc")
std::string TlHdf5Utils::getChildName(const std::string& path,
                                      std::string* pRest) {
  std::size_t offset = 0;
  if ((path.length() > 0) && (path[0] == TlHdf5Utils::delim_)) {
    offset = 1;
  }

  std::string answer = "";
  std::string rest = "";
  const std::size_t next = path.find_first_of(TlHdf5Utils::delim_, offset);
  if (next == std::string::npos) {
    answer = path.substr(offset);
    rest = "";
  } else {
    answer = path.substr(offset, next - offset);
    rest = path.substr(next + 1);
  }

  if (pRest != NULL) {
    *pRest = rest;
  }

  // std::cerr << "getChildName: " << path << "> " << answer << ":" << rest <<
  // std::endl;

  return answer;
}

// ---------------------------------------------------------------------
std::string TlHdf5Utils::getDirName(const std::string& path) {
  std::string dirName, baseName;
  TlHdf5Utils::getDirBaseName(path, &dirName, &baseName);

  return dirName;
}

std::string TlHdf5Utils::getBaseName(const std::string& path) {
  std::string dirName, baseName;
  TlHdf5Utils::getDirBaseName(path, &dirName, &baseName);

  return baseName;
}

void TlHdf5Utils::getDirBaseName(const std::string& path, std::string* pDirName,
                                 std::string* pBaseName) {
  assert(pDirName != NULL);
  assert(pBaseName != NULL);

  const std::size_t length = path.length();
  std::size_t offset = length - 1;
  if ((path.length() > 0) && (path[length - 1] == TlHdf5Utils::delim_)) {
    offset = length - 2;
  }

  const std::size_t next = path.find_last_of(TlHdf5Utils::delim_, offset);
  if (next == std::string::npos) {
    *pBaseName = path;
    *pDirName = "";
  } else {
    *pBaseName = path.substr(next + 1);
    *pDirName = path.substr(0, next);
  }
}
