#ifndef TLHDF5UTILS_H
#define TLHDF5UTILS_H

#include <string>
#include <vector>

#include <H5Cpp.h>

class TlHdf5Utils {
 public:
  typedef hsize_t size_type;

 public:
  TlHdf5Utils(const std::string& path);
  ~TlHdf5Utils();

 public:
  void flush() const;

 public:
  /// 指定されたパスのグループを作成する
  void createGroup(const std::string& path);

 public:
  /// 指定されたパスにグループがある場合はtrueを返す
  bool hasGroup(const std::string& path);

 public:
  /// 指定されたパスにDataSetがある場合はtrueを返す
  bool hasDataSet(const std::string& path);

  // -----------------------------------------------------------------
  // create dataset
  // -----------------------------------------------------------------
 public:
  /// create dataset which contains int array
  void createDataSet_int(const std::string& path, const std::size_t size);

  /// create dataset which contains double array
  void createDataSet_double(const std::string& path, const std::size_t size);

 protected:
  void createDataSet_common(const std::string& path,
                            const H5::PredType& predType,
                            const std::size_t size);

  // -----------------------------------------------------------------
  // write
  // -----------------------------------------------------------------
 public:
  void write(const std::string& name, const char value);
  void write(const std::string& name, const unsigned char value);
  void write(const std::string& name, const int value);
  void write(const std::string& name, const unsigned int value);
  void write(const std::string& name, const long value);
  void write(const std::string& name, const unsigned long value);
  void write(const std::string& name, const float value);
  void write(const std::string& name, const double value);
  void write(const std::string& name, const std::string& values);

  void write(const std::string& name, const int* pValues, std::size_t size);
  void write(const std::string& name, const unsigned int* pValues,
             std::size_t size);
  void write(const std::string& name, const long* pValues, std::size_t size);
  void write(const std::string& name, const unsigned long* pValues,
             std::size_t size);
  void write(const std::string& name, const float* pValues, std::size_t size);
  void write(const std::string& name, const double* pValues, std::size_t size);

  void write(const std::string& name, const std::vector<char>& values);
  void write(const std::string& name, const std::vector<unsigned char>& values);
  void write(const std::string& name, const std::vector<int>& values);
  void write(const std::string& name, const std::vector<unsigned int>& values);
  void write(const std::string& name, const std::vector<long>& values);
  void write(const std::string& name, const std::vector<unsigned long>& values);
  void write(const std::string& name, const std::vector<float>& values);
  void write(const std::string& name, const std::vector<double>& values);
  void write(const std::string& name, const std::vector<std::string>& values);

 protected:
  template <typename T>
  void write_common(const std::string& path, const T value,
                    const H5::PredType& predType);
  template <typename T>
  void write_common(const std::string& path, const T* pValues,
                    const H5::PredType& predType, const std::size_t size);
  template <typename T>
  void write_common(const std::string& path, const std::vector<T>& values,
                    const H5::PredType& predType);
  void prepareToWrite(const std::string& path);

  // -----------------------------------------------------------------
  // get data from dataset
  // -----------------------------------------------------------------
 public:
  void get(const std::string& path, char* pOut) const;
  void get(const std::string& path, unsigned char* pOut) const;
  void get(const std::string& path, int* pOut) const;
  void get(const std::string& path, unsigned int* pOut) const;
  void get(const std::string& path, long* pOut) const;
  void get(const std::string& path, unsigned long* pOut) const;
  void get(const std::string& path, float* pOut) const;
  void get(const std::string& path, double* pOut) const;
  void get(const std::string& path, std::string* pOut) const;

  void get(const std::string& path, double* pOut, const std::size_t size) const;

  void get(const std::string& path, std::vector<char>* pOut) const;
  void get(const std::string& path, std::vector<unsigned char>* pOut) const;
  void get(const std::string& path, std::vector<int>* pOut) const;
  void get(const std::string& path, std::vector<unsigned int>* pOut) const;
  void get(const std::string& path, std::vector<long>* pOut) const;
  void get(const std::string& path, std::vector<unsigned long>* pOut) const;
  void get(const std::string& path, std::vector<float>* pOut) const;
  void get(const std::string& path, std::vector<double>* pOut) const;
  void get(const std::string& path, std::vector<std::string>* pOut) const;

 protected:
  template <typename T>
  void get(const std::string& path, const H5::PredType& predType,
           T* value) const;
  template <typename T>
  void get_common(const std::string& path, const H5::PredType& predType,
                  T* pOut, const std::size_t size) const;
  template <typename T>
  void get(const std::string& path, const H5::PredType& predType,
           std::vector<T>* pOut) const;

 public:
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const std::vector<char>& data);
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const std::vector<unsigned char>& data);
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const std::vector<int>& data);
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const std::vector<unsigned int>& data);
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const std::vector<long>& data);
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const std::vector<unsigned long>& data);
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const std::vector<float>& data);
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const std::vector<double>& data);

  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           std::vector<char>* pOut);
  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           std::vector<unsigned char>* pOut);
  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           std::vector<int>* pOut);
  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           std::vector<unsigned int>* pOut);
  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           std::vector<long>* pOut);
  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           std::vector<unsigned long>* pOut);
  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           std::vector<float>* pOut);
  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           std::vector<double>* pOut);

 protected:
  template <typename T>
  void setSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const H5::PredType& predType,
                           const std::vector<T>& data);
  template <typename T>
  void getSelectedElements(const std::string& path,
                           const std::vector<hsize_t>& coord,
                           const H5::PredType& predType, std::vector<T>* pOut);

  // -----------------------------------------------------------------
  // attribute
  // -----------------------------------------------------------------
 public:
  void setAttr(const std::string& path, const std::string& attrName,
               const char value);
  void setAttr(const std::string& path, const std::string& attrName,
               const unsigned char value);
  void setAttr(const std::string& path, const std::string& attrName,
               const int value);
  void setAttr(const std::string& path, const std::string& attrName,
               const unsigned int value);
  void setAttr(const std::string& path, const std::string& attrName,
               const long value);
  void setAttr(const std::string& path, const std::string& attrName,
               const unsigned long value);
  void setAttr(const std::string& path, const std::string& attrName,
               const float value);
  void setAttr(const std::string& path, const std::string& attrName,
               const double value);
  void getAttr(const std::string& path, const std::string& attrName,
               char* pValue);
  void getAttr(const std::string& path, const std::string& attrName,
               unsigned char* pValue);
  void getAttr(const std::string& path, const std::string& attrName,
               int* pValue);
  void getAttr(const std::string& path, const std::string& attrName,
               unsigned int* pValue);
  void getAttr(const std::string& path, const std::string& attrName,
               long* pValue);
  void getAttr(const std::string& path, const std::string& attrName,
               unsigned long* pValue);
  void getAttr(const std::string& path, const std::string& attrName,
               float* pValue);
  void getAttr(const std::string& path, const std::string& attrName,
               double* pValue);

 protected:
  template <typename T>
  void setAttr(const std::string& path, const std::string& attrName,
               const H5::PredType& predType, const T pValue);
  template <typename T>
  void getAttr(const std::string& path, const std::string& attrName,
               const H5::PredType& predType, T* pValue);

  // -----------------------------------------------------------------
 protected:
  bool hasGroup(const H5::CommonFG* pParent, const std::string& name);

 protected:
  bool hasDataSet(const H5::CommonFG* pParent, const std::string& name);

 protected:
  bool hasObject(const H5::CommonFG* pParent, const std::string& name,
                 const H5O_type_t requestType);
  bool hasChildObject(const H5::CommonFG* pParent, const std::string& name,
                      const H5O_type_t requestType);

  H5::Group getGroup(const H5::CommonFG* pParent, const std::string& name);

  void removeObject(H5::CommonFG* pParent, const std::string& name);

  static std::string getChildName(const std::string& path, std::string* pRest);

 private:
  static std::string getDirName(const std::string& path);
  static std::string getBaseName(const std::string& path);
  static void getDirBaseName(const std::string& path, std::string* pDirName,
                             std::string* pBaseName);

 private:
  static const char delim_;
  H5::H5File file_;
};

#endif  // TLHDF5UTILS_H
