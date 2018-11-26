#include <cmath>

#include "tl_dense_vector_object.h"
#include "TlUtils.h"
#include "tl_dense_vector_impl_object.h"
#include "tl_vector_utils.h"

#ifdef HAVE_HDF5
#include "TlHdf5Utils.h"
#endif  // HAVE_HDF5

// ---------------------------------------------------------------------------
// constructor & destructor
// ---------------------------------------------------------------------------
TlDenseVectorObject::TlDenseVectorObject()
    : pImpl_(NULL), log_(TlLogging::getInstance()) {}

TlDenseVectorObject::~TlDenseVectorObject() {}

// ---------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------
TlDenseVectorObject::size_type TlDenseVectorObject::getSize() const {
  return this->pImpl_->getSize();
}

void TlDenseVectorObject::resize(const TlDenseVectorObject::size_type newSize) {
  this->pImpl_->resize(newSize);
}

double TlDenseVectorObject::get(const index_type index) const {
  return this->pImpl_->get(index);
}

void TlDenseVectorObject::set(const index_type index, const double value) {
  this->pImpl_->set(index, value);
}

void TlDenseVectorObject::add(const index_type index, const double value) {
  this->pImpl_->add(index, value);
}

void TlDenseVectorObject::mul(const index_type index, const double value) {
  this->pImpl_->mul(index, value);
}

// ---------------------------------------------------------------------------
// operations
// ---------------------------------------------------------------------------
double TlDenseVectorObject::getMaxAbsoluteElement(
    TlDenseVectorObject::index_type* index) const {
  double maxVal = 0.0;
  TlDenseVectorObject::index_type argmax = -1;
  const TlDenseVectorObject::index_type size = this->getSize();
  for (TlDenseVectorObject::index_type i = 0; i < size; ++i) {
    const double v = std::fabs(this->get(i));
    if (maxVal < v) {
      maxVal = v;
      argmax = i;
    }
  }

  if (index != NULL) {
    *index = argmax;
  }

  return maxVal;
}

double TlDenseVectorObject::sum() const { return this->pImpl_->sum(); }

double TlDenseVectorObject::norm() const { return this->pImpl_->norm(); }

double TlDenseVectorObject::norm2() const { return this->pImpl_->norm2(); }

TlDenseVectorObject::index_type TlDenseVectorObject::argmax(
    const TlDenseVectorObject::index_type& begin,
    const TlDenseVectorObject::index_type& end) const {
  return this->pImpl_->argmax(begin, end);
}

void TlDenseVectorObject::sortByGreater() {
  return this->pImpl_->sortByGreater();
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
bool TlDenseVectorObject::load(const std::string& filePath) {
  bool answer = false;

  TlDenseVectorObject::index_type size = 0;
  const TlVectorUtils::FileSize headerSize =
      TlVectorUtils::getHeaderInfo(filePath, &size);
  if (headerSize > 0) {
    this->resize(size);

    std::fstream fs;
    fs.open(filePath.c_str(), std::ios::in | std::ios::binary);
    if (!fs.fail()) {
      fs.seekg(headerSize);

      double v;
      for (TlDenseVectorObject::index_type i = 0; i < size; ++i) {
        fs.read(reinterpret_cast<char*>(&v), sizeof(double));
        this->set(i, v);
      }
      answer = true;
    }

    fs.close();
  } else {
    this->log_.critical(
        TlUtils::format("load failed.: %d@%s", __FILE__, __LINE__));
    answer = false;
  }

  return answer;
}

bool TlDenseVectorObject::save(const std::string& filePath) const {
  bool answer = false;

  std::fstream fs;
  fs.open(filePath.c_str(), std::ios::out | std::ios::binary);

  if (!fs.fail()) {
    const TlDenseVectorObject::size_type size = this->getSize();
    fs.write(reinterpret_cast<const char*>(&size),
             sizeof(TlDenseVectorObject::size_type));
    for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
      const double v = this->get(i);
      fs.write(reinterpret_cast<const char*>(&v), sizeof(double));
    }
    answer = true;
  }

  fs.close();

  return answer;
}

bool TlDenseVectorObject::loadText(const std::string& filePath) {
  bool answer = false;

  std::ifstream ifs;
  ifs.open(filePath.c_str(), std::ios::in);
  if (!ifs.fail()) {
    // load contents
    std::string line = "";
    std::getline(ifs, line);  // read 1st line

    if (line == "TEXT") {
      std::string tmp = "";
      ifs >> tmp;
      const int size = std::atoi(tmp.c_str());
      ifs >> tmp;  // equal to 'size'
      ifs >> tmp;  // equal to '0'
      this->resize(size);
      for (int i = 0; i < size; ++i) {
        ifs >> tmp;
        const double v = std::atof(tmp.c_str());
        this->set(i, v);
      }
      answer = true;
    } else {
      this->log_.critical(TlUtils::format(
          "illegal format: %s (%d@%s)", filePath.c_str(), __LINE__, __FILE__));
    }

  } else {
    this->log_.critical(TlUtils::format("could not open file: %s (%d@%s)",
                                        filePath.c_str(), __LINE__, __FILE__));
  }
  ifs.close();

  return answer;
}

void TlDenseVectorObject::outputText(std::ostream& os) const {
  const TlDenseVectorObject::size_type nSize = this->getSize();

  os << "TEXT\n";
  os << nSize << "\n";
  os << nSize << "\n";
  os << "0\n";

  for (TlDenseVectorObject::size_type j = 0; j < nSize; j += 10) {
    for (TlDenseVectorObject::size_type i = j; ((i < j + 10) && (i < nSize));
         ++i) {
      os << TlUtils::format("  %10.4lf", this->get(i));
    }
    os << std::endl;
  }
  os << std::endl;
}

#ifdef HAVE_HDF5
bool TlDenseVectorObject::saveHdf5(const std::string& filepath,
                                   const std::string& h5path) const {
  TlHdf5Utils h5(filepath);

  const TlDenseVectorObject::size_type size = this->getSize();
  std::vector<double> buf(size);
  for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
    buf[i] = this->get(i);
  }
  h5.write(h5path, &(buf[0]), size);
  h5.setAttr(h5path, "size", size);

  return true;
}

bool TlDenseVectorObject::loadHdf5(const std::string& filepath,
                                   const std::string& h5path) {
  TlHdf5Utils h5(filepath);

  TlDenseVectorObject::size_type size;
  h5.getAttr(h5path, "size", &size);

  this->resize(size);

  std::vector<double> buf(size);
  h5.get(h5path, &(buf[0]), size);
  for (TlDenseVectorObject::size_type i = 0; i < size; ++i) {
    this->set(i, buf[i]);
  }

  return true;
}

#endif  // HAVE_HDF5

// -----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream, const TlDenseVectorObject& mat) {
  const TlDenseVectorObject::size_type size = mat.getSize();
  for (TlDenseVectorObject::size_type ord = 0; ord < size; ord += 10) {
    stream << "\n";
    for (TlDenseVectorObject::size_type j = ord; (j < ord + 10) && (j < size);
         ++j) {
      stream << TlUtils::format("   %5d th", j + 1);
    }
    stream << "\n";

    for (TlDenseVectorObject::size_type j = ord; (j < ord + 10) && (j < size);
         ++j) {
      stream << "-----------";
    }
    stream << "----\n\n";

    for (TlDenseVectorObject::size_type j = ord; (j < ord + 10) && (j < size);
         ++j) {
      stream << TlUtils::format(" %10.6lf", mat.get(j));
    }
    stream << "\n\n";
  }

  return stream;
}
