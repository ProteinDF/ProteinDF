#ifndef TLMATRIX_MMAP_OBJECT_H
#define TLMATRIX_MMAP_OBJECT_H

#include <string>
#include <vector>

#include "tl_matrix_object.h"

class TlDenseMatrixMmapObject : public TlMatrixObject {
public:
    TlDenseMatrixMmapObject(const TlMatrixObject::MatrixType matrixType,
                            const std::string& filePath, const index_type row,
                            const index_type col);
    TlDenseMatrixMmapObject(const TlMatrixObject::MatrixType matrixType,
                            const std::string& filePath);
    virtual ~TlDenseMatrixMmapObject();

    // void resize(const index_type newRow, const index_type newCol);

public:
    virtual std::size_t getMemSize() const;

    virtual double get(index_type row, index_type col) const;
    virtual void set(index_type row, index_type col, double value);
    virtual void add(index_type row, index_type col, double value);

    virtual void setRowVector(const index_type row,
                              const std::vector<double>& v);
    virtual void setRowVector(const index_type row,
                              const std::valarray<double>& v);
    virtual void setColVector(const index_type col,
                              const std::vector<double>& v);
    virtual void setColVector(const index_type col,
                              const std::valarray<double>& v);

    virtual std::vector<double> getRowVector(const index_type row) const;
    virtual std::vector<double> getColVector(const index_type col) const;
    virtual std::size_t getRowVector(const index_type row, double* pBuf,
                                     std::size_t maxCount) const;
    virtual std::size_t getColVector(const index_type col, double* pBuf,
                                     std::size_t maxCount) const;

public:
    virtual bool saveText(const std::string& filePath) const;
    virtual void saveText(std::ostream& os) const;

    virtual bool saveCsv(const std::string& filePath) const;
    virtual void saveCsv(std::ostream& os) const;

protected:
    // virtual TlMatrixObject::size_type getNumOfElements() const;
    // virtual size_type getIndex(const index_type row,
    //                            const index_type col) const = 0;
    // virtual TlDenseMatrixMmapObject* copy(const std::string& path) const = 0;

private:
    virtual bool load(const std::string& path);
    virtual bool save(const std::string& path) const;

protected:
    void createNewFile();
    void openFile();

protected:
    void deleteMmap();

private:
    void getHeaderInfo();
    void newMmap();
    void syncMmap();

protected:
    std::string filePath_;

    char* mmapBegin_;
    double* dataBegin_;

    std::size_t headerSize_;
    std::size_t fileSize_;
};

#endif  // TLMATRIX_MMAP_OBJECT_H
