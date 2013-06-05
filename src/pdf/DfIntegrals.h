#ifndef DFINTEGRALS_H
#define DFINTEGRALS_H

#include "DfObject.h"
#include "TlSerializeData.h"

class DfHpq;
class DfOverlap;
class DfEri;
class DfXMatrix;
class DfInvMatrix;
class DfCD;
class DfGridFreeXC;
class DfGenerateGrid;
class TlSymmetricMatrix;

/// 1 電子ハミルトニアン, 2 中心積分、1 中心積分(Na)を計算するクラス
class DfIntegrals : public DfObject {
protected:
    enum CalcState {
        Hpq = 1,
        Spq = 2,
        Sab2 = 4,
        Sgd = 8,
        Na = 16,
        Sab = 32,
        CD = 64,
        X = 128,
        INV = 256,
        CHOLESKY_VECTORS_XC = 512,
        GRID_FREE = 1024,
        GRID = 2048
    };

public:
    explicit DfIntegrals(TlSerializeData* param = NULL,
                         const std::string& saveParamPath = "");
    virtual ~DfIntegrals();

public:
    void main();

protected:
    virtual void saveParam();

protected:
    virtual void createHpqMatrix();
    virtual void createOverlapMatrix();
    virtual void createERIMatrix();
    void createCholeskyVectors();

    void createInverseMatrixes();
    void createXMatrix();
    void createCholeskyVectors_XC();
    void prepareGridFree();
    void createGrids();

protected:
    virtual DfXMatrix* getDfXMatrixObject();
    virtual DfInvMatrix* getDfInvMatrixObject();
    virtual DfCD* getDfCDObject();
    virtual DfGridFreeXC* getDfGridFreeXCObject();
    virtual DfGenerateGrid* getDfGenerateGridObject();

protected:
    virtual void saveInvSquareVMatrix(const TlSymmetricMatrix& v);

protected:
    virtual void outputStartTitle(const std::string& stepName, const char lineChar = '-');
    virtual void outputEndTitle(const std::string& stepName ="", const char lineChar = '-');

protected:
    std::string saveParamPath_;
};

#endif // DFINTEGRALS_H

