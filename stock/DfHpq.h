#ifndef DFHPQ_H
#define DFHPQ_H

#include <vector>

#include "DfObject.h"
#include "TlMatrixObject.h"
#include "TlSymmetricMatrix.h"
#include "TlPosition.h"
#include "TlOrbitalInfo.h"

/// 1電子ハミルトニアンの計算を行う
/// 得られた密度行列、交換相関ポテンシャルをもとに一電子部分（運動エネルギー項と原子核引力項）の微分計算、
/// およびPulay 補正項の計算を行う
class DfHpq : public DfObject {
public:
    DfHpq(TlSerializeData* pPdfParam);
    virtual ~DfHpq();

public:
    virtual void getHpq(TlSymmetricMatrix* pHpq,
                        TlSymmetricMatrix* pHpq2);

    std::vector<double> getESP(const TlMatrixObject* pPpq,
                               const std::vector<TlPosition>& grids);

protected:
    void getHpq_core(TlMatrixObject* pHpq, TlMatrixObject* pHpq2);

    std::vector<double> getESP_core(const std::vector<IJShellPair>& IJShellList,
                                    const TlMatrixObject* pPpq,
                                    const std::vector<TlPosition>& grids);

    void getMemory();

    virtual void makeTable();

    int auxSet();    // aux set
    void fmtSet();   // FMT set

    void getNuclearAttractionIntegrals(const std::size_t orbA,
                                       const std::size_t orbB,
                                       int nqA, int nqB, int npA, int npB,
                                       int izA, int izB, int icA, int icB,
                                       const TlPosition& posC,
                                       double* pOutput);

    void kinSS(const int npA, const int npB,
               const double* Za, const double* Zb, const double* Ca, const double* Cb,
               const double* A, const double* B, double* E);
    void kinPS(const int npA, const int npB,
               const double* Za, const double* Zb, const double* Ca, const double* Cb,
               const double* A, const double* B, double* E);
    void kinPP(const int npA, const int npB,
               const double* Za, const double* Zb, const double* Ca, const double* Cb,
               const double* A, const double* B, double* E);
    void kinDS(const int npA, const int npB,
               const double* Za, const double* Zb, const double* Ca, const double* Cb,
               const double* A, const double* B, double* E);
    void kinDP(const int npA, const int npB,
               const double* Za, const double* Zb, const double* Ca, const double* Cb,
               const double* A, const double* B, double* E);
    void kinDD(const int npA, const int npB,
               const double* Za, const double* Zb, const double* Ca, const double* Cb,
               const double* A, const double* B, double* E);

    // FMT calculations by backward recursive method
    void fmtRecursive(int,int,int);

    // main routine of one-electron Hamiltonian
    void hpqcalc(const std::vector<IJShellPair>& IJShellList,
                 TlMatrixObject* Hpq, TlMatrixObject* Hpq2 = NULL);

protected:
    virtual std::vector<IJShellPair> getQueue(int maxNumOfPairs, bool initialize = false);
    bool isCutoffUsingOverlap(std::size_t orbA, std::size_t orbB,
                              double cutvalue) const;

    void resetCounter();
    void countProgress(const int p_type, const int q_type);
    void countCutoff(const int p_type, const int q_type);
    int getProgress() const;

//     void testcode();
    
protected:
    const static double FPAI;        // 2*PAI**(5/2)
    const static double D33;         // 1/0.03
    const static double SQR3I;       // 1/sqrt(3)
    const static double SQR3I2;      // 2/sqrt(3)

    static const std::size_t TF_SIZE;
    static const std::size_t RMI_SIZE;
    static const std::size_t GA_SIZE;
    static const std::size_t EDAT_SIZE;
    static const std::size_t SDAT_SIZE;
    static const std::size_t ADAT_SIZE;

    TlOrbitalInfo* pOrbitalInfo_;

    int m_nNumOfTotalAtoms; // dummy atomを含む全原子数
    int m_nNumOfRealAtoms;  // dummy atomを除いた全原子数

    double cutvalue;    // cut value

    // for auxSet & fmtSet
    double* TF;
    double* RMI;
    double* GA;
    double* EDAT;
    double* SDAT;
    double** FDAT;
    double* ADAT;

    /// 計算ブロックのサイズ
    int blockSize_;

    int m_nChargeExtrapolateNumber;

    std::vector<std::vector<std::size_t> > shellList_; // shellList_[shell type][index] = AO index

    std::vector<std::size_t> totalCounter_;
    std::vector<std::size_t> progressCounter_;
    std::vector<std::size_t> cutoffCounter_;
};

#endif // DFHPQ_H


