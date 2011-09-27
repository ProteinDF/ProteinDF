#ifndef DFTOTALENERGY_PARALLEL_H
#define DFTOTALENERGY_PARALLEL_H

#include "DfTotalEnergy.h"
#include "TlDistributeVector.h"

class DfTotalEnergy_Parallel : public DfTotalEnergy {
public:
    DfTotalEnergy_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy_Parallel();

public:
    virtual void exec();

    /// ダミー原子を除いた時の全エネルギー計算を行う
    virtual void calculate_real_energy();

protected:
    void logger(const std::string& str) const;

protected:
    void exec_LAPACK();
    void exec_ScaLAPACK();

    // double calculate_one_electron_part(const TlDistributeSymmetricMatrix& D);

protected:
    virtual void output();

protected:
    bool m_bUseDistributeMatrix;
};

#endif // DFTOTALENERGY_PARALLEL_H
