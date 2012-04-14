#ifndef DFINTEGRALS_PARALLEL_H
#define DFINTEGRALS_PARALLEL_H

#include "DfIntegrals.h"

class DfHpq;
class DfOverlap;

class DfIntegrals_Parallel : public DfIntegrals {
public:
    explicit DfIntegrals_Parallel(TlSerializeData* param = NULL);
    virtual ~DfIntegrals_Parallel();

protected:
    virtual void logger(const std::string& str) const;

protected:
    virtual DfCD* getDfCDObject();
    virtual DfXMatrix* getDfXMatrixObject();
    virtual DfInvMatrix* getDfInvMatrixObject();
    virtual DfGenerateGrid* getDfGenerateGridObject();

protected:
    virtual void createHpqMatrix();
    void createHpqMatrix_LAPACK();
    void createHpqMatrix_ScaLAPACK();

    virtual void createOverlapMatrix();
    void createOverlapMatrix_LAPACK();
    void createOverlapMatrix_ScaLAPACK();

    virtual void createERIMatrix();
    void createERIMatrix_LAPACK();
    void createERIMatrix_ScaLAPACK();

protected:
    virtual void outputStartTitle(const std::string& stepName, const char lineChar = '-');
    virtual void outputEndTitle(const std::string& stepName ="", const char lineChar = '-');
};

#endif // DFINTEGRALS_PARALLEL_H
