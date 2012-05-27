#ifndef DFSCF_PARALLEL_H
#define DFSCF_PARALLEL_H

#include "DfScf.h"
#include "Fl_GlobalinputX.h"

class DfScf_Parallel : public DfScf {
public:
    DfScf_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfScf_Parallel();

protected:
    virtual void logger(const std::string& str) const;

protected:
    virtual DfXCFunctional* getDfXCFunctional();

protected:
    virtual void saveParam() const;
    virtual void setScfParam();

protected:
    virtual void diffDensityMatrix();

    virtual DfDensityFittingObject* getDfDensityFittingObject();

    virtual void doXCIntegral();

    virtual void doThreeIndexIntegral();

    virtual DfJMatrix* getDfJMatrixObject();
    
    virtual DfKMatrix* getDfKMatrixObject();

    virtual DfFockMatrix* getDfFockMatrixObject();

    virtual DfTransFmatrix* getDfTransFmatrixObject(bool isExecDiis);

    virtual void doLevelShift();

    virtual DfDiagonal* getDfDiagonalObject();

    virtual DfTransatob* getDfTransatobObject();

    virtual DfDmatrix* getDfDmatrixObject();

    virtual DfTotalEnergy* getDfTotalEnergyObject();

    virtual DfPopulation* getDfPopulationObject();

    virtual DfSummary* getDfSummaryObject();

    virtual bool judge();
    virtual int checkConverge();

    virtual DfConverge* getDfConverge();

    virtual void cleanup();

    virtual bool checkMaxIteration();
};

#endif // DFSCF_PARALLEL_H
