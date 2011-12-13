#ifndef PROTEINDF_PARALLEL_H
#define PROTEINDF_PARALLEL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "ProteinDF.h"

class ProteinDF_Parallel : public ProteinDF {
public:
    ProteinDF_Parallel();
    virtual ~ProteinDF_Parallel();

protected:
    virtual void loadParam(const std::string& requestFilePath = "");
    virtual void saveParam() const;
    
    virtual void setupGlobalCondition();

protected:
    virtual DfIntegrals* getDfIntegralsObject();
    virtual DfInitialGuess* getDfInitialGuessObject();
    virtual DfForce* getDfForceObject();
    
    virtual void inputData();

    virtual DfScf* createDfScfInstance();

    virtual void startlogo();
    virtual void endlogo();
    virtual void stepStartTitle(const std::string& stepName);
    virtual void stepEndTitle();
};

#endif // PROTEINDF_H

