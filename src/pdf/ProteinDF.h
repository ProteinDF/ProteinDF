#ifndef PROTEINDF_H
#define PROTEINDF_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <string>
#include "TlSerializeData.h"

class DfIntegrals;
class DfInitialGuess;
class DfForce;

class ProteinDF {
public:
    ProteinDF();
    virtual ~ProteinDF();

public:
    void run();
    void restart(const std::string& restartParamFilePath);

protected:
    void exec();
    
protected:
    virtual void logger(const std::string& str) const;

    virtual void saveParam() const;
    virtual void loadParam(const std::string& requestFilePath = "");

protected:
    virtual void setupGlobalCondition();
    virtual void manageMemory();
    
protected:
    virtual void inputData();
    void stepCreate();
    void stepIntegral();
    virtual void stepGuess();
    virtual void stepScf();
    virtual void stepForce();
    
protected:
    virtual DfIntegrals* getDfIntegralsObject();
    virtual DfInitialGuess* getDfInitialGuessObject();
    virtual DfForce* getDfForceObject();

protected:
    virtual void startlogo();
    virtual void endlogo();
    virtual void stepStartTitle(const std::string& stepName);
    virtual void stepEndTitle();

protected:
    TlSerializeData pdfParam_;

    bool showCacheReport_;
};

#endif // PROTEINDF_H

