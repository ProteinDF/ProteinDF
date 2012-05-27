#ifndef PROTEINDF_H
#define PROTEINDF_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <string>
#include "TlSerializeData.h"
#include "TlLogging.h"

class DfIntegrals;
class DfInitialGuess;
class DfForce;
class DfScf;

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
    virtual void saveParam() const;
    virtual void loadParam(const std::string& requestFilePath = "");

protected:
    void setupGlobalCondition();
    virtual void setupGlobalCondition_extra();
    virtual void manageMemory();
    
protected:
    virtual void startlogo();
    void startlogo(const std::string& version,
                   const std::string& info = "");

    virtual void inputData();
    void stepCreate();
    void stepIntegral();
    virtual void stepGuess();
    virtual DfScf* createDfScfInstance();
    void stepScf();
    virtual void stepForce();

    virtual void endlogo();
    void endlogo(const std::string& reports);
    
protected:
    virtual DfIntegrals* getDfIntegralsObject();
    virtual DfInitialGuess* getDfInitialGuessObject();
    virtual DfForce* getDfForceObject();

protected:
    virtual void stepStartTitle(const std::string& stepName);
    virtual void stepEndTitle();

protected:
    TlLogging& log_;
    TlSerializeData pdfParam_;

    bool showCacheReport_;
};

#endif // PROTEINDF_H

