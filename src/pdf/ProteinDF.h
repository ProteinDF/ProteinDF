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
    void restart();

protected:
    void exec();
    
protected:
    virtual void logger(const std::string& str) const;
    virtual void save_Fl_Globalinput() const;

    virtual void saveParam() const;
    virtual void loadParam();

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
    std::string pdfParamPath_;
    TlSerializeData pdfParam_;

    bool showCacheReport_;
};

#endif // PROTEINDF_H

