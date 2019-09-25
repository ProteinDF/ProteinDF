// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef PROTEINDF_H
#define PROTEINDF_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <string>
#include "TlLogging.h"
#include "TlSerializeData.h"

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
    void startlogo(const std::string& version, const std::string& info = "");
    void checkEnvironment();

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

#endif  // PROTEINDF_H
