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
    
    virtual void setupGlobalCondition_extra();

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

