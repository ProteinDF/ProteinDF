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

#ifndef DFSCF_H
#define DFSCF_H

#include <fstream>
#include <map>

#include "DfConverge.h"
#include "DfFockMatrix.h"
#include "DfGridFreeXC.h"
#include "DfObject.h"
#include "DfXCFunctional.h"

#include "TlTime.h"

class DfDensityFittingObject;
class DfJMatrix;
class DfKMatrix;
class DfTransFmatrix;
class DfDiagonal;
class DfTransatob;
class DfDmatrix;
class DfTotalEnergy;
class DfPopulation;
class DfSummary;

/// SCF 繰り返し計算を行うクラス
class DfScf : public DfObject {
   protected:
    enum ScfAccelerationType {
        SCF_ACCELERATION_SIMPLE = 0,
        SCF_ACCELERATION_ANDERSON,
        SCF_ACCELERATION_DIIS,
        SCF_ACCELERATION_MIX
    };

    enum DampObjectType {
        DAMP_NONE = 0,
        DAMP_DENSITY_MATRIX,
        DAMP_DENSITY,
        DAMP_FOCK
    };

    enum SYMMETRIC_MATRIX_NAME {
        DIFF_DENSITY_MATRIX,
        DIFF_DENSITY_MATRIX_ALPHA,
        DIFF_DENSITY_MATRIX_BETA,
    };

   public:
    DfScf(TlSerializeData* pPdfParam);
    virtual ~DfScf();

   public:
    int dfScfMain();

   protected:
    /// メンバ変数へ値を設定する
    virtual void setScfParam();

   protected:
    int execScfLoop();
    virtual void setScfRestartPoint(const std::string& str);

    // update法に伴い、差電子密度行列を計算する
    virtual void diffDensityMatrix();

    void doDensityFitting();
    virtual DfDensityFittingObject* getDfDensityFittingObject();

    virtual void doXCIntegral();

    // virtual void execScfLoop_XcEneFit();

    virtual void doThreeIndexIntegral();

    void buildXcMatrix();
    virtual DfGridFreeXC* getDfGridFreeXcObject();
    virtual DfXCFunctional* getDfXCFunctional();

    void buildJMatrix();
    virtual DfJMatrix* getDfJMatrixObject();

    void buildKMatrix();
    virtual DfKMatrix* getDfKMatrixObject();

    void buildFock();
    virtual DfFockMatrix* getDfFockMatrixObject();

    void transformFock();
    virtual DfTransFmatrix* getDfTransFmatrixObject(bool isExecDiis);

    virtual void doLevelShift();

    void diagonal();
    virtual DfDiagonal* getDfDiagonalObject();

    void execScfLoop_EndFock_TransC();
    virtual DfTransatob* getDfTransatobObject();

    virtual void calcDensityMatrix();
    virtual DfDmatrix* getDfDmatrixObject();

    // Total Energy -----------------------------------------------------------
    virtual void calcTotalEnergy();
    virtual void calcTotalRealEnergy();

    template <class DfTotalEnergyClass>
    void calcTotalEnergy_tmpl() {
        TlTime timer;
        this->loggerStartTitle("Total Energy");

        DfTotalEnergyClass dfTotalEnergy(this->pPdfParam_);
        dfTotalEnergy.run();

        this->loggerEndTitle();
        (*this->pPdfParam_)["stat"]["elapsed_time"]["total_energy"]
                           [this->m_nIteration] = timer.getElapseTime();
    };

    template <class DfTotalEnergyClass>
    void calcTotalRealEnergy_tmpl() {
        this->loggerStartTitle("Total Energy derived from point charges");

        DfTotalEnergyClass dfTotalEnergy(this->pPdfParam_);
        dfTotalEnergy.outputForDummy();

        this->loggerEndTitle();
    }

    // virtual DfTotalEnergy* getDfTotalEnergyObject();

    // ------------------------------------------------------------------------
    virtual DfPopulation* getDfPopulationObject();
    void calcPopulation();

    virtual DfSummary* getDfSummaryObject();
    void summarize();

    virtual bool judge();
    virtual bool checkConverge();

    void converge();
    virtual DfConverge* getDfConverge();

    virtual void cleanup();

    virtual bool checkMaxIteration();

   protected:
    virtual void saveParam() const;

   protected:
    DampObjectType m_nDampObject;
    ScfAccelerationType m_nScfAcceleration;

    int m_nConvergenceCounter;

    bool isUseNewEngine_;
};

// ----------------------------------------------------------------------------
// template
// ----------------------------------------------------------------------------

// calculate total energy
// template <class T>
// void DfScf::calcTotalEnergy_tmpl<T>() {
//     TlTime timer;
//     this->loggerStartTitle("Total Energy");

//     T dfTotalEnergy(this->pPdfParam_);
//     dfTotalEnergy.run();

//     this->loggerEndTitle();
//     (*this->pPdfParam_)["stat"]["elapsed_time"]["total_energy"]
//                        [this->m_nIteration] = timer.getElapseTime();
// }

// template <class DfTotalEnergyClass>
// void DfScf::calcTotalRealEnergy<DfTotalEnergyClass>() {
//     this->loggerStartTitle("Total Energy derived from point charges");

//     DfTotalEnergyClass dfTotalEnergy(this->pPdfParam_);
//     dfTotalEnergy.calculate_real_energy();

//     this->loggerEndTitle();
// }

#endif  // DFSCF_H
