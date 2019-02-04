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

#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#include <cstdlib>
#include <iostream>
#include "CnError.h"
#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlLogging.h"
#include "TlSystem.h"
#include "TlTime.h"
#include "TlUtils.h"

const std::string DfObject::m_sWorkDirPath = "fl_Work";
const std::string DfObject::m_sTempDirPath = "fl_Work";

int DfObject::objectCount_ = 0;
TlMatrixCache DfObject::matrixCache_;
int DfObject::rank_ = 0;

const std::string DfObject::m_sRunTypeSuffix[DfObject::RUN_MAXINDEX] = {
    "undefined",   "rks",       "uks_alpha",  "uks_beta", "roks",
    "roks_closed", "roks_open", "roks_alpha", "roks_beta"};

DfObject::DfObject(TlSerializeData* pPdfParam)
    : pPdfParam_(pPdfParam), log_(TlLogging::getInstance()) {
  this->setParam(*pPdfParam);
  ++DfObject::objectCount_;
}

DfObject::~DfObject() {
  --DfObject::objectCount_;
  if (DfObject::objectCount_ == 0) {
    // DfObjectの個数が0になったので、行列cacheを空にする
    this->matrixCache_.flush();
  }
}

void DfObject::setParam(const TlSerializeData& data) {
  this->numOfThreads_ = 1;
#ifdef _OPENMP
  this->numOfThreads_ = omp_get_max_threads();
#endif  // _OPENMP

  // computational resource
  // this->procMaxMemSize_ = 1024UL * 1024UL * 1024UL;
  // {
  //   const std::string memSizeStr =
  //       TlUtils::toUpper(data["memory_size"].getStr());
  //   if (memSizeStr.empty() == false) {
  //     std::size_t value = std::atol(memSizeStr.c_str());
  //     if (memSizeStr.rfind("MB") != std::string::npos) {
  //       value *= (1024UL * 1024UL);
  //     } else if (memSizeStr.rfind("GB") != std::string::npos) {
  //       value *= (1024UL * 1024UL * 1024UL);
  //     }
  //
  //     this->procMaxMemSize_ = value;
  //   }
  // }

  // this->isEnableMmap_ = data["use_mapfile"].getBoolean();
  // this->isWorkOnDisk_ = data["work_on_disk"].getBoolean();
  this->localTempDir_ = data["local_temp_dir"].getStr();
  if (this->localTempDir_ == "") {
    this->localTempDir_ = "/tmp/";
  }

  this->isRestart_ = false;
  if (TlUtils::toUpper(data["restart"].getStr()) == "YES") {
    this->isRestart_ = true;
  }

  this->isUseNewEngine_ = data["new_engine"].getBoolean();

  // SCF type
  const std::string methodType = TlUtils::toUpper(data["method"].getStr());
  if ((methodType == "RKS") || (methodType == "NSP")) {
    this->m_nMethodType = METHOD_RKS;
  } else if ((methodType == "UKS") || (methodType == "SP")) {
    this->m_nMethodType = METHOD_UKS;
  } else if (methodType == "ROKS") {
    this->m_nMethodType = METHOD_ROKS;
  } else {
    this->log_.warn("no method is specified. use rks method.");
    this->m_nMethodType = METHOD_RKS;
  }

  // model
  this->m_nNumOfAtoms = data["num_of_atoms"].getInt();
  {
    const Fl_Geometry geom(data["coordinates"]);
    this->m_nNumOfDummyAtoms = geom.getNumOfDummyAtoms();
  }
  // this->numOfRealAtoms_ = this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;

  this->m_nIteration = data["num_of_iterations"].getInt();
  this->m_nNumOfAOs = data["num_of_AOs"].getInt();
  this->m_nNumOfMOs = data["num_of_MOs"].getInt();
  this->m_nNumOfAux = data["num_of_auxCDs"].getInt();
  this->numOfAuxXC_ = data["num_of_auxXCs"].getInt();

  this->m_nNumOfElectrons = data["method/rks/electrons"].getInt();
  this->m_nNumOfAlphaElectrons = data["method/uks/alpha_electrons"].getInt();
  this->m_nNumOfBetaElectrons = data["method/uks/beta_electrons"].getInt();
  this->numOfClosedShellElectrons_ =
      data["method/roks/closed_electrons"].getInt();
  this->numOfOpenShellElectrons_ = data["method/roks/open_electrons"].getInt();

  // guess
  this->initialGuessType_ = GUESS_UNKNOWN;
  {
    const std::string guess = TlUtils::toUpper(data["guess"].getStr());
    if ((guess == "RHO") || (guess == "ATOM_RHO")) {
      this->initialGuessType_ = GUESS_RHO;
    } else if (guess == "FILE_RHO") {
      this->initialGuessType_ = GUESS_FILE_RHO;
    } else if ((guess == "LCAO") || (guess == "FILE_LCAO")) {
      this->initialGuessType_ = GUESS_LCAO;
    } else if ((guess == "DENSITY") || (guess == "DENSITY_MATRIX")) {
      this->initialGuessType_ = GUESS_DENSITY;
    } else if (guess == "HUCKEL") {
      this->initialGuessType_ = GUESS_HUCKEL;
    } else if (guess == "CORE") {
      this->initialGuessType_ = GUESS_CORE;
    } else if (guess == "HARRIS") {
      this->initialGuessType_ = GUESS_HARRIS;
    } else {
      if (guess.empty() == true) {
        std::cerr << "initial guess parameter is not configured." << std::endl;
      } else {
        std::cerr << "unknown initial guess parameter: " << guess << std::endl;
      }
    }
  }

  // calculaton properties ===================================================
  this->chargeExtrapolateNumber_ = data["charge_extrapolate_number"].getInt();
  // disk utilization(no == DIRECT, yes == DISK)
  this->m_bDiskUtilization =
      (TlUtils::toUpper(data["disk-utilization"].getStr()) == "YES") ? true
                                                                     : false;
  this->m_bMemorySave =
      (TlUtils::toUpper(data["scf-memory-saving"].getStr()) == "YES") ? true
                                                                      : false;

  this->isUpdateMethod_ = data["update_method"].getBoolean();

  // J
  {
    this->J_engine_ = J_ENGINE_RI_J;
    const std::string J_engine = TlUtils::toUpper(data["J_engine"].getStr());
    if (J_engine == "RI_J") {
      this->J_engine_ = J_ENGINE_RI_J;
    } else if (J_engine == "CD") {
      this->J_engine_ = J_ENGINE_CD;
    } else if (J_engine == "CONVENTIONAL") {
      this->J_engine_ = J_ENGINE_CONVENTIONAL;
    }
  }

  // K
  {
    this->K_engine_ = K_ENGINE_CONVENTIONAL;
    const std::string K_engine = TlUtils::toUpper(data["K_engine"].getStr());
    if (K_engine == "RI_K") {
      this->K_engine_ = K_ENGINE_RI_K;
    } else if (K_engine == "CD") {
      this->K_engine_ = K_ENGINE_CD;
    } else if (K_engine == "FASTCDK") {
      this->K_engine_ = K_ENGINE_FASTCDK;
    } else if (K_engine == "CONVENTIONAL") {
      this->K_engine_ = K_ENGINE_CONVENTIONAL;
    } else {
      this->log_.warn(
          TlUtils::format("unknown parameter: K_engine=%s", K_engine.c_str()));
      this->log_.warn("use conventional engine for K");
      this->K_engine_ = K_ENGINE_CONVENTIONAL;
    }
  }

  // XC functional
  {
    this->XC_engine_ = XC_ENGINE_GRID;
    const std::string XC_engine = TlUtils::toUpper(data["XC_engine"].getStr());
    if (XC_engine == "GRID") {
      this->XC_engine_ = XC_ENGINE_GRID;
    } else if (XC_engine == "GRIDFREE") {
      this->XC_engine_ = XC_ENGINE_GRIDFREE;
    } else if (XC_engine == "GRIDFREE_CD") {
      this->XC_engine_ = XC_ENGINE_GRIDFREE_CD;
    } else {
      this->log_.warn(TlUtils::format("unknown parameter: XC_engine=%s",
                                      XC_engine.c_str()));
      this->log_.warn("use grid engine for XC");
      this->XC_engine_ = XC_ENGINE_GRID;
    }
  }

  this->m_sXCFunctional = TlUtils::toUpper(data["xc_functional"].getStr());
  {
    const char nLastChar =
        this->m_sXCFunctional[this->m_sXCFunctional.length() - 1];
    this->m_bIsXCFitting = (nLastChar == '~') ? true : false;
  }
  {
    this->isDFT_ = true;
    if (this->m_sXCFunctional == "HF") {
      this->isDFT_ = false;
    }
  }
  this->m_bIsUpdateXC =
      (TlUtils::toUpper(data["xc-update"].getStr()) == "NO") ? false : true;
  this->isDedicatedBasisForGridFree_ =
      data["gridfree/dedicated_basis"].getBoolean();

  // Grimme empirical dispersion check
  {
    this->enableGrimmeDispersion_ = false;
    int len = this->m_sXCFunctional.size();
    if (len > 2) {
      const std::string suffix = this->m_sXCFunctional.substr(len - 2, 2);
      if (suffix == "-D") {
        this->enableGrimmeDispersion_ = true;
      }
    }
  }

  // RI-K
  this->isRI_K_ = data["RI_K"].getBoolean();

  // matrix operation
  this->linearAlgebraPackage_ = DfObject::LAP_LAPACK;
  this->updateLinearAlgebraPackageParam(
      data["linear_algebra_package"].getStr());

  this->m_bUsingSCALAPACK = false;
#ifdef HAVE_SCALAPACK
  this->m_bUsingSCALAPACK =
      (TlUtils::toUpper(data["linear_algebra_package"].getStr()) == "SCALAPACK")
          ? true
          : false;
#endif  // HAVE_SCALAPACK

  this->isSaveDistributedMatrixToLocalDisk_ = false;
#ifdef HAVE_SCALAPACK
  if (data.hasKey("save_distributed_matrix_to_local_disk") == true) {
    this->isSaveDistributedMatrixToLocalDisk_ =
        data["save_distributed_matrix_to_local_disk"].getBoolean();
  }
#endif  // HAVE_SCALAPACK

  this->localDiskPath_ = "/tmp";
  if (data.hasKey("local_disk_path") == true) {
    this->localDiskPath_ = data["local_disk_path"].getStr();
  }

  // for HPC ============================================================
  this->grainSize_ = 100;
  if (data.hasKey("omp_grain_size") == true) {
    this->grainSize_ = data["omp_grain_size"].getInt();
  }
  this->grainSize_ *= this->numOfThreads_;

  this->isMasterSlave_ = false;
  {
    const std::string parallelProcessingType =
        TlUtils::toUpper(data["parallel_processing_type"].getStr());
    if ((parallelProcessingType == "MASTER_SLAVE") ||
        (parallelProcessingType == "MASTER-SLAVE") ||
        (parallelProcessingType == "MS")) {
      this->isMasterSlave_ = true;
    }
  }

  // for DEBUG ===============================================================
  this->isFileWarning = data["debug/file_warning"].getBoolean();
  this->isSaveJMatrix_ = data["debug/save_J"].getBoolean();
  this->enableExperimentalCode_ = data["experimental_code"].getBoolean();

  // for memory ==============================================================
  // this->isUseCache_ = (*(this->pPdfParam_))["use_matrix_cache"].getBoolean();
  // if (this->isUseCache_ == true) {
  //   this->matrixCache_.setMaxMemSize(this->procMaxMemSize_);
  // } else {
  this->matrixCache_.setMaxMemSize(0);
  // }
  const bool isForceLoadingFromDisk =
      (*(this->pPdfParam_))["force_loading_from_disk"].getBoolean();
  this->matrixCache_.forceLoadingFromDisk(isForceLoadingFromDisk);

  // setup
  TlSerializeData& paramFileBaseName =
      (*(this->pPdfParam_))["control"]["file_base_name"];
  paramFileBaseName["Hpq_matrix"] = "Hpq.mat";
  paramFileBaseName["Hpq2_matrix"] = "Hpq2.mat";
  paramFileBaseName["Spq_matrix"] = "Spq.mat";
  paramFileBaseName["Sab_matrix"] = "Sab.mat";
  paramFileBaseName["Nalpha_vtr"] = "Nalpha.vtr";
  paramFileBaseName["Sab2_matrix"] = "Sab2.mat";
  paramFileBaseName["Sgd_matrix"] = "Sgd.mat";
  if (paramFileBaseName["SabInv_matrix"].getStr().empty() == true) {
    paramFileBaseName["SabInv_matrix"] = "Sabinv.mat";
  }
  if (paramFileBaseName["SgdInv_matrix"].getStr().empty() == true) {
    paramFileBaseName["SgdInv_matrix"] = "Sgdinv.mat";
  }
  if (paramFileBaseName["I2PQ_vtr"].getStr().empty() == true) {
    paramFileBaseName["I2PQ_vtr"] = "I2PQ.vtr";
  }
  if (paramFileBaseName["I2PR_vtr"].getStr().empty() == true) {
    paramFileBaseName["I2PR_vtr"] = "I2PR.vtr";
  }
  if (paramFileBaseName["I2PQ_XC_vtr"].getStr().empty() == true) {
    paramFileBaseName["I2PQ_XC_vtr"] = "I2PQ.XC.vtr";
  }
  if (paramFileBaseName["Ljk_matrix"].getStr().empty() == true) {
    paramFileBaseName["Ljk_matrix"] = "Ljk.mat";
  }
  if (paramFileBaseName["Lk_matrix"].getStr().empty() == true) {
    paramFileBaseName["Lk_matrix"] = "Lk.mat";
  }
  if (paramFileBaseName["Lxc_matrix"].getStr().empty() == true) {
    paramFileBaseName["Lxc_matrix"] = "Lxc.mat";
  }
  if (paramFileBaseName["X_matrix"].getStr().empty() == true) {
    paramFileBaseName["X_matrix"] = "X.mat";
  }
  if (paramFileBaseName["Xinv_matrix"].getStr().empty() == true) {
    paramFileBaseName["Xinv_matrix"] = "Xinv.mat";
  }
  if (paramFileBaseName["XEigval_vtr"].getStr().empty() == true) {
    paramFileBaseName["XEigval_vtr"] = "XEigval.vtr";
  }
  if (paramFileBaseName["diff_density_matrix"].getStr().empty() == true) {
    paramFileBaseName["diff_density_matrix"] = "dP.%s.mat";
  }
  if (paramFileBaseName["spin_density_matrix"].getStr().empty() == true) {
    paramFileBaseName["spin_density_matrix"] = "spin_density.%s.mat";
  }

  if (paramFileBaseName["occupation_vtr"].getStr().empty() == true) {
    paramFileBaseName["occupation_vtr"] = "occupation.%s.vtr";
  }
  if (paramFileBaseName["Ppq_matrix"].getStr().empty() == true) {
    paramFileBaseName["Ppq_matrix"] = "Ppq.%s.mat";
  }
  paramFileBaseName["P1pq_matrix"] = "P1pq.%s.mat";
  paramFileBaseName["P2pq_matrix"] = "P2pq.%s.mat";
  paramFileBaseName["HFx_matrix"] = "HFx.%s.mat";
  paramFileBaseName["Fpq_matrix"] = "Fpq.%s.mat";
  paramFileBaseName["Fprime_matrix"] = "Fprime.%s.mat";
  paramFileBaseName["Fxc_matrix"] = "Fxc.%s.mat";
  paramFileBaseName["Exc_matrix"] = "Exc.%s.mat";
  paramFileBaseName["FxcPure_matrix"] = "FxcPure.%s.mat";
  paramFileBaseName["J_matrix"] = "J.%s.mat";
  paramFileBaseName["C_matrix"] = "C.%s.mat";
  paramFileBaseName["Cprime_matrix"] = "Cprime.%s.mat";
  paramFileBaseName["grid_matrix"] = "grid.mat";
  paramFileBaseName["T_alpha"] = "T_alpha.%s.vtr";

  // for GridFree
  if (paramFileBaseName["GF_S_matrix"].getStr().empty() == true) {
    paramFileBaseName["GF_S_matrix"] = "GF_S.mat";
  }
  if (paramFileBaseName["GF_Stilde_matrix"].getStr().empty() == true) {
    paramFileBaseName["GF_Stilde_matrix"] = "GF_S.mat";
  }
  if (paramFileBaseName["GF_omega_matrix"].getStr().empty() == true) {
    paramFileBaseName["GF_omega_matrix"] = "GF_omega.mat";
  }
  if (paramFileBaseName["GF_V_matrix"].getStr().empty() == true) {
    paramFileBaseName["GF_V_matrix"] = "GF_V.mat";
  }
  if (paramFileBaseName["GF_VEigval_vtr"].getStr().empty() == true) {
    paramFileBaseName["GF_VEigval_vtr"] = "GF_VEigval.vtr";
  }

  if (paramFileBaseName["dipoleVelocityIntegrals_x"].getStr().empty() == true) {
    paramFileBaseName["dipoleVelocityIntegrals_x"] = "dSdx.mat";
  }
  if (paramFileBaseName["dipoleVelocityIntegrals_y"].getStr().empty() == true) {
    paramFileBaseName["dipoleVelocityIntegrals_y"] = "dSdy.mat";
  }
  if (paramFileBaseName["dipoleVelocityIntegrals_z"].getStr().empty() == true) {
    paramFileBaseName["dipoleVelocityIntegrals_z"] = "dSdz.mat";
  }

  // vectors
  if (paramFileBaseName["rho_vector"].getStr().empty() == true) {
    paramFileBaseName["rho_vector"] = "rho.%s.vtr";
  }
  if (paramFileBaseName["myu_vector"].getStr().empty() == true) {
    paramFileBaseName["myu_vector"] = "myu.%s.vtr";
  }
  if (paramFileBaseName["nyu_vector"].getStr().empty() == true) {
    paramFileBaseName["nyu_vector"] = "nyu.%s.vtr";
  }
  if (paramFileBaseName["eigenvalues"].getStr().empty() == true) {
    paramFileBaseName["eigenvalues"] = "eigenvalues.%s.vtr";
  }

  // for Population
  if (paramFileBaseName["pop/gross/orb"].getStr().empty() == true) {
    paramFileBaseName["pop/gross/orb"] = "pop.gross.orb.%s.vtr";
  }
  if (paramFileBaseName["pop/gross/atom"].getStr().empty() == true) {
    paramFileBaseName["pop/gross/atom"] = "pop.gross.atom.%s.vtr";
  }
  if (paramFileBaseName["pop/mulliken"].getStr().empty() == true) {
    paramFileBaseName["pop/mulliken"] = "pop.mulliken.%s.vtr";
  }

  // for lo
  paramFileBaseName["Clo_matrix"] = "Clo.%s.mat";
}

void DfObject::updateLinearAlgebraPackageParam(const std::string& keyword) {
  const std::string linearAlgebraPackage = TlUtils::toUpper(keyword);
#ifdef HAVE_EIGEN
  if (linearAlgebraPackage == "EIGEN") {
    this->linearAlgebraPackage_ = DfObject::LAP_EIGEN;
  }
#endif  // HAVE_EIGEN
#ifdef HAVE_LAPACK
  if (linearAlgebraPackage == "LAPACK") {
    this->linearAlgebraPackage_ = DfObject::LAP_LAPACK;
  }
#endif  // HAVE_LAPACK
#ifdef HAVE_VIENNACL
  if (linearAlgebraPackage == "VIENNACL") {
    this->linearAlgebraPackage_ = DfObject::LAP_VIENNACL;
  }
#endif  // HAVE_VIENNACL
}

void DfObject::logger(const std::string& str) const {
  // TlLogX& log = TlLogX::getInstance();
  // log << str;
  this->log_.info(str);
}

void DfObject::loggerTime(const std::string& str) const {
  // std::string out = str;
  // int size = out.size();
  // if (size > 0) {
  //     if (out[size -1] == '\n') {
  //         out.erase(size -1, 1);
  //     }

  //     const std::string timeStr = "[" + TlTime::getNow() + "]";
  //     TlUtils::pad(out, (72 - timeStr.length()), ' ');
  //     out += (timeStr + "\n");
  //     this->logger(out);
  // }
  this->log_.info(str);
}

void DfObject::loggerStartTitle(const std::string& stepName,
                                const char lineChar) const {
  // std::string line = "";
  // TlUtils::pad(line, 72, lineChar);

  // const std::string timeString = TlUtils::format("[%s %s]",
  //                                                TlTime::getNowDate().c_str(),
  //                                                TlTime::getNowTime().c_str());

  // std::string title = ">>>> " + stepName + " ";
  // TlUtils::pad(title, (72 - timeString.length()), ' ');
  // title += timeString;

  // // 出力
  // const std::string str = line + "\n" + title + "\n";

  std::string line = "";
  TlUtils::pad(line, 50, lineChar);
  const std::string str = line + "\n" + ">>>> " + stepName + "\n";
  this->logger(str);
}

void DfObject::loggerEndTitle(const std::string& stepName,
                              const char lineChar) const {
  // const std::string timeString = TlUtils::format("[%s %s]",
  //                                                TlTime::getNowDate().c_str(),
  //                                                TlTime::getNowTime().c_str());

  // std::string title = "<<<< " + stepName + " ";
  // TlUtils::pad(title, (72 - timeString.length()), ' ');
  // title += timeString;

  // // 出力
  // const std::string str = "\n" + title + "\n";

  const std::string str = "<<<<\n";
  this->logger(str);
}

// =====================================================================
std::string DfObject::makeFilePath(const std::string& baseFileName,
                                   const std::string& suffix) const {
  std::string base =
      (*(this->pPdfParam_))["control"]["file_base_name"][baseFileName].getStr();
  if (suffix.empty() != true) {
    base = TlUtils::format(base.c_str(), suffix.c_str());
  }

  std::string path;
  if (this->isSaveDistributedMatrixToLocalDisk_ == true) {
    path = this->localDiskPath_ + "/" + base + "." + TlUtils::xtos(this->rank_);
  } else {
    path = DfObject::m_sWorkDirPath + "/" + base;
  }

  return path;
}

std::string DfObject::getHpqMatrixPath() {
  return this->makeFilePath("Hpq_matrix");
}

std::string DfObject::getHpq2MatrixPath() {
  return this->makeFilePath("Hpq2_matrix");
}

std::string DfObject::getSpqMatrixPath() {
  return this->makeFilePath("Spq_matrix");
}

std::string DfObject::getSabMatrixPath() {
  return this->makeFilePath("Sab_matrix");
}

std::string DfObject::getSab2MatrixPath() {
  return this->makeFilePath("Sab2_matrix");
}

std::string DfObject::getSgdMatrixPath() {
  return this->makeFilePath("Sgd_matrix");
}

std::string DfObject::getSabInvMatrixPath() {
  return this->makeFilePath("SabInv_matrix");
}

std::string DfObject::getSgdInvMatrixPath() {
  return this->makeFilePath("SgdInv_matrix");
}

std::string DfObject::getI2pqVtrPath() {
  return this->makeFilePath("I2PQ_vtr");
}

std::string DfObject::getI2prVtrPath() {
  return this->makeFilePath("I2PR_vtr");
}

std::string DfObject::getI2pqVtrXCPath() {
  return this->makeFilePath("I2PQ_XC_vtr");
}

std::string DfObject::getLjkMatrixPath() {
  return this->makeFilePath("Ljk_matrix");
}

std::string DfObject::getLkMatrixPath() {
  return this->makeFilePath("Lk_matrix");
}

std::string DfObject::getLxcMatrixPath() {
  return this->makeFilePath("Lxc_matrix");
}

std::string DfObject::getXMatrixPath() {
  return this->makeFilePath("X_matrix");
}

std::string DfObject::getXInvMatrixPath() {
  return this->makeFilePath("Xinv_matrix");
}

std::string DfObject::getXEigvalVtrPath() {
  return this->makeFilePath("XEigval_vtr");
}

std::string DfObject::getNalphaPath() {
  return this->makeFilePath("Nalpha_vtr");
}

std::string DfObject::getOccupationPath(const RUN_TYPE runType) {
  return this->makeFilePath("occupation_vtr",
                            DfObject::m_sRunTypeSuffix[runType]);
}

std::string DfObject::getGridDataFilePath() const {
  return DfObject::m_sWorkDirPath + "/grids.dat";
}

std::string DfObject::getGridMatrixPath(const int iteration) const {
  return this->makeFilePath("grid_matrix", TlUtils::xtos(iteration));
}

std::string DfObject::getDiffDensityMatrixPath(const RUN_TYPE runType,
                                               const int iteration) const {
  return this->makeFilePath(
      "diff_density_matrix",
      DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration));
}

std::string DfObject::getSpinDensityMatrixPath(const RUN_TYPE runType,
                                               const int iteration) const {
  return this->makeFilePath(
      "spin_density_matrix",
      DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration));
}

std::string DfObject::getPpqMatrixPath(const RUN_TYPE nRunType,
                                       const int nIteration) const {
  return this->makeFilePath("Ppq_matrix", DfObject::m_sRunTypeSuffix[nRunType] +
                                              TlUtils::xtos(nIteration));
}

std::string DfObject::getP1pqMatrixPath(const int iteration) {
  return this->makeFilePath("P1pq_matrix", TlUtils::xtos(iteration));
}

std::string DfObject::getP2pqMatrixPath(const int iteration) {
  return this->makeFilePath("P2pq_matrix", TlUtils::xtos(iteration));
}

std::string DfObject::getHFxMatrixPath(const RUN_TYPE nRunType,
                                       const int nIteration) {
  return this->makeFilePath("HFx_matrix", DfObject::m_sRunTypeSuffix[nRunType] +
                                              TlUtils::xtos(nIteration));
}

std::string DfObject::getFpqMatrixPath(const RUN_TYPE nRunType,
                                       const int nIteration) const {
  return this->makeFilePath("Fpq_matrix", DfObject::m_sRunTypeSuffix[nRunType] +
                                              TlUtils::xtos(nIteration));
}

std::string DfObject::getFprimeMatrixPath(const RUN_TYPE runType,
                                          const int iteration,
                                          const std::string& fragment) {
  std::string suffix = "";
  if (fragment.empty() != true) {
    suffix = fragment + ".";
  }
  suffix += DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration);
  return this->makeFilePath("Fprime_matrix", suffix);
}

std::string DfObject::getFxcMatrixPath(const RUN_TYPE nRunType,
                                       const int nIteration) {
  return this->makeFilePath("Fxc_matrix", DfObject::m_sRunTypeSuffix[nRunType] +
                                              TlUtils::xtos(nIteration));
}

std::string DfObject::getExcMatrixPath(const RUN_TYPE runType,
                                       const int iteration) {
  return this->makeFilePath("Exc_matrix", DfObject::m_sRunTypeSuffix[runType] +
                                              TlUtils::xtos(iteration));
}

std::string DfObject::getFxcPureMatrixPath(const RUN_TYPE nRunType,
                                           const int nIteration) {
  return this->makeFilePath(
      "FxcPure_matrix",
      DfObject::m_sRunTypeSuffix[nRunType] + TlUtils::xtos(nIteration));
}

std::string DfObject::getJMatrixPath(const int iteration) {
  return this->makeFilePath("J_matrix", TlUtils::xtos(iteration));
}

std::string DfObject::getCMatrixPath(const RUN_TYPE runType, int iteration,
                                     const std::string& fragment) const {
  std::string suffix = "";
  if (fragment.empty() != true) {
    suffix = fragment + ".";
  }
  suffix += DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration);
  return this->makeFilePath("C_matrix", suffix);
}

std::string DfObject::getCprimeMatrixPath(const RUN_TYPE runType, int iteration,
                                          const std::string& fragment) {
  std::string suffix = "";
  if (fragment.empty() != true) {
    suffix = fragment + ".";
  }
  suffix += DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration);
  return this->makeFilePath("Cprime_matrix", suffix);
}

std::string DfObject::getGfSMatrixPath() const {
  return this->makeFilePath("GF_S_matrix");
}

std::string DfObject::getGfStildeMatrixPath() const {
  return this->makeFilePath("GF_Stilde_matrix");
}

std::string DfObject::getGfOmegaMatrixPath() const {
  return this->makeFilePath("GF_omega_matrix");
}

std::string DfObject::getGfVMatrixPath() const {
  return this->makeFilePath("GF_V_matrix");
}

std::string DfObject::getGfVEigvalVtrPath() const {
  return this->makeFilePath("GF_VEigval_vtr");
}

std::string DfObject::getRhoPath(const RUN_TYPE nRunType,
                                 const int nIteration) const {
  return this->makeFilePath("rho_vector", DfObject::m_sRunTypeSuffix[nRunType] +
                                              TlUtils::xtos(nIteration));
}

std::string DfObject::getMyuPath(const RUN_TYPE nRunType,
                                 const int nIteration) const {
  return this->makeFilePath("myu_vector", DfObject::m_sRunTypeSuffix[nRunType] +
                                              TlUtils::xtos(nIteration));
}

std::string DfObject::getNyuPath(const RUN_TYPE nRunType,
                                 const int nIteration) const {
  return this->makeFilePath("myu_vector", DfObject::m_sRunTypeSuffix[nRunType] +
                                              TlUtils::xtos(nIteration));
}

std::string DfObject::getTalphaPath(const RUN_TYPE runType,
                                    const int iteration) const {
  return this->makeFilePath("T_alpha", DfObject::m_sRunTypeSuffix[runType] +
                                           TlUtils::xtos(iteration));
}

std::string DfObject::getEigenvaluesPath(const RUN_TYPE runType,
                                         const int iteration) const {
  return this->makeFilePath("eigenvalues", DfObject::m_sRunTypeSuffix[runType] +
                                               TlUtils::xtos(iteration));
}

std::string DfObject::getDipoleVelocityIntegralsXPath() const {
  return this->makeFilePath("dipoleVelocityIntegrals_x");
}

std::string DfObject::getDipoleVelocityIntegralsYPath() const {
  return this->makeFilePath("dipoleVelocityIntegrals_y");
}

std::string DfObject::getDipoleVelocityIntegralsZPath() const {
  return this->makeFilePath("dipoleVelocityIntegrals_z");
}

std::string DfObject::getPopGrossOrbPath(const RUN_TYPE runType,
                                         const int iteration) const {
  return this->makeFilePath(
      "pop/gross/orb",
      DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration));
}

std::string DfObject::getPopGrossAtomPath(const RUN_TYPE runType,
                                          const int iteration) const {
  return this->makeFilePath(
      "pop/gross/atom",
      DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration));
}

std::string DfObject::getPopMullikenPath(const RUN_TYPE runType,
                                         const int iteration) const {
  return this->makeFilePath(
      "pop/mulliken",
      DfObject::m_sRunTypeSuffix[runType] + TlUtils::xtos(iteration));
}

std::string DfObject::getCloMatrixPath(const RUN_TYPE runType,
                                       const int iteration) const {
  return this->makeFilePath("Clo_matrix", DfObject::m_sRunTypeSuffix[runType] +
                                              TlUtils::xtos(iteration));
}

// ----------------------------------------------------------------------------
// properties
// ----------------------------------------------------------------------------
int DfObject::getNumOfAtoms() const { return this->m_nNumOfAtoms; }

int DfObject::iteration() const { return this->m_nIteration; }
