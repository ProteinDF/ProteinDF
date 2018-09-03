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

#include "DfTaskCtrl.h"
#include "TlFmt.h"
#include "TlMath.h"
#include "TlOrbitalInfo.h"

DfTaskCtrl::DfTaskCtrl(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {
  TlOrbitalInfo orbitalInfo((*pPdfParam)["coordinates"],
                            (*pPdfParam)["basis_set"]);
  this->maxShellType_ = orbitalInfo.getMaxShellType();

  this->lengthScaleParameter_ = 3.0;
  if ((*pPdfParam)["length_scale_parameter"].getStr() != "") {
    this->lengthScaleParameter_ =
        (*pPdfParam)["length_scale_parameter"].getDouble();
  }

  this->cutoffThreshold_ = 1.0E-10;
  if ((*pPdfParam)["cut_value"].getStr().empty() != true) {
    this->cutoffThreshold_ = (*pPdfParam)["cut_value"].getDouble();
  }

  this->cutoffEpsilon_density_ = this->cutoffThreshold_;
  if ((*pPdfParam)["cutoff_density"].getStr().empty() != true) {
    this->cutoffEpsilon_density_ = (*pPdfParam)["cutoff_density"].getDouble();
  }

  this->cutoffEpsilon_distribution_ = this->cutoffThreshold_;
  if ((*pPdfParam)["cutoff_distribution"].getStr().empty() != true) {
    this->cutoffEpsilon_distribution_ =
        (*pPdfParam)["cutoff_distribution"].getDouble();
  }

  // this->cutoffEpsilon_primitive_ = this->cutoffThreshold_ * 0.01;
  // if ((*pPdfParam)["cutoff_primitive"].getStr().empty() != true) {
  //     this->cutoffEpsilon_primitive_ =
  //     (*pPdfParam)["cutoff_primitive"].getDouble();
  // }
}

DfTaskCtrl::~DfTaskCtrl() {}

// -----------------------------------------------------------------------------

void DfTaskCtrl::setCutoffThreshold(const double value) {
  this->cutoffThreshold_ = value;
}

double DfTaskCtrl::getCutoffThreshold() const { return this->cutoffThreshold_; }

void DfTaskCtrl::setCutoffEpsilon_distribution(const double value) {
  this->cutoffEpsilon_distribution_ = value;
}

double DfTaskCtrl::getCutoffEpsilon_distribution() const {
  return this->cutoffEpsilon_distribution_;
}

void DfTaskCtrl::setCutoffEpsilon_density(const double value) {
  this->cutoffEpsilon_density_ = value;
}

double DfTaskCtrl::getCutoffEpsilon_density() const {
  return this->cutoffEpsilon_density_;
}

// void DfTaskCtrl::setCutoffEpsilon_primitive(const double value)
// {
//     this->cutoffEpsilon_primitive_ = value;
// }

// double DfTaskCtrl::getCutoffEpsilon_primitive() const
// {
//     return this->cutoffEpsilon_primitive_;
// }

// -----------------------------------------------------------------------------

void DfTaskCtrl::clearCutoffStats(const int maxShellType) {
  const int numOfShellPairType = maxShellType * maxShellType;
  const int numOfShellQuartetType = numOfShellPairType * numOfShellPairType;

  this->cutoffAll_density_.clear();
  this->cutoffAlive_density_.clear();
  this->cutoffAll_density_.resize(numOfShellPairType, 0);
  this->cutoffAlive_density_.resize(numOfShellPairType, 0);

  this->cutoffAll_distribution_.clear();
  this->cutoffAlive_distribution_.clear();
  this->cutoffAll_distribution_.resize(numOfShellPairType, 0);
  this->cutoffAlive_distribution_.resize(numOfShellPairType, 0);

  this->cutoffAll_schwarz_.clear();
  this->cutoffAlive_schwarz_.clear();
  this->cutoffAll_schwarz_.resize(numOfShellQuartetType, 0);
  this->cutoffAlive_schwarz_.resize(numOfShellQuartetType, 0);
}

DfTaskCtrl::ShellArrayTable DfTaskCtrl::makeShellArrayTable(
    const TlOrbitalInfoObject& orbitalInfo) {
  const int maxShellType = orbitalInfo.getMaxShellType();
  ShellArrayTable shellArrayTable(maxShellType);
  const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();
  for (int i = 0; i < maxShellType; ++i) {
    shellArrayTable[i].reserve(maxShellIndex);
  }

  index_type shellIndex = 0;
  while (shellIndex < maxShellIndex) {
    // shellType: 0=s, 1=p, 2=d
    const int shellType = orbitalInfo.getShellType(shellIndex);
    const int steps = 2 * shellType + 1;

    shellArrayTable[shellType].push_back(shellIndex);
    shellIndex += steps;
  }

  return shellArrayTable;
}

// ----------------------------------------------------------------------------
// queue
// ----------------------------------------------------------------------------
bool DfTaskCtrl::getQueue(const std::size_t maxIndeces,
                          const std::size_t maxGrainSize,
                          std::vector<std::size_t>* pTasks,
                          bool initialize) {
  assert(pTasks != NULL);
  static std::size_t currentIndex = 0;

  if (initialize == true) {
    currentIndex = 0;
    
    return true;
  }

  pTasks->clear();
  pTasks->reserve(maxGrainSize);
  std::size_t numOfTasks = 0;
  for (; currentIndex < maxIndeces; ++currentIndex) {
    pTasks->push_back(currentIndex);
    ++numOfTasks;

    if (numOfTasks > maxGrainSize) {
      break;
    }
  }

  return (pTasks->empty() != true);
}

bool DfTaskCtrl::getQueue(const TlOrbitalInfoObject& orbitalInfo,
                          const int maxGrainSize, std::vector<Task>* pTaskList,
                          bool initialize) {
  assert(pTaskList != NULL);

  const int maxShellType = orbitalInfo.getMaxShellType();
  static ShellArrayTable shellArrayTable;
  static int shellTypeP = maxShellType - 1;
  static std::size_t shellArrayIndexP = 0;

  pTaskList->clear();
  pTaskList->reserve(maxGrainSize);

  if (initialize == true) {
    this->clearCutoffStats(maxShellType);
    shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    shellTypeP = maxShellType - 1;
    shellArrayIndexP = 0;

    return true;
  }

  int grainSize = 0;
  Task task;
  for (; shellTypeP >= 0;) {
    const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
    const size_t shellArraySizeP = shellArrayP.size();

    for (; shellArrayIndexP < shellArraySizeP;) {
      const index_type shellIndexP = shellArrayP[shellArrayIndexP];
      task.shellIndex1 = shellIndexP;

      pTaskList->push_back(task);
      ++shellArrayIndexP;
      ++grainSize;

      if (grainSize >= maxGrainSize) {
        return true;
      }
    }
    shellArrayIndexP = 0;
    --shellTypeP;
  }

  return (pTaskList->empty() != true);
}

bool DfTaskCtrl::getQueue2(const TlOrbitalInfoObject& orbitalInfo,
                           const bool isCutoffByDistribution,
                           const int maxGrainSize,
                           std::vector<Task2>* pTaskList, bool initialize) {
  assert(pTaskList != NULL);

  const int maxShellType = orbitalInfo.getMaxShellType();
  static ShellArrayTable shellArrayTable;
  static int shellTypeP = maxShellType - 1;
  static int shellTypeQ = maxShellType - 1;
  static std::size_t shellArrayIndexP = 0;
  static std::size_t shellArrayIndexQ = 0;
  static DistributedCutoffTable dct;
  static std::size_t progress = 0;
  static std::size_t progress_all = 0;
  static int progress_stage = 0;
  static std::size_t division = 0;

  pTaskList->clear();
  pTaskList->reserve(maxGrainSize);

  if (initialize == true) {
    this->clearCutoffStats(maxShellType);
    shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    shellTypeP = maxShellType - 1;
    shellTypeQ = maxShellType - 1;
    shellArrayIndexP = 0;
    shellArrayIndexQ = 0;

    if (isCutoffByDistribution == true) {
      dct = this->makeDistributedCutoffTable(orbitalInfo);

      // setup progress status
      progress = 0;
      progress_all =
          this->getTotalCalcAmount2(orbitalInfo, shellArrayTable, dct);
      progress_stage = 0;
      division = progress_all * 0.1;
      // pre-screening report
      this->prescreeningReport();
    } else {
      progress = 0;
      progress_all = this->getTotalCalcAmount2(orbitalInfo, shellArrayTable);
      progress_stage = 0;
      division = progress_all * 0.1;
    }

    this->log_.info(TlUtils::format("progress_all: %ld", progress_all));

    return true;
  }

  int grainSize = 0;
  Task2 task;
  for (; shellTypeP >= 0;) {
    const int maxStepsP = 2 * shellTypeP + 1;
    const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
    const size_t shellArraySizeP = shellArrayP.size();

    for (; shellTypeQ >= 0;) {
      const int maxStepsQ = 2 * shellTypeQ + 1;

      for (; shellArrayIndexP < shellArraySizeP;) {
        const index_type shellIndexP = shellArrayP[shellArrayIndexP];
        task.shellIndex1 = shellIndexP;

        ShellArray shellArrayQ;
        if (isCutoffByDistribution == true) {
          shellArrayQ = dct[shellIndexP][shellTypeQ];
        } else {
          shellArrayQ = shellArrayTable[shellTypeQ];
        }
        ShellArray::iterator qItEnd = std::upper_bound(
            shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
        const std::size_t shellArraySizeQ =
            std::distance(shellArrayQ.begin(), qItEnd);

        for (; shellArrayIndexQ < shellArraySizeQ;) {
          const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
          task.shellIndex2 = shellIndexQ;

          pTaskList->push_back(task);
          ++shellArrayIndexQ;
          ++grainSize;
          progress += maxStepsP * maxStepsQ;

          if (grainSize >= maxGrainSize) {
            if (progress >= (progress_stage * division)) {
              const double progressLevel =
                  (double)progress / (double)progress_all * 100.0;
              this->log_.info(
                  TlUtils::format("progress: %3.0f%% done.", progressLevel));

              progress_stage = progress / division + 1;
            }
            return true;
          }
        }
        shellArrayIndexQ = 0;
        ++shellArrayIndexP;
      }
      shellArrayIndexP = 0;
      --shellTypeQ;
    }
    --shellTypeP;
    shellTypeQ = maxShellType - 1;
  }

  if (pTaskList->empty() != true) {
    const double progressLevel =
        (double)progress / (double)progress_all * 100.0;
    this->log_.info(TlUtils::format("progress: %3.0f%% done.", progressLevel));
    return true;
  }

  // this->cutoffReport();
  return false;
}

bool DfTaskCtrl::getQueue2(const TlOrbitalInfoObject& orbitalInfo1,
                           const TlOrbitalInfoObject& orbitalInfo2,
                           const bool isCutoffByDistribution,
                           const int maxGrainSize,
                           std::vector<Task2>* pTaskList, bool initialize) {
  assert(pTaskList != NULL);

  const int maxShellType1 = orbitalInfo1.getMaxShellType();
  const int maxShellType2 = orbitalInfo2.getMaxShellType();
  static ShellArrayTable shellArrayTable1;
  static ShellArrayTable shellArrayTable2;
  static int shellTypeP = maxShellType1 - 1;
  static int shellTypeQ = maxShellType2 - 1;
  static std::size_t shellArrayIndexP = 0;
  static std::size_t shellArrayIndexQ = 0;
  static DistributedCutoffTable dct;
  static std::size_t progress = 0;
  static std::size_t progress_all = 0;
  static int progress_stage = 0;
  static std::size_t division = 0;

  pTaskList->clear();
  pTaskList->reserve(maxGrainSize);

  if (initialize == true) {
    this->clearCutoffStats(std::max(maxShellType1, maxShellType2));
    shellArrayTable1 = this->makeShellArrayTable(orbitalInfo1);
    shellArrayTable2 = this->makeShellArrayTable(orbitalInfo2);
    shellTypeP = maxShellType1 - 1;
    shellTypeQ = maxShellType2 - 1;
    shellArrayIndexP = 0;
    shellArrayIndexQ = 0;

    if (isCutoffByDistribution == true) {
      dct = this->makeDistributedCutoffTable(orbitalInfo1, orbitalInfo2);

      // setup progress status
      progress = 0;
      progress_all = this->getTotalCalcAmount2(
          orbitalInfo1, orbitalInfo2, shellArrayTable1, shellArrayTable2, dct);
      progress_stage = 0;
      division = progress_all * 0.1;
      // pre-screening report
      this->prescreeningReport();
    } else {
      progress = 0;
      progress_all = this->getTotalCalcAmount2(
          orbitalInfo1, orbitalInfo2, shellArrayTable1, shellArrayTable2);
      progress_stage = 0;
      division = progress_all * 0.1;
    }

    return true;
  }

  int grainSize = 0;
  Task2 task;
  for (; shellTypeP >= 0;) {
    const int maxStepsP = 2 * shellTypeP + 1;
    const ShellArray& shellArrayP = shellArrayTable1[shellTypeP];
    const size_t shellArraySizeP = shellArrayP.size();

    for (; shellTypeQ >= 0;) {
      const int maxStepsQ = 2 * shellTypeQ + 1;

      for (; shellArrayIndexP < shellArraySizeP;) {
        const index_type shellIndexP = shellArrayP[shellArrayIndexP];
        task.shellIndex1 = shellIndexP;

        ShellArray shellArrayQ;
        if (isCutoffByDistribution == true) {
          shellArrayQ = dct[shellIndexP][shellTypeQ];
        } else {
          shellArrayQ = shellArrayTable2[shellTypeQ];
        }
        const std::size_t shellArraySizeQ = shellArrayQ.size();

        for (; shellArrayIndexQ < shellArraySizeQ;) {
          const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
          task.shellIndex2 = shellIndexQ;

          pTaskList->push_back(task);
          ++shellArrayIndexQ;
          ++grainSize;
          progress += maxStepsP * maxStepsQ;

          if (grainSize >= maxGrainSize) {
            if (progress >= (progress_stage * division)) {
              const double progressLevel =
                  (double)progress / (double)progress_all * 100.0;
              this->log_.info(
                  TlUtils::format("progress: %3.0f%% done.", progressLevel));

              progress_stage = progress / division + 1;
            }
            return true;
          }
        }
        shellArrayIndexQ = 0;
        ++shellArrayIndexP;
      }
      shellArrayIndexP = 0;
      --shellTypeQ;
    }
    --shellTypeP;
    shellTypeQ = maxShellType2 - 1;
  }

  if (pTaskList->empty() != true) {
    const double progressLevel =
        (double)progress / (double)progress_all * 100.0;
    this->log_.info(TlUtils::format("progress: %3.0f%% done.", progressLevel));
    return true;
  }

  this->cutoffReport();
  return false;
}

bool DfTaskCtrl::getQueue4(const TlOrbitalInfoObject& orbitalInfo,
                           const TlSparseSymmetricMatrix& schwarzTable,
                           const int maxGrainSize,
                           std::vector<Task4>* pTaskList, bool initialize) {
  const int maxShellType = orbitalInfo.getMaxShellType();
  static ShellArrayTable shellArrayTable;
  static ShellPairArrayTable shellPairArrayTable;
  static int shellTypeP = maxShellType - 1;
  static int shellTypeQ = maxShellType - 1;
  static int shellTypeR = maxShellType - 1;
  static int shellTypeS = maxShellType - 1;
  static std::size_t prIndex = 0;
  static std::size_t shellArrayIndexQ = 0;
  static std::size_t shellArrayIndexS = 0;
  static DistributedCutoffTable dct;
  static std::size_t progress = 0;
  static std::size_t progress_all = 0;
  static int progress_stage = 0;
  static std::size_t division = 0;

  if (initialize == true) {
    this->clearCutoffStats(maxShellType);
    shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    shellPairArrayTable =
        this->getShellPairArrayTable(orbitalInfo, shellArrayTable);
    shellPairArrayTable = this->selectShellPairArrayTableByDensity(
        shellPairArrayTable, orbitalInfo);
    shellTypeP = maxShellType - 1;
    shellTypeQ = maxShellType - 1;
    shellTypeR = maxShellType - 1;
    shellTypeS = maxShellType - 1;
    prIndex = 0;
    shellArrayIndexQ = 0;
    shellArrayIndexS = 0;

    // make DistributedCutoffTable
    dct = this->makeDistributedCutoffTable(orbitalInfo);

    // setup progress status
    progress = 0;
    progress_all =
        this->getTotalCalcAmount4(orbitalInfo, shellPairArrayTable, dct);
    progress_stage = 0;
    division = progress_all * 0.1;

    // pre-screening report
    this->prescreeningReport();

    return true;
  }

  assert(pTaskList != NULL);
  pTaskList->clear();
  pTaskList->reserve(maxGrainSize);

  int grainSize = 0;
  DfTaskCtrl::Task4 task;
  for (; shellTypeP >= 0;) {
    const int maxStepsP = 2 * shellTypeP + 1;
    for (; shellTypeR >= 0;) {
      const int maxStepsR = 2 * shellTypeR + 1;
      const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
      const ShellPairArray& shellPairArray_PR =
          shellPairArrayTable[shellPairType_PR];
      const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();

      for (; prIndex < numOfShellPairArray_PR;) {
        const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
        const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;
        task.shellIndex1 = shellIndexP;
        task.shellIndex3 = shellIndexR;

        for (; shellTypeQ >= 0;) {
          const int maxStepsQ = 2 * shellTypeQ + 1;
          const ShellArray shellArrayQ = dct[shellIndexP][shellTypeQ];
          ShellArray::const_iterator qItEnd = std::upper_bound(
              shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
          const std::size_t shellArraySizeQ =
              std::distance(shellArrayQ.begin(), qItEnd);

          for (; shellTypeS >= 0;) {
            const int maxStepsS = 2 * shellTypeS + 1;
            const ShellArray shellArrayS = dct[shellIndexR][shellTypeS];

            for (; shellArrayIndexQ < shellArraySizeQ;) {
              const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
              task.shellIndex2 = shellIndexQ;

              const index_type maxShellIndexS =
                  (shellIndexP == shellIndexR) ? shellIndexQ : shellIndexR;
              ShellArray::const_iterator sItEnd = std::upper_bound(
                  shellArrayS.begin(), shellArrayS.end(), maxShellIndexS);
              const std::size_t shellArraySizeS =
                  std::distance(shellArrayS.begin(), sItEnd);
              for (; shellArrayIndexS < shellArraySizeS;) {
                const index_type shellIndexS = shellArrayS[shellArrayIndexS];
                task.shellIndex4 = shellIndexS;

                // schwarz cutoff
                const int shellQuartetType = this->getShellQuartetType(
                    orbitalInfo, orbitalInfo, shellTypeP, shellTypeQ,
                    shellTypeR, shellTypeS);
                const bool isAlive = this->isAliveBySchwarzCutoff(
                    shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
                    shellQuartetType, schwarzTable, this->cutoffThreshold_);
                ++shellArrayIndexS;
                progress += maxStepsP * maxStepsQ * maxStepsR * maxStepsS;
                if (isAlive == true) {
                  pTaskList->push_back(task);
                  ++grainSize;

                  if (grainSize >= maxGrainSize) {
                    if (progress >= (progress_stage * division)) {
                      const double progressLevel =
                          (double)progress / (double)progress_all * 100.0;
                      this->log_.info(TlUtils::format("progress: %3.0f%% done.",
                                                      progressLevel));

                      progress_stage = progress / division + 1;
                    }
                    return true;
                  }
                }
              }
              shellArrayIndexS = 0;
              ++shellArrayIndexQ;
            }
            shellArrayIndexQ = 0;
            --shellTypeS;
          }
          shellTypeS = maxShellType - 1;
          --shellTypeQ;
        }
        shellTypeQ = maxShellType - 1;
        ++prIndex;
      }
      prIndex = 0;
      --shellTypeR;
    }
    shellTypeR = maxShellType - 1;
    --shellTypeP;
  }

  if (pTaskList->empty() != true) {
    const double progressLevel =
        (double)progress / (double)progress_all * 100.0;
    this->log_.info(TlUtils::format("progress: %3.0f%% done.", progressLevel));
    return true;
  }

  // this->cutoffReport();
  return false;
}

bool DfTaskCtrl::getQueue4(const TlOrbitalInfoObject& orbitalInfo_PQ,
                           const TlOrbitalInfoObject& orbitalInfo_RS,
                           const TlSparseSymmetricMatrix& schwarzTable_PQ,
                           const TlSparseSymmetricMatrix& schwarzTable_RS,
                           const int maxGrainSize,
                           std::vector<Task4>* pTaskList, bool initialize) {
  const int maxShellType_PQ = orbitalInfo_PQ.getMaxShellType();
  const int maxShellType_RS = orbitalInfo_RS.getMaxShellType();
  static ShellArrayTable shellArrayTable_PQ;
  static ShellArrayTable shellArrayTable_RS;
  static ShellPairArrayTable shellPairArrayTable;
  static int shellTypeP = maxShellType_PQ - 1;
  static int shellTypeQ = maxShellType_PQ - 1;
  static int shellTypeR = maxShellType_RS - 1;
  static int shellTypeS = maxShellType_RS - 1;
  static std::size_t prIndex = 0;
  static std::size_t shellArrayIndexQ = 0;
  static std::size_t shellArrayIndexS = 0;
  static DistributedCutoffTable dct_PQ;
  static DistributedCutoffTable dct_RS;
  static std::size_t progress = 0;
  static std::size_t progress_all = 0;
  static int progress_stage = 0;
  static std::size_t division = 0;

  if (initialize == true) {
    this->clearCutoffStats(std::max(maxShellType_PQ, maxShellType_RS));
    shellArrayTable_PQ = this->makeShellArrayTable(orbitalInfo_PQ);
    shellArrayTable_RS = this->makeShellArrayTable(orbitalInfo_RS);

    // P-R
    // shellPairArrayTable = this->getShellPairArrayTable(orbitalInfo_PQ,
    //                                                    shellArrayTable_PQ,
    //                                                    orbitalInfo_RS,
    //                                                    shellArrayTableRS);
    {
      shellPairArrayTable.clear();
      shellPairArrayTable.resize(maxShellType_PQ * maxShellType_RS);
      for (int shellTypeP = maxShellType_PQ - 1; shellTypeP >= 0;
           --shellTypeP) {
        const ShellArray& shellArrayP = shellArrayTable_PQ[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeR = maxShellType_RS - 1; shellTypeR >= 0;
             --shellTypeR) {
          const ShellArray& shellArrayR = shellArrayTable_RS[shellTypeR];
          ShellArray::const_iterator rItEnd = shellArrayR.end();

          const int shellPairType_PR =
              shellTypeP * maxShellType_PQ + shellTypeR;
          for (ShellArray::const_iterator pIt = shellArrayP.begin();
               pIt != pItEnd; ++pIt) {
            const index_type indexP = *pIt;

            for (ShellArray::const_iterator rIt = shellArrayR.begin();
                 rIt != rItEnd; ++rIt) {
              const index_type indexR = *rIt;

              ShellPair shellPair(indexP, indexR);
              shellPairArrayTable[shellPairType_PR].push_back(shellPair);
            }
          }
        }
      }
    }
    shellPairArrayTable = this->selectShellPairArrayTableByDensity(
        shellPairArrayTable, orbitalInfo_PQ, orbitalInfo_RS);
    shellTypeP = maxShellType_PQ - 1;
    shellTypeQ = maxShellType_PQ - 1;
    shellTypeR = maxShellType_RS - 1;
    shellTypeS = maxShellType_RS - 1;
    prIndex = 0;
    shellArrayIndexQ = 0;
    shellArrayIndexS = 0;

    // make DistributedCutoffTable
    dct_PQ = this->makeDistributedCutoffTable(orbitalInfo_PQ);
    dct_RS = this->makeDistributedCutoffTable(orbitalInfo_RS);

    // setup progress status
    progress = 0;
    progress_all = this->getTotalCalcAmount4(
        orbitalInfo_PQ, orbitalInfo_RS, shellPairArrayTable, dct_PQ, dct_RS);
    progress_stage = 0;
    division = progress_all * 0.1;

    // pre-screening report
    this->prescreeningReport();

    return true;
  }

  assert(pTaskList != NULL);
  pTaskList->clear();
  pTaskList->reserve(maxGrainSize);

  int grainSize = 0;
  DfTaskCtrl::Task4 task;
  for (; shellTypeP >= 0;) {
    const int maxStepsP = 2 * shellTypeP + 1;
    for (; shellTypeR >= 0;) {
      const int maxStepsR = 2 * shellTypeR + 1;
      const int shellPairType_PR = shellTypeP * maxShellType_PQ + shellTypeR;
      const ShellPairArray& shellPairArray_PR =
          shellPairArrayTable[shellPairType_PR];
      const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();

      for (; prIndex < numOfShellPairArray_PR;) {
        const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
        const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;
        task.shellIndex1 = shellIndexP;
        task.shellIndex3 = shellIndexR;

        for (; shellTypeQ >= 0;) {
          const int maxStepsQ = 2 * shellTypeQ + 1;
          const ShellArray shellArrayQ = dct_PQ[shellIndexP][shellTypeQ];
          ShellArray::const_iterator qItEnd = std::upper_bound(
              shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
          const std::size_t shellArraySizeQ =
              std::distance(shellArrayQ.begin(), qItEnd);

          for (; shellTypeS >= 0;) {
            const int maxStepsS = 2 * shellTypeS + 1;
            const ShellArray shellArrayS = dct_RS[shellIndexR][shellTypeS];

            for (; shellArrayIndexQ < shellArraySizeQ;) {
              const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
              task.shellIndex2 = shellIndexQ;

              const index_type maxShellIndexS = shellIndexR;
              ShellArray::const_iterator sItEnd = std::upper_bound(
                  shellArrayS.begin(), shellArrayS.end(), maxShellIndexS);
              const std::size_t shellArraySizeS =
                  std::distance(shellArrayS.begin(), sItEnd);
              for (; shellArrayIndexS < shellArraySizeS;) {
                const index_type shellIndexS = shellArrayS[shellArrayIndexS];
                task.shellIndex4 = shellIndexS;

                // schwarz cutoff
                const int shellQuartetType = this->getShellQuartetType(
                    orbitalInfo_PQ, orbitalInfo_RS, shellTypeP, shellTypeQ,
                    shellTypeR, shellTypeS);
                const bool isAlive = this->isAliveBySchwarzCutoff(
                    shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
                    shellQuartetType, schwarzTable_PQ, schwarzTable_RS,
                    this->cutoffThreshold_);
                ++shellArrayIndexS;
                progress += maxStepsP * maxStepsQ * maxStepsR * maxStepsS;
                if (isAlive == true) {
                  pTaskList->push_back(task);
                  ++grainSize;

                  if (grainSize >= maxGrainSize) {
                    if (progress >= (progress_stage * division)) {
                      const double progressLevel =
                          (double)progress / (double)progress_all * 100.0;
                      this->log_.info(TlUtils::format("progress: %3.0f%% done.",
                                                      progressLevel));

                      progress_stage = progress / division + 1;
                    }
                    return true;
                  }
                }
              }
              shellArrayIndexS = 0;
              ++shellArrayIndexQ;
            }
            shellArrayIndexQ = 0;
            --shellTypeS;
          }
          shellTypeS = maxShellType_RS - 1;
          --shellTypeQ;
        }
        shellTypeQ = maxShellType_PQ - 1;
        ++prIndex;
      }
      prIndex = 0;
      --shellTypeR;
    }
    shellTypeR = maxShellType_RS - 1;
    --shellTypeP;
  }

  if (pTaskList->empty() != true) {
    const double progressLevel =
        (double)progress / (double)progress_all * 100.0;
    this->log_.info(TlUtils::format("progress: %3.0f%% done.", progressLevel));
    return true;
  }

  // this->cutoffReport();
  return false;
}

bool DfTaskCtrl::getQueue4_K(
    const TlOrbitalInfoObject& orbitalInfo,
    const TlSparseSymmetricMatrix& schwarzTable, const TlMatrixObject& P,
    const std::vector<TlMatrixObject::index_type>& rowIndexes,
    const std::vector<TlMatrixObject::index_type>& colIndexes,
    const int maxGrainSize, std::vector<Task4>* pTaskList, bool initialize) {
  assert(pTaskList != NULL);

  const int maxShellType = orbitalInfo.getMaxShellType();
  static ShellArrayTable shellArrayTable;
  static ShellPairArrayTable shellPairArrayTable;
  // static int shellTypeP = maxShellType -1;
  // static int shellTypeQ = maxShellType -1;
  // static int shellTypeR = maxShellType -1;
  // static int shellTypeS = maxShellType -1;
  // static std::size_t prIndex = 0;
  // static std::size_t shellArrayIndexQ = 0;
  // static std::size_t shellArrayIndexS = 0;

  pTaskList->clear();
  pTaskList->reserve(maxGrainSize);

  if (initialize == true) {
    this->clearCutoffStats(maxShellType);
    shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    shellPairArrayTable =
        this->getShellPairArrayTable(orbitalInfo, shellArrayTable);
    shellPairArrayTable = this->selectShellPairArrayTableByDensity(
        shellPairArrayTable, orbitalInfo);
    // shellTypeP = maxShellType -1;
    // shellTypeQ = maxShellType -1;
    // shellTypeR = maxShellType -1;
    // shellTypeS = maxShellType -1;
    // prIndex = 0;
    // shellArrayIndexQ = 0;
    // shellArrayIndexS = 0;

    return true;
  }

  // int grainSize = 0;
  DfTaskCtrl::Task4 task;

  const index_type numOfRowIndexes = rowIndexes.size();
  const index_type numOfColIndexes = colIndexes.size();
  for (index_type index_p = 0; index_p < numOfRowIndexes; ++index_p) {
    const index_type shellIndexP = rowIndexes[index_p];
    const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
    const int maxStepsP = shellTypeP * 2 + 1;

    task.shellIndex1 = shellIndexP;

    for (index_type index_r = 0; index_r < numOfColIndexes; ++index_r) {
      const index_type shellIndexR = colIndexes[index_r];
      const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
      const int maxStepsR = shellTypeR * 2 + 1;

      const double max_P_pr = std::max<double>(
          1.0E-20,
          this->getMaxValue(P, shellIndexP, maxStepsP, shellIndexR, maxStepsR));
      task.shellIndex3 = shellIndexR;

      for (index_type index_q = 0; index_q < numOfRowIndexes; ++index_q) {
        const index_type shellIndexQ = rowIndexes[index_q];
        if (shellIndexQ > shellIndexP) {
          continue;
        }
        const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
        const int maxStepsQ = shellTypeQ * 2 + 1;

        const double max_P_qr = this->getMaxValue(P, shellIndexQ, maxStepsQ,
                                                  shellIndexR, maxStepsR);
        const double max_P_pqr = std::max(max_P_pr, max_P_qr);
        task.shellIndex2 = shellIndexQ;

        const index_type maxShellIndexS =
            (shellIndexP == shellIndexR) ? shellIndexQ : shellIndexR;
        for (index_type index_s = 0; index_s < numOfColIndexes; ++index_s) {
          const index_type shellIndexS = colIndexes[index_s];
          if (shellIndexS > maxShellIndexS) {
            continue;
          }
          const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
          const int maxStepsS = shellTypeS * 2 + 1;

          const double max_P_ps = this->getMaxValue(P, shellIndexP, maxStepsP,
                                                    shellIndexS, maxStepsS);
          const double max_P_qs = this->getMaxValue(P, shellIndexQ, maxStepsQ,
                                                    shellIndexS, maxStepsS);
          // schwarz cutoff
          const double max_P_pqs = std::max(max_P_ps, max_P_qs);
          const double max_P_pqrs = std::max(max_P_pqr, max_P_pqs);
          const double cutoffThreshold = this->cutoffThreshold_ / max_P_pqrs;
          const int shellQuartetType =
              this->getShellQuartetType(orbitalInfo, orbitalInfo, shellTypeP,
                                        shellTypeQ, shellTypeR, shellTypeS);
          const bool isAlive = this->isAliveBySchwarzCutoff(
              shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
              shellQuartetType, schwarzTable, cutoffThreshold);
          if (isAlive == true) {
            task.shellIndex4 = shellIndexS;
            pTaskList->push_back(task);
          }
        }
      }
    }
  }

  if (pTaskList->empty() != true) {
    return true;
  }

  return false;
}

bool DfTaskCtrl::getQueue4_K0(
    const TlOrbitalInfoObject& orbitalInfo,
    const TlSparseSymmetricMatrix& schwarzTable, const TlMatrixObject& P,
    const std::vector<TlMatrixObject::index_type>& rowIndexes,
    const std::vector<TlMatrixObject::index_type>& colIndexes,
    const int maxGrainSize, std::vector<Task4>* pTaskList, bool initialize) {
  assert(pTaskList != NULL);

  const int maxShellType = orbitalInfo.getMaxShellType();
  static ShellArrayTable shellArrayTable;
  static ShellPairArrayTable shellPairArrayTable;
  static int shellTypeP = maxShellType - 1;
  static int shellTypeQ = maxShellType - 1;
  static int shellTypeR = maxShellType - 1;
  static int shellTypeS = maxShellType - 1;
  static std::size_t prIndex = 0;
  static std::size_t shellArrayIndexQ = 0;
  static std::size_t shellArrayIndexS = 0;

  pTaskList->clear();
  pTaskList->reserve(maxGrainSize);

  if (initialize == true) {
    this->clearCutoffStats(maxShellType);
    shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    shellPairArrayTable =
        this->getShellPairArrayTable(orbitalInfo, shellArrayTable);
    shellPairArrayTable = this->selectShellPairArrayTableByDensity(
        shellPairArrayTable, orbitalInfo);
    shellTypeP = maxShellType - 1;
    shellTypeQ = maxShellType - 1;
    shellTypeR = maxShellType - 1;
    shellTypeS = maxShellType - 1;
    prIndex = 0;
    shellArrayIndexQ = 0;
    shellArrayIndexS = 0;

    return true;
  }

  int grainSize = 0;
  DfTaskCtrl::Task4 task;
  for (; shellTypeP >= 0;) {
    const int maxStepsP = 2 * shellTypeP + 1;

    for (; shellTypeR >= 0;) {
      const int maxStepsR = 2 * shellTypeR + 1;

      const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
      const ShellPairArray& shellPairArray_PR =
          shellPairArrayTable[shellPairType_PR];
      const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();

      for (; prIndex < numOfShellPairArray_PR;) {
        const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
        const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;
        task.shellIndex1 = shellIndexP;
        task.shellIndex3 = shellIndexR;

        // if ((std::binary_search(rowIndexes.begin(), rowIndexes.end(),
        // shellIndexP) != true) ||
        //     (std::binary_search(colIndexes.begin(), colIndexes.end(),
        //     shellIndexR) != true)) { continue;
        // }
        const double max_P_pr = std::max<double>(
            1.0E-20, this->getMaxValue(P, shellIndexP, maxStepsP, shellIndexR,
                                       maxStepsR));

        for (; shellTypeQ >= 0;) {
          const int maxStepsQ = 2 * shellTypeQ + 1;
          const ShellArray shellArrayQ = this->selectShellArrayByDistribution(
              shellArrayTable[shellTypeQ], shellIndexP, orbitalInfo);
          ShellArray::const_iterator qItEnd = std::upper_bound(
              shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
          const std::size_t shellArraySizeQ =
              std::distance(shellArrayQ.begin(), qItEnd);

          for (; shellTypeS >= 0;) {
            const int maxStepsS = 2 * shellTypeS + 1;
            const ShellArray shellArrayS = this->selectShellArrayByDistribution(
                shellArrayTable[shellTypeS], shellIndexR, orbitalInfo);

            for (; shellArrayIndexQ < shellArraySizeQ;) {
              const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
              task.shellIndex2 = shellIndexQ;

              // if (std::binary_search(rowIndexes.begin(), rowIndexes.end(),
              // shellIndexQ) != true) {
              //     continue;
              // }
              const double max_P_qr = this->getMaxValue(
                  P, shellIndexQ, maxStepsQ, shellIndexR, maxStepsR);
              const double max_P_pqr = std::max(max_P_pr, max_P_qr);

              const index_type maxShellIndexS =
                  (shellIndexP == shellIndexR) ? shellIndexQ : shellIndexR;
              ShellArray::const_iterator sItEnd = std::upper_bound(
                  shellArrayS.begin(), shellArrayS.end(), maxShellIndexS);
              const std::size_t shellArraySizeS =
                  std::distance(shellArrayS.begin(), sItEnd);
              for (; shellArrayIndexS < shellArraySizeS;) {
                const index_type shellIndexS = shellArrayS[shellArrayIndexS];
                task.shellIndex4 = shellIndexS;

                // if (std::binary_search(colIndexes.begin(), colIndexes.end(),
                // shellIndexS) != true) {
                //     continue;
                // }
                const double max_P_ps = this->getMaxValue(
                    P, shellIndexP, maxStepsP, shellIndexS, maxStepsS);
                const double max_P_qs = this->getMaxValue(
                    P, shellIndexQ, maxStepsQ, shellIndexS, maxStepsS);

                // schwarz cutoff
                const double max_P_pqs = std::max(max_P_ps, max_P_qs);
                const double max_P_pqrs = std::max(max_P_pqr, max_P_pqs);
                const double cutoffThreshold =
                    this->cutoffThreshold_ / max_P_pqrs;
                const int shellQuartetType = this->getShellQuartetType(
                    orbitalInfo, orbitalInfo, shellTypeP, shellTypeQ,
                    shellTypeR, shellTypeS);
                const bool isAlive = this->isAliveBySchwarzCutoff(
                    shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
                    shellQuartetType, schwarzTable, cutoffThreshold);
                ++shellArrayIndexS;
                if (isAlive == true) {
                  pTaskList->push_back(task);
                  ++grainSize;

                  if (grainSize >= maxGrainSize) {
                    return true;
                  }
                }
              }
              shellArrayIndexS = 0;
              ++shellArrayIndexQ;
            }
            shellArrayIndexQ = 0;
            --shellTypeS;
          }
          shellTypeS = maxShellType - 1;
          --shellTypeQ;
        }
        shellTypeQ = maxShellType - 1;
        ++prIndex;
      }
      prIndex = 0;
      --shellTypeR;
    }
    shellTypeR = maxShellType - 1;
    --shellTypeP;
  }

  if (pTaskList->empty() != true) {
    return true;
  }

  return false;
}

double DfTaskCtrl::getMaxValue(const TlMatrixObject& P,
                               const TlMatrixObject::index_type row,
                               const TlMatrixObject::index_type d_row,
                               const TlMatrixObject::index_type col,
                               const TlMatrixObject::index_type d_col) {
  const TlMatrixObject::index_type maxRow = row + d_row;
  const TlMatrixObject::index_type maxCol = col + d_col;
  double answer = 0.0;
  for (TlMatrixObject::index_type r = row; r < maxRow; ++r) {
    for (TlMatrixObject::index_type c = col; c < maxCol; ++c) {
      answer = std::max(answer, std::fabs(P.getLocal(r, c)));
    }
  }

  return answer;
}

bool DfTaskCtrl::getQueue_Force4(const TlOrbitalInfoObject& orbitalInfo,
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const int maxGrainSize,
                                 std::vector<Task4>* pTaskList,
                                 bool initialize) {
  assert(pTaskList != NULL);

  const int maxShellType = orbitalInfo.getMaxShellType();
  static ShellArrayTable shellArrayTable;
  static ShellPairArrayTable shellPairArrayTable;
  static int shellTypeP = maxShellType - 1;
  static int shellTypeQ = maxShellType - 1;
  static int shellTypeR = maxShellType - 1;
  static int shellTypeS = maxShellType - 1;
  static std::size_t prIndex = 0;
  static index_type shellArrayIndexQ = 0;
  static index_type shellArrayIndexS = 0;

  pTaskList->clear();
  pTaskList->reserve(maxGrainSize);

  if (initialize == true) {
    this->clearCutoffStats(maxShellType);
    shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    shellPairArrayTable =
        this->getShellPairArrayTable(orbitalInfo, shellArrayTable);
    shellPairArrayTable = this->selectShellPairArrayTableByDensity(
        shellPairArrayTable, orbitalInfo);
    shellTypeP = maxShellType - 1;
    shellTypeQ = maxShellType - 1;
    shellTypeR = maxShellType - 1;
    shellTypeS = maxShellType - 1;
    prIndex = 0;
    shellArrayIndexQ = 0;
    shellArrayIndexS = 0;

    return true;
  }

  int grainSize = 0;
  DfTaskCtrl::Task4 task;
  for (; shellTypeP >= 0;) {
    for (; shellTypeR >= 0;) {
      const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
      const ShellPairArray& shellPairArray_PR =
          shellPairArrayTable[shellPairType_PR];
      const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();

      for (; prIndex < numOfShellPairArray_PR;) {
        const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
        const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;
        task.shellIndex1 = shellIndexP;
        task.shellIndex3 = shellIndexR;

        for (; shellTypeQ >= 0;) {
          // const int maxStepsQ = 2 * shellTypeQ + 1;
          const ShellArray shellArrayQ = this->selectShellArrayByDistribution(
              shellArrayTable[shellTypeQ], shellIndexP, orbitalInfo);
          ShellArray::const_iterator qItEnd = std::upper_bound(
              shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
          const index_type shellArraySizeQ =
              std::distance(shellArrayQ.begin(), qItEnd);

          for (; shellTypeS >= 0;) {
            // const int maxStepsS = 2 * shellTypeS + 1;
            const ShellArray shellArrayS = this->selectShellArrayByDistribution(
                shellArrayTable[shellTypeS], shellIndexR, orbitalInfo);
            ShellArray::const_iterator sItEnd = std::upper_bound(
                shellArrayS.begin(), shellArrayS.end(), shellIndexR);
            const index_type shellArraySizeS =
                std::distance(shellArrayS.begin(), sItEnd);

            for (; shellArrayIndexQ < shellArraySizeQ;) {
              const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];
              task.shellIndex2 = shellIndexQ;

              for (; shellArrayIndexS < shellArraySizeS;) {
                const index_type shellIndexS = shellArrayS[shellArrayIndexS];
                task.shellIndex4 = shellIndexS;

                // schwarz cutoff
                const int shellQuartetType = this->getShellQuartetType(
                    orbitalInfo, orbitalInfo, shellTypeP, shellTypeQ,
                    shellTypeR, shellTypeS);
                const bool isAlive = this->isAliveBySchwarzCutoff(
                    shellIndexP, shellIndexQ, shellIndexP, shellIndexQ,
                    shellQuartetType, schwarzTable, this->cutoffThreshold_);
                ++shellArrayIndexS;
                if (isAlive == true) {
                  pTaskList->push_back(task);
                  ++grainSize;

                  if (grainSize >= maxGrainSize) {
                    return true;
                  }
                }
              }
              shellArrayIndexS = 0;
              ++shellArrayIndexQ;
            }
            shellArrayIndexQ = 0;
            --shellTypeS;
          }
          shellTypeS = maxShellType - 1;
          --shellTypeQ;
        }
        shellTypeQ = maxShellType - 1;
        ++prIndex;
      }
      prIndex = 0;
      --shellTypeR;
    }
    shellTypeR = maxShellType - 1;
    --shellTypeP;
  }

  if (pTaskList->empty() != true) {
    return true;
  }

  return false;
}

// J. Chem. Phys.,105,2726 (1996)
// eq.32
DfTaskCtrl::ShellArray DfTaskCtrl::selectShellArrayByDistribution(
    const ShellArray& inShellArray, const index_type companionShellIndex,
    const TlOrbitalInfoObject& orbitalInfo) {
  ShellArray answer;
  answer.reserve(inShellArray.size());

  const TlPosition posB = orbitalInfo.getPosition(companionShellIndex);
  const int shellTypeB = orbitalInfo.getShellType(companionShellIndex);
  // orbitalInfoのPGTOリストは指数が小さい順にソートされているため、
  // 最初(index=0)の指数のみをチェックすれば良い。
  const double exponentB = orbitalInfo.getExponent(companionShellIndex, 0);

  // check
  const int maxShellType = orbitalInfo.getMaxShellType();
  static const double INV_EQ32_COEF =
      1.0 / (std::pow(2.0 * TlMath::PI(), 0.25) * TlMath::PI());
  const double threshold = this->cutoffEpsilon_distribution_ * INV_EQ32_COEF;
  ShellArray::const_iterator itEnd = inShellArray.end();
  for (ShellArray::const_iterator it = inShellArray.begin(); it != itEnd;
       ++it) {
    const int shellPairType =
        orbitalInfo.getShellType(*it) * maxShellType + shellTypeB;
    const double distance2 =
        posB.squareDistanceFrom(orbitalInfo.getPosition(*it));
    const double exponentA = orbitalInfo.getExponent(*it, 0);

    const double zetaP = exponentA + exponentB;
    const double zeta = exponentA * exponentB / zetaP;

    const double exponent = -zeta * distance2;
    const double coef = 1.0 / (std::pow(zetaP, 1.25));

    if (coef * std::exp(exponent) >= threshold) {
      answer.push_back(*it);

#pragma omp critical(DfTaskCtrl__selectShellArrayByDistribution_alive)
      { ++(this->cutoffAlive_distribution_[shellPairType]); }
    }

#pragma omp critical(DfTaskCtrl__selectShellArrayByDistribution_all)
    { ++(this->cutoffAll_distribution_[shellPairType]); }
  }

  // swap technique
  ShellArray(answer).swap(answer);

  return answer;
}

DfTaskCtrl::ShellArray DfTaskCtrl::selectShellArrayByDistribution(
    const ShellArray& inShellArray, const index_type companionShellIndex) {
  ShellArray::const_iterator it = std::upper_bound(
      inShellArray.begin(), inShellArray.end(), companionShellIndex);
  ShellArray answer(std::distance(inShellArray.begin(), it));
  std::copy(inShellArray.begin(), it, answer.begin());

  return answer;
}

// J. Chem. Phys.,105,2726 (1996)
// eq.32
DfTaskCtrl::DistributedCutoffTable DfTaskCtrl::makeDistributedCutoffTable(
    const TlOrbitalInfoObject& orbitalInfo) {
  static const double INV_EQ32_COEF =
      1.0 / (std::pow(2.0 * TlMath::PI(), 0.25) * TlMath::PI());
  const double threshold = this->cutoffEpsilon_distribution_ * INV_EQ32_COEF;
  const int maxShellType = orbitalInfo.getMaxShellType();
  // const index_type numOfAOs = orbitalInfo.getNumOfOrbitals();

  const std::vector<index_type> orbList =
      orbitalInfo.getStartIndexArrayOfShellGroup();
  const int orbListSize = orbList.size();

  DistributedCutoffTable answer;

  for (int i = 0; i < orbListSize; ++i) {
    const index_type indexI = orbList[i];
    answer[indexI].clear();
    answer[indexI].resize(maxShellType + 1);
    for (int j = 0; j < maxShellType + 1; ++j) {
      answer[indexI][j].reserve(orbListSize);
    }
  }

  for (int i = 0; i < orbListSize; ++i) {
    const index_type indexI = orbList[i];
    const int shellTypeI = orbitalInfo.getShellType(indexI);
    const TlPosition posI = orbitalInfo.getPosition(indexI);
    // orbitalInfoのPGTOリストは指数が小さい順にソートされているため、
    // 最初(index=0)の指数のみをチェックすれば良い。
    const double exponentI = orbitalInfo.getExponent(indexI, 0);

    for (int j = 0; j < i; ++j) {
      const index_type indexJ = orbList[j];
      const int shellTypeJ = orbitalInfo.getShellType(indexJ);

      // judge
      const int shellPairType = shellTypeI * maxShellType + shellTypeJ;

      const TlPosition posJ = orbitalInfo.getPosition(indexJ);
      const double distance2 = posJ.squareDistanceFrom(posI);

      const double exponentJ = orbitalInfo.getExponent(indexJ, 0);
      const double zetaIJ = exponentI + exponentJ;
      const double zeta = exponentI * exponentJ / zetaIJ;
      const double exponent = -zeta * distance2;

      const double coef = 1.0 / (std::pow(zetaIJ, 1.25));

      if (coef * std::exp(exponent) >= threshold) {
        answer[indexI][shellTypeJ].push_back(indexJ);
        answer[indexJ][shellTypeI].push_back(indexI);

#pragma omp critical(DfTaskCtrl__makeDistributedCutoffTable_alive)
        { ++(this->cutoffAlive_distribution_[shellPairType]); }
      }

#pragma omp critical(DfTaskCtrl__makeDistributedCutoffTable_all)
      { ++(this->cutoffAll_distribution_[shellPairType]); }
    }
    answer[indexI][shellTypeI].push_back(indexI);
  }

  // sort and optimize
  for (int i = 0; i < orbListSize; ++i) {
    const index_type indexI = orbList[i];

    for (int j = 0; j < maxShellType + 1; ++j) {
      // swap technique
      ShellArray(answer[indexI][j]).swap(answer[indexI][j]);

      // sort
      std::sort(answer[indexI][j].begin(), answer[indexI][j].end());
    }
  }

  return answer;
}

// J. Chem. Phys.,105,2726 (1996)
// eq.32
DfTaskCtrl::DistributedCutoffTable DfTaskCtrl::makeDistributedCutoffTable(
    const TlOrbitalInfoObject& orbitalInfo1,
    const TlOrbitalInfoObject& orbitalInfo2) {
  static const double INV_EQ32_COEF =
      1.0 / (std::pow(2.0 * TlMath::PI(), 0.25) * TlMath::PI());
  const double threshold = this->cutoffEpsilon_distribution_ * INV_EQ32_COEF;
  // const int maxShellType1 = orbitalInfo1.getMaxShellType();
  const int maxShellType2 = orbitalInfo2.getMaxShellType();
  // const index_type numOfOrbitals1 = orbitalInfo1.getNumOfOrbitals();
  // const index_type numOfOrbitals2 = orbitalInfo2.getNumOfOrbitals();

  const std::vector<index_type> orbList1 =
      orbitalInfo1.getStartIndexArrayOfShellGroup();
  const int orbListSize1 = orbList1.size();
  const std::vector<index_type> orbList2 =
      orbitalInfo2.getStartIndexArrayOfShellGroup();
  const int orbListSize2 = orbList2.size();

  DistributedCutoffTable answer;

  for (int i = 0; i < orbListSize1; ++i) {
    const index_type indexI = orbList1[i];
    answer[indexI].clear();
    answer[indexI].resize(maxShellType2 + 1);
    for (int j = 0; j < maxShellType2 + 1; ++j) {
      answer[indexI][j].reserve(orbListSize2);
    }
  }

  for (int i = 0; i < orbListSize1; ++i) {
    const index_type indexI = orbList1[i];
    const int shellTypeI = orbitalInfo1.getShellType(indexI);
    const TlPosition posI = orbitalInfo1.getPosition(indexI);
    // orbitalInfoのPGTOリストは指数が小さい順にソートされているため、
    // 最初(index=0)の指数のみをチェックすれば良い。
    const double exponentI = orbitalInfo1.getExponent(indexI, 0);

    for (int j = 0; j < orbListSize2; ++j) {
      const index_type indexJ = orbList2[j];
      const int shellTypeJ = orbitalInfo2.getShellType(indexJ);

      // judge
      const int shellPairType = shellTypeI * maxShellType2 + shellTypeJ;

      const TlPosition posJ = orbitalInfo2.getPosition(indexJ);
      const double distance2 = posJ.squareDistanceFrom(posI);

      const double exponentJ = orbitalInfo2.getExponent(indexJ, 0);
      const double zetaIJ = exponentI + exponentJ;
      const double zeta = exponentI * exponentJ / zetaIJ;
      const double exponent = -zeta * distance2;

      const double coef = 1.0 / (std::pow(zetaIJ, 1.25));

      if (coef * std::exp(exponent) >= threshold) {
        answer[indexI][shellTypeJ].push_back(indexJ);

#pragma omp critical(DfTaskCtrl__makeDistributedCutoffTable_alive)
        { ++(this->cutoffAlive_distribution_[shellPairType]); }
      }

#pragma omp critical(DfTaskCtrl__makeDistributedCutoffTable_all)
      { ++(this->cutoffAll_distribution_[shellPairType]); }
    }
  }

  // sort and optimize
  for (int i = 0; i < orbListSize1; ++i) {
    const index_type indexI = orbList1[i];

    for (int j = 0; j < maxShellType2 + 1; ++j) {
      // swap technique
      ShellArray(answer[indexI][j]).swap(answer[indexI][j]);

      // sort
      std::sort(answer[indexI][j].begin(), answer[indexI][j].end());
    }
  }

  return answer;
}

DfTaskCtrl::ShellPairArrayTable DfTaskCtrl::getShellPairArrayTable(
    const TlOrbitalInfoObject& orbitalInfo,
    const ShellArrayTable& shellArrayTable) {
  const int maxShellType = orbitalInfo.getMaxShellType();
  ShellPairArrayTable shellPairArrayTable(maxShellType * maxShellType);

  for (int shellTypeP = maxShellType - 1; shellTypeP >= 0; --shellTypeP) {
    const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
    ShellArray::const_iterator pItEnd = shellArrayP.end();

    for (int shellTypeR = maxShellType - 1; shellTypeR >= 0; --shellTypeR) {
      const ShellArray& shellArrayR = shellArrayTable[shellTypeR];
      ShellArray::const_iterator rItEnd = shellArrayR.end();

      const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
      for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd;
           ++pIt) {
        const index_type indexP = *pIt;

        for (ShellArray::const_iterator rIt = shellArrayR.begin();
             rIt != rItEnd; ++rIt) {
          const index_type indexR = *rIt;

          if (indexP >= indexR) {
            ShellPair shellPair(indexP, indexR);
            shellPairArrayTable[shellPairType_PR].push_back(shellPair);
          }
        }
      }
    }
  }

  return shellPairArrayTable;
}

DfTaskCtrl::ShellPairArrayTable DfTaskCtrl::getShellPairArrayTable(
    const TlOrbitalInfoObject& orbitalInfo1,
    const ShellArrayTable& shellArrayTable1,
    const TlOrbitalInfoObject& orbitalInfo2,
    const ShellArrayTable& shellArrayTable2) {
  const int maxShellType1 = orbitalInfo1.getMaxShellType();
  const int maxShellType2 = orbitalInfo2.getMaxShellType();
  ShellPairArrayTable shellPairArrayTable(maxShellType1 * maxShellType2);

  for (int shellTypeP = maxShellType1 - 1; shellTypeP >= 0; --shellTypeP) {
    const ShellArray& shellArrayP = shellArrayTable1[shellTypeP];
    ShellArray::const_iterator pItEnd = shellArrayP.end();

    for (int shellTypeR = maxShellType2 - 1; shellTypeR >= 0; --shellTypeR) {
      const ShellArray& shellArrayR = shellArrayTable2[shellTypeR];
      ShellArray::const_iterator rItEnd = shellArrayR.end();

      const int shellPairType_PR = shellTypeP * maxShellType1 + shellTypeR;
      for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd;
           ++pIt) {
        const index_type indexP = *pIt;

        for (ShellArray::const_iterator rIt = shellArrayR.begin();
             rIt != rItEnd; ++rIt) {
          const index_type indexR = *rIt;

          if (indexP >= indexR) {
            ShellPair shellPair(indexP, indexR);
            shellPairArrayTable[shellPairType_PR].push_back(shellPair);
          }
        }
      }
    }
  }

  return shellPairArrayTable;
}

// J. Chem. Phys.,105,2726 (1996)
// eq.31
// 1/r cutoff
DfTaskCtrl::ShellPairArrayTable DfTaskCtrl::selectShellPairArrayTableByDensity(
    const ShellPairArrayTable& inShellPairArrayTable,
    const TlOrbitalInfoObject& orbitalInfo) {
  const int maxShellType = this->maxShellType_;
  assert(inShellPairArrayTable.size() ==
         static_cast<std::size_t>(maxShellType * maxShellType));

  const double cutoffThreshold = this->cutoffEpsilon_density_;
  const double CONTRIBUTE_COEF = 2.0 * std::pow(TlMath::PI(), 2.5);

  TlFmt& FmT = TlFmt::getInstance();
  static const double coef[6] = {-0.017450254, 0.132520568,  -0.047915444,
                                 0.792267596,  -0.583721015, 0.697593555};
  // l = 1.0のときのgamma値
  static const double gamma1[6] = {0.03,        0.081722098, 0.222616709,
                                   0.606423483, 1.651939972, 4.5};
  std::vector<double> gamma(6);
  {
    const double inv_lprime = 1.0 / this->lengthScaleParameter_;
    const double ll = inv_lprime * inv_lprime;
    for (int i = 0; i < 6; ++i) {
      gamma[i] = gamma1[i] * ll;
    }
  }

  const int maxShellPairType = maxShellType * maxShellType;
  ShellPairArrayTable answer(maxShellPairType);
  for (int shellPairType = 0; shellPairType < maxShellPairType;
       ++shellPairType) {
    const ShellPairArray& shellPairArray = inShellPairArrayTable[shellPairType];
    const std::size_t shellPairArraySize = shellPairArray.size();
    ShellPairArray tmp;
    tmp.reserve(shellPairArraySize);

    for (std::size_t shellPairIndex = 0; shellPairIndex < shellPairArraySize;
         ++shellPairIndex) {
      const index_type shellIndexA = shellPairArray[shellPairIndex].shellIndex1;
      const index_type shellIndexB = shellPairArray[shellPairIndex].shellIndex2;

      const TlPosition posA = orbitalInfo.getPosition(shellIndexA);
      const int numOfContractionsA =
          orbitalInfo.getCgtoContraction(shellIndexA);
      const TlPosition posB = orbitalInfo.getPosition(shellIndexB);
      const int numOfContractionsB =
          orbitalInfo.getCgtoContraction(shellIndexB);
      const double AB2 = posB.squareDistanceFrom(posA);

      double judge = 0.0;
      for (int pgtoIndexA = 0; pgtoIndexA < numOfContractionsA; ++pgtoIndexA) {
        const double coefA =
            orbitalInfo.getCoefficient(shellIndexA, pgtoIndexA);
        const double zetaA = orbitalInfo.getExponent(shellIndexA, pgtoIndexA);

        for (int pgtoIndexB = 0; pgtoIndexB < numOfContractionsB;
             ++pgtoIndexB) {
          const double coefB =
              orbitalInfo.getCoefficient(shellIndexB, pgtoIndexB);
          const double zetaB = orbitalInfo.getExponent(shellIndexB, pgtoIndexB);

          const double coefAB = coefA * coefB;
          const double zetaAB = zetaA * zetaB;
          const double zetaA_B = zetaA + zetaB;
          for (int i = 0; i < 6; ++i) {
            const double zetaAgamma = zetaA * gamma[i];
            const double zetaBgamma = zetaB * gamma[i];
            const double param = zetaAB + zetaAgamma + zetaBgamma;
            const double term1 = CONTRIBUTE_COEF / (std::sqrt(zetaA_B) * param);
            const double term2 = std::exp(-zetaAB * gamma[i] * AB2 / param);
            const double T = zetaAB * zetaAB * AB2 / (zetaA_B * param);
            double term3 = 0.0;
            FmT.getFmT(0, T, &term3);
            judge += coefAB * coef[i] * term1 * term2 * term3;
          }
        }
      }

      if (std::fabs(judge) > cutoffThreshold) {
        tmp.push_back(shellPairArray[shellPairIndex]);

#pragma omp critical(DfTaskCtrl__selectShellPairArrayTableByDensity_add_alive)
        { ++(this->cutoffAlive_density_[shellPairType]); }
      }

#pragma omp critical(DfTaskCtrl__selectShellPairArrayTableByDensity_add_all)
      { ++(this->cutoffAll_density_[shellPairType]); }
    }

    // swap technique
    ShellPairArray(tmp).swap(tmp);

    answer[shellPairType] = tmp;
  }

  return answer;
}

DfTaskCtrl::ShellPairArrayTable DfTaskCtrl::selectShellPairArrayTableByDensity(
    const ShellPairArrayTable& inShellPairArrayTable,
    const TlOrbitalInfoObject& orbitalInfo1,
    const TlOrbitalInfoObject& orbitalInfo2) {
  const int maxShellType1 = orbitalInfo1.getMaxShellType();
  const int maxShellType2 = orbitalInfo2.getMaxShellType();
  assert(inShellPairArrayTable.size() ==
         static_cast<std::size_t>(maxShellType1 * maxShellType2));

  const double cutoffThreshold = this->cutoffEpsilon_density_;
  const double CONTRIBUTE_COEF = 2.0 * std::pow(TlMath::PI(), 2.5);

  TlFmt& FmT = TlFmt::getInstance();
  static const double coef[6] = {-0.017450254, 0.132520568,  -0.047915444,
                                 0.792267596,  -0.583721015, 0.697593555};
  // l = 1.0のときのgamma値
  static const double gamma1[6] = {0.03,        0.081722098, 0.222616709,
                                   0.606423483, 1.651939972, 4.5};
  std::vector<double> gamma(6);
  {
    const double inv_lprime = 1.0 / this->lengthScaleParameter_;
    const double ll = inv_lprime * inv_lprime;
    for (int i = 0; i < 6; ++i) {
      gamma[i] = gamma1[i] * ll;
    }
  }

  const int maxShellPairType = maxShellType1 * maxShellType2;
  ShellPairArrayTable answer(maxShellPairType);
  for (int shellPairType = 0; shellPairType < maxShellPairType;
       ++shellPairType) {
    const ShellPairArray& shellPairArray = inShellPairArrayTable[shellPairType];
    const std::size_t shellPairArraySize = shellPairArray.size();
    ShellPairArray tmp;
    tmp.reserve(shellPairArraySize);

    for (std::size_t shellPairIndex = 0; shellPairIndex < shellPairArraySize;
         ++shellPairIndex) {
      const index_type shellIndexA = shellPairArray[shellPairIndex].shellIndex1;
      const index_type shellIndexB = shellPairArray[shellPairIndex].shellIndex2;

      const TlPosition posA = orbitalInfo1.getPosition(shellIndexA);
      const int numOfContractionsA =
          orbitalInfo1.getCgtoContraction(shellIndexA);
      const TlPosition posB = orbitalInfo2.getPosition(shellIndexB);
      const int numOfContractionsB =
          orbitalInfo2.getCgtoContraction(shellIndexB);
      const double AB2 = posB.squareDistanceFrom(posA);

      double judge = 0.0;
      for (int pgtoIndexA = 0; pgtoIndexA < numOfContractionsA; ++pgtoIndexA) {
        const double coefA =
            orbitalInfo1.getCoefficient(shellIndexA, pgtoIndexA);
        const double zetaA = orbitalInfo1.getExponent(shellIndexA, pgtoIndexA);

        for (int pgtoIndexB = 0; pgtoIndexB < numOfContractionsB;
             ++pgtoIndexB) {
          const double coefB =
              orbitalInfo2.getCoefficient(shellIndexB, pgtoIndexB);
          const double zetaB =
              orbitalInfo2.getExponent(shellIndexB, pgtoIndexB);

          const double coefAB = coefA * coefB;
          const double zetaAB = zetaA * zetaB;
          const double zetaA_B = zetaA + zetaB;
          for (int i = 0; i < 6; ++i) {
            const double zetaAgamma = zetaA * gamma[i];
            const double zetaBgamma = zetaB * gamma[i];
            const double param = zetaAB + zetaAgamma + zetaBgamma;
            const double term1 = CONTRIBUTE_COEF / (std::sqrt(zetaA_B) * param);
            const double term2 = std::exp(-zetaAB * gamma[i] * AB2 / param);
            const double T = zetaAB * zetaAB * AB2 / (zetaA_B * param);
            double term3 = 0.0;
            FmT.getFmT(0, T, &term3);
            judge += coefAB * coef[i] * term1 * term2 * term3;
          }
        }
      }

      if (std::fabs(judge) > cutoffThreshold) {
        tmp.push_back(shellPairArray[shellPairIndex]);

#pragma omp critical(DfTaskCtrl__selectShellPairArrayTableByDensity_add_alive)
        { ++(this->cutoffAlive_density_[shellPairType]); }
      }

#pragma omp critical(DfTaskCtrl__selectShellPairArrayTableByDensity_add_all)
      { ++(this->cutoffAll_density_[shellPairType]); }
    }

    // swap technique
    ShellPairArray(tmp).swap(tmp);

    answer[shellPairType] = tmp;
  }

  return answer;
}

TlSparseSymmetricMatrix DfTaskCtrl::makeSchwarzTable(
    const TlOrbitalInfoObject& orbitalInfo, DfEriEngine* pEngine) {
  assert(pEngine != NULL);

  const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();
  TlSparseSymmetricMatrix schwarz(maxShellIndex);

  DfEriEngine engine;
  engine.setPrimitiveLevelThreshold(0.0);

  for (index_type shellIndexP = 0; shellIndexP < maxShellIndex;) {
    const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
    const int maxStepsP = 2 * shellTypeP + 1;

    for (index_type shellIndexQ = 0; shellIndexQ < maxShellIndex;) {
      const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
      const int maxStepsQ = 2 * shellTypeQ + 1;

      const DfEriEngine::AngularMomentum2 queryPQ(0, 0, shellTypeP, shellTypeQ);
      // const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
      //                                                       shellIndexP,
      //                                                       shellIndexQ,
      //                                                       0.0);
      // pEngine->calc(queryPQ, queryPQ, PQ, PQ);
      pEngine->calc(0, orbitalInfo, shellIndexP, 0, orbitalInfo, shellIndexQ, 0,
                    orbitalInfo, shellIndexP, 0, orbitalInfo, shellIndexQ);

      double maxValue = 0.0;
      const int maxIndex = maxStepsP * maxStepsQ;
      for (int index = 0; index < maxIndex; ++index) {
        maxValue = std::max(maxValue, std::fabs(engine.WORK[index]));
      }
      schwarz.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));

      shellIndexQ += maxStepsQ;
    }
    shellIndexP += maxStepsP;
  }

  return schwarz;
}

// int DfTaskCtrl::getShellQuartetType(const TlOrbitalInfoObject& orbitalInfo,
//                                     const int shellTypeP,
//                                     const int shellTypeQ,
//                                     const int shellTypeR,
//                                     const int shellTypeS)
// {
//     const int maxShellType = orbitalInfo.getMaxShellType();
//     // schwarz cutoff
//     const int shellQuartetType =
//         ((shellTypeP * maxShellType + shellTypeQ) * maxShellType +
//         shellTypeR) * maxShellType + shellTypeS;

//     return shellQuartetType;
// }

int DfTaskCtrl::getShellQuartetType(const TlOrbitalInfoObject& orbitalInfo_PQ,
                                    const TlOrbitalInfoObject& orbitalInfo_RS,
                                    const int shellTypeP, const int shellTypeQ,
                                    const int shellTypeR,
                                    const int shellTypeS) {
  const int maxShellType_PQ = orbitalInfo_PQ.getMaxShellType();
  const int maxShellType_RS = orbitalInfo_RS.getMaxShellType();
  // schwarz cutoff
  const int shellQuartetType =
      ((shellTypeP * maxShellType_PQ + shellTypeQ) * maxShellType_PQ +
       shellTypeR) *
          maxShellType_RS +
      shellTypeS;

  return shellQuartetType;
}

bool DfTaskCtrl::isAliveBySchwarzCutoff(
    const index_type shellIndexP, const index_type shellIndexQ,
    const index_type shellIndexR, const index_type shellIndexS,
    const int shellQuartetType, const TlSparseSymmetricMatrix& schwarzTable,
    const double threshold) {
  bool answer = false;

  const double sqrt_pqpq = schwarzTable.get(shellIndexP, shellIndexQ);
  const double sqrt_rsrs = schwarzTable.get(shellIndexR, shellIndexS);

  if ((sqrt_pqpq * sqrt_rsrs) >= threshold) {
    answer = true;

#pragma omp critical(DfTaskCtrl__isAliveBySchwarzCutoff_alive)
    { ++(this->cutoffAlive_schwarz_[shellQuartetType]); }
  }

#pragma omp critical(DfTaskCtrl__isAliveBySchwarzCutoff_all)
  { ++(this->cutoffAll_schwarz_[shellQuartetType]); }

  return answer;
}

bool DfTaskCtrl::isAliveBySchwarzCutoff(
    const index_type shellIndexP, const index_type shellIndexQ,
    const index_type shellIndexR, const index_type shellIndexS,
    const int shellQuartetType, const TlSparseSymmetricMatrix& schwarzTable_PQ,
    const TlSparseSymmetricMatrix& schwarzTable_RS, const double threshold) {
  bool answer = false;

  const double sqrt_pqpq = schwarzTable_PQ.get(shellIndexP, shellIndexQ);
  const double sqrt_rsrs = schwarzTable_RS.get(shellIndexR, shellIndexS);

  if ((sqrt_pqpq * sqrt_rsrs) >= threshold) {
    answer = true;

#pragma omp critical(DfTaskCtrl__isAliveBySchwarzCutoff_alive)
    { ++(this->cutoffAlive_schwarz_[shellQuartetType]); }
  }

#pragma omp critical(DfTaskCtrl__isAliveBySchwarzCutoff_all)
  { ++(this->cutoffAll_schwarz_[shellQuartetType]); }

  return answer;
}

void DfTaskCtrl::prescreeningReport() {
  const int maxShellType = this->maxShellType_;
  static const char typeStr2[][3] = {"SS", "SP", "SD", "SF", "SG", "PS", "PP",
                                     "PD", "PF", "PG", "DS", "DP", "DD", "DF",
                                     "DG", "FS", "FP", "FD", "FF", "FG", "GS",
                                     "GP", "GD", "GF", "GG"};

  // cutoff report for Epsilon1
  bool hasCutoff_density = false;
  for (int shellTypeA = 0;
       ((hasCutoff_density == false) && (shellTypeA < maxShellType));
       ++shellTypeA) {
    for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
      const int shellPairType = shellTypeA * maxShellType + shellTypeB;
      if (this->cutoffAll_density_[shellPairType] != 0) {
        hasCutoff_density = true;
        break;
      }
    }
  }
  if (hasCutoff_density == true) {
    this->log_.info("pre-screening: density cutoff report");
    this->log_.info(
        TlUtils::format("epsilon(density): %e", this->cutoffEpsilon_density_));
    this->log_.info(TlUtils::format("length scale parameter: %f",
                                    this->lengthScaleParameter_));
    this->log_.info("type: alive / all (ratio)");
    for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
      for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
        const int shellPairType = shellTypeA * maxShellType + shellTypeB;

        if (this->cutoffAll_density_[shellPairType] > 0) {
          const double ratio =
              (double)this->cutoffAlive_density_[shellPairType] /
              (double)this->cutoffAll_density_[shellPairType] * 100.0;
          this->log_.info(TlUtils::format(
              " %2s: %12ld / %12ld (%6.2f%%)", typeStr2[shellPairType],
              this->cutoffAlive_density_[shellPairType],
              this->cutoffAll_density_[shellPairType], ratio));
        }
      }
    }
    this->log_.info("\n");
  }

  // cutoff report for distribition
  bool hasCutoff_distribution = false;
  for (int shellTypeA = 0;
       ((hasCutoff_distribution == false) && (shellTypeA < maxShellType));
       ++shellTypeA) {
    for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
      const int shellPairType = shellTypeA * maxShellType + shellTypeB;
      if (this->cutoffAll_distribution_[shellPairType] != 0) {
        hasCutoff_distribution = true;
        break;
      }
    }
  }
  if (hasCutoff_distribution == true) {
    this->log_.info("prescreening: distribution cutoff report");
    this->log_.info(TlUtils::format("epsilon(distribution): %e",
                                    this->cutoffEpsilon_distribution_));
    this->log_.info("type: alive / all (ratio)");
    for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
      for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
        const int shellPairType = shellTypeA * maxShellType + shellTypeB;

        if (this->cutoffAll_distribution_[shellPairType] > 0) {
          const double ratio =
              (double)this->cutoffAlive_distribution_[shellPairType] /
              (double)this->cutoffAll_distribution_[shellPairType] * 100.0;
          this->log_.info(TlUtils::format(
              " %2s: %12ld / %12ld (%6.2f%%)", typeStr2[shellPairType],
              this->cutoffAlive_distribution_[shellPairType],
              this->cutoffAll_distribution_[shellPairType], ratio));
        }
      }
    }
    this->log_.info("\n");
  }
}

void DfTaskCtrl::cutoffReport() {
  const int maxShellType = this->maxShellType_;
  std::vector<std::string> typeStr4(maxShellType * maxShellType * maxShellType *
                                    maxShellType);
  {
    static const char typeChar[] = "SPDFG";
    std::string tmp(4, 'X');
    int index = 0;
    for (int i = 0; i < maxShellType; ++i) {
      tmp[0] = typeChar[i];
      for (int j = 0; j < maxShellType; ++j) {
        tmp[1] = typeChar[j];
        for (int k = 0; k < maxShellType; ++k) {
          tmp[2] = typeChar[k];
          for (int l = 0; l < maxShellType; ++l) {
            tmp[3] = typeChar[l];
            typeStr4[index] = tmp;
            ++index;
          }
        }
      }
    }
  }

  // cutoff report for schwarz
  bool hasCutoffSchwarz = false;
  for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
    for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
      const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
      for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
        const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
        for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
          const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
          if (this->cutoffAll_schwarz_[shellTypeABCD] != 0) {
            hasCutoffSchwarz = true;
            break;
          }
        }
      }
    }
  }
  if (hasCutoffSchwarz == true) {
    this->log_.info("schwarz cutoff report");
    this->log_.info(TlUtils::format("threshold: %e", this->cutoffThreshold_));
    this->log_.info("type: alive / all (ratio)");
    for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
      for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
        const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
        for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
          const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
          for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
            const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;

            if (this->cutoffAll_schwarz_[shellTypeABCD] > 0) {
              const double ratio =
                  (double)this->cutoffAlive_schwarz_[shellTypeABCD] /
                  (double)this->cutoffAll_schwarz_[shellTypeABCD] * 100.0;
              this->log_.info(TlUtils::format(
                  " %4s: %12ld / %12ld (%6.2f%%)",
                  typeStr4[shellTypeABCD].c_str(),
                  this->cutoffAlive_schwarz_[shellTypeABCD],
                  this->cutoffAll_schwarz_[shellTypeABCD], ratio));
            }
          }
        }
      }
    }
  }
}

std::size_t DfTaskCtrl::getTotalCalcAmount2(
    const TlOrbitalInfoObject& orbitalInfo,
    const ShellArrayTable& shellArrayTable) {
  std::size_t answer = 0;
  const int maxShellType = orbitalInfo.getMaxShellType();
  for (int shellTypeP = 0; shellTypeP < maxShellType; ++shellTypeP) {
    const int maxStepsP = 2 * shellTypeP + 1;
    const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
    const size_t shellArraySizeP = shellArrayP.size();

    for (int shellTypeQ = 0; shellTypeQ < maxShellType; ++shellTypeQ) {
      const int maxStepsQ = 2 * shellTypeQ + 1;

      for (std::size_t shellArrayIndexP = 0; shellArrayIndexP < shellArraySizeP;
           ++shellArrayIndexP) {
        const index_type shellIndexP = shellArrayP[shellArrayIndexP];

        const ShellArray& shellArrayQ = shellArrayTable[shellTypeQ];
        ShellArray::const_iterator qItEnd = std::upper_bound(
            shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
        const std::size_t shellArraySizeQ =
            std::distance(shellArrayQ.begin(), qItEnd);

        answer += shellArraySizeQ * maxStepsP * maxStepsQ;
      }
    }
  }

  return answer;
}

std::size_t DfTaskCtrl::getTotalCalcAmount2(
    const TlOrbitalInfoObject& orbitalInfo1,
    const TlOrbitalInfoObject& orbitalInfo2,
    const ShellArrayTable& shellArrayTable1,
    const ShellArrayTable& shellArrayTable2) {
  std::size_t answer = 0;
  const int maxShellType1 = orbitalInfo1.getMaxShellType();
  const int maxShellType2 = orbitalInfo2.getMaxShellType();
  for (int shellTypeP = 0; shellTypeP < maxShellType1; ++shellTypeP) {
    const int maxStepsP = 2 * shellTypeP + 1;
    const ShellArray& shellArrayP = shellArrayTable1[shellTypeP];
    const size_t shellArraySizeP = shellArrayP.size();

    for (int shellTypeQ = 0; shellTypeQ < maxShellType2; ++shellTypeQ) {
      const int maxStepsQ = 2 * shellTypeQ + 1;

      for (std::size_t shellArrayIndexP = 0; shellArrayIndexP < shellArraySizeP;
           ++shellArrayIndexP) {
        // const index_type shellIndexP = shellArrayP[shellArrayIndexP];

        const ShellArray& shellArrayQ = shellArrayTable2[shellTypeQ];
        const std::size_t shellArraySizeQ = shellArrayQ.size();

        answer += shellArraySizeQ * maxStepsP * maxStepsQ;
      }
    }
  }

  return answer;
}

std::size_t DfTaskCtrl::getTotalCalcAmount2(
    const TlOrbitalInfoObject& orbitalInfo,
    const ShellArrayTable& shellArrayTable, const DistributedCutoffTable& dct) {
  std::size_t answer = 0;
  const int maxShellType = orbitalInfo.getMaxShellType();
  for (int shellTypeP = 0; shellTypeP < maxShellType; ++shellTypeP) {
    const int maxStepsP = 2 * shellTypeP + 1;
    const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
    const size_t shellArraySizeP = shellArrayP.size();

    for (int shellTypeQ = 0; shellTypeQ < maxShellType; ++shellTypeQ) {
      const int maxStepsQ = 2 * shellTypeQ + 1;

      for (std::size_t shellArrayIndexP = 0; shellArrayIndexP < shellArraySizeP;
           ++shellArrayIndexP) {
        const index_type shellIndexP = shellArrayP[shellArrayIndexP];

        DistributedCutoffTable::const_iterator dctItQ = dct.find(shellIndexP);
        assert(dctItQ != dct.end());
        const ShellArray shellArrayQ = dctItQ->second[shellTypeQ];

        ShellArray::const_iterator qItEnd = std::upper_bound(
            shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
        const std::size_t shellArraySizeQ =
            std::distance(shellArrayQ.begin(), qItEnd);

        answer += shellArraySizeQ * maxStepsP * maxStepsQ;
      }
    }
  }

  return answer;
}

std::size_t DfTaskCtrl::getTotalCalcAmount2(
    const TlOrbitalInfoObject& orbitalInfo1,
    const TlOrbitalInfoObject& orbitalInfo2,
    const ShellArrayTable& shellArrayTable1,
    const ShellArrayTable& shellArrayTable2,
    const DistributedCutoffTable& dct) {
  std::size_t answer = 0;
  const int maxShellType1 = orbitalInfo1.getMaxShellType();
  const int maxShellType2 = orbitalInfo2.getMaxShellType();
  for (int shellTypeP = 0; shellTypeP < maxShellType1; ++shellTypeP) {
    const int maxStepsP = 2 * shellTypeP + 1;
    const ShellArray& shellArrayP = shellArrayTable1[shellTypeP];
    const size_t shellArraySizeP = shellArrayP.size();

    for (int shellTypeQ = 0; shellTypeQ < maxShellType2; ++shellTypeQ) {
      const int maxStepsQ = 2 * shellTypeQ + 1;

      for (std::size_t shellArrayIndexP = 0; shellArrayIndexP < shellArraySizeP;
           ++shellArrayIndexP) {
        const index_type shellIndexP = shellArrayP[shellArrayIndexP];

        DistributedCutoffTable::const_iterator dctItQ = dct.find(shellIndexP);
        assert(dctItQ != dct.end());
        const ShellArray shellArrayQ = dctItQ->second[shellTypeQ];
        const std::size_t shellArraySizeQ = shellArrayQ.size();

        answer += shellArraySizeQ * maxStepsP * maxStepsQ;
      }
    }
  }

  return answer;
}

std::size_t DfTaskCtrl::getTotalCalcAmount4(
    const TlOrbitalInfoObject& orbitalInfo,
    const ShellPairArrayTable& shellPairArrayTable,
    const DistributedCutoffTable& dct) {
  std::size_t answer = 0;

  // calc prescreened size
  const int maxShellType = orbitalInfo.getMaxShellType();
  for (int shellTypeP = 0; shellTypeP < maxShellType; ++shellTypeP) {
    const int maxStepsP = 2 * shellTypeP + 1;

    for (int shellTypeR = 0; shellTypeR < maxShellType; ++shellTypeR) {
      const int maxStepsR = 2 * shellTypeR + 1;

      const int shellPairType_PR = shellTypeP * maxShellType + shellTypeR;
      const ShellPairArray& shellPairArray_PR =
          shellPairArrayTable[shellPairType_PR];
      const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();

      for (std::size_t prIndex = 0; prIndex < numOfShellPairArray_PR;
           ++prIndex) {
        const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
        const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;

        for (int shellTypeQ = 0; shellTypeQ < maxShellType; ++shellTypeQ) {
          const int maxStepsQ = 2 * shellTypeQ + 1;

          // const ShellArray shellArrayQ = dct[shellIndexP][shellTypeQ];
          DistributedCutoffTable::const_iterator dctItP = dct.find(shellIndexP);
          assert(dctItP != dct.end());
          const ShellArray shellArrayQ = dctItP->second[shellTypeQ];

          ShellArray::const_iterator qItEnd = std::upper_bound(
              shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
          const std::size_t shellArraySizeQ =
              std::distance(shellArrayQ.begin(), qItEnd);

          for (int shellTypeS = 0; shellTypeS < maxShellType; ++shellTypeS) {
            const int maxStepsS = 2 * shellTypeS + 1;

            // const ShellArray shellArrayS = dct[shellIndexR][shellTypeS];
            DistributedCutoffTable::const_iterator dctItR =
                dct.find(shellIndexR);
            assert(dctItR != dct.end());
            const ShellArray shellArrayS = dctItR->second[shellTypeS];

            for (std::size_t shellArrayIndexQ = 0;
                 shellArrayIndexQ < shellArraySizeQ; ++shellArrayIndexQ) {
              const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];

              const index_type maxShellIndexS =
                  (shellIndexP == shellIndexR) ? shellIndexQ : shellIndexR;
              ShellArray::const_iterator sItEnd = std::upper_bound(
                  shellArrayS.begin(), shellArrayS.end(), maxShellIndexS);
              const std::size_t shellArraySizeS =
                  std::distance(shellArrayS.begin(), sItEnd);
              answer += shellArraySizeS * maxStepsP * maxStepsQ * maxStepsR *
                        maxStepsS;
            }
          }
        }
      }
    }
  }

  return answer;
}

std::size_t DfTaskCtrl::getTotalCalcAmount4(
    const TlOrbitalInfoObject& orbitalInfo_PQ,
    const TlOrbitalInfoObject& orbitalInfo_RS,
    const ShellPairArrayTable& shellPairArrayTable,
    const DistributedCutoffTable& dct_PQ,
    const DistributedCutoffTable& dct_RS) {
  std::size_t answer = 0;

  // calc prescreened size
  const int maxShellType_PQ = orbitalInfo_PQ.getMaxShellType();
  const int maxShellType_RS = orbitalInfo_RS.getMaxShellType();
  for (int shellTypeP = 0; shellTypeP < maxShellType_PQ; ++shellTypeP) {
    const int maxStepsP = 2 * shellTypeP + 1;

    for (int shellTypeR = 0; shellTypeR < maxShellType_RS; ++shellTypeR) {
      const int maxStepsR = 2 * shellTypeR + 1;

      const int shellPairType_PR = shellTypeP * maxShellType_PQ + shellTypeR;
      const ShellPairArray& shellPairArray_PR =
          shellPairArrayTable[shellPairType_PR];
      const std::size_t numOfShellPairArray_PR = shellPairArray_PR.size();

      for (std::size_t prIndex = 0; prIndex < numOfShellPairArray_PR;
           ++prIndex) {
        const index_type shellIndexP = shellPairArray_PR[prIndex].shellIndex1;
        const index_type shellIndexR = shellPairArray_PR[prIndex].shellIndex2;

        for (int shellTypeQ = 0; shellTypeQ < maxShellType_PQ; ++shellTypeQ) {
          const int maxStepsQ = 2 * shellTypeQ + 1;

          // const ShellArray shellArrayQ = dct[shellIndexP][shellTypeQ];
          DistributedCutoffTable::const_iterator dctItP =
              dct_PQ.find(shellIndexP);
          assert(dctItP != dct_PQ.end());
          const ShellArray shellArrayQ = dctItP->second[shellTypeQ];

          ShellArray::const_iterator qItEnd = std::upper_bound(
              shellArrayQ.begin(), shellArrayQ.end(), shellIndexP);
          const std::size_t shellArraySizeQ =
              std::distance(shellArrayQ.begin(), qItEnd);

          for (int shellTypeS = 0; shellTypeS < maxShellType_RS; ++shellTypeS) {
            const int maxStepsS = 2 * shellTypeS + 1;

            // const ShellArray shellArrayS = dct[shellIndexR][shellTypeS];
            DistributedCutoffTable::const_iterator dctItR =
                dct_RS.find(shellIndexR);
            assert(dctItR != dct_RS.end());
            const ShellArray shellArrayS = dctItR->second[shellTypeS];

            for (std::size_t shellArrayIndexQ = 0;
                 shellArrayIndexQ < shellArraySizeQ; ++shellArrayIndexQ) {
              const index_type shellIndexQ = shellArrayQ[shellArrayIndexQ];

              const index_type maxShellIndexS =
                  (shellIndexP == shellIndexR) ? shellIndexQ : shellIndexR;
              ShellArray::const_iterator sItEnd = std::upper_bound(
                  shellArrayS.begin(), shellArrayS.end(), maxShellIndexS);
              const std::size_t shellArraySizeS =
                  std::distance(shellArrayS.begin(), sItEnd);
              answer += shellArraySizeS * maxStepsP * maxStepsQ * maxStepsR *
                        maxStepsS;
            }
          }
        }
      }
    }
  }

  return answer;
}
