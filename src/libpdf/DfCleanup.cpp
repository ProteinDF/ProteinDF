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

#include <string>
#include "DfCleanup.h"
#include "TlFile.h"

DfCleanup::DfCleanup(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
}


DfCleanup::~DfCleanup()
{
}


void DfCleanup::cleanup()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->cleanup(RUN_RKS, this->m_nIteration);
        break;

    case METHOD_UKS:
        this->cleanup(RUN_UKS_ALPHA, this->m_nIteration);
        this->cleanup(RUN_UKS_BETA, this->m_nIteration);
        break;

    default:
        break;
    }
}


void DfCleanup::cleanup(const RUN_TYPE runType, const int iteration)
{
    this->cleanupFxc(runType, iteration);
    this->cleanupFpq(runType, iteration);
    this->cleanupFprime(runType, iteration);
    this->cleanupCprime(runType, iteration);
    this->cleanupC(runType, iteration);
    this->cleanupP(runType, iteration);
}


void DfCleanup::cleanupFxc(const RUN_TYPE runType, const int iteration)
{
    std::string path = "";
    if (this->m_bIsUpdateXC == true) {
        // current iterationの1つ前を削除
        path = DfObject::getFxcMatrixPath(runType, iteration -1);
    } else {
        // current iterationを削除
        path = DfObject::getFxcMatrixPath(runType, iteration);
    }

    this->log_.info(TlUtils::format("remove: %s", path.c_str()));
    this->matrixCache_.erase(path);
    TlFile::remove(path);
}


void DfCleanup::cleanupHFx(const RUN_TYPE runType, const int iteration)
{
    std::string path = "";
    if (this->m_bIsUpdateXC == true) {
        // current iterationの1つ前を削除
        path = DfObject::getHFxMatrixPath(runType, iteration -1);
    } else {
        // current iterationを削除
        path = DfObject::getHFxMatrixPath(runType, iteration);
    }

    this->log_.info(TlUtils::format("remove: %s", path.c_str()));
    this->matrixCache_.erase(path);
    TlFile::remove(path);
}

void DfCleanup::cleanupFpq(const RUN_TYPE runType, const int iteration)
{
    int itr = iteration;
    if (this->m_bDiskUtilization == false) {
        --itr;
    }

    const std::string path = DfObject::getFpqMatrixPath(runType, itr);

    this->log_.info(TlUtils::format("remove: %s", path.c_str()));
    this->matrixCache_.erase(path);
    TlFile::remove(path);
}


void DfCleanup::cleanupFprime(const RUN_TYPE runType, const int iteration)
{
    // current iterationを削除
    const std::string path = DfObject::getFprimeMatrixPath(runType, iteration);

    this->log_.info(TlUtils::format("remove: %s", path.c_str()));
    this->matrixCache_.erase(path);
    TlFile::remove(path);
}


void DfCleanup::cleanupCprime(const RUN_TYPE runType, const int iteration)
{
    // 1つ前はDfDmatrixのチェックで利用
    const std::string path = DfObject::getCprimeMatrixPath(runType, iteration -1);

    this->log_.info(TlUtils::format("remove: %s", path.c_str()));
    this->matrixCache_.erase(path);
    TlFile::remove(path);
}


void DfCleanup::cleanupC(const RUN_TYPE runType, const int iteration)
{
    // current iterationの1つ前を削除(中断からの継続のため)
    const std::string path = DfObject::getCMatrixPath(runType, iteration -1);

    this->log_.info(TlUtils::format("remove: %s", path.c_str()));
    this->matrixCache_.erase(path);
    TlFile::remove(path);
}


void DfCleanup::cleanupP(const RUN_TYPE runType, const int iteration)
{
    int itr = iteration -1;
    if (this->m_bDiskUtilization == false) {
        // update法
        --itr;
    }
    // current iterationの1つ前を削除(current iterationは次回のF作成に使用)
    const std::string PpqPath = DfObject::getPpqMatrixPath(runType, itr);
    const std::string P1pqPath = DfObject::getP1pqMatrixPath(itr);
    const std::string P2pqPath = DfObject::getP2pqMatrixPath(itr);

    this->log_.info(TlUtils::format("remove: %s", PpqPath.c_str()));
    this->matrixCache_.erase(PpqPath);
    TlFile::remove(PpqPath);

    this->log_.info(TlUtils::format("remove: %s", P1pqPath.c_str()));
    this->matrixCache_.erase(P1pqPath);
    TlFile::remove(P1pqPath);

    this->log_.info(TlUtils::format("remove: %s", P2pqPath.c_str()));
    this->matrixCache_.erase(P2pqPath);
    TlFile::remove(P2pqPath);
}
