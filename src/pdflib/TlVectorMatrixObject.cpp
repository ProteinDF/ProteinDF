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

#include <iostream>
#include "TlVectorMatrixObject.h"
#include "TlMemManager.h"
#include "TlLogging.h"
#include "TlUtils.h"

TlVectorMatrixObject::TlVectorMatrixObject(const index_type numOfVectors,
                                           const index_type sizeOfVector,
                                           const int numOfSubunits,
                                           const int subunitID,
                                           const bool isUsingMemManager)
    : numOfVectors_(0),
      sizeOfVector_(0),
      numOfSubunits_(numOfSubunits),
      subunitID_(subunitID),
      numOfLocalVectors_(0),
      reservedVectorSize_(0),
      isUsingMemManager_(isUsingMemManager) {
    this->resize(numOfVectors, sizeOfVector);
}


TlVectorMatrixObject::~TlVectorMatrixObject()
{
    const index_type numOfLocalVectors = this->numOfLocalVectors_;

    if (this->isUsingMemManager_ == true) {
        TlMemManager& rMemManager = TlMemManager::getInstance();
        const index_type reservedVectorSize = this->reservedVectorSize_;
        for (index_type i = 0; i < numOfLocalVectors; ++i) {
            rMemManager.deallocate((char*)this->data_[i], sizeof(double) * reservedVectorSize);
            this->data_[i] = NULL;
        }
    } else {
        for (index_type i = 0; i < numOfLocalVectors; ++i) {
            delete[] this->data_[i];
            this->data_[i] = NULL;
        }
    }
    this->data_.clear();
}


void TlVectorMatrixObject::resize(const index_type newNumOfVectors,
                                  const index_type newSizeOfVector)
                                  
{
    index_type prevNumOfLocalVectors = this->numOfLocalVectors_;
    const div_t turns = std::div(newNumOfVectors, this->numOfSubunits_);
    index_type newNumOfLocalVectors = turns.quot;
    if (this->subunitID_ < turns.rem) {
        newNumOfLocalVectors += 1;
    }
    this->numOfLocalVectors_ = newNumOfLocalVectors;

    if (newNumOfLocalVectors > prevNumOfLocalVectors) {
        // ベクトル数を拡大
        this->data_.resize(newNumOfLocalVectors, NULL);
        const index_type reservedVectorSize = this->reservedVectorSize_;
        for (index_type i = prevNumOfLocalVectors; i < newNumOfLocalVectors; ++i) {
            double* pNew = NULL;
            if (this->isUsingMemManager_ == true) {
                TlMemManager& rMemManager = TlMemManager::getInstance();
                pNew = reinterpret_cast<double*>(rMemManager.allocate(sizeof(double) * reservedVectorSize));
            } else {
                pNew = new double[reservedVectorSize];
                this->data_[i] = pNew;
            }
        }
    } else if (newNumOfLocalVectors < prevNumOfLocalVectors) {
        // ベクトル数を縮小
        if (this->isUsingMemManager_ == true) {
            TlMemManager& rMemManager = TlMemManager::getInstance();
            const index_type reservedVectorSize = this->reservedVectorSize_;
            for (index_type i = newNumOfLocalVectors; i < prevNumOfLocalVectors; ++i) {
                rMemManager.deallocate((char*)this->data_[i], sizeof(double) * reservedVectorSize);
                this->data_[i] = NULL;
            }
        } else {
            for (index_type i = newNumOfLocalVectors; i < prevNumOfLocalVectors; ++i) {
                delete[] this->data_[i];
                this->data_[i] = NULL;
            }
        }
        this->data_.resize(newNumOfLocalVectors);
    }
    
    this->numOfVectors_ = newNumOfVectors;
    this->reserveVectorSize(newSizeOfVector);
    this->sizeOfVector_ = newSizeOfVector;
}


void TlVectorMatrixObject::reserveVectorSize(index_type newReservedVectorSize) {
    newReservedVectorSize = std::max(this->sizeOfVector_, newReservedVectorSize);
    const index_type prevReservedVectorSize = this->reservedVectorSize_;

    if (prevReservedVectorSize < newReservedVectorSize) {
        this->reservedVectorSize_ = newReservedVectorSize;
        const index_type numOfLocalVectors = this->numOfLocalVectors_;
        for (index_type i = 0; i < numOfLocalVectors; ++i) {
            double* pNew = NULL;
            if (this->isUsingMemManager_ == true) {
                TlMemManager& rMemManager = TlMemManager::getInstance();
                pNew = reinterpret_cast<double*>(rMemManager.allocate(sizeof(double) * newReservedVectorSize));
            } else {
                pNew = new double[newReservedVectorSize];
            }

            for (index_type j = 0; j < newReservedVectorSize; ++j) {
                pNew[j] = 0.0;
            }
            
            if (this->data_[i] != NULL) {
                for (index_type j = 0; j < prevReservedVectorSize; ++j) {
                    pNew[j] = this->data_[i][j];
                }
                if (this->isUsingMemManager_ == true) {
                    TlMemManager& rMemManager = TlMemManager::getInstance();
                    rMemManager.deallocate((char*)this->data_[i], sizeof(double)*prevReservedVectorSize);
                } else {
                    delete[] this->data_[i];
                }
                this->data_[i] = NULL;
            }

            this->data_[i] = pNew;
        }
    }

    this->reservedVectorSize_ = newReservedVectorSize;
}


void TlVectorMatrixObject::set(const index_type vectorIndex,
                               const index_type index,
                               const double value)
{
    const div_t turns = std::div(vectorIndex, this->numOfSubunits_);
    if (turns.rem == this->subunitID_) {
        const index_type localVectorIndex = turns.quot;
        this->data_[localVectorIndex][index] = value;
    }
}


double TlVectorMatrixObject::get(const index_type vectorIndex,
                                 const index_type index) const
{
    assert((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_));
    assert((0 <= index) && (index < this->vectorSize_));

    double answer = 0.0;
    const div_t turns = std::div(vectorIndex, this->numOfSubunits_);
    if (turns.rem == this->subunitID_) {
        const index_type localVectorIndex = turns.quot;
        answer = this->data_[localVectorIndex][index];
    }

    return answer;
}

std::vector<double> TlVectorMatrixObject::getVector(const index_type vectorIndex) const
{
    const index_type vectorSize = this->sizeOfVector_;
    std::vector<double> answer;
    const div_t turns = std::div(vectorIndex, this->numOfSubunits_);
    if (turns.rem == this->subunitID_) {
        const index_type localVectorIndex = turns.quot;
        answer.resize(vectorSize);
        std::copy(this->data_[localVectorIndex],
                  this->data_[localVectorIndex] + vectorSize,
                  answer.begin());
    }

    return answer;
}


int TlVectorMatrixObject::getSubunitID(const index_type vectorIndex) const
{
    assert((0 <= vectorIndex) && (vectorIndex < this->numOfVectors_));
    const div_t turns = std::div(vectorIndex, this->numOfSubunits_);
    return turns.rem;
}


void TlVectorMatrixObject::save(const std::string& basename) const 
{
    const index_type numOfVectors = this->numOfVectors_;
    const index_type sizeOfVector = this->sizeOfVector_;

    std::ofstream ofs;
    const std::string path = TlUtils::format("%s.part%d.mat",
                                             basename.c_str(),
                                             this->subunitID_);
    ofs.open(path.c_str(), std::ofstream::out | std::ofstream::binary);

    // header
    ofs.write(reinterpret_cast<const char*>(&numOfVectors), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&sizeOfVector), sizeof(index_type));
    ofs.write(reinterpret_cast<const char*>(&(this->numOfSubunits_)), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&(this->subunitID_)), sizeof(int));
    // ofs.write((const char*)&numOfVectors, sizeof(index_type));
    // ofs.write((const char*)&sizeOfVector, sizeof(index_type));
    // ofs.write((const char*)&(this->numOfSubunits_), sizeof(int));
    // ofs.write((const char*)&(this->subunitID_), sizeof(int));

    // data
    const index_type numOfLocalVectors = this->data_.size();
    assert(numOfLocalVectors = this->numOfLocalVectors_);
    for (index_type i = 0; i < numOfLocalVectors; ++i) {
        ofs.write(reinterpret_cast<const char*>(&(this->data_[i][0])), sizeof(double) * sizeOfVector);
    }
    
    ofs.close();
}


void TlVectorMatrixObject::load(const std::string& basename)
{
    TlLogging& log = TlLogging::getInstance();

    std::ifstream ifs;
    const std::string path = TlUtils::format("%s.part%d.mat",
                                             basename.c_str(),
                                             this->subunitID_);
    ifs.open(path.c_str(), std::ifstream::in);
    if (ifs.good()) {
        // header
        index_type numOfVectors = 0;
        index_type sizeOfVector = 0;
        int allProcs = 0;
        int rank = 0;
        ifs.read((char*)&numOfVectors, sizeof(index_type));
        ifs.read((char*)&sizeOfVector, sizeof(index_type));
        ifs.read((char*)&allProcs, sizeof(int));
        ifs.read((char*)&rank, sizeof(int));
        this->resize(numOfVectors, sizeOfVector);
        assert((allProcs == this->allProcs_) && (rank == this->rank_));
        
        // data
        const div_t turns = std::div(numOfVectors, allProcs);
        index_type numOfLocalVectors = turns.quot;
        if (rank < turns.rem) {
            ++numOfLocalVectors;
        }
        
        for (int i = 0; i < numOfLocalVectors; ++i) {
            ifs.read((char*)&(this->data_[i][0]), sizeof(double) * sizeOfVector);
        }
    } else {
        log.error(TlUtils::format("cannot open file: %s", path.c_str()));
    }

    ifs.close();
}




