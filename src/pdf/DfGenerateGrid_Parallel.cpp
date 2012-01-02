#include <cassert>
#include "DfGenerateGrid_Parallel.h"
#include "DfXCFunctional.h"
#include "TlCommunicate.h"
#include "GridDataManager.h"
#include "TlUtils.h"
#include "TlFileMatrix.h"

DfGenerateGrid_Parallel::DfGenerateGrid_Parallel(TlSerializeData* pPdfParam)
    : DfGenerateGrid(pPdfParam)
{
}

DfGenerateGrid_Parallel::~DfGenerateGrid_Parallel()
{
}

void DfGenerateGrid_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfGenerateGrid::logger(str);
    }
}

void DfGenerateGrid_Parallel::makeTable()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfGenerateGrid::makeTable();
    }

    rComm.broadcast(this->maxRadii_);
}


void DfGenerateGrid_Parallel::generateGrid()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    
    // if (this->isMasterSlave_ == true) {
    //     this->generateGrid_MS();
    // } else {
        this->generateGrid_DC();
    // }

    // gather
    index_type numOfRowsOfGlobalGridMatrix = this->grdMat_.getNumOfRows();
    rComm.allReduce_SUM(numOfRowsOfGlobalGridMatrix);
    if (rComm.isMaster() == true) {
        const index_type numOfColsOfGlobalGridMatrix = this->grdMat_.getNumOfCols();
        TlFileMatrix grdMat(this->getGridMatrixPath(0),
                            numOfRowsOfGlobalGridMatrix,
                            numOfColsOfGlobalGridMatrix);

        index_type currentNumOfRows = 0;
        // from 0
        grdMat.setBlockMatrix(currentNumOfRows, 0,
                              this->grdMat_);
        currentNumOfRows += this->grdMat_.getNumOfRows();
        
        // from the others
        for (int i = 1; i < numOfProcs; ++i) {
            TlMatrix tmpGrdMat;
            rComm.receiveData(tmpGrdMat, i);
            grdMat.setBlockMatrix(currentNumOfRows, 0,
                                  tmpGrdMat);
            currentNumOfRows += tmpGrdMat.getNumOfRows();
        }
    } else {
        rComm.sendData(this->grdMat_, 0);
    }
}

void DfGenerateGrid_Parallel::generateGrid_DC()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int nEndAtomNumber = this->m_nNumOfAtoms;

    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();
    const int nRange = nEndAtomNumber;
    const int nInterval = (nRange + (nProc -1)) / nProc; // +(nProc-1) は余り対策
    const int nLocalStart = nRank * nInterval; // nProc = 0, 1, 2, ...
    const int nLocalEnd   = std::min((nLocalStart + nInterval), nEndAtomNumber);

    // std::map<int, std::vector<double> > atomCoordX;
    // std::map<int, std::vector<double> > atomCoordY;
    // std::map<int, std::vector<double> > atomCoordZ;
    // std::map<int, std::vector<double> > atomWeight;
    std::size_t numOfGrids = 0;
    for (int atom = nLocalStart; atom < nLocalEnd; ++atom) {
        std::vector<double> coordX;
        std::vector<double> coordY;
        std::vector<double> coordZ;
        std::vector<double> weight;

        if (this->m_gridType == SG_1) {
            DfGenerateGrid::generateGrid_SG1(atom, &coordX, &coordY, &coordZ, &weight);
        } else {
            DfGenerateGrid::generateGrid(atom, &coordX, &coordY, &coordZ, &weight);
        }

        // atomCoordX[atom] = coordX;
        // atomCoordY[atom] = coordY;
        // atomCoordZ[atom] = coordZ;
        // atomWeight[atom] = weight;

        // store grid matrix
        const std::size_t numOfAtomGrids = weight.size();
        this->grdMat_.resize(numOfGrids + numOfAtomGrids, this->numOfColsOfGrdMat_);
        for (std::size_t i = 0; i < numOfAtomGrids; ++i) {
            this->grdMat_.set(numOfGrids, 0, coordX[i]);
            this->grdMat_.set(numOfGrids, 1, coordY[i]);
            this->grdMat_.set(numOfGrids, 2, coordZ[i]);
            this->grdMat_.set(numOfGrids, 3, weight[i]);
            this->grdMat_.set(numOfGrids, 4, atom);
            ++numOfGrids;
        }
    }

}

void DfGenerateGrid_Parallel::generateGrid_MS()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    //int startAtom = 0;
    int endAtom = this->m_nNumOfAtoms - this->m_nNumOfDummyAtoms;
    const int rangeAtom = 1;

    enum {
        REQUEST_JOB,
        SUBMIT_RESULTS_X,
        SUBMIT_RESULTS_Y,
        SUBMIT_RESULTS_Z,
        SUBMIT_RESULTS_W
    };


    if (rComm.isMaster() == true) {
        GridDataManager gdm(this->gridDataFilePath_);
        int finishMsgCount = 0;
        const int numOfSlaves = rComm.getNumOfProc() -1;

        int currentAtom = 0;
        while (finishMsgCount < numOfSlaves) {
            int msg = 0;
            int src = 0;
            rComm.receiveDataFromAnySource(msg, &src, TAG_GENGRID_MSG_TO_ROOT);

            switch (msg) {
            case REQUEST_JOB: {
                if (currentAtom < endAtom) {
                    const int lastAtom = std::min((currentAtom + rangeAtom), endAtom);
                    const int range = lastAtom - currentAtom;

                    std::vector<int> atomList(range);
                    for (int i = 0 ; i < range; ++i) {
                        atomList[i] = currentAtom + i;
                    }
                    currentAtom += range;

                    rComm.sendData(range, src, TAG_GENGRID_SEND_RANGE);
                    rComm.sendData(atomList, src, TAG_GENGRID_SEND_ATOMLIST);
                } else {
                    // end of calc.
                    const int range = 0;
                    rComm.sendData(range, src, TAG_GENGRID_SEND_RANGE);
                    ++finishMsgCount;
                }
            }
            break;

            case SUBMIT_RESULTS_X: {
                int atom = 0;
                std::vector<double> coordX;
                rComm.receiveData(atom, src, TAG_GENGRID_SEND_ATOM);
                rComm.receiveData(coordX, src, TAG_GENGRID_SEND_DATA);
                gdm.setData(atom, GridDataManager::COORD_X, coordX);
            }
            break;

            case SUBMIT_RESULTS_Y: {
                int atom = 0;
                std::vector<double> coordY;
                rComm.receiveData(atom, src, TAG_GENGRID_SEND_ATOM);
                rComm.receiveData(coordY, src, TAG_GENGRID_SEND_DATA);
                gdm.setData(atom, GridDataManager::COORD_Y, coordY);
            }
            break;

            case SUBMIT_RESULTS_Z: {
                int atom = 0;
                std::vector<double> coordZ;
                rComm.receiveData(atom, src, TAG_GENGRID_SEND_ATOM);
                rComm.receiveData(coordZ, src, TAG_GENGRID_SEND_DATA);
                gdm.setData(atom, GridDataManager::COORD_Z, coordZ);
            }
            break;

            case SUBMIT_RESULTS_W: {
                int atom = 0;
                std::vector<double> weight;
                rComm.receiveData(atom, src, TAG_GENGRID_SEND_ATOM);
                rComm.receiveData(weight, src, TAG_GENGRID_SEND_DATA);
                gdm.setData(atom, GridDataManager::GRID_WEIGHT, weight);
            }
            break;

            default:
                this->logger(" program error: DfGenerateGrid_Parallel::generateGrid_MS()\n");
                break;
            }
        }

    } else {
        const int root = 0;
        while (true) {
            int msg = REQUEST_JOB;
            rComm.sendData(msg, root, TAG_GENGRID_MSG_TO_ROOT);

            // recieve works from Master
            rComm.receiveData(msg, root, TAG_GENGRID_SEND_RANGE);

            if (msg == 0) {
                break;
            } else {
                std::vector<int> atomList;
                rComm.receiveData(atomList, root, TAG_GENGRID_SEND_ATOMLIST);

                // do work
                std::map<int, std::vector<double> > atomCoordX;
                std::map<int, std::vector<double> > atomCoordY;
                std::map<int, std::vector<double> > atomCoordZ;
                std::map<int, std::vector<double> > atomWeight;
                for (std::vector<int>::const_iterator p = atomList.begin(); p != atomList.end(); ++p) {
                    std::vector<double> coordX;
                    std::vector<double> coordY;
                    std::vector<double> coordZ;
                    std::vector<double> weight;

                    if (this->m_gridType == SG_1) {
                        DfGenerateGrid::generateGrid_SG1(*p, &coordX, &coordY, &coordZ, &weight);
                    } else {
                        DfGenerateGrid::generateGrid(*p, &coordX, &coordY, &coordZ, &weight);
                    }

                    atomCoordX[*p] = coordX;
                    atomCoordY[*p] = coordY;
                    atomCoordZ[*p] = coordZ;
                    atomWeight[*p] = weight;
                }

                // send resuls for Master
                for (std::map<int, std::vector<double> >::const_iterator p = atomCoordX.begin(); p != atomCoordX.end(); ++p) {
                    int atom = p->first;
                    std::vector<double> data = p->second;

                    rComm.sendData(SUBMIT_RESULTS_X, root, TAG_GENGRID_MSG_TO_ROOT);
                    rComm.sendData(atom, root, TAG_GENGRID_SEND_ATOM);
                    rComm.sendData(data, root, TAG_GENGRID_SEND_DATA);
                }

                for (std::map<int, std::vector<double> >::const_iterator p = atomCoordY.begin(); p != atomCoordY.end(); ++p) {
                    int atom = p->first;
                    std::vector<double> data = p->second;

                    rComm.sendData(SUBMIT_RESULTS_Y, root, TAG_GENGRID_MSG_TO_ROOT);
                    rComm.sendData(atom, root, TAG_GENGRID_SEND_ATOM);
                    rComm.sendData(data, root, TAG_GENGRID_SEND_DATA);
                }

                for (std::map<int, std::vector<double> >::const_iterator p = atomCoordZ.begin(); p != atomCoordZ.end(); ++p) {
                    int atom = p->first;
                    std::vector<double> data = p->second;

                    rComm.sendData(SUBMIT_RESULTS_Z, root, TAG_GENGRID_MSG_TO_ROOT);
                    rComm.sendData(atom, root, TAG_GENGRID_SEND_ATOM);
                    rComm.sendData(data, root, TAG_GENGRID_SEND_DATA);
                }

                for (std::map<int, std::vector<double> >::const_iterator p = atomWeight.begin(); p != atomWeight.end(); ++p) {
                    int atom = p->first;
                    std::vector<double> data = p->second;

                    rComm.sendData(SUBMIT_RESULTS_W, root, TAG_GENGRID_MSG_TO_ROOT);
                    rComm.sendData(atom, root, TAG_GENGRID_SEND_ATOM);
                    rComm.sendData(data, root, TAG_GENGRID_SEND_DATA);
                }
            }
        }
    }

    rComm.barrier();
}

