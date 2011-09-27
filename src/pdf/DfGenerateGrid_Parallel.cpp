#include <cassert>
#include "DfGenerateGrid_Parallel.h"
#include "DfXCFunctional.h"
#include "TlCommunicate.h"
#include "GridDataManager.h"
#include "TlUtils.h"

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

    if (this->isMasterSlave_ == true) {
        this->generateGrid_MS();
    } else {
        this->generateGrid_DC();
    }

    // grid matrix
    // 最適化すること
    if (rComm.isMaster() == true) {
        int numOfCols = 4; // x, y, z, weight
        const int coef = (this->m_nMethodType == METHOD_RKS) ? 1 : 2;
        {
            DfXCFunctional dfXcFunctional(this->pPdfParam_);
            if (dfXcFunctional.getFunctionalType() == DfXCFunctional::LDA) {
                numOfCols += coef * 1; // rho only
            } else if (dfXcFunctional.getFunctionalType() == DfXCFunctional::GGA) {
                numOfCols += coef * 4; // rho, gradRhoX, gradRhoY, gradRhoZ
            }
        }
        
        TlMatrix gridMtx(1, numOfCols);
        GridDataManager gdm(this->gridDataFilePath_);
        const int endAtom = this->numOfRealAtoms_;
        std::size_t numOfGrids = 0;
        for (int atom = 0; atom < endAtom; ++atom) {
            std::vector<double> coordX = gdm.getData(atom, GridDataManager::COORD_X);
            std::vector<double> coordY = gdm.getData(atom, GridDataManager::COORD_Y);
            std::vector<double> coordZ = gdm.getData(atom, GridDataManager::COORD_Z);
            std::vector<double> weight = gdm.getData(atom, GridDataManager::GRID_WEIGHT);

            const std::size_t gridSize = coordX.size();
            gridMtx.resize(numOfGrids + gridSize, numOfCols);
            for (std::size_t i = 0; i < gridSize; ++i) {
                gridMtx.set(numOfGrids + i, 0, coordX[i]);
                gridMtx.set(numOfGrids + i, 1, coordY[i]);
                gridMtx.set(numOfGrids + i, 2, coordZ[i]);
                gridMtx.set(numOfGrids + i, 3, weight[i]);
            }
            numOfGrids += gridSize;
        }
        this->saveGridMatrix(gridMtx);
    }
    //
    
    rComm.barrier();
}

void DfGenerateGrid_Parallel::generateGrid_DC()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    //const int nStartAtomNumber = 0;
    const int nEndAtomNumber = this->m_nNumOfAtoms;

    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();
    const int nRange = nEndAtomNumber;
    const int nInterval = (nRange + (nProc -1)) / nProc; // +(nProc-1) は余り対策
    const int nLocalStart = nRank * nInterval; // nProc = 0, 1, 2, ...
    const int nLocalEnd   = std::min((nLocalStart + nInterval), nEndAtomNumber);

    // A set of grid points is generated around each nucleus in the system.
    // Loop for the number of atom
    rComm.barrier();

    std::map<int, std::vector<double> > atomCoordX;
    std::map<int, std::vector<double> > atomCoordY;
    std::map<int, std::vector<double> > atomCoordZ;
    std::map<int, std::vector<double> > atomWeight;
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

        atomCoordX[atom] = coordX;
        atomCoordY[atom] = coordY;
        atomCoordZ[atom] = coordZ;
        atomWeight[atom] = weight;
//     std::cerr << TlUtils::format("DfGenerateGrid_Parallel::generateGrid() atom=%d, size=%d",
//               atom, coordX.size())
//        << std::endl;
    }

    // all reduce
    if (rComm.isMaster() == true) {
        GridDataManager gdm(this->gridDataFilePath_);
        // for master
        {
            for (std::map<int, std::vector<double> >::const_iterator p = atomCoordX.begin();
                    p != atomCoordX.end(); ++p) {
                gdm.setData(p->first, GridDataManager::COORD_X, p->second);
            }
            for (std::map<int, std::vector<double> >::const_iterator p = atomCoordY.begin();
                    p != atomCoordY.end(); ++p) {
                gdm.setData(p->first, GridDataManager::COORD_Y, p->second);
            }
            for (std::map<int, std::vector<double> >::const_iterator p = atomCoordZ.begin();
                    p != atomCoordZ.end(); ++p) {
                gdm.setData(p->first, GridDataManager::COORD_Z, p->second);
            }
            for (std::map<int, std::vector<double> >::const_iterator p = atomWeight.begin();
                    p != atomWeight.end(); ++p) {
                gdm.setData(p->first, GridDataManager::GRID_WEIGHT, p->second);
            }
        }
        // for slave
        for (int proc = 1; proc < nProc; ++proc) {
            int numOfAtoms = 0;
            rComm.receiveData(numOfAtoms, proc);

            for (int i = 0; i < numOfAtoms; ++i) {
                int atom = 0;
                std::vector<double> coordX;

                rComm.receiveData(atom, proc);
                rComm.receiveData(coordX, proc);
                gdm.setData(atom, GridDataManager::COORD_X, coordX);
            }

            for (int i = 0; i < numOfAtoms; ++i) {
                int atom = 0;
                std::vector<double> coordY;

                rComm.receiveData(atom, proc);
                rComm.receiveData(coordY, proc);
                gdm.setData(atom, GridDataManager::COORD_Y, coordY);
            }

            for (int i = 0; i < numOfAtoms; ++i) {
                int atom = 0;
                std::vector<double> coordZ;

                rComm.receiveData(atom, proc);
                rComm.receiveData(coordZ, proc);
                gdm.setData(atom, GridDataManager::COORD_Z, coordZ);
            }

            for (int i = 0; i < numOfAtoms; ++i) {
                int atom = 0;
                std::vector<double> weight;

                rComm.receiveData(atom, proc);
                rComm.receiveData(weight, proc);
                gdm.setData(atom, GridDataManager::GRID_WEIGHT, weight);
            }

            rComm.barrier();
        }
    } else {
        for (int proc = 1; proc < nProc; ++proc) {
            if (proc == nRank) {
                const int numOfAtoms = atomCoordX.size();
                assert(static_cast<std::size_t>(numOfAtoms) == atomCoordY.size());
                assert(static_cast<std::size_t>(numOfAtoms) == atomCoordZ.size());
                assert(static_cast<std::size_t>(numOfAtoms) == atomWeight.size());
                rComm.sendData(numOfAtoms);

                for (std::map<int, std::vector<double> >::const_iterator p = atomCoordX.begin();
                        p != atomCoordX.end(); ++p) {
                    rComm.sendData(p->first);
                    rComm.sendData(p->second);
                }
                for (std::map<int, std::vector<double> >::const_iterator p = atomCoordY.begin();
                        p != atomCoordY.end(); ++p) {
                    rComm.sendData(p->first);
                    rComm.sendData(p->second);
                }
                for (std::map<int, std::vector<double> >::const_iterator p = atomCoordZ.begin();
                        p != atomCoordZ.end(); ++p) {
                    rComm.sendData(p->first);
                    rComm.sendData(p->second);
                }
                for (std::map<int, std::vector<double> >::const_iterator p = atomWeight.begin();
                        p != atomWeight.end(); ++p) {
                    rComm.sendData(p->first);
                    rComm.sendData(p->second);
                }
            }

            rComm.barrier();
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

