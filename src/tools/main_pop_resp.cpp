#include <iostream>
#include <cassert>
#include <cstdlib>

#include "TlGetopt.h"
#include "TlUtils.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlPrdctbl.h"
#include "TlLebedevGrid.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "TlPosition.h"
#include "TlMath.h"
#include "TlEspField.h"

#include "DfObject.h"
#include "Fl_Geometry.h"

#define AU2ANG 0.5291772108
#define ANG2AU 1.889762

// std::vector<TlPosition> getMerzKollmanGrids(const TlSerializeData& param);
// std::vector<TlPosition> getMKGridsOnAtom(const double vdwr);
// void makeMat(const TlSerializeData& param,
//              const std::vector<TlPosition>& grids,
//              TlMatrix* pA, TlVector* pB);

void help(const std::string& progName) {
    std::cout << TlUtils::format("%s [options] args...", progName.c_str()) << std::endl;
    std::cout << "calc population of CEEM."
              << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p PATH      set ProteinDF parameter file. default = pdfparam.mpac" << std::endl;
    std::cout << " -d PATH      set density matrix file. default is presumed by parameter file." << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


// radii: atomic unit
std::vector<TlPosition> getMKGridsOnAtom(const TlPosition& center, const double radii)
{
    TlLebedevGrid lebGrid;
    const std::vector<int> gridList = lebGrid.getSupportedGridNumber();
    
    std::vector<TlPosition> grids;
    const double r = radii * AU2ANG; // to angstroam unit
        
    const double area = 4.0 * TlMath::PI() * r * r;
    std::vector<int>::const_iterator it = std::upper_bound(gridList.begin(),
                                                           gridList.end(),
                                                           static_cast<int>(area));
    const int numOfGrids = *it;
    std::cout << TlUtils::format("area=%8.3f ANG^2 grids=%d", area, numOfGrids) << std::endl;
        
    std::vector<TlPosition> layerGrids;
    std::vector<double> layerWeights;
    lebGrid.getGrids(numOfGrids, &layerGrids, &layerWeights);
    for (int grid = 0; grid < numOfGrids; ++grid) {
        layerGrids[grid] *= radii;
        layerGrids[grid].shiftBy(center);
    }
        
    grids.insert(grids.end(), layerGrids.begin(), layerGrids.end());

    return grids;
}


bool isInMolecule(const TlPosition& p, double coef, const Fl_Geometry& flGeom)
{
    bool answer = false;

    const int numOfAtoms = flGeom.getNumOfAtoms();
    for (int i = 0; i < numOfAtoms; ++i) {
        const TlAtom atom = flGeom.getAtom(i);
        const std::string symbol = atom.getSymbol();
        if (symbol != "X") {
            const double distance = p.distanceFrom(atom.getPosition());
            const double vdwr = TlPrdctbl::getVdwRadii(TlPrdctbl::getAtomicNumber(symbol)) * ANG2AU;
            if (distance < vdwr * coef) {
                answer = true;
                break;
            }
        }
    }

    return answer;
}


std::vector<TlPosition> getMerzKollmanGrids(const TlSerializeData& param)
{
    static const double layers[] = {1.4, 1.6, 1.8, 2.0};
    static const int numOfLayers = sizeof(layers)/sizeof(layers[0]);

    Fl_Geometry flGeom(param["coordinates"]); // 単位はa.u.

    std::vector<TlPosition> allGrids;
    std::size_t numOfGrids = 0;
    std::size_t screened = 0;
    const int numOfAtoms = flGeom.getNumOfAtoms();
    for (int i = 0; i < numOfAtoms; ++i) {
        const TlAtom atom = flGeom.getAtom(i);
        const std::string symbol = atom.getSymbol();
        if (symbol != "X") {
            const double vdwr = TlPrdctbl::getVdwRadii(TlPrdctbl::getAtomicNumber(symbol)) * ANG2AU;

            for (int layer_index = 0; layer_index < numOfLayers; ++layer_index) {
                const double coef = layers[layer_index];
                
                // define the number of grids
                std::vector<TlPosition> grids = getMKGridsOnAtom(atom.getPosition(), coef * vdwr);
                
                // check in molecule
                std::vector<TlPosition>::iterator itEnd = grids.end();
                for (std::vector<TlPosition>::iterator it = grids.begin(); it != itEnd; ++it) {
                    if (isInMolecule(*it, coef, flGeom) != true) {
                        allGrids.push_back(*it);
                    } else {
                        ++screened;
                    }
                    ++numOfGrids;
                }
            }            
        }
    }

    std::cerr << TlUtils::format("screened %ld/%ld", screened, numOfGrids) << std::endl;
    return allGrids;
}


void makeMat(const TlSerializeData& param,
             const std::vector<TlPosition>& grids,
             const TlVector& esps,
             TlSymmetricMatrix* pA, TlVector* pB)
{
    const std::size_t numOfGrids = grids.size();
    assert(esps.getSize() == numOfGrids);
    
    const Fl_Geometry flGeom(param["coordinates"]); // 単位はa.u.
    const int numOfAllAtoms = flGeom.getNumOfAtoms();
    const int numOfRealAtoms = numOfAllAtoms - flGeom.getNumOfDummyAtoms();

    std::vector<int> realAtoms(numOfRealAtoms);
    {
        int realAtomIndex = 0;
        for (int i = 0; i < numOfAllAtoms; ++i) {
            const TlAtom atom = flGeom.getAtom(i);
            if (atom.getSymbol() != "X") {
                realAtoms[realAtomIndex] = i;
                ++realAtomIndex;
            }
        }
        assert(realAtomIndex == numOfRealAtoms);
    }
    
    // make 1/r distance table
    TlMatrix d(numOfRealAtoms, numOfGrids);
    std::cerr << TlUtils::format("# of atoms: %d", numOfRealAtoms) << std::endl;
    std::cerr << TlUtils::format("# of grids: %d", numOfGrids) << std::endl;
    for (int a = 0; a < numOfRealAtoms; ++a) {
        const TlPosition posA = flGeom.getCoordinate(realAtoms[a]);
        for (std::size_t i = 0; i < numOfGrids; ++i) {
            const double r_ai = posA.distanceFrom(grids[i]);

            d.set(a, i, 1.0 / r_ai);
        }
    }
    d.save("d.mat");
    
    // make M & b
    assert(pA != NULL);
    assert(pB != NULL);
    pA->resize(numOfRealAtoms +1);
    pB->resize(numOfRealAtoms +1);
    for (int a = 0; a < numOfRealAtoms; ++a) {
        const TlVector r_a = d.getRowVector(a);
        assert(r_a.getSize() == numOfGrids);
        
        // a == b
        {
            TlVector r_a2 = r_a;
            r_a2.dot(r_a);
            pA->set(a, a, r_a2.sum());
        }

        // a != b
        for (int b = 0; b < a; ++b) {
            TlVector r_b = d.getRowVector(b);
            r_b.dot(r_a);
            pA->set(a, b, r_b.sum());
        }

        // for Lagurange
        pA->set(a, numOfRealAtoms, 1.0);

        //
        {
            TlVector r = r_a;
            (*pB)[a] = r.dot(esps).sum();
        }
    }
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "p:d:hv");
    
    const bool verbose = (opt["v"] == "defined");
    if ((opt["h"] == "defined")) {
        help(opt[0]);
        return EXIT_SUCCESS;
    }

    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }
    std::string PMatrixFilePath = "";
    if (opt["d"].empty() != true) {
        PMatrixFilePath = opt["d"];
    }

    // パラメータファイルの読み込み
    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }

    std::vector<TlPosition> grids = getMerzKollmanGrids(param);
    const std::size_t numOfGrids = grids.size();
    if (verbose) {
        std::cerr << TlUtils::format("# grids: %ld", numOfGrids) << std::endl;
    }
    
    // 密度行列の読み込み
    TlSymmetricMatrix P;
    if (PMatrixFilePath == "") {
        const int iteration = param["num_of_iterations"].getInt();
        DfObject::RUN_TYPE runType = DfObject::RUN_RKS;
        DfObject dfObject(&param);
        PMatrixFilePath = dfObject.getPpqMatrixPath(runType, iteration);
    }
    if (verbose) {
        std::cerr << "loading: " << PMatrixFilePath << std::endl;
    }
    P.load(PMatrixFilePath);

    // ESP計算
    TlEspField espFld(param);
    const TlVector esp = espFld.makeEspFld(P, grids);
    
    //
    TlSymmetricMatrix A;
    TlVector B;
    makeMat(param, grids, esp, &A, &B);
    // A.save("A.mat");
    // B.save("B.vtr");

    // solve
    A.inverse();
    // A.save("Ainv.mat");

    TlVector x = A * B;
    x.save("resp_x.vtr");
    x.print(std::cout);
}


