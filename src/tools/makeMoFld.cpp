#include <cassert>
#include <cstdlib>
#include <vector>
#include <map>
#include <fstream>
#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlAtom.h"
#include "TlMoField.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "makeField_common.h"

void help(const std::string& progName)
{
    std::cout << TlUtils::format("%s [options] <MO index> [<MO index> ...]", progName.c_str()) << std::endl;
    std::cout << std::endl;
    std::cout << "output volume data of molecular orbital" << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p <path>    set ProteinDF parameter file(default: pdfparam.mpac)" << std::endl;
    std::cout << " -l <path>    set LCAO matrix file" << std::endl;
    std::cout << " -f <path>    save AVS field file" << std::endl;
    std::cout << " -m <path>    save message pack file" << std::endl;
    std::cout << " -c <path>    save cube file" << std::endl;
    std::cout << " -h           show help message (this)" << std::endl;
    std::cout << " -v           verbose mode" << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "c:f:hl:m:p:v");

    const bool verbose = (opt["v"] == "defined");
    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }
    std::string lcaoMatrixPath = "";
    if (opt["l"].empty() != true) {
        lcaoMatrixPath = opt["l"];
    }
    std::string mpacFilePath = "";
    if (opt["m"].empty() != true) {
        mpacFilePath = opt["m"];
    }
    std::string fieldFilePath = "";
    if (opt["f"].empty() != true) {
        fieldFilePath = opt["f"];
    }
    std::string cubeFilePath = "";
    if (opt["c"].empty() != true) {
        cubeFilePath = opt["c"];
    }

    if ((lcaoMatrixPath.empty() && mpacFilePath.empty() && cubeFilePath.empty()) ||
        (opt["h"] == "defined")) {
        help(opt[0]);
        return EXIT_SUCCESS;
    }
    
    // パラメータファイルの読み込み
    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }

    // LCAO行列の読み込み
    TlMatrix C;
    if (lcaoMatrixPath.empty() == true) {
        const int iteration = param["num_of_iterations"].getInt();
        DfObject::RUN_TYPE runType = DfObject::RUN_RKS;
        DfObject dfObject(&param);
        lcaoMatrixPath = dfObject.getCMatrixPath(runType, iteration);
    }
    std::cerr << "C matrix path: " << lcaoMatrixPath << std::endl;
    C.load(lcaoMatrixPath);

    // 計算サイズの決定
    TlPosition startPos;
    TlPosition endPos;
    getDefaultSize(param, &startPos, &endPos);
    compensateRange(&startPos, &endPos);

    const TlPosition gridPitch(GRID_PITCH, GRID_PITCH, GRID_PITCH);
    int numOfGridX, numOfGridY, numOfGridZ;
    std::vector<TlPosition> grids = makeGrids(startPos, endPos, gridPitch,
                                              &numOfGridX, &numOfGridY, &numOfGridZ);
    const std::size_t numOfGrids = grids.size();

    if (verbose == true) {
        std::cerr << TlUtils::format("start point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                     startPos.x(), startPos.y(), startPos.z())
                  << std::endl;
        std::cerr << TlUtils::format("end point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                     endPos.x(), endPos.y(), endPos.z())
                  << std::endl;
    }
    
    // MO計算
    TlMoField moFld(param);
    std::map<std::string, std::vector<double> > storeData;
    // i = 1 から始まる理由はopt[0]がプログラム名のため
    for (int i = 1; i < opt.getCount(); ++i) {
        const int MO_index = std::atoi(opt[i].c_str());
        const std::vector<double>values = moFld.makeMoFld(C.getColVector(MO_index), grids);

        const std::string key = TlUtils::format("MO_%d", MO_index);
        storeData[key] = values;
    }

    // convert grid unit to angstrom
    for (std::size_t i = 0; i < numOfGrids; ++i) {
        grids[i] *= ANG_PER_AU;
    }
    
    // save to mpac
    if (mpacFilePath.empty() != true) {
        std::cerr << "save message pack file: " << mpacFilePath << std::endl;

        TlSerializeData output;
        output["version"] = "2010.0";
        output["num_of_grids_x"] = numOfGridX;
        output["num_of_grids_y"] = numOfGridY;
        output["num_of_grids_z"] = numOfGridZ;
        for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            TlSerializeData pos;
            pos.pushBack(grids[gridIndex].x());
            pos.pushBack(grids[gridIndex].y());
            pos.pushBack(grids[gridIndex].z());
            
            output["coord"].pushBack(pos);
        }

        for (std::map<std::string, std::vector<double> >::const_iterator p = storeData.begin();
             p != storeData.end(); ++p) {
            const std::string key = p->first;
            const std::vector<double>& value = p->second;
            std::cerr << key << std::endl;

            for (std::size_t i = 0; i < numOfGrids; ++i) {
                output[key].setAt(i, value[i]);
            }
        }

        TlMsgPack mpac(output);
        mpac.save(mpacFilePath);
     }

    // save to fld
    if (fieldFilePath.empty() != true) {
        std::cerr << "save AVS field file: " << fieldFilePath << std::endl;

        std::string label = "";
        std::vector<std::vector<double> > data;
        for (std::map<std::string, std::vector<double> >::const_iterator p = storeData.begin();
             p != storeData.end(); ++p) {
            const std::string key = p->first;
            const std::vector<double>& value = p->second;

            label = TlUtils::format("%s %s", label.c_str(), key.c_str());
            data.push_back(value);
        }

        saveFieldData(numOfGridX, numOfGridY, numOfGridZ, grids, data, label, fieldFilePath);
    }

    // save to cube
    if (cubeFilePath.empty() != true) {
        const Fl_Geometry flGeom(param["coordinates"]);
        const int numOfAtoms = flGeom.getNumOfAtoms();

        std::string output = "";
        output += "comment line: \n";
        output += "comment line: \n";
        output += TlUtils::format("%5d % 12.6f % 12.6f % 12.6f\n",
                                  numOfAtoms, startPos.x(), startPos.y(), startPos.z());

        output += TlUtils::format("%5d % 12.6f % 12.6f % 12.6f\n", numOfGridX, gridPitch.x(), 0.0, 0.0);
        output += TlUtils::format("%5d % 12.6f % 12.6f % 12.6f\n", numOfGridY, 0.0, gridPitch.y(), 0.0);
        output += TlUtils::format("%5d % 12.6f % 12.6f % 12.6f\n", numOfGridZ, 0.0, 0.0, gridPitch.z());

        for (int i = 0; i < numOfAtoms; ++i) {
            const int atomic_number = TlAtom::getElementNumber(flGeom.getAtom(i));
            const double charge = flGeom.getCharge(i);
            const TlPosition p = flGeom.getCoordinate(i);
            output += TlUtils::format("%5d % 12.6f % 12.6f % 12.6f % 12.6f\n",
                                      atomic_number, charge, p.x(), p.y(), p.z());
        }

        for (std::map<std::string, std::vector<double> >::const_iterator p = storeData.begin();
             p != storeData.end(); ++p) {

            const std::string key = p->first;
            const std::vector<double>& value = p->second;

            std::string dat_str = "";
            //int counter = 0;
            for (int x = 0; x < numOfGridX; ++x) {
                for (int y = 0; y < numOfGridY; ++y) {
                    for (int z = 0; z < numOfGridZ; ++z) {
                        const int index = (z*numOfGridY +y)*numOfGridX +x;
                        dat_str += TlUtils::format("% 12.5E ", value[index]);
                        if (z % 6 == 5) {
                            dat_str += "\n";
                        }
                        // ++counter;
                    }
                    dat_str += "\n";
                }
            }

            const std::string path = TlUtils::format("%s_%s.cube",
                                                     cubeFilePath.c_str(),
                                                     key.c_str());
            std::cerr << "save CUBE file: " << path << std::endl;
            std::ofstream ofs(path.c_str());
            ofs << output << dat_str << std::endl;
            ofs.close();
        }
    }
    
    return EXIT_SUCCESS;
}


