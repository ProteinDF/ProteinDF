#include <cassert>
#include <cstdlib>

#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlDensField.h"
#include "TlEspField.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "makeField_common.h"

void help(const std::string& progName)
{
    std::cout << TlUtils::format("%s [options] args...", progName.c_str()) << std::endl;
    std::cout << "output volume data of density and ESP."
              << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p FILE      set ProteinDF parameter file.(default: pdfparam.mpac)" << std::endl;
    std::cout << " -c FILE      set coordinate data" << std::endl;
    std::cout << " -d FILE      set density matrix file. default is presumed by parameter file." << std::endl;
    std::cout << " -f FILE      save AVS field file." << std::endl;
    std::cout << " -m FILE      save message pack file." << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "d:f:hm:p:v");

    const bool verbose = (opt["v"] == "defined");
    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }
    std::string crdFilePath = "";
    if (opt["c"].empty() != true) {
        crdFilePath = opt["c"];
    }
    std::string PMatrixFilePath = "";
    if (opt["d"].empty() != true) {
        PMatrixFilePath = opt["d"];
    }
    std::string mpacFilePath = "";
    if (opt["m"].empty() != true) {
        mpacFilePath = opt["m"];
    }
    std::string fieldFilePath = (mpacFilePath != "") ? "" : "densESP.fld";
    if (opt["f"].empty() != true) {
        fieldFilePath = opt["f"];
    }

    if ((opt["h"] == "defined")) {
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

    // 密度行列の読み込み
    TlSymmetricMatrix P;
    if (PMatrixFilePath == "") {
        const int iteration = param["model"]["iterations"].getInt();
        DfObject::RUN_TYPE runType = DfObject::RUN_RKS;
        DfObject dfObject(&param);
        PMatrixFilePath = dfObject.getPpqMatrixPath(runType, iteration);
    }
    P.load(PMatrixFilePath);

    // 計算サイズの決定
    std::vector<TlPosition> grids;
    int numOfGridX = 0;
    int numOfGridY = 0;
    int numOfGridZ = 0;
    if (crdFilePath.empty() == true) {
        TlPosition startPos;
        TlPosition endPos;
        getDefaultSize(param, &startPos, &endPos);
        compensateRange(&startPos, &endPos);
        
        const TlPosition gridPitch(GRID_PITCH, GRID_PITCH, GRID_PITCH);
        grids = makeGrids(startPos, endPos, gridPitch,
                          &numOfGridX, &numOfGridY, &numOfGridZ);

        if (verbose == true) {
            std::cerr << TlUtils::format("start point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                         startPos.x(), startPos.y(), startPos.z())
                      << std::endl;
            std::cerr << TlUtils::format("end point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                         endPos.x(), endPos.y(), endPos.z())
                      << std::endl;
        }
    } else {
        TlMsgPack crdFile;
        if (verbose == true) {
            std::cerr << TlUtils::format("reading crd file: %s", crdFilePath.c_str())
                      << std::endl;
        }
        crdFile.load(crdFilePath);
        const TlSerializeData crdData = crdFile.getSerializeData();
        numOfGridX = crdData["num_of_grids_x"].getInt();
        numOfGridY = crdData["num_of_grids_y"].getInt();
        numOfGridZ = crdData["num_of_grids_z"].getInt();
        const std::size_t numOfGrids = std::size_t(numOfGridX) * std::size_t(numOfGridY) * std::size_t(numOfGridZ);
        grids.resize(numOfGrids);
        for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            const TlSerializeData p = crdData["coord"].getAt(gridIndex);
            grids[gridIndex] = TlPosition(p.getAt(0).getDouble(),
                                          p.getAt(1).getDouble(),
                                          p.getAt(2).getDouble());
        }
    }

    // 電子密度計算
    TlDensField densFld(param);
    const std::vector<double> values_den = densFld.makeDensFld(P, grids);

    // ESP計算
    TlEspField espFld(param);
    const std::vector<double> values_esp = espFld.makeEspFld(P, grids);


    // save to mpac
    if (mpacFilePath.empty() != true) {
        std::cerr << "save message pack file: " << mpacFilePath << std::endl;

        TlSerializeData output;
        output["version"] = "2010.0";
        output["num_of_grids_x"] = numOfGridX;
        output["num_of_grids_y"] = numOfGridY;
        output["num_of_grids_z"] = numOfGridZ;

        const std::size_t numOfGrids = std::size_t(numOfGridX) * std::size_t(numOfGridY) * std::size_t(numOfGridZ);
        output["coord"].resize(numOfGrids);
        for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            TlSerializeData pos;
            pos.pushBack(grids[gridIndex].x());
            pos.pushBack(grids[gridIndex].y());
            pos.pushBack(grids[gridIndex].z());
            
            output["coord"].setAt(gridIndex, pos);
        }

        output["density"].resize(numOfGrids);
        for (std::size_t i = 0; i < numOfGrids; ++i) {
            output["density"].setAt(i, values_den[i]);
        }

        output["ESP"].resize(numOfGrids);
        for (std::size_t i = 0; i < numOfGrids; ++i) {
            output["ESP"].setAt(i, values_esp[i]);
        }

        TlMsgPack mpac(output);
        mpac.save(mpacFilePath);
     }

    // save to fld
    if (fieldFilePath.empty() != true) {
        if (verbose == true) {
            std::cerr << "save AVS field file: " << fieldFilePath << std::endl;
        }

        // convert grid unit to angstrom
        const std::size_t numOfGrids = grids.size();
        for (std::size_t i = 0; i < numOfGrids; ++i) {
            grids[i] *= ANG_PER_AU;
        }

        std::string label = "density ESP";
        std::vector<std::vector<double> > data;
        data.push_back(values_den);
        data.push_back(values_esp);
        
        saveFieldData(numOfGridX, numOfGridY, numOfGridZ, grids, data, label, fieldFilePath);
    }

    return EXIT_SUCCESS;
}


