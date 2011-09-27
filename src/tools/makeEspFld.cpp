#include <cassert>
#include <cstdlib>
#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlParameter.h"
#include "TlEspField.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "makeField_common.h"

void help(const std::string& progName)
{
    std::cout << TlUtils::format("%s [options] args...", progName.c_str()) << std::endl;
    std::cout << "output volume data of ESP."
              << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p <path>    set ProteinDF parameter file. default = pdfparam.mpac" << std::endl;
    std::cout << " -d <path>    set density matrix file. default is presumed by parameter file." << std::endl;
    std::cout << " -f <path>    save AVS field file." << std::endl;
    std::cout << " -m <path>    save message pack file." << std::endl;
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
    std::string PMatrixFilePath = "";
    if (opt["d"].empty() != true) {
        PMatrixFilePath = opt["d"];
    }
    std::string mpacFilePath = "";
    if (opt["m"].empty() != true) {
        mpacFilePath = opt["m"];
    }
    std::string fieldFilePath = (mpacFilePath != "") ? "" : "ESP.fld";
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

    // ESP計算
    TlEspField espFld(param);
    const std::vector<double> values = espFld.makeEspFld(P, grids);

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

        for (std::size_t i = 0; i < numOfGrids; ++i) {
            output["ESP"].setAt(i, values[i]);
        }

        TlMsgPack mpac(output);
        mpac.save(mpacFilePath);
     }

    // save to fld
    if (fieldFilePath.empty() != true) {
        std::cerr << "save AVS field file: " << fieldFilePath << std::endl;

        std::string label = "ESP";
        std::vector<std::vector<double> > data;
        data.push_back(values);
        
        saveFieldData(numOfGridX, numOfGridY, numOfGridZ, grids, data, label, fieldFilePath);
    }

    return EXIT_SUCCESS;
}


