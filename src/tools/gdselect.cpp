#include <iostream>
#include <cstdlib>

#include "GridDataManager.h"
#include "TlVector.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "f:hvs:");

    const bool isVerbose = (opt["v"] == "defined");

    std::string gridDataFilePath = "./fl_Work/grids.dat";
    if (opt["f"] != "") {
        gridDataFilePath = opt["f"];
    }

    std::string savePath = "";
    if (opt["s"].empty() != true) {
        savePath = opt["s"];
    }
    
    int atomIndex = 0;
    if (opt[1] != "") {
        atomIndex = std::atoi(opt[1].c_str());
    }

    int chunkType = 0;
    if (opt[2] != "") {
        chunkType = std::atoi(opt[2].c_str());
    }
    if ((chunkType < 1) || (chunkType > 16)) {
        chunkType = 1;
    }

    static const std::string chunkTypeStr[] = {
        "UNDEFINED",
        "COORD_X",
        "COORD_Y",
        "COORD_Z",
        "GRID_WEIGHT",
        "DENSITY",
        "DENSITY_ALPHA",
        "DENSITY_BETA",
        "GRAD_DENSITY_X",
        "GRAD_DENSITY_Y",
        "GRAD_DENSITY_Z",
        "GRAD_DENSITY_X_ALPHA",
        "GRAD_DENSITY_Y_ALPHA",
        "GRAD_DENSITY_Z_ALPHA",
        "GRAD_DENSITY_X_BETA",
        "GRAD_DENSITY_Y_BETA",
        "GRAD_DENSITY_Z_BETA"
    };

    if (isVerbose == true) {
        std::cerr << "GridData: " << gridDataFilePath << std::endl;
        std::cerr << "atom index: " << atomIndex << std::endl;
        std::cerr << "chunk type: " << chunkTypeStr[chunkType] << std::endl;
    }

    GridDataManager gdm(gridDataFilePath);

    const std::vector<double> tmp = gdm.getData(atomIndex, GridDataManager::ChunkType(chunkType));
    const TlVector output(tmp);

    if (savePath.empty() == true) {
        output.print(std::cout);
    } else {
        std::cerr << "save: " << savePath << std::endl;
        output.save(savePath);
    }

    return EXIT_SUCCESS;
}


