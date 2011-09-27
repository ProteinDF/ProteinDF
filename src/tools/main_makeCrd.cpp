#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <cassert>
#include <cstdlib>

#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlDensField.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "makeField_common.h"

#define DEFAULT_OUTPUT_FILE_NAME "crd.mpac"


void help(const std::string& progName)
{
    const std::string outputFileName = DEFAULT_OUTPUT_FILE_NAME;
    std::cout << TlUtils::format("Usage: %s [options] [file]", progName.c_str()) << std::endl;
    std::cout << TlUtils::format("output grid coordinates msgpack file(default: %s) to calculate physical values.",
                                 outputFileName.c_str())
              << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p <path>    set ProteinDF parameter file(default: pdfparam.mpac)." << std::endl;
    std::cout << " -m <path>    save message pack file." << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hm:p:v");

    const bool isVerboseMode = (opt["v"] == "defined");
    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }
    if ((opt["h"] == "defined")) {
        help(opt[0]);
        return EXIT_SUCCESS;
    }

    std::string mpacFilePath = DEFAULT_OUTPUT_FILE_NAME;
    if (opt.getCount() > 1) {
        mpacFilePath = opt[1];
    }
    
    // パラメータファイルの読み込み
    TlSerializeData param;
    {
        TlMsgPack mpac;
        if (isVerboseMode == true) {
            std::cerr << TlUtils::format("loading ProteinDF parameter: %s", pdfParamPath.c_str())
                      << std::endl;
        }
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }

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

    if (isVerboseMode == true) {
        std::cerr << TlUtils::format("start point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                     startPos.x(), startPos.y(), startPos.z())
                  << std::endl;
        std::cerr << TlUtils::format("end point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                     endPos.x(), endPos.y(), endPos.z())
                  << std::endl;
    }

    // save to mpac
    if (isVerboseMode == true) {
        std::cerr << "save message pack file: " << mpacFilePath << std::endl;
    }

    TlSerializeData output;
    output["version"] = VERSION; // set in "config.h"
    output["num_of_grids_x"] = numOfGridX;
    output["num_of_grids_y"] = numOfGridY;
    output["num_of_grids_z"] = numOfGridZ;

    output["coord"].resize(numOfGrids);
    for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
        TlSerializeData pos;
        pos.pushBack(grids[gridIndex].x());
        pos.pushBack(grids[gridIndex].y());
        pos.pushBack(grids[gridIndex].z());
        
        output["coord"].setAt(gridIndex, pos);
    }
    
    TlMsgPack mpac(output);
    mpac.save(mpacFilePath);

    return EXIT_SUCCESS;
}


