#include <cstdlib>
#include <iostream>
#include <vector>

#include "TlLebedevGrid.h"
#include "TlGetopt.h"
#include "TlPosition.h"
#include "TlUtils.h"

// output Lebedev grid information

void showHelp()
{
    std::cout << "pdf-info-grid [OPTIONS] [grid_size]" << std::endl;
    std::cout << std::endl;
    std::cout << " display gird information." << std::endl;
    std::cout << " -d:           debug mode" << std::endl;
    std::cout << " -h:           help" << std::endl;
};

void debugMode()
{
    TlLebedevGrid lebGrid;
    const std::vector<int> supportedGrids = lebGrid.getSupportedGridNumber();
    for (std::vector<int>::const_iterator it = supportedGrids.begin();
         it != supportedGrids.end(); ++it) {

        std::vector<TlPosition> grids;
        std::vector<double> weights;

        const int numOfGrids = *it;
        lebGrid.getGrids(numOfGrids, &grids, &weights);
        for (int i = 0; i < numOfGrids; ++i) {
            std::cout << numOfGrids << "\t"
                      << i << "\t"
                      << grids[i].x() << "\t"
                      << grids[i].y() << "\t"
                      << grids[i].z() << "\t"
                      << weights[i]
                      << std::endl;
        }
    }
};

int main(int argc, char* argv[])
{
    // parse options
    TlGetopt opt(argc, argv, "dh");
    const bool isShowHelp = (opt["h"] == "defined");
    const bool isDebugMode = (opt["d"] == "defined");

    if (isShowHelp) {
        showHelp();
        std::exit(0);
    }

    if (isDebugMode) {
        debugMode();
        std::exit(0);
    }

    TlLebedevGrid lebGrid;
    const std::vector<int> supportedGrids = lebGrid.getSupportedGridNumber();

    const int numOfArgs = opt.getCount() -1; // opt[0] means the name of this program.
    for (int i = 0; i < numOfArgs; ++i) {
        const int gridSize = std::atoi(opt[i +1].c_str());

        std::vector<int>::const_iterator it = std::find(supportedGrids.begin(),
                                                        supportedGrids.end(),
                                                        gridSize);
        if (it != supportedGrids.end()) {
            std::vector<TlPosition> grids;
            std::vector<double> weights;
            lebGrid.getGrids(gridSize, &grids, &weights);

            std::cout << TlUtils::format("# grid = %d", gridSize) << std::endl;
            for (int i = 0; i < gridSize; ++i) {
                std::cout << TlUtils::format("point(% 12.8e, % 12.8e, % 12.8e) weight=% 12.8e",
                                             grids[i].x(), grids[i].y(), grids[i].z(), weights[i]) << std::endl;
            }
        } else {
            std::cerr << TlUtils::format("# grid = %d, is not defined.", gridSize) << std::endl;
        }
    }

    

    return EXIT_SUCCESS;
};

