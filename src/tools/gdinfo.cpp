#include <iostream>
#include <cstdlib>

#include "GridDataManager.h"
#include "TlVector.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "f:hv");

    const bool isVerbose = (opt["v"] == "defined");
    //const bool bGuessMode = (opt["f"] == "defined");

    std::string gridDataFilePath = "./fl_Work/grids.dat";
    if (opt["f"] != "") {
        gridDataFilePath = opt["f"];
    }

    if (isVerbose == true) {
        std::cerr << "GridData: " << gridDataFilePath << std::endl;
    }

    GridDataManager gdm(gridDataFilePath);
    std::cout << gdm.str() << std::endl;

    return EXIT_SUCCESS;
}


