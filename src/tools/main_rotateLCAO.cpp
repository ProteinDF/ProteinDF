#include <cstdlib>
#include <string>
#include <iostream>
#include "TlGetopt.h"
#include "TlRotateLCAO.h"
#include "TlMsgPack.h"

void showHelp()
{
    std::cout << "rotateLCAO <TARGET_PDFPARAM>" << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  -i <PATH>: define input matrix path. default: rotate_lcao.input" << std::endl;
    std::cout << "  -o <PATH>: define output matrix path. default: rotate_lcao.output" << std::endl;
    std::cout << "  -m <PATH>: define rotate matrix." << std::endl;
    std::cout << "  -h: show help(this)" << std::endl;
    std::cout << "  -v: verbose" << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hi:o:m:v");
    if (opt["h"] == "defined") {
        showHelp();
        return EXIT_SUCCESS;
    }
    const bool verbose = (opt["v"] == "defined");

    //const std::string refPdfParamPath = opt[1];
    const std::string newPdfParamPath = opt[1];
    std::string inLcaoPath = "rotate_lcao.input";
    if (opt["i"] != "") {
        inLcaoPath = opt["i"];
    }
    std::string outLcaoPath = "rotate_lcao.output";
    if (opt["o"] != "") {
        outLcaoPath = opt["o"];
    }
    std::string rotateMatrixPath = "rotate.matrix";
    if (opt["m"] != "") {
        rotateMatrixPath = opt["m"];
    }
    
    
//     TlMsgPack inMsgPack;
//     if (verbose == true) {
//         std::cerr << "reference parameter: " << refPdfParamPath << std::endl;
//     }
//     inMsgPack.load(refPdfParamPath);
//     const TlSerializeData inParam = inMsgPack.getSerializeData();

    TlMsgPack outMsgPack;
    if (verbose == true) {
        std::cerr << "target parameter: " << newPdfParamPath << std::endl;
    }
    outMsgPack.load(newPdfParamPath);
    const TlSerializeData outParam = outMsgPack.getSerializeData();
    
    TlMatrix lcao;
    lcao.load(inLcaoPath);

    TlMatrix rot;
    rot.load(rotateMatrixPath);

    const TlOrbitalInfo orbInfo(outParam["coordinates"], outParam["basis_sets"]);
    //orbInfo.printCGTOs(std::cerr);
    
    TlRotateLCAO rotateLCAO(orbInfo);
    const TlMatrix newLCAO = rotateLCAO.exec(lcao, rot);

    newLCAO.save(outLcaoPath);
}
