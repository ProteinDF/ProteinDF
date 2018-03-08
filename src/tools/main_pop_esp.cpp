#include <iostream>
#include <cassert>
#include <cstdlib>

#include "TlGetopt.h"
#include "TlUtils.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"

#include "TlEspPop.h"


void help(const std::string& progName) {
    std::cout << TlUtils::format("%s [options] args...", progName.c_str()) << std::endl;
    std::cout << "calc population of ESP (MK) charge."
              << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p PATH      set ProteinDF parameter file. default = pdfparam.mpac" << std::endl;
    std::cout << " -d PATH      set density matrix file. default is presumed by parameter file." << std::endl;
    std::cout << " -m PATH      save messagepack file for grids and ESPs" << std::endl;
    std::cout << " -A PATH      save design matrix" << std::endl;
    std::cout << " -y PATH      save target vector" << std::endl;
    std::cout << " -x PATH      save model coef vector" << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


// ---------------------------------------------------------------------
int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "p:d:m:x:y:hv");

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

    std::string mpacFilePath = "grid-esp.mpac";
    if (opt["m"].empty() != true) {
        mpacFilePath = opt["m"];
    }

    std::string designMatrixPath = "esp_design.mat";
    if (opt["A"].empty() != true) {
        designMatrixPath = opt["A"];
    }

    std::string targetVectorPath = "esp_target.vtr";
    if (! opt["y"].empty()) {
        targetVectorPath = opt["y"];
    }

    std::string modelCoefPath = "esp_model.vtr";
    if (opt["x"].empty() != true) {
        modelCoefPath = opt["x"];
    }

    // パラメータファイルの読み込み
    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }

    //
    TlEspPop espPop(param);
    espPop.setRespRestriction(TlEspPop::REST_NONE);
    espPop.verbose(verbose);
    espPop.saveMpacFilePath(mpacFilePath);
    espPop.saveDesignMatrixPath(designMatrixPath);
    espPop.savePredictedVectorPath(targetVectorPath);
    espPop.saveModelCoefVectorPath(modelCoefPath);

    espPop.exec(PMatrixFilePath);

    return EXIT_SUCCESS;
}
