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
    std::cout << "calc population of RESP(ESP)."
              << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p PATH      set ProteinDF parameter file. default = pdfparam.mpac" << std::endl;
    std::cout << " -d PATH      set density matrix file. default is presumed by parameter file." << std::endl;
    std::cout << " -m PATH      save messagepack file for grids and ESPs" << std::endl;
    std::cout << " -A PATH      save design matrix" << std::endl;
    std::cout << " -y PATH      save predicted vector" << std::endl;
    std::cout << " -x PATH      save model coef vector" << std::endl;
    std::cout << " -r [Q, H]    RESP restriction function" << std::endl;
    std::cout << " -a value     RESP parameter a" << std::endl;
    std::cout << " -b value     RESP parameter b" << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


// ---------------------------------------------------------------------
int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "A:a:b:d:hm:p:r:vx:y:");
    
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

    std::string designMatrixPath = "MK_design.mat";
    if (opt["A"].empty() != true) {
        designMatrixPath = opt["A"];
    }

    std::string predictedPath = "MK_predicted.vtr";
    if (opt["y"].empty() != true) {
        predictedPath = opt["y"];
    }
    std::string modelCoefPath = "MK_model.vtr";
    if (opt["x"].empty() != true) {
        modelCoefPath = opt["x"];
    }
    
    double param_a = 0.0005;
    if (opt["a"].empty() != true) {
        param_a = std::atof(opt["a"].c_str());
    }
    
    double param_b = 0.1;
    if (opt["b"].empty() != true) {
        param_a = std::atof(opt["b"].c_str());
    }

    TlEspPop::RESP_RESTRICTION rest = TlEspPop::REST_NONE;
    if (opt["r"].empty() != true) {
        std::string value = opt["r"];
        value = TlUtils::toUpper(value);
        switch (value[0]) {
        case 'Q':
            rest = TlEspPop::REST_QUADRIC;
            break;

        case 'H':
            rest = TlEspPop::REST_HYPERBOLIC;
            break;

        default:
            rest = TlEspPop::REST_NONE;
            break;
        }
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
    espPop.setRespRestriction(rest);
    espPop.setRestrictionParameterA(param_a);
    espPop.setRestrictionParameterB(param_b);
    espPop.verbose(verbose);
    espPop.saveMpacFilePath(mpacFilePath);
    espPop.saveDesignMatrixPath(designMatrixPath);
    espPop.savePredictedVectorPath(predictedPath);
    espPop.saveModelCoefVectorPath(modelCoefPath);

    const double totalCharge = 0.0;
    espPop.exec(totalCharge, PMatrixFilePath);
}


