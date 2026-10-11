#include <climits>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "PdfKeyword.h"
#include "PdfUserInput.h"
#include "TlFile.h"
#include "TlSerializeData.h"
#include "gtest/gtest.h"

namespace {

    void ensurePdfHome() {
        const char* env = std::getenv("PDF_HOME");
        std::string pdfHome = (env != nullptr) ? env : "";
        if (!pdfHome.empty() &&
            TlFile::isExistFile(pdfHome + "/data/pdfkwd.mpac") &&
            TlFile::isExistFile(pdfHome + "/data/alias.mpac")) {
            return;
        }

        const std::vector<std::string> searchBases = {
            pdfHome,
            ".",
            "..",
            "../..",
            "build-check/install",
            "build-regress/install",
            "build/install"};

        const std::vector<std::string> subDirs = {
            "/install",
            "/build-check/install",
            "/build-regress/install",
            ""};

        for (size_t i = 0; i < searchBases.size(); ++i) {
            if (searchBases[i].empty())
                continue;
            for (size_t j = 0; j < subDirs.size(); ++j) {
                std::string candidate = searchBases[i] + subDirs[j];
                if (TlFile::isExistFile(candidate + "/data/pdfkwd.mpac") &&
                    TlFile::isExistFile(candidate + "/data/alias.mpac")) {
                    char resolved[PATH_MAX];
                    if (realpath(candidate.c_str(), resolved) != nullptr) {
                        setenv("PDF_HOME", resolved, 1);
                    } else {
                        setenv("PDF_HOME", candidate.c_str(), 1);
                    }
                    return;
                }
            }
        }
    }

}  // namespace

TEST(PdfXCFunctionalDefault, OmittedXCDefaultB3LYP) {
    ensurePdfHome();

    const std::string filename = "test_userInput_omitted.txt";
    {
        std::ofstream ofs(filename.c_str());
        ofs << ">>>>MAIN\nstep-control = [create]\n>>>>SCF\nmethod = rks\n";
        ofs.close();
    }

    PdfKeyword pdfKwd;
    TlSerializeData param = pdfKwd.getDefault();

    PdfUserInput userInput(filename);
    userInput.load();
    TlSerializeData inputParam = userInput.getSerializeData();
    pdfKwd.convertAlias(&inputParam);

    // Verify omitted detection
    EXPECT_FALSE(inputParam.hasKey("xc_functional"));

    param.merge(inputParam);
    std::remove(filename.c_str());

    EXPECT_EQ(std::string("b3lyp"), param["xc_functional"].getStr());
}

TEST(PdfXCFunctionalDefault, ExplicitXCPotentialPreserved) {
    ensurePdfHome();

    const std::string filename = "test_userInput_svwn.txt";
    {
        std::ofstream ofs(filename.c_str());
        ofs << ">>>>MAIN\nstep-control = [create]\n>>>>SCF\nmethod = rks\nxc-potential = SVWN\n";
        ofs.close();
    }

    PdfKeyword pdfKwd;
    TlSerializeData param = pdfKwd.getDefault();

    PdfUserInput userInput(filename);
    userInput.load();
    TlSerializeData inputParam = userInput.getSerializeData();
    pdfKwd.convertAlias(&inputParam);

    EXPECT_TRUE(inputParam.hasKey("xc_functional"));

    param.merge(inputParam);
    std::remove(filename.c_str());

    EXPECT_EQ(std::string("SVWN"), param["xc_functional"].getStr());
}

TEST(PdfXCFunctionalDefault, ExplicitXCFunctionalPreserved) {
    ensurePdfHome();

    const std::string filename = "test_userInput_hf.txt";
    {
        std::ofstream ofs(filename.c_str());
        ofs << ">>>>MAIN\nstep-control = [create]\n>>>>SCF\nmethod = rks\nxc_functional = HF\n";
        ofs.close();
    }

    PdfKeyword pdfKwd;
    TlSerializeData param = pdfKwd.getDefault();

    PdfUserInput userInput(filename);
    userInput.load();
    TlSerializeData inputParam = userInput.getSerializeData();
    pdfKwd.convertAlias(&inputParam);

    EXPECT_TRUE(inputParam.hasKey("xc_functional"));

    param.merge(inputParam);
    std::remove(filename.c_str());

    EXPECT_EQ(std::string("HF"), param["xc_functional"].getStr());
}
