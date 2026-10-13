#include <sys/stat.h>
#include <unistd.h>

#include <cstdlib>
#include <fstream>
#include <string>

#include "DfInputdata.h"
#include "PdfKeyword.h"
#include "PdfUserInput.h"
#include "TlFile.h"
#include "TlSerializeData.h"
#include "gtest/gtest.h"

namespace {

    class ScopedPdfHome {
    public:
        ScopedPdfHome() {
            const char* current = std::getenv("PDF_HOME");
            if (current != nullptr) {
                hasOriginal_ = true;
                original_ = current;
            } else {
                hasOriginal_ = false;
            }

#if defined(PDFKWD_MPAC_PATH) && defined(ALIAS_MPAC_PATH)
            char tempTemplate[] = "/tmp/test_pdfhome_XXXXXX";
            char* dir = mkdtemp(tempTemplate);
            if (dir != nullptr) {
                tempDir_ = dir;
                std::string dataDir = tempDir_ + "/data";
                mkdir(dataDir.c_str(), 0755);

                if (TlFile::isExistFile(PDFKWD_MPAC_PATH)) {
                    symlink(PDFKWD_MPAC_PATH, (dataDir + "/pdfkwd.mpac").c_str());
                }
                if (TlFile::isExistFile(ALIAS_MPAC_PATH)) {
                    symlink(ALIAS_MPAC_PATH, (dataDir + "/alias.mpac").c_str());
                }

                setenv("PDF_HOME", tempDir_.c_str(), 1);
            }
#endif
        }

        ~ScopedPdfHome() {
            if (hasOriginal_) {
                setenv("PDF_HOME", original_.c_str(), 1);
            } else {
                unsetenv("PDF_HOME");
            }

            if (!tempDir_.empty()) {
                std::string dataDir = tempDir_ + "/data";
                unlink((dataDir + "/pdfkwd.mpac").c_str());
                unlink((dataDir + "/alias.mpac").c_str());
                rmdir(dataDir.c_str());
                rmdir(tempDir_.c_str());
            }
        }

    private:
        bool hasOriginal_ = false;
        std::string original_;
        std::string tempDir_;
    };

}  // namespace

TEST(PdfXCFunctionalDefault, OmittedXCDefaultB3LYP) {
    ScopedPdfHome scopedPdfHome;

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
    ScopedPdfHome scopedPdfHome;

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
    ScopedPdfHome scopedPdfHome;

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

// Tests for scf-memory-saving override:
// ~なし・省略(既定b3lyp)はyes、~付きはnoのまま、~付きyes指定はyes
TEST(DfInputdataCheck, ScfMemorySavingOverriddenToYesWhenOmitted) {
    ScopedPdfHome scopedPdfHome;

    const std::string filename = "test_userInput_memsave_omitted.txt";
    {
        std::ofstream ofs(filename.c_str());
        ofs << ">>>>MAIN\nstep-control = [create]\n>>>>SCF\nmethod = rks\nscf-memory-saving = no\n";
        ofs.close();
    }

    PdfKeyword pdfKwd;
    TlSerializeData param = pdfKwd.getDefault();

    PdfUserInput userInput(filename);
    userInput.load();
    TlSerializeData inputParam = userInput.getSerializeData();
    pdfKwd.convertAlias(&inputParam);

    param.merge(inputParam);
    pdfKwd.convertAlias(&param);
    DfInputdata::fixScfMemorySaving(&param);
    std::remove(filename.c_str());

    EXPECT_EQ(std::string("yes"), param["scf-memory-saving"].getStr());
}

TEST(DfInputdataCheck, ScfMemorySavingOverriddenToYesWithNormalXC) {
    ScopedPdfHome scopedPdfHome;

    const std::string filename = "test_userInput_memsave_normal.txt";
    {
        std::ofstream ofs(filename.c_str());
        ofs << ">>>>MAIN\nstep-control = [create]\n>>>>SCF\nmethod = rks\nscf-memory-saving = no\nxc-potential = b3lyp\n";
        ofs.close();
    }

    PdfKeyword pdfKwd;
    TlSerializeData param = pdfKwd.getDefault();

    PdfUserInput userInput(filename);
    userInput.load();
    TlSerializeData inputParam = userInput.getSerializeData();
    pdfKwd.convertAlias(&inputParam);

    param.merge(inputParam);
    pdfKwd.convertAlias(&param);
    DfInputdata::fixScfMemorySaving(&param);
    std::remove(filename.c_str());

    EXPECT_EQ(std::string("yes"), param["scf-memory-saving"].getStr());
}

TEST(DfInputdataCheck, ScfMemorySavingPreservedNoWithTildeXC) {
    ScopedPdfHome scopedPdfHome;

    const std::string filename = "test_userInput_memsave_tilde_no.txt";
    {
        std::ofstream ofs(filename.c_str());
        ofs << ">>>>MAIN\nstep-control = [create]\n>>>>SCF\nmethod = rks\nscf-memory-saving = no\nxc-potential = svwn~\n";
        ofs.close();
    }

    PdfKeyword pdfKwd;
    TlSerializeData param = pdfKwd.getDefault();

    PdfUserInput userInput(filename);
    userInput.load();
    TlSerializeData inputParam = userInput.getSerializeData();
    pdfKwd.convertAlias(&inputParam);

    param.merge(inputParam);
    pdfKwd.convertAlias(&param);
    DfInputdata::fixScfMemorySaving(&param);
    std::remove(filename.c_str());

    EXPECT_EQ(std::string("no"), param["scf-memory-saving"].getStr());
}

TEST(DfInputdataCheck, ScfMemorySavingPreservedYesWithTildeXC) {
    ScopedPdfHome scopedPdfHome;

    const std::string filename = "test_userInput_memsave_tilde_yes.txt";
    {
        std::ofstream ofs(filename.c_str());
        ofs << ">>>>MAIN\nstep-control = [create]\n>>>>SCF\nmethod = rks\nscf-memory-saving = yes\nxc-potential = svwn~\n";
        ofs.close();
    }

    PdfKeyword pdfKwd;
    TlSerializeData param = pdfKwd.getDefault();

    PdfUserInput userInput(filename);
    userInput.load();
    TlSerializeData inputParam = userInput.getSerializeData();
    pdfKwd.convertAlias(&inputParam);

    param.merge(inputParam);
    pdfKwd.convertAlias(&param);
    DfInputdata::fixScfMemorySaving(&param);
    std::remove(filename.c_str());

    EXPECT_EQ(std::string("yes"), param["scf-memory-saving"].getStr());
}
