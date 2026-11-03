#include <dirent.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include "DfInputdata.h"
#include "PdfKeyword.h"
#include "PdfUserInput.h"
#include "ProteinDF.h"
#include "TlFile.h"
#include "TlLogging.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlUtils.h"
#include "gtest/gtest.h"

namespace {

    void removeDirectoryRecursive(const std::string& dirPath) {
        DIR* dir = opendir(dirPath.c_str());
        if (dir == NULL) {
            return;
        }

        struct dirent* entry = NULL;
        while ((entry = readdir(dir)) != NULL) {
            if ((std::strcmp(entry->d_name, ".") == 0) ||
                (std::strcmp(entry->d_name, "..") == 0)) {
                continue;
            }

            const std::string childPath = dirPath + "/" + entry->d_name;
            struct stat st;
            if (lstat(childPath.c_str(), &st) == 0) {
                if (S_ISDIR(st.st_mode)) {
                    removeDirectoryRecursive(childPath);
                } else {
                    unlink(childPath.c_str());
                }
            }
        }
        closedir(dir);
        rmdir(dirPath.c_str());
    }

    class ScopedPdfHomeWithBasis {
    public:
        ScopedPdfHomeWithBasis() {
            const char* current = std::getenv("PDF_HOME");
            if (current != NULL) {
                hasOriginal_ = true;
                original_ = current;
            }

#if defined(PDFKWD_MPAC_PATH) && defined(ALIAS_MPAC_PATH) && defined(BASIS2_PATH)
            char tempTemplate[] = "/tmp/test_pdfkwd_home_XXXXXX";
            char* dir = mkdtemp(tempTemplate);
            if (dir != NULL) {
                tempDir_ = dir;
                const std::string dataDir = tempDir_ + "/data";
                mkdir(dataDir.c_str(), 0755);

                symlink(PDFKWD_MPAC_PATH, (dataDir + "/pdfkwd.mpac").c_str());
                symlink(ALIAS_MPAC_PATH, (dataDir + "/alias.mpac").c_str());
                symlink(BASIS2_PATH, (dataDir + "/basis2").c_str());

                setenv("PDF_HOME", tempDir_.c_str(), 1);
            }
#endif
        }

        ~ScopedPdfHomeWithBasis() {
            if (hasOriginal_) {
                setenv("PDF_HOME", original_.c_str(), 1);
            } else {
                unsetenv("PDF_HOME");
            }

            if (!tempDir_.empty()) {
                removeDirectoryRecursive(tempDir_);
            }
        }

    private:
        bool hasOriginal_ = false;
        std::string original_;
        std::string tempDir_;
    };

    class ScopedTempWorkDir {
    public:
        ScopedTempWorkDir() {
            char cwdBuf[PATH_MAX];
            if (getcwd(cwdBuf, sizeof(cwdBuf)) != NULL) {
                originalCwd_ = cwdBuf;
            }

            char tempTemplate[] = "/tmp/test_pdfkwd_work_XXXXXX";
            char* dir = mkdtemp(tempTemplate);
            if (dir != NULL) {
                tempDir_ = dir;
                mkdir((tempDir_ + "/fl_Work").c_str(), 0755);
                chdir(tempDir_.c_str());
                TlLogging::getInstance().setFilePath(tempDir_ + "/fl_Out_Std");
            }
        }

        ~ScopedTempWorkDir() {
            TlLogging::getInstance().setFilePath("output.log");
            if (!originalCwd_.empty()) {
                chdir(originalCwd_.c_str());
            }
            if (!tempDir_.empty()) {
                removeDirectoryRecursive(tempDir_);
            }
        }

        const std::string& path() const {
            return tempDir_;
        }

    private:
        std::string originalCwd_;
        std::string tempDir_;
    };

    TlSerializeData loadPdfParam(const std::string& path) {
        TlMsgPack mpac;
        mpac.load(path);
        return mpac.getSerializeData();
    }

    std::string makeH2Input(const std::string& dampingType,
                            const std::string& extraKeywords = "") {
        std::string text;
        text += ">>>>MAIN\n";
        text += "step-control = [create integral guess scf]\n";
        text += ">>>>INTEGRAL\n";
        text += "cut-value = 1.0e-10\n";
        text += ">>>>GUESS\n";
        text += ">>>>SCF\n";
        text += "charge-extrapolate-number = 0\n";
        text += "scf-start-guess = core\n";
        text += "orbital_independence_threshold = 0.00\n";
        text += "orbital_independence_threshold/canonical = 0.00\n";
        text += "orbital_independence_threshold/lowdin = 0.00\n";
        text += "orbital-overlap-correspondence = off\n";
        text += "orbital-overlap-correspondence-first = off\n";
        text += "max-iteration = 50\n";
        text += "method = rks\n";
        text += "method/rks/electron-number = 2\n";
        text += "method/rks/occlevel = [ 1 ]\n";
        text += "convergence/type = density\n";
        text += "convergence/threshold = 1e-6\n";
        text += "convergence/threshold-energy = 1e-8\n";
        text += "scf-acceleration = damping\n";
        text += "scf-acceleration/damping/damping-factor = 0.67\n";
        text += "scf-acceleration/damping/damping-type = " + dampingType + "\n";
        text += "xc-potential = SVWN\n";
        text += "J_engine = CD\n";
        text += "K_engine = CD\n";
        text += "XC_engine = grid\n";
        text += "CDAM_tau = 1.0E-10\n";
        text += "CD_epsilon = 1.0E-6\n";
        if (!extraKeywords.empty()) {
            text += extraKeywords;
        }
        text += ">>>>MOLECULE\n";
        text += "geometry/cartesian/unit = angstrom\n";
        text += "geometry/cartesian/input = {\n";
        text += "  H 0.00 0.00 0.00\n";
        text += "  H 0.74 0.00 0.00\n";
        text += "}end\n";
        text += "basis-set/orbital = {\n";
        text += "  H = \"O-STO-3G.H\"\n";
        text += "}end\n";
        text += "basis-set/density-auxiliary = {\n";
        text += "  H = \"A-DZVP.H\"\n";
        text += "}end\n";
        text += "basis-set/exchange-auxiliary = {\n";
        text += "  H = \"A-DZVP.H\"\n";
        text += "}end\n";
        return text;
    }

}  // namespace

// (1) 古い綴り3件が別名で新しい綴りになり、辞書の既定値がコードの読むキーに入ること
TEST(PdfKeywordSpelling, DictionaryDefaultsAndAliasConversion) {
    ScopedPdfHomeWithBasis scopedPdfHome;

    PdfKeyword pdfKwd;
    TlSerializeData defParam = pdfKwd.getDefault();

    // 辞書の既定値がコードの読むキーに入ること
    EXPECT_TRUE(defParam.hasKey("charge_extrapolate_number"));
    EXPECT_EQ("1", defParam["charge_extrapolate_number"].getStr());
    EXPECT_EQ(1, defParam["charge_extrapolate_number"].getInt());

    EXPECT_TRUE(defParam.hasKey("scf_acceleration/damping/start"));
    EXPECT_EQ("0", defParam["scf_acceleration/damping/start"].getStr());
    EXPECT_EQ(0, defParam["scf_acceleration/damping/start"].getInt());

    // 1. charge-extrapolate-number -> charge_extrapolate_number
    {
        TlSerializeData input;
        input["charge-extrapolate-number"] = "2";
        pdfKwd.convertAlias(&input);
        EXPECT_FALSE(input.hasKey("charge-extrapolate-number"));
        EXPECT_TRUE(input.hasKey("charge_extrapolate_number"));
        EXPECT_EQ("2", input["charge_extrapolate_number"].getStr());

        TlSerializeData param = defParam;
        param.merge(input);
        EXPECT_EQ(2, param["charge_extrapolate_number"].getInt());
    }

    // 2. scf_acceleration/daming/start -> scf_acceleration/damping/start
    {
        TlSerializeData input;
        input["scf_acceleration/daming/start"] = "3";
        pdfKwd.convertAlias(&input);
        EXPECT_FALSE(input.hasKey("scf_acceleration/daming/start"));
        EXPECT_TRUE(input.hasKey("scf_acceleration/damping/start"));
        EXPECT_EQ("3", input["scf_acceleration/damping/start"].getStr());

        TlSerializeData param = defParam;
        param.merge(input);
        EXPECT_EQ(3, param["scf_acceleration/damping/start"].getInt());
    }

    // 3. scf-acceleration/damping/start -> scf_acceleration/damping/start
    {
        TlSerializeData input;
        input["scf-acceleration/damping/start"] = "4";
        pdfKwd.convertAlias(&input);
        EXPECT_FALSE(input.hasKey("scf-acceleration/damping/start"));
        EXPECT_TRUE(input.hasKey("scf_acceleration/damping/start"));
        EXPECT_EQ("4", input["scf_acceleration/damping/start"].getStr());

        TlSerializeData param = defParam;
        param.merge(input);
        EXPECT_EQ(4, param["scf_acceleration/damping/start"].getInt());
    }
}

// (2) ダンピングの種類 density_matrix・density・fock で、開始0の既定のまま、修正前と同じ収束値になること
TEST(PdfKeywordSpelling, DampingDensityMatrixMatchesPreFixEnergy) {
    ScopedPdfHomeWithBasis scopedPdfHome;
    ScopedTempWorkDir workDir;
    ASSERT_FALSE(workDir.path().empty());

    {
        std::ofstream ofs("fl_Userinput");
        ofs << makeH2Input("density_matrix");
    }

    ProteinDF pdf;
    pdf.run();

    TlSerializeData param = loadPdfParam("pdfparam.mpac");
    ASSERT_TRUE(param["control"]["scf_converged"].getBoolean());

    const int lastItr = param["num_of_iterations"].getInt();
    const std::string lastItrKey = TlUtils::xtos(lastItr);
    const double convergedTE = param["TEs"][lastItrKey].getDouble();

    const double expectedTE = -1.12120912970981;
    EXPECT_NEAR(convergedTE, expectedTE, 1e-8);
}

TEST(PdfKeywordSpelling, DampingDensityMatchesPreFixEnergy) {
    ScopedPdfHomeWithBasis scopedPdfHome;
    ScopedTempWorkDir workDir;
    ASSERT_FALSE(workDir.path().empty());

    {
        std::ofstream ofs("fl_Userinput");
        ofs << makeH2Input("density");
    }

    ProteinDF pdf;
    pdf.run();

    TlSerializeData param = loadPdfParam("pdfparam.mpac");
    ASSERT_TRUE(param["control"]["scf_converged"].getBoolean());

    const int lastItr = param["num_of_iterations"].getInt();
    const std::string lastItrKey = TlUtils::xtos(lastItr);
    const double convergedTE = param["TEs"][lastItrKey].getDouble();

    const double expectedTE = -1.12120912970981;
    EXPECT_NEAR(convergedTE, expectedTE, 1e-8);
}

TEST(PdfKeywordSpelling, DampingFockMatchesPreFixEnergy) {
    ScopedPdfHomeWithBasis scopedPdfHome;
    ScopedTempWorkDir workDir;
    ASSERT_FALSE(workDir.path().empty());

    {
        std::ofstream ofs("fl_Userinput");
        ofs << makeH2Input("fock");
    }

    ProteinDF pdf;
    pdf.run();

    TlSerializeData param = loadPdfParam("pdfparam.mpac");
    ASSERT_TRUE(param["control"]["scf_converged"].getBoolean());

    const int lastItr = param["num_of_iterations"].getInt();
    const std::string lastItrKey = TlUtils::xtos(lastItr);
    const double convergedTE = param["TEs"][lastItrKey].getDouble();

    const double expectedTE = -1.12120912970981;
    EXPECT_NEAR(convergedTE, expectedTE, 1e-8);
}
