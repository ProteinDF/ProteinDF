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

#include "DfObject.h"
#include "ProteinDF.h"
#include "TlFile.h"
#include "TlLogging.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlUtils.h"
#include "df_total_energy_lapack.h"
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

            char tempTemplate[] = "/tmp/test_scf_edens_home_XXXXXX";
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

            char tempTemplate[] = "/tmp/test_scf_edens_work_XXXXXX";
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

    std::string makeH2OInput(const bool isGridFree, const bool useDensityGuess,
                             const int maxIteration) {
        std::string text;
        text += ">>>>MAIN\n";
        text += "step-control = [create integral guess scf]\n";
        text += ">>>>INTEGRAL\n";
        text += "cut-value = 1.0e-10\n";
        text += ">>>>GUESS\n";
        text += ">>>>SCF\n";
        text += "charge-extrapolate-number = 0\n";
        if (useDensityGuess) {
            text += "scf-start-guess = density\n";
            text += "guess/normalize_density_matrix = no\n";
        } else {
            text += "scf-start-guess = core\n";
        }
        text += "orbital_independence_threshold = 0.00\n";
        text += "orbital_independence_threshold/canonical = 0.00\n";
        text += "orbital_independence_threshold/lowdin = 0.00\n";
        text += "orbital-overlap-correspondence = off\n";
        text += "orbital-overlap-correspondence-first = off\n";
        text += TlUtils::format("max-iteration = %d\n", maxIteration);
        text += "method = rks\n";
        text += "method/rks/electron-number = 10\n";
        text += "method/rks/occlevel = [ 1 - 5 ]\n";
        text += "convergence/type = density\n";
        text += "convergence/threshold = 1e-7\n";
        text += "convergence/threshold-energy = 1e-9\n";
        text += "scf-acceleration = damping\n";
        text += "scf-acceleration/damping/damping-factor = 0.50\n";
        text += "scf-acceleration/damping/damping-type = density_matrix\n";
        text += "xc-potential = BLYP\n";
        text += "J_engine = CD\n";
        text += "K_engine = CD\n";
        if (isGridFree) {
            text += "XC_engine = gridfree_CD\n";
            text += "gridfree/dedicated_basis = true\n";
        } else {
            text += "XC_engine = grid\n";
        }
        text += "CDAM_tau = 1.0E-10\n";
        text += "CD_epsilon = 1.0E-6\n";
        text += "cleanup = no\n";
        text += ">>>>MOLECULE\n";
        text += "geometry/cartesian/unit = angstrom\n";
        text += "geometry/cartesian/input = {\n";
        text += "  O 0.00000000 0.00000000 0.00000000\n";
        text += "  H 1.00000000 0.00000000 0.00000000\n";
        text += "  H 0.00000000 1.00000000 0.00000000\n";
        text += "}end\n";
        text += "basis-set/orbital = {\n";
        text += "  H = \"O-6-31G.H\"\n";
        text += "  O = \"O-6-31G.O\"\n";
        text += "}end\n";
        text += "basis-set/density-auxiliary = {\n";
        text += "  H = \"A-DZVP2.H\"\n";
        text += "  O = \"A-DZVP2.O\"\n";
        text += "}end\n";
        text += "basis-set/exchange-auxiliary = {\n";
        text += "  H = \"A-DZVP2.H\"\n";
        text += "  O = \"A-DZVP2.O\"\n";
        text += "}end\n";
        if (isGridFree) {
            text += "basis-set/gridfree = {\n";
            text += "  H = \"O-cc-pVDZ-sp.H\"\n";
            text += "  O = \"O-cc-pVDZ-sp.O\"\n";
            text += "}end\n";
        }
        return text;
    }

    TlSerializeData loadPdfParam(const std::string& path) {
        TlMsgPack mpac;
        mpac.load(path);
        return mpac.getSerializeData();
    }

    class TestDfObject : public DfObject {
    public:
        explicit TestDfObject(TlSerializeData* pPdfParam)
            : DfObject(pPdfParam) {}

        using DfObject::getExcMatrix;
        using DfObject::getHpqMatrix;
        using DfObject::getJMatrix;
        using DfObject::getPInMatrix;
    };

    void verifyFirstIterationFromConvergedDensity(const bool isGridFree) {
        ScopedPdfHomeWithBasis scopedPdfHome;

        double convergedTE = 0.0;
        std::string convergedDensityCopyPath;

        char savedDensityTemplate[] = "/tmp/test_scf_edens_dens_XXXXXX";
        char* savedDir = mkdtemp(savedDensityTemplate);
        ASSERT_TRUE(savedDir != NULL);
        const std::string savedDirStr(savedDir);
        convergedDensityCopyPath = savedDirStr + "/guess.density.rks.mat";

        {
            ScopedTempWorkDir workDir;
            ASSERT_FALSE(workDir.path().empty());

            {
                std::ofstream ofs("fl_Userinput");
                ofs << makeH2OInput(isGridFree, false, 150);
            }

            ProteinDF pdf;
            pdf.run();

            TlSerializeData param = loadPdfParam("pdfparam.mpac");
            ASSERT_TRUE(param["control"]["scf_converged"].getBoolean());

            const int lastItr = param["num_of_iterations"].getInt();
            ASSERT_GT(lastItr, 1);

            const std::string lastItrKey = TlUtils::xtos(lastItr);
            convergedTE = param["TEs"][lastItrKey].getDouble();

            TestDfObject dfObj(&param);
            const std::string pOutPath =
                dfObj.getPOutMatrixPath(DfObject::RUN_RKS, lastItr);
            ASSERT_TRUE(TlFile::isExistFile(pOutPath));
            ASSERT_EQ(0, TlFile::copy(pOutPath, convergedDensityCopyPath));

            // At iteration 1 (core guess), P_in(1) != P_out(1); verify energy terms use P_in(1)
            {
                const TlDenseSymmetricMatrix_Lapack Hpq =
                    dfObj.getHpqMatrix<TlDenseSymmetricMatrix_Lapack>();
                const TlDenseSymmetricMatrix_Lapack PIn1 =
                    dfObj.getPInMatrix<TlDenseSymmetricMatrix_Lapack>(
                        DfObject::RUN_RKS, 1);
                const TlDenseSymmetricMatrix_Lapack J1 =
                    dfObj.getJMatrix<TlDenseSymmetricMatrix_Lapack>(1);

                TlDenseSymmetricMatrix_Lapack Eh1 = Hpq;
                Eh1.dotInPlace(PIn1);
                EXPECT_NEAR(Eh1.sum(),
                            param["total_energy"]["1"]["h"].getDouble(),
                            1.0e-12);

                TlDenseSymmetricMatrix_Lapack Ej1 = J1;
                Ej1.dotInPlace(PIn1);
                EXPECT_NEAR(0.5 * Ej1.sum(),
                            param["total_energy"]["1"]["J"].getDouble(),
                            1.0e-12);

                if (isGridFree) {
                    TlDenseSymmetricMatrix_Lapack Exc1 =
                        dfObj.getExcMatrix<TlDenseSymmetricMatrix_Lapack>(
                            DfObject::RUN_RKS, 1);
                    Exc1.dotInPlace(PIn1);
                    EXPECT_NEAR(Exc1.sum(),
                                param["total_energy"]["1"]["XC"].getDouble(),
                                1.0e-12);
                }
            }

            // Verify post-SCF energy component recalculation (used by pdf-info-ecomp)
            TlSerializeData recalcParam = param;
            {
                DfTotalEnergy_Lapack recalc(&recalcParam);
                recalc.run();
                EXPECT_NEAR(convergedTE,
                            recalcParam["TEs"][lastItrKey].getDouble(), 1.0e-12);
            }
        }

        double firstIterTE = 0.0;
        {
            ScopedTempWorkDir workDir1;
            ASSERT_FALSE(workDir1.path().empty());

            ASSERT_EQ(0,
                      TlFile::copy(convergedDensityCopyPath, "./guess.density.rks.mat"));
            {
                std::ofstream ofs("fl_Userinput");
                ofs << makeH2OInput(isGridFree, true, 1);
            }

            ProteinDF pdf1;
            pdf1.run();

            const TlSerializeData param1 = loadPdfParam("pdfparam.mpac");
            ASSERT_EQ(1, param1["num_of_iterations"].getInt());
            firstIterTE = param1["TEs"]["1"].getDouble();
        }

        removeDirectoryRecursive(savedDirStr);

        EXPECT_NEAR(convergedTE, firstIterTE, 1.0e-8);
    }

}  // namespace

TEST(ScfEnergyDensityConsistency, H2OGridFirstIterMatchesConvergedEnergy) {
    verifyFirstIterationFromConvergedDensity(false);
}

TEST(ScfEnergyDensityConsistency, H2OGridFreeFirstIterMatchesConvergedEnergy) {
    verifyFirstIterationFromConvergedDensity(true);
}
