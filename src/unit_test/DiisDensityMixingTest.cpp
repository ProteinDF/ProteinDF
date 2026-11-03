// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include <dirent.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

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

    // DfTotalEnergy_tmpl caches E_nuc_ across runs in a static variable;
    // reset it before each test so a prior calculation with different geometry
    // (e.g. H2) does not corrupt the nuclear repulsion energy of H2O.
    struct DfTotalEnergyResetter : public DfTotalEnergy_Lapack {
        static void reset() {
            E_nuc_ = 0.0;
            E_nuc_X_ = 0.0;
        }
    };

    void removeDirectoryRecursive(const std::string& dirPath) {
        DIR* dir = opendir(dirPath.c_str());
        if (dir == NULL) {
            return;
        }

        struct dirent* entry = NULL;
        while ((entry = readdir(dir)) != NULL) {
            if ((std::strcmp(entry->d_name, ".") == 0) || (std::strcmp(entry->d_name, "..") == 0)) {
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

    // Limit OpenMP threads for small molecular unit tests to avoid heavy
    // initialization overhead of thread-local ERI engines.
    class ScopedOmpNumThreads {
    public:
        explicit ScopedOmpNumThreads(int threads) {
#ifdef _OPENMP
            original_ = omp_get_max_threads();
            omp_set_num_threads(threads);
#else
            (void)threads;
#endif
        }

        ~ScopedOmpNumThreads() {
#ifdef _OPENMP
            omp_set_num_threads(original_);
#endif
        }

    private:
#ifdef _OPENMP
        int original_ = 1;
#endif
    };

    class ScopedPdfHomeWithBasis {
    public:
        ScopedPdfHomeWithBasis() {
            const char* current = std::getenv("PDF_HOME");
            if (current != NULL) {
                hasOriginal_ = true;
                original_ = current;
            }

            char tempTemplate[] = "/tmp/test_diis_mixing_home_XXXXXX";
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

            char tempTemplate[] = "/tmp/test_diis_mixing_work_XXXXXX";
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

    std::string makeH2OHfCdInput(const std::string& scfAcc, const std::string& dampType) {
        std::string text;
        text += ">>>>MAIN\n";
        text += "step-control = [create integral guess scf]\n";
        text += ">>>>INTEGRAL\n";
        text += "cut-value = 1.0e-9\n";
        text += ">>>>GUESS\n";
        text += ">>>>SCF\n";
        text += "charge-extrapolate-number = 0\n";
        text += "scf-start-guess = core\n";
        text += "orbital_independence_threshold = 0.00\n";
        text += "orbital_independence_threshold/canonical = 0.00\n";
        text += "orbital_independence_threshold/lowdin = 0.00\n";
        text += "orbital-overlap-correspondence = off\n";
        text += "orbital-overlap-correspondence-first = off\n";
        text += "max-iteration = 30\n";
        text += "method = rks\n";
        text += "method/rks/electron-number = 10\n";
        text += "method/rks/occlevel = [ 1 - 5 ]\n";
        text += "convergence/type = density\n";
        text += "convergence/threshold = 1.0e-3\n";
        text += "convergence/threshold-energy = 1.0e-5\n";
        text += "scf-acceleration = " + scfAcc + "\n";
        text += "scf-acceleration/damping/damping-factor = 0.50\n";
        text += "scf-acceleration/damping/damping-type = " + dampType + "\n";
        text += "xc_functional = hf\n";
        text += "J_engine = CD\n";
        text += "K_engine = CD\n";
        text += "CDAM_tau = 1.0E-9\n";
        text += "CD_epsilon = 1.0E-5\n";
        text += "cleanup = no\n";
        text += ">>>>MOLECULE\n";
        text += "geometry/cartesian/unit = angstrom\n";
        text += "geometry/cartesian/input = {\n";
        text += "  O 0.00000000 0.00000000 0.00000000\n";
        text += "  H 1.00000000 0.00000000 0.00000000\n";
        text += "  H 0.00000000 1.00000000 0.00000000\n";
        text += "}end\n";
        text += "basis-set/orbital = {\n";
        text += "  H = \"O-STO-3G.H\"\n";
        text += "  O = \"O-STO-3G.O\"\n";
        text += "}end\n";
        text += "basis-set/density-auxiliary = {\n";
        text += "  H = \"A-DZVP2.H\"\n";
        text += "  O = \"A-DZVP2.O\"\n";
        text += "}end\n";
        text += "basis-set/exchange-auxiliary = {\n";
        text += "  H = \"A-DZVP2.H\"\n";
        text += "  O = \"A-DZVP2.O\"\n";
        text += "}end\n";
        return text;
    }

    std::string makeH2ORiJDftInput(const std::string& scfAcc, const std::string& dampType) {
        std::string text;
        text += ">>>>MAIN\n";
        text += "step-control = [create integral guess scf]\n";
        text += ">>>>INTEGRAL\n";
        text += "cut-value = 1.0e-8\n";
        text += ">>>>GUESS\n";
        text += ">>>>SCF\n";
        text += "charge-extrapolate-number = 0\n";
        text += "scf-start-guess = core\n";
        text += "orbital_independence_threshold = 0.00\n";
        text += "orbital_independence_threshold/canonical = 0.00\n";
        text += "orbital_independence_threshold/lowdin = 0.00\n";
        text += "orbital-overlap-correspondence = off\n";
        text += "orbital-overlap-correspondence-first = off\n";
        text += "max-iteration = 25\n";
        text += "method = rks\n";
        text += "method/rks/electron-number = 10\n";
        text += "method/rks/occlevel = [ 1 - 5 ]\n";
        text += "convergence/type = density\n";
        text += "convergence/threshold = 1.0e-3\n";
        text += "convergence/threshold-energy = 1.0e-5\n";
        text += "scf-acceleration = " + scfAcc + "\n";
        text += "scf-acceleration/damping/damping-factor = 0.50\n";
        text += "scf-acceleration/damping/damping-type = " + dampType + "\n";
        text += "scf_acceleration/diis/start = 2\n";
        text += "xc_functional = hf\n";
        text += "J_engine = RI_J\n";
        text += "K_engine = CD\n";
        text += "CDAM_tau = 1.0E-7\n";
        text += "CD_epsilon = 1.0E-4\n";
        text += "cleanup = no\n";
        text += ">>>>MOLECULE\n";
        text += "geometry/cartesian/unit = angstrom\n";
        text += "geometry/cartesian/input = {\n";
        text += "  O 0.00000000 0.00000000 0.00000000\n";
        text += "  H 1.00000000 0.00000000 0.00000000\n";
        text += "  H 0.00000000 1.00000000 0.00000000\n";
        text += "}end\n";
        text += "basis-set/orbital = {\n";
        text += "  H = \"O-STO-3G.H\"\n";
        text += "  O = \"O-STO-3G.O\"\n";
        text += "}end\n";
        text += "basis-set/density-auxiliary = {\n";
        text += "  H = \"A-DZVP.H\"\n";
        text += "  O = \"A-DZVP.O\"\n";
        text += "}end\n";
        text += "basis-set/exchange-auxiliary = {\n";
        text += "  H = \"A-DZVP.H\"\n";
        text += "  O = \"A-DZVP.O\"\n";
        text += "}end\n";
        return text;
    }

}  // namespace

TEST(DiisDensityMixingTest, H2OHfCdDensityMatrixDiisConverges) {
    ScopedOmpNumThreads scopedOmp(1);
    DfTotalEnergyResetter::reset();
    ScopedPdfHomeWithBasis scopedPdfHome;
    ScopedTempWorkDir workDir;
    ASSERT_FALSE(workDir.path().empty());

    {
        std::ofstream ofs("fl_Userinput");
        ofs << makeH2OHfCdInput("diis", "density_matrix");
    }

    ProteinDF pdf;
    pdf.run();

    const TlSerializeData param = loadPdfParam("pdfparam.mpac");
    EXPECT_TRUE(param["control"]["scf_converged"].getBoolean());

    // Verify DIIS was actually activated
    EXPECT_TRUE(param.hasKey("internal/dfconverge_damping_diis/e_max_pass_itr"));
    EXPECT_NE(-1, param["internal/dfconverge_damping_diis/e_max_pass_itr"].getInt());

    const int lastItr = param["num_of_iterations"].getInt();
    // DIIS converges faster than damping-only (20 iterations)
    EXPECT_LT(lastItr, 20);
    const std::string lastItrKey = TlUtils::format("%d", lastItr);
    const double convergedTE = param["TEs"][lastItrKey].getDouble();

    // Baseline damping energy: -74.9611654993 Ha (actual diff with DIIS is 6.0e-6 Ha)
    EXPECT_NEAR(convergedTE, -74.9611654993, 1.0e-5);
}

TEST(DiisDensityMixingTest, H2ORiJDftDensityDiisConverges) {
    ScopedOmpNumThreads scopedOmp(1);
    DfTotalEnergyResetter::reset();
    ScopedPdfHomeWithBasis scopedPdfHome;
    ScopedTempWorkDir workDir;
    ASSERT_FALSE(workDir.path().empty());

    {
        std::ofstream ofs("fl_Userinput");
        ofs << makeH2ORiJDftInput("diis", "density");
    }

    ProteinDF pdf;
    pdf.run();

    const TlSerializeData param = loadPdfParam("pdfparam.mpac");
    EXPECT_TRUE(param["control"]["scf_converged"].getBoolean());

    // Verify DIIS was actually activated for rho
    EXPECT_TRUE(param.hasKey("internal/dfconverge_damping_diis/rho_e_max_pass_itr"));
    EXPECT_NE(-1, param["internal/dfconverge_damping_diis/rho_e_max_pass_itr"].getInt());

    const int lastItr = param["num_of_iterations"].getInt();
    // DIIS converges faster than damping-only (17 iterations)
    EXPECT_LT(lastItr, 17);
    const std::string lastItrKey = TlUtils::format("%d", lastItr);
    const double convergedTE = param["TEs"][lastItrKey].getDouble();

    // Baseline damping energy: -74.9684817181 Ha (actual diff with DIIS is 4.4e-7 Ha)
    EXPECT_NEAR(convergedTE, -74.9684817181, 1.0e-6);
}
