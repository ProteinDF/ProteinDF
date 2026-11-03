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
        text += "max-iteration = 150\n";
        text += "method = rks\n";
        text += "method/rks/electron-number = 10\n";
        text += "method/rks/occlevel = [ 1 - 5 ]\n";
        text += "convergence/type = density\n";
        text += "convergence/threshold = 1e-7\n";
        text += "convergence/threshold-energy = 1e-9\n";
        text += "scf-acceleration = " + scfAcc + "\n";
        text += "scf-acceleration/damping/damping-factor = 0.50\n";
        text += "scf-acceleration/damping/damping-type = " + dampType + "\n";
        text += "xc_functional = hf\n";
        text += "J_engine = CD\n";
        text += "K_engine = CD\n";
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
        text += "  H = \"O-cc-pVDZ-sp.H\"\n";
        text += "  O = \"O-cc-pVDZ-sp.O\"\n";
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
        text += "max-iteration = 150\n";
        text += "method = rks\n";
        text += "method/rks/electron-number = 10\n";
        text += "method/rks/occlevel = [ 1 - 5 ]\n";
        text += "convergence/type = density\n";
        text += "convergence/threshold = 2.0e-4\n";
        text += "convergence/threshold-energy = 1.0e-7\n";
        text += "scf-acceleration = " + scfAcc + "\n";
        text += "scf-acceleration/damping/damping-factor = 0.50\n";
        text += "scf-acceleration/damping/damping-type = " + dampType + "\n";
        text += "xc_functional = BLYP\n";
        text += "J_engine = RI_J\n";
        text += "XC_engine = grid\n";
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
        return text;
    }

}  // namespace

TEST(DiisDensityMixingTest, H2OHfCdDensityMatrixDiisConverges) {
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

    const int lastItr = param["num_of_iterations"].getInt();
    const std::string lastItrKey = TlUtils::format("%d", lastItr);
    const double convergedTE = param["TEs"][lastItrKey].getDouble();

    // Target energy: -76.0195066806 Ha (within 1e-6) and within 1e-8 of baseline damping (-76.0195062160 Ha)
    EXPECT_NEAR(convergedTE, -76.0195066806, 1.0e-6);
    EXPECT_NEAR(convergedTE, -76.0195062160, 1.0e-8);
}

TEST(DiisDensityMixingTest, H2ORiJDftDensityDiisConverges) {
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

    const int lastItr = param["num_of_iterations"].getInt();
    const std::string lastItrKey = TlUtils::format("%d", lastItr);
    const double convergedTE = param["TEs"][lastItrKey].getDouble();

    // Baseline damping energy: -76.3603238893 Ha (within 1e-8: actual diff 7.6e-9)
    EXPECT_NEAR(convergedTE, -76.3603238893, 1.0e-8);
}
