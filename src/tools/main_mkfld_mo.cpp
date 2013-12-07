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

#include <cassert>
#include <cstdlib>
#include <vector>
#include <map>
#include <fstream>
#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlAtom.h"
#include "TlMoField.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "mkfld_common.h"

void help(const std::string& progName)
{
    std::cout << TlUtils::format("%s [options] <MO index> [<MO index> ...]", progName.c_str()) << std::endl;
    std::cout << std::endl;
    std::cout << "output volume data of molecular orbital" << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p <path>    set ProteinDF parameter file(default: pdfparam.mpac)" << std::endl;
    std::cout << " -l <path>    set LCAO matrix file" << std::endl;
    std::cout << " -f           save AVS field (.fld) file." << std::endl;
    std::cout << " -c           save cube (.cube) file" << std::endl;
    std::cout << " -m           save message pack (.mpac) file." << std::endl;
    std::cout << " -o PREFIX    output prefix (default: MO)" << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "cfhl:mo:p:v");

    const bool verbose = (opt["v"] == "defined");
    if ((opt["h"] == "defined")) {
        help(opt[0]);
        return EXIT_SUCCESS;
    }
    
    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }
    std::string lcaoMatrixPath = "";
    if (opt["l"].empty() != true) {
        lcaoMatrixPath = opt["l"];
    }

    bool isSaveAvsFieldFile = false;
    if (opt["f"] == "defined") {
        isSaveAvsFieldFile = true;
    }

    bool isSaveCubeFile = false;
    if (opt["c"] == "defined") {
        isSaveCubeFile = true;
    }
    
    bool isSaveMpacFile = false;
    if (opt["m"] == "defined") {
        isSaveMpacFile = true;
    }

    if ((!isSaveAvsFieldFile) && (!isSaveCubeFile) && (!isSaveMpacFile)) {
        isSaveMpacFile = true;
    }

    std::string outputPrefix = "MO";
    if (! opt["o"].empty()) {
        outputPrefix = opt["o"];
    }

    // パラメータファイルの読み込み
    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }

    // LCAO行列の読み込み
    TlMatrix C;
    if (lcaoMatrixPath.empty()) {
        const int iteration = param["num_of_iterations"].getInt();
        DfObject::RUN_TYPE runType = DfObject::RUN_RKS;
        DfObject dfObject(&param);
        lcaoMatrixPath = dfObject.getCMatrixPath(runType, iteration);
    }
    std::cerr << "C matrix path: " << lcaoMatrixPath << std::endl;
    C.load(lcaoMatrixPath);

    // 計算サイズの決定
    TlPosition startPos;
    TlPosition endPos;
    getDefaultSize(param, &startPos, &endPos); // a.u.で取得
    compensateRange(&startPos, &endPos);

    const TlPosition gridPitch(GRID_PITCH, GRID_PITCH, GRID_PITCH);
    int numOfGridX, numOfGridY, numOfGridZ;
    std::vector<TlPosition> grids = makeGrids(startPos, endPos, gridPitch,
                                              &numOfGridX, &numOfGridY, &numOfGridZ);
    const std::size_t numOfGrids = grids.size();

    if (verbose == true) {
        std::cerr << TlUtils::format("start point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                     startPos.x(), startPos.y(), startPos.z())
                  << std::endl;
        std::cerr << TlUtils::format("end point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                     endPos.x(), endPos.y(), endPos.z())
                  << std::endl;
    }
    
    // MO計算
    TlMoField moFld(param);
    std::map<std::string, std::vector<double> > storeData;
    // i = 1 から始まる理由はopt[0]がプログラム名のため
    for (int i = 1; i < opt.getCount(); ++i) {
        const int MO_index = std::atoi(opt[i].c_str());
        const std::vector<double>values = moFld.makeMoFld(C.getColVector(MO_index), grids);

        const std::string key = TlUtils::format("%d", MO_index);
        storeData[key] = values;
    }

    // convert grid unit to angstrom
    // for (std::size_t i = 0; i < numOfGrids; ++i) {
    //     grids[i] *= ANG_PER_AU;
    // }
    
    // save to mpac
    if (isSaveMpacFile) {
        const std::string mpacFilePath = outputPrefix + ".mpac";
        std::cerr << "save message pack file: " << mpacFilePath << std::endl;

        TlSerializeData output;
        output["version"] = "2013.0";
        output["num_of_grids"] = numOfGrids;
        for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            TlSerializeData pos;
            pos.pushBack(grids[gridIndex].x() * ANG_PER_AU);
            pos.pushBack(grids[gridIndex].y() * ANG_PER_AU);
            pos.pushBack(grids[gridIndex].z() * ANG_PER_AU);
            
            output["grids"].pushBack(pos);
        }
        output["grid_unit"] = "angstrom";

        for (std::map<std::string, std::vector<double> >::const_iterator p = storeData.begin();
             p != storeData.end(); ++p) {
            const std::string key = p->first;
            const std::vector<double>& value = p->second;

            for (std::size_t i = 0; i < numOfGrids; ++i) {
                output["data"][key].setAt(i, value[i]);
            }
        }

        TlMsgPack mpac(output);
        mpac.save(mpacFilePath);
     }

    // save to fld
    if (isSaveAvsFieldFile) {
        const std::string fieldFilePath = outputPrefix + ".mpac";
        std::cerr << "save AVS field file: " << fieldFilePath << std::endl;

        std::string label = "";
        std::vector<std::vector<double> > data;
        for (std::map<std::string, std::vector<double> >::const_iterator p = storeData.begin();
             p != storeData.end(); ++p) {
            const std::string key = p->first;
            const std::vector<double>& value = p->second;

            label = TlUtils::format("%s %s", label.c_str(), key.c_str());
            data.push_back(value);
        }

        saveFieldData(numOfGridX, numOfGridY, numOfGridZ, grids, data, label, fieldFilePath);
    }

    // save to cube
    if (isSaveCubeFile) {
        const Fl_Geometry flGeom(param["coordinates"]);
        const int numOfAtoms = flGeom.getNumOfAtoms();
        std::vector<TlAtom> atoms(numOfAtoms);
        for (int i = 0; i < numOfAtoms; ++i) {
            TlAtom atom(flGeom.getAtom(i));
            atom.setCharge(flGeom.getCharge(i));
            atom.moveTo(flGeom.getCoordinate(i));
            //atom.moveTo(flGeom.getCoordinate(i) * ANG_PER_AU); // angstroamに変換
            atoms[i] = atom;
        }

        for (std::map<std::string, std::vector<double> >::const_iterator p = storeData.begin();
             p != storeData.end(); ++p) {

            const std::string key = p->first;
            const std::vector<double>& values = p->second;

            const std::string label = TlUtils::format("MO_%s", key.c_str());
            std::vector<double> values_au = values;

            std::string path = TlUtils::format("%s%s.cube",
                                               outputPrefix.c_str(),
                                               key.c_str());
            std::cerr << "save CUBE file: " << path << std::endl;
            saveCubeData(atoms,
                         numOfGridX, numOfGridY, numOfGridZ,
                         startPos, gridPitch,
                         values, label, path);
        }
    }
    
    return EXIT_SUCCESS;
}


