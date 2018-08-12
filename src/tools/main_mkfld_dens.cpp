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
#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlDensField.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "mkfld_common.h"
#include "tl_dense_symmetric_matrix_lapack.h"

void help(const std::string& progName) {
  std::cout << TlUtils::format("%s [options] args...", progName.c_str())
            << std::endl;
  std::cout << "output volume data of density." << std::endl;
  std::cout << "OPTIONS:" << std::endl;
  std::cout
      << " -p PATH      set ProteinDF parameter file. default = pdfparam.mpac"
      << std::endl;
  std::cout << " -d PATH      set density matrix file. default is presumed by "
               "parameter file."
            << std::endl;
  std::cout << " -f           save AVS field (.fld) file." << std::endl;
  std::cout << " -c           save cube (.cube) file" << std::endl;
  std::cout << " -m           save message pack (.mpac) file." << std::endl;
  std::cout << " -o PREFIX    output prefix (default: density)" << std::endl;
  std::cout << " -h           show help message (this)." << std::endl;
  std::cout << " -v           show message verbosely." << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "cd:fhmo:p:v");

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

  std::string outputPrefix = "density";
  if (!opt["o"].empty()) {
    outputPrefix = opt["o"];
  }

  // パラメータファイルの読み込み
  TlSerializeData param;
  {
    TlMsgPack mpac;
    mpac.load(pdfParamPath);
    param = mpac.getSerializeData();
  }

  // 密度行列の読み込み
  TlDenseSymmetricMatrix_Lapack P;
  if (PMatrixFilePath == "") {
    const int iteration = param["num_of_iterations"].getInt();
    DfObject::RUN_TYPE runType = DfObject::RUN_RKS;
    DfObject dfObject(&param);
    PMatrixFilePath = dfObject.getPpqMatrixPath(runType, iteration);
  }
  if (verbose) {
    std::cerr << "loading: " << PMatrixFilePath << std::endl;
  }
  P.load(PMatrixFilePath);

  // 計算サイズの決定
  TlPosition startPos;
  TlPosition endPos;
  getDefaultSize(param, &startPos, &endPos);
  compensateRange(&startPos, &endPos);

  TlPosition gridPitch(GRID_PITCH, GRID_PITCH, GRID_PITCH);
  int numOfGridX, numOfGridY, numOfGridZ;
  std::vector<TlPosition> grids = makeGrids(
      startPos, endPos, gridPitch, &numOfGridX, &numOfGridY, &numOfGridZ);
  const std::size_t numOfGrids = grids.size();

  if (verbose == true) {
    std::cerr << TlUtils::format("start point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                 startPos.x(), startPos.y(), startPos.z())
              << std::endl;
    std::cerr << TlUtils::format("end point(a.u.) = (% 8.3f, % 8.3f, % 8.3f)",
                                 endPos.x(), endPos.y(), endPos.z())
              << std::endl;
  }

  // 電子密度計算
  TlDensField densFld(param);
  const std::vector<double> values = densFld.makeDensFld(P, grids);

  // convert grid unit to angstrom
  // for (std::size_t i = 0; i < numOfGrids; ++i) {
  //     grids[i] *= ANG_PER_AU;
  // }
  // startPos *= ANG_PER_AU;
  // gridPitch *= ANG_PER_AU;

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

    for (std::size_t i = 0; i < numOfGrids; ++i) {
      output["data"]["density"].setAt(i, values[i]);
    }

    TlMsgPack mpac(output);
    mpac.save(mpacFilePath);
  }

  // save to fld
  if (isSaveAvsFieldFile) {
    const std::string fieldFilePath = outputPrefix + ".fld";
    std::cerr << "save AVS field file: " << fieldFilePath << std::endl;

    std::string label = "density";
    std::vector<std::vector<double> > data;
    data.push_back(values);

    saveFieldData(numOfGridX, numOfGridY, numOfGridZ, grids, data, label,
                  fieldFilePath);
  }

  // save to cube file
  if (isSaveCubeFile) {
    const std::string cubeFilePath = outputPrefix + ".cube";
    std::cerr << "save CUBE file: " << cubeFilePath << std::endl;

    const Fl_Geometry flGeom(param["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();
    std::vector<TlAtom> atoms(numOfAtoms);
    for (int i = 0; i < numOfAtoms; ++i) {
      TlAtom atom(flGeom.getAtom(i));
      atom.setCharge(flGeom.getCharge(i));
      atom.moveTo(flGeom.getCoordinate(i));
      // atom.moveTo(flGeom.getCoordinate(i) * ANG_PER_AU); // angstromに変換
      atoms[i] = atom;
    }

    const std::string label = "density";
    saveCubeData(atoms, numOfGridX, numOfGridY, numOfGridZ, startPos, gridPitch,
                 values, label, cubeFilePath);
  }

  return EXIT_SUCCESS;
}
