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

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "Fl_Db_Basis.h"
#include "Fl_Gto.h"
#include "TlGetopt.h"
#include "TlSerializeData.h"
#include "TlUtils.h"

void showHelp(const std::string& progname) {
  std::cout << TlUtils::format("%s [OPTIONS] basisset_name ...",
                               progname.c_str())
            << std::endl;
  std::cout << std::endl;
  std::cout << " display basisset information." << std::endl;
  std::cout << " -v:           verbose" << std::endl;
}

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "hv");
  // const bool isVerbose = (opt["v"] == "defined");
  const bool isShowHelp = (opt["h"] == "defined");

  if (isShowHelp) {
    showHelp(opt[0]);
    return EXIT_SUCCESS;
  }

  if (opt.getCount() <= 1) {
    std::cerr << "please input basisset name" << std::endl;
    std::cerr << std::endl;
    showHelp(opt[0]);
    return EXIT_FAILURE;
  }

  std::string basisName = opt[1];
  Fl_Db_Basis flDbBasis(basisName);

  TlSerializeData data;
  data["name"] = basisName;
  const int numOfCGTOs = flDbBasis.getTotalcgto();
  for (int i = 0; i < numOfCGTOs; ++i) {
    TlSerializeData cGTO;
    cGTO["shell_type"] = TlUtils::format("%c", flDbBasis.getShell(i));
    cGTO["scale_factor"] = flDbBasis.getScalefactor(i);
    const int contractions = flDbBasis.getContraction(i);

    for (int j = 0; j < contractions; ++j) {
      TlSerializeData pGTO;
      pGTO["exp"] = flDbBasis.getExpornent(i, j);
      pGTO["coef"] = flDbBasis.getCoefficient(i, j);
      cGTO["pGTOs"].pushBack(pGTO);
    }

    data["cGTOs"].pushBack(cGTO);
  }

  TlSerializeData gto_data;
  gto_data["x"] = data;

  Fl_Gto gto(gto_data);
  gto.show();

  return EXIT_SUCCESS;
}
