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
#include <iostream>
#include <string>

#include "DfInputdata.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlUtils.h"

void usage(const std::string& name);

int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "hvo:");

  if (opt["h"] == "defined") {
    usage(opt[0]);
    return EXIT_SUCCESS;
  }

  bool verbose = false;
  if (opt["v"] == "defined") {
    verbose = true;
  }

  std::string outputPath = "init_pdfparam.mpac";
  if (!opt["o"].empty()) {
    outputPath = opt["o"];
  }

  DfInputdata dfInput;
  const bool isReadUserInput = false;
  TlSerializeData param = dfInput.main(isReadUserInput);

  if (verbose) {
    std::cout << TlUtils::format("save: %s", outputPath.c_str()) << std::endl;
  }
  TlMsgPack mpac(param);
  mpac.save(outputPath);

  return EXIT_SUCCESS;
}

void usage(const std::string& name) {
  std::cout << "%s [OPTION]" << std::endl;
  std::cout << "initialize ProteinDF parameters" << std::endl;
  std::cout << std::endl;
  std::cout << "  -o FILE: output MsgPack file" << std::endl;
}
