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
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TlGetopt.h"
#include "TlStringTokenizer.h"
#include "TlUtils.h"

#define EPS 1.0E-5

std::map<std::string, std::string> getAvsFieldDataHeader(std::ifstream* ifs) {
    assert(ifs != NULL);
    std::map<std::string, std::string> param;

    bool isFF = false;
    std::stringstream os;
    do {
        char c;
        ifs->read(&c, sizeof(char));

        if (c == '\f') {  // ''
            if (isFF == true) {
                break;
            } else {
                isFF = true;
                continue;
            }
        }

        if ((c == '\n') || (c == '\r')) {
            std::string line = os.str();
            os.str("");

            if ((line.empty() == true) || (line[0] == '#')) {
                continue;
            }

            TlStringTokenizer token(line);
            std::string w1 = token.nextToken();
            std::string w2 = token.nextToken();
            std::string w3 = token.nextToken();

            if (w2 == "=") {
                param[w1] = w3;
            }
        }

        os << c;
        isFF = false;
    } while (ifs->eof() == false);

    return param;
}

std::size_t getDataSize(const std::map<std::string, std::string>& param) {
    std::map<std::string, std::string> tmp = param;

    int ndim = std::atoi(tmp["ndim"].c_str());
    std::size_t size = 0;

    if (ndim > 0) {
        size = 1;
        for (int i = 0; i < ndim; ++i) {
            const std::string key = TlUtils::format("dim%d", i + 1);
            size *= std::atol(tmp[key].c_str());
        }
    }

    return size;
}

std::string makeHeader(const std::map<std::string, std::string>& param) {
    std::map<std::string, std::string> tmp = param;

    std::string str = "# AVS field file\n";
    str += "ndim = " + tmp["ndim"] + "\n";
    const int ndim = std::atoi(tmp["ndim"].c_str());
    for (int i = 0; i < ndim; ++i) {
        const std::string key = TlUtils::format("dim%d", i + 1);
        str += TlUtils::format("%s = %s\n", key.c_str(), tmp[key].c_str());
    }

    str += "nspace = " + tmp["nspace"] + "\n";
    str += "veclen = " + tmp["veclen"] + "\n";
    str += "data = " + tmp["data"] + "\n";
    str += "field = " + tmp["field"] + "\n";
    str += "\f\f";

    return str;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "d:v");
    // const bool isVerbose = (opt["v"] == "defined");

    const std::string file1 = opt[1];
    const std::string file2 = opt[2];
    const std::string file3 = opt[3];

    std::ifstream ifs1(file1.c_str(), std::ios::in | std::ios::binary);
    //   std::ifstream ifs2(file1.c_str(), std::ios::in | std::ios::binary);
    std::map<std::string, std::string> param1 = getAvsFieldDataHeader(&ifs1);
    //   std::map<std::string, std::string> param2 =
    //   getAvsFieldDataHeader(&ifs2);

    // debug
    std::cerr << "1st FLD:" << file1 << std::endl;
    for (std::map<std::string, std::string>::const_iterator p = param1.begin();
         p != param1.end(); ++p) {
        std::cerr << p->first << " = " << p->second << std::endl;
    }
    std::cerr << std::endl;
    //   std::cerr << "2nd FLD:" << file2 << std::endl;
    //   for (std::map<std::string, std::string>::const_iterator p =
    //   param2.begin();
    //        p != param2.end(); ++p) {
    //     std::cerr << p->first << " = " << p->second << std::endl;
    //   }

    const std::size_t size1 = getDataSize(param1);
    //   const std::size_t size2 = getDataSize(param2);
    //   if (size1 != size2) {
    //     std::cerr << TlUtils::format("data size mismatch. (%ld, %ld)", size1,
    //     size2)
    //        << std::endl;
    //     return EXIT_FAILURE;
    //   }
    //   std::cerr << "number of data = " << size1 << std::endl;

    const int veclen1 = std::atoi(param1["veclen"].c_str());
    // const int veclen2 = std::atoi(param2["veclen"].c_str());
    // const int veclen3 = veclen1 + veclen2;
    const int veclen3 = 2;

    std::ofstream ofs(file3.c_str(), std::ios::out | std::ios::binary);
    std::map<std::string, std::string> param3 = param1;
    param3["veclen"] = TlUtils::xtos(veclen3);
    const std::string header = makeHeader(param3);
    ofs.write(header.c_str(), header.size());

    std::ifstream ifs2(file2.c_str(), std::ios::in);

    const int dataSize = sizeof(float);
    // data
    {
        char* pBuf = new char[dataSize];
        for (std::size_t i = 0; i < size1; ++i) {
            for (int j = 0; j < veclen1; ++j) {
                ifs1.read(pBuf, dataSize);
                ofs.write(pBuf, dataSize);
            }
            //       for (int j = 0; j < veclen2; ++j) {
            //  ifs2.read(pBuf, dataSize);
            //  ofs.write(pBuf, dataSize);
            //       }
            {
                std::string line = "";
                std::getline(ifs2, line);
                float v = std::atof(line.c_str());
                std::cerr << v << std::endl;
                ofs.write((char*)&v, sizeof(float));
            }
        }
        delete[] pBuf;
        pBuf = NULL;
    }

    // coord
    {
        std::size_t numOfItems = size1 * 3;
        float v1 = 0.0;
        // float v2 = 0.0;
        // bool isWarned = false;
        for (std::size_t i = 0; i < numOfItems; ++i) {
            ifs1.read((char*)&v1, sizeof(float));
            // ifs2.read((char*)&v2, sizeof(float));
            //       if ((isWarned == false) && (std::fabs(v1 - v2) > EPS)) {
            //  std::cerr << "warning!: 1st file coord block is different from
            //  2nd one.:"
            //        << "(" << v1 << ", " << v2 << ")"
            //        << std::endl;
            //  //isWarned = true;
            //       }
            ofs.write((char*)&v1, sizeof(float));
        }
    }

    return EXIT_SUCCESS;
}
