#include <cassert>
#include <iostream>
#include <cstdlib>
#include <string>

#include "PdfKeyword.h"
#include "TlGetopt.h"
#include "TlSerializeData.h"
#include "TlMsgPack.h"

void usage();

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "acjm:");
    const bool showAll = (opt["a"] == "defined");
    const bool isCSV = (opt["c"] == "defined");
    const bool isJapanese = (opt["j"] == "defined");
    const std::string mpacFilePath = opt["m"];

    PdfKeyword kwd;

    if (isCSV == true) {
        if (isJapanese == true) {
            std::string csv = kwd.getCSV_J(showAll);
            std::cout << csv << std::endl;
        } else {
            std::string csv = kwd.getCSV(showAll);
            std::cout << csv << std::endl;
        }
    } else {
        kwd.printDefault(std::cout);
    }

    if (mpacFilePath.empty() != true) {
        const TlSerializeData data = kwd.getSerializeData();
        TlMsgPack mpac(data);
        mpac.save(mpacFilePath);
    }
    
    return EXIT_SUCCESS;
}


void usage()
{
    std::cout << "calculation Mulliken Population." << std::endl;
    std::cout << "usage: " << std::endl;
    std::cout << "prog [options] P_matrix_path S_matrix_path" << std::endl;
}

