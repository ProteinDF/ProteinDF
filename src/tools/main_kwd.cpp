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
    const std::string mpacFilePath = opt["m"];

    PdfKeyword kwd;

    if (isCSV == true) {
        std::string csv = kwd.getCSV(showAll);
        std::cout << csv << std::endl;
    } else if (mpacFilePath.empty() != true) {
        const TlSerializeData data = kwd.getSerializeData();
        TlMsgPack mpac(data);
        mpac.save(mpacFilePath);
    } else {
        kwd.printDefault(std::cout);
    }

    return EXIT_SUCCESS;
}


void usage()
{
    std::cout << "output ProteinDF keyword" << std::endl;
    std::cout << std::endl;
    std::cout << "  -a:      output hidden parameters" << std::endl;
    std::cout << "  -c:      output CSV format" << std::endl;
    std::cout << "  -m FILE: output MsgPack file" << std::endl;
}
