#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <cstdlib>
#include <ios>
#include <iostream>
#include <fstream>

#include "ProteinDF.h"
#include "TlGetopt.h"
#include "TlSerializeData.h"
#include "TlLogging.h"

#ifdef __FUJITSU
#define PDF_MAIN MAIN__
#else
#define PDF_MAIN main
#endif // __FUJITSU

int PDF_MAIN(int argc, char *argv[])
{
    // setup parameters
    TlGetopt opt(argc, argv, "dro:");

    bool isRestart = false;
    if (opt["r"] == "defined") {
        isRestart = true;
    }

    TlLogging& log = TlLogging::getInstance();
    std::string output = "fl_Out_Std";
    if (opt["o"].empty() != true) {
        output = opt["o"];
    }
    log.setFilePath(output);

    if (opt["d"] == "defined") {
        log.setLevel(TlLogging::DEBUG);
    }

    // do ProteinDF
    ProteinDF PDF;
    if (isRestart == true) {
        PDF.restart("pdfparam.mpac");
    } else {
        PDF.run();
    }

    return EXIT_SUCCESS;
}

