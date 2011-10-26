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

#ifdef __FUJITSU
#define PDF_MAIN __MAIN
#else
#define PDF_MAIN main
#endif // __FUJITSU

int PDF_MAIN(int argc, char *argv[])
{
    TlGetopt opt(argc, argv, "r");

    // setup parameters
    bool isRestart = false;
    if (opt["r"] == "defined") {
        isRestart = true;
    }
    
    ProteinDF PDF;
    if (isRestart == true) {
        PDF.restart("pdfparam.mpac");
    } else {
        PDF.run();
    }

    return EXIT_SUCCESS;
}

