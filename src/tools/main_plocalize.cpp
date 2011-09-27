#include <cstdlib>
#include <iostream>
#include <string>
#include "DfLocalize_Parallel.h"
#include "TlCommunicate.h"
#include "TlMsgPack.h"

int main(int argc, char *argv[])
{
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

    TlSerializeData param;
    if (rComm.isMaster() == true) {
        TlMsgPack mpac;
        mpac.load("pdfparam.mpac");
        param = mpac.getSerializeData();
    }
    rComm.broadcast(param);
    
    DfLocalize_Parallel lo(&param);
    lo.localize();

    rComm.finalize();
    return EXIT_SUCCESS;
}




