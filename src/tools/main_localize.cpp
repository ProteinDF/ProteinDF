#include <iostream>
#include <string>
#include "DfLocalize.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"

int main(int argc, char *argv[])
{
    TlMsgPack mpac;
    mpac.load("pdfparam.mpac");
    TlSerializeData param = mpac.getSerializeData();
    
    DfLocalize lo(&param);

    lo.localize();

    return 0;
}




