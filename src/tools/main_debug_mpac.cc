#include <cstdlib>
#include <iostream>

#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlUtils.h"

void showHelp(const std::string& name) {
    std::cout << TlUtils::format("%s <MSGPACK_FILE_PATH>", name.c_str()) << std::endl;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hv");

    // parameters - common
    if ((opt.getCount() < 2) || (opt["h"] == "defined")) {
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }

    const bool isVerbose = (opt["v"] == "defined");

    const std::string mpacPath = opt[1];
    std::cout << "msgpack file path: " << mpacPath << std::endl;

    TlMsgPack mpac;
    mpac.load(mpacPath);
    const TlSerializeData data = mpac.getSerializeData();
    std::cout << data.str() << std::endl;

    return EXIT_SUCCESS;
}
