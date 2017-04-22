#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "dhv");

    const bool bVerbose = (opt["v"] == "defined");
    // const bool bDiagonalMode = (opt["d"] == "defined");

    std::string filepath = opt[1];
    if (bVerbose) {
        std::cerr << "loading... " << filepath << std::endl;
    }

    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
    if (ifs.fail()) {
        std::cerr << "could not found: " << filepath << std::endl;
        abort();
    }

    std::size_t size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));
    std::cout << "size = " << size << std::endl;
    
    int shellIndex1 = 0;
    int shellIndex2 = 0;
    for (std::size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&shellIndex1), sizeof(int));
        ifs.read(reinterpret_cast<char*>(&shellIndex2), sizeof(int));
        std::cout << shellIndex1 << ", " << shellIndex2 << std::endl;
    }

    ifs.close();
    
    return EXIT_SUCCESS;
}
