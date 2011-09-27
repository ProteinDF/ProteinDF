#include <unistd.h>
#include <cstdlib>
#include <string>
#include <iostream>

#include "TlCommunicate.h"

#define HOSTNAME_LEN 256

int main(int argc, char *argv[])
{
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

    char* pHostName = new char[HOSTNAME_LEN];
    (void)gethostname(pHostName, HOSTNAME_LEN);
    const std::string hostname(pHostName);
    delete []pHostName;
    pHostName = NULL;

    for (int i = 0; i < rComm.getNumOfProc(); ++i) {
        if (i == rComm.getRank()) {
            std::cout << "[" << i << "] "
                      << hostname
                      << std::endl;
        }
    }

    rComm.finalize();
    return EXIT_SUCCESS;
}

