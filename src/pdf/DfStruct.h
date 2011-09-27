#ifndef DFSTRUCT_H
#define DFSTRUCT_H

struct IJShellPair {
public:
    explicit IJShellPair(int i = 0, int j = 0) : nIShell(i), nJShell(j) {
    };

    IJShellPair(const IJShellPair& rhs) : nIShell(rhs.nIShell), nJShell(rhs.nJShell) {
    };

public:
    int nIShell;
    int nJShell;
};

#endif // DFSTRUCT_H
