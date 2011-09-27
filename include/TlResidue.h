#ifndef CLASS_TLRESIDUE
#define CLASS_TLRESIDUE

#include <string>

class TlResidue {
public:
    TlResidue();
    ~TlResidue() {
    }

    // return number of amino acid date on information table
    int number();

    // return p th amino acid label on information table
    std::string label(int p);

    static bool isResidue(const std::string& str);

private:
    static const char* residueTable[];
};

#endif

