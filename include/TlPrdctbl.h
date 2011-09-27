#include <string>

#ifndef TLPERIODICTABLE_H
#define TLPERIODICTABLE_H

#include <string>

class TlPrdctbl {
public:
    TlPrdctbl();
    ~TlPrdctbl();

public:
    static int getAtomicNumber(const std::string& str);
    static std::string getSymbol(int n);

private:
    static const char* periodicTable[];
};

#endif // TLPERIODICTABLE_H
