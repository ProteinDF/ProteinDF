#ifndef TLINTVECTOR_H
#define TLINTVECTOR_H

#include <vector>

class TlIntVector {
public:
    TlIntVector(const std::vector<int>& v);
    ~TlIntVector();

public:
    void save(const std::string& path);
    void load(const std::string& path);
};

#endif // TLINTVECTOR_H
